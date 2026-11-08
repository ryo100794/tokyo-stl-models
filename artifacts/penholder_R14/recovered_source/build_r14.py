#!/usr/bin/env python3
"""Offline R14 rebuild from saved R13 paths, ports and the R12 landscape.

Run stages sequentially: rails, base, assemble, export. Caches are outside the
deliverable. All design computation uses ENU millimetres and Mesh64. STL is
quantized by its format to float32 and is independently reloaded before use.
"""
from __future__ import annotations
import argparse, gc, hashlib, json, os, sys, time
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import manifold3d as mf
import trimesh
from geometry_r14 import (read_glb, closed_group, save_npz, load_npz, unite,
                         arrays, capsule_path, object_group, write_glb)
from rail_components_r14 import build_line, color_for


def dump(p,d):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()


def read_arrays(obj):
    m=obj.to_mesh64()
    return np.asarray(m.vert_properties)[:,:3].copy(),np.asarray(m.tri_verts).copy()


def topology(mesh):
    edges=mesh.edges_sorted
    _,cnt=np.unique(edges,axis=0,return_counts=True)
    # The link of every vertex of a closed 2-manifold is a single cycle.
    # Edge closedness alone misses point-contact singularities.
    f=np.asarray(mesh.faces)
    vf=mesh.vertex_faces
    bad_link=0
    for v,row in enumerate(vf):
        row=row[row>=0]
        if not len(row):continue
        tris=f[row];other=tris[tris!=v].reshape(-1,2)
        ns,ic,degree=np.unique(other,return_inverse=True,return_counts=True)
        if np.any(degree!=2):bad_link+=1;continue
        pairs=ic.reshape(-1,2);nbr={int(x):[] for x in range(len(ns))}
        for a,b in pairs:nbr[int(a)].append(int(b));nbr[int(b)].append(int(a))
        seen=set();todo=[0]
        while todo:
            q=todo.pop()
            if q in seen:continue
            seen.add(q);todo.extend(nbr[q])
        if len(seen)!=len(ns):bad_link+=1
    area=np.linalg.norm(np.cross(mesh.triangles[:,1]-mesh.triangles[:,0],
                                 mesh.triangles[:,2]-mesh.triangles[:,0]),axis=1)/2
    return dict(vertices=len(mesh.vertices),triangles=len(mesh.faces),
        watertight=bool(mesh.is_watertight),winding_consistent=bool(mesh.is_winding_consistent),
        boundary_edges=int((cnt==1).sum()),nonmanifold_edges=int((cnt>2).sum()),
        nonmanifold_vertex_links=bad_link,zero_area_triangles=int((area==0).sum()),
        duplicate_faces=int(len(f)-len(np.unique(np.sort(f,axis=1),axis=0))),
        volume_mm3=float(mesh.volume),bounds_mm=mesh.bounds.tolist(),dimensions_mm=mesh.extents.tolist())


def strict_basic(m):
    return bool(m.is_watertight and m.is_winding_consistent and m.volume>0
                and np.all(m.area_faces>0) and np.all(np.isfinite(m.vertices)))


def export_stl(obj,path,transform,tolerances=(.002,.005,.01,.02),repair_numeric=False):
    """Bounded simplification followed by real STL float32 serialization.

    No holes are filled. Failed candidates are never placed in models/.
    Tiny contact defects are removed by Manifold's bounded simplification;
    final coordinate-welded edge and vertex-link tests are independent.
    """
    trials=[]
    for tolerance in tolerances:
        solid=obj.simplify(tolerance)
        v,f=read_arrays(solid);v=v@transform[:3,:3].T+transform[:3,3]
        m=trimesh.Trimesh(v,f,process=False)
        raw=trimesh.exchange.stl.export_stl(m)
        import io
        saved=trimesh.load_mesh(io.BytesIO(raw),file_type='stl',process=True)
        before=len(saved.faces)
        saved.update_faces(saved.nondegenerate_faces(height=1e-10))
        saved.update_faces(saved.unique_faces());saved.remove_unreferenced_vertices()
        repair_used=False
        intersects=0
        if repair_numeric and strict_basic(saved):
            import pymeshlab
            check=pymeshlab.MeshSet();check.add_mesh(pymeshlab.Mesh(saved.vertices,saved.faces))
            check.compute_selection_by_self_intersections_per_face()
            intersects=int(check.current_mesh().selected_face_number())
            del check
        if repair_numeric and (not strict_basic(saved) or intersects):
            import pymeshfix
            fix=pymeshfix.MeshFix(saved.vertices,saved.faces,verbose=False)
            fix.repair(joincomp=False,remove_smallest_components=False)
            saved=trimesh.Trimesh(fix.points,fix.faces,process=True)
            repair_used=True
        trials.append(dict(tolerance_mm=tolerance,triangles=len(saved.faces),
            removed_faces=before-len(saved.faces),basic_passed=strict_basic(saved)))
        if not trials[-1]['basic_passed']:
            del saved,m,v,f,solid;gc.collect();continue
        # Re-write cleaned coordinates and reload again; the delivered bytes
        # must pass, rather than only the pre-serialization mesh.
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(trimesh.exchange.stl.export_stl(saved))
        del saved,m,v,f,solid;gc.collect()
        saved=trimesh.load_mesh(path,process=True)
        info=topology(saved)
        if info['nonmanifold_vertex_links'] or not strict_basic(saved):
            path.unlink();del saved;gc.collect();continue
        info.update(file=path.name,sha256=sha(path),bytes=path.stat().st_size,
                    simplification_tolerance_mm=tolerance,attempts=trials,
                    hole_filling_performed=repair_used,
                    numeric_repair='MeshFix; joincomp=False; remove_smallest_components=False' if repair_used else None,
                    self_intersecting_faces_before_repair=intersects,
                    format='binary STL',units='mm')
        return info
    raise RuntimeError(f'No acceptable STL candidate for {path.name}: {trials}')


def rails(root,cache,data,ports):
    report=[]
    for line in data['lines']:
        outer,core=build_line(line,ports,cache)
        # Persist the independent cavities with float64 coordinates.
        report.append(dict(line_id=line['line_id'],name_ja=line['name_ja'],
            outer_triangles=outer.num_tri(),core_triangles=core.num_tri(),
            core_volume_mm3=core.volume(),core_status=str(core.status())))
        del outer,core;gc.collect()
        dump(root/'reports/rail_build_R14.json',report)


def base(root,cache):
    groups,_,reference_sha=read_glb(root/'reference/subway_R12_lite.glb')
    records=[]
    for i,g in enumerate(groups):
        if 5<=i<=18:continue
        p=cache/f'base_{i:02d}.npz'
        if p.exists():obj=load_npz(p)
        else:
            print('base',i,g['name'],flush=True)
            margin=.006 if i in [0,1,3,4] else 0.
            obj=closed_group(g,contact_margin_mm=margin).simplify(.002)
            save_npz(obj,p)
        records.append(dict(index=i,name=g['name'],triangles=obj.num_tri(),volume_mm3=obj.volume(),
            inherited_tangent_primitive_vertex_normal_expansion_mm=.006 if i in [0,1,3,4] else 0.))
        del obj;gc.collect()
        dump(root/'reports/base_build_R14.json',dict(reference_sha256=reference_sha,groups=records))


def assemble(root,cache,data):
    # Work in two batches, saving float64 checkpoints to keep memory bounded.
    all_core=unite([load_npz(cache/f"{l['line_id']}_core.npz") for l in data['lines']])
    save_npz(all_core,cache/'all_cores.npz')
    del all_core;gc.collect()
    bp=cache/'base_union.npz'
    if not bp.exists():
        print('union base',flush=True)
        base_obj=unite([load_npz(p) for p in sorted(cache.glob('base_[0-9][0-9].npz'))])
        save_npz(base_obj,bp);del base_obj;gc.collect()
    pp=cache/'positive.npz'
    if not pp.exists():
        print('union rails with base',flush=True)
        positive=unite([load_npz(bp)]+[load_npz(cache/f"{l['line_id']}_outer.npz") for l in data['lines']])
        save_npz(positive,pp);del positive;gc.collect()
    print('subtract all cavities',flush=True)
    final=load_npz(pp)-load_npz(cache/'all_cores.npz')
    final=final.simplify(.002)
    assert final.status()==mf.Error.NoError
    save_npz(final,cache/'full_hollow.npz')
    dump(root/'reports/assembly_R14.json',dict(triangles=final.num_tri(),
        volume_mm3=final.volume(),status=str(final.status()),float64_checkpoints=True))


def export(root,cache,data,ports):
    obj=load_npz(cache/'full_hollow.npz')
    v,_=read_arrays(obj);highest=float(v[:,2].max());del v;gc.collect()
    transform=np.diag([1.,-1.,-1.,1.]);transform[2,3]=highest
    coords=dict(revision='R14',design_units='mm',design_axes=['east','north','up'],
        print_units='mm',rotation_about_X_degrees=180.,linear_determinant=1.,
        geographic_ENU_mm_to_print_XYZ_mm_row_major=transform.tolist(),
        print_transform_baked_into_vertex_coordinates=True,
        horizontal_scale=120000,vertical_scale=1800,height_multiplier_reapplied=False,
        route_source='R13 adjusted ENU-mm paths; no new survey or geographic rescaling',
        shared_coordinates_module='source/model_coordinates.py',
        shared_coordinates_applied_as_identity_to_already_scaled_input=True)
    from model_coordinates import ModelCoordinates
    identity=ModelCoordinates(0.,0.,mm_per_m=1000.,crs='LOCAL:inherited-model-metres')
    test=np.array([[1.23,-4.56,7.89]])
    assert np.allclose(identity.points(test/1000),test,atol=1e-12)
    dump(root/'data/coordinate_systems_R14.json',coords)
    reports=[]
    print('export whole body',flush=True)
    info=export_stl(obj,root/'models/subway_R14_hollow_print.stl',transform,tolerances=(.02,.03,.04))
    reports.append(info);dump(root/'reports/saved_meshes_R14.json',reports)
    del obj;gc.collect()
    # The optional indexed 3MF is built from the same reloaded STL, so they
    # share exactly the same float32 physical surface and print orientation.
    mesh=trimesh.load_mesh(root/'models/subway_R14_hollow_print.stl',process=True)
    mesh.units='mm'
    path=root/'models/subway_R14_hollow_print.3mf'
    path.write_bytes(trimesh.exchange.threemf.export_3MF(mesh,compresslevel=6))
    reports.append(dict(file=path.name,sha256=sha(path),bytes=path.stat().st_size,
                        format='3MF',same_surface_as_STL=True))
    del mesh;gc.collect()
    dump(root/'reports/saved_meshes_R14.json',reports)
    for line in data['lines']:
        lid=line['line_id'];print('export line',lid,flush=True)
        # A separate line part contains its own lumen. Foreign line cutters
        # are only needed in the assembled model, not in a stand-alone line.
        part=load_npz(cache/f'{lid}_outer.npz')-load_npz(cache/f'{lid}_core.npz')
        save_npz(part,cache/f'{lid}_hollow.npz')
        info=export_stl(part,root/'models'/f'rail_R14_{lid}.stl',transform)
        info.update(line_id=lid,name_ja=line['name_ja'],standalone_line_reference=True,
            installed_lumen_validation_uses_whole_body=True)
        reports.append(info);dump(root/'reports/saved_meshes_R14.json',reports)
        del part;gc.collect()
    print('export landscape part',flush=True)
    landscape=load_npz(cache/'base_union.npz')-load_npz(cache/'all_cores.npz')
    reports.append(export_stl(landscape,root/'models/landscape_R14.stl',transform))
    dump(root/'reports/saved_meshes_R14.json',reports)


def export_parts(root,cache,data):
    with np.load(cache/'full_hollow.npz') as saved:highest=float(saved['vertices'][:,2].max())
    transform=np.diag([1.,-1.,-1.,1.]);transform[2,3]=highest
    reports=[]
    for line in data['lines']:
        lid=line['line_id'];print('export separate line',lid,flush=True)
        part=load_npz(cache/f'{lid}_outer.npz')-load_npz(cache/f'{lid}_core.npz')
        save_npz(part,cache/f'{lid}_hollow.npz')
        info=export_stl(part,root/'models'/f'rail_R14_{lid}.stl',transform,tolerances=(.01,.02),repair_numeric=True)
        info.update(line_id=lid,name_ja=line['name_ja'],standalone_line_reference=True,
            installed_lumen_validation_uses_whole_body=True)
        reports.append(info);dump(root/'reports/individual_rails_R14.json',reports)
        del part;gc.collect()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--stage',choices=['rails','base','assemble','export','parts'],required=True)
    a=p.parse_args();root=a.root.resolve();cache=a.cache.resolve();cache.mkdir(parents=True,exist_ok=True)
    data=json.loads((root/'data/rail_paths_R14.json').read_text())
    ports=json.loads((root/'data/ports_R14.json').read_text())['ports']
    assert len(data['lines'])==14 and len(ports)==69
    started=time.time()
    if a.stage=='rails':rails(root,cache,data,ports)
    elif a.stage=='base':base(root,cache)
    elif a.stage=='assemble':assemble(root,cache,data)
    elif a.stage=='export':export(root,cache,data,ports)
    elif a.stage=='parts':export_parts(root,cache,data)
    print('COMPLETE',a.stage,time.time()-started,flush=True)


if __name__=='__main__':main()
