#!/usr/bin/env python3
"""Build R14 float64 geometry and serialize standalone parts in a fresh output root.

Whole-body output is owned by reproduce_body_r14.py and its explicit saved recipe.
This entrypoint never qualifies printing or physical flow.
"""
from __future__ import annotations
import argparse, gc, hashlib, json, os, time
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import manifold3d as mf
import trimesh
from geometry_r14 import read_glb, closed_group, save_npz, load_npz, unite
from rail_components_r14 import build_line

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


def export_parts(root, cache, data):
    # Keep the saved ENU-to-print transform; do not infer a new pose from a rebuild.
    from export_repaired_parts import export_one
    transform=np.asarray(json.loads((root/'data/coordinate_systems_R14.json').read_text())[
        'geographic_ENU_mm_to_print_XYZ_mm_row_major'],dtype=float)
    reports=[]
    for line in data['lines']:
        lid=line['line_id']
        part=load_npz(cache/f'{lid}_outer.npz')-load_npz(cache/f'{lid}_core.npz')
        row=export_one(part,root/'models'/f'rail_R14_{lid}.stl',transform)
        row.update(line_id=lid,name_ja=line['name_ja'],standalone_line_reference=True,
                   installed_lumen_validation_uses_whole_body=True)
        reports.append(row)
        dump(root/'reports/repaired_part_exports_R14.json',reports)
        del part;gc.collect()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True,help='R14 working output root; use a fresh copy of inputs')
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--stage',choices=['rails','base','assemble','parts'],required=True,
                   help='No whole-body export stage: use reproduce_body_r14.py with the saved repair recipe')
    a=p.parse_args();root=a.root.resolve();cache=a.cache.resolve()
    cache.mkdir(parents=True,exist_ok=True);(root/'reports').mkdir(parents=True,exist_ok=True)
    data=json.loads((root/'data/rail_paths_R14.json').read_text())
    ports=json.loads((root/'data/ports_R14.json').read_text())['ports']
    if len(data['lines'])!=14 or len(ports)!=69:
        raise ValueError('R14 input contract requires 14 lines and 69 ports')
    started=time.time()
    if a.stage=='rails':rails(root,cache,data,ports)
    elif a.stage=='base':base(root,cache)
    elif a.stage=='assemble':assemble(root,cache,data)
    else:export_parts(root,cache,data)
    print('COMPLETE',a.stage,time.time()-started,flush=True)


if __name__=='__main__':main()
