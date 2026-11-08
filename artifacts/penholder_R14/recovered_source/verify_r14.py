#!/usr/bin/env python3
"""Independent tests over the delivered STL surface and saved fluid cores.

Self-intersections use MeshLab/VCGLib's non-adjacent-triangle selection.
Finite passage rays and radial thickness samples are not global proofs of
printability or a substitute for resin washing, leakage and pigment tests.
"""
import argparse, gc, hashlib, json, os, time
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import trimesh
import manifold3d as mf
import vtk
from mesh_ray_checks_r13 import locator_for,first_hit,basis,ray_winding,probe_rays


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
def load_core(p):
    with np.load(p) as d:
        return mf.Manifold(mf.Mesh64(d['vertices'].astype(np.float64),d['faces'].astype(np.uint64)))


def design_mesh(root):
    path=root/'models/subway_R14_hollow_print.stl'
    m=trimesh.load_mesh(path,process=True)
    matrix=np.asarray(json.loads((root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    m.apply_transform(np.linalg.inv(matrix))
    return path,m


def self_intersections(root):
    import pymeshlab
    rows=[]
    for p in sorted((root/'models').glob('*.stl')):
        t=time.time();ms=pymeshlab.MeshSet()
        welded=trimesh.load_mesh(p,process=True)
        ms.add_mesh(pymeshlab.Mesh(welded.vertices,welded.faces))
        ms.compute_selection_by_self_intersections_per_face()
        count=int(ms.current_mesh().selected_face_number())
        row=dict(file=p.name,sha256=digest(p),triangles=ms.current_mesh().face_number(),
            self_intersecting_faces=count,passed=count==0,seconds=time.time()-t,
            method='MeshLab/VCGLib self-intersecting non-adjacent triangle selection; no geometry modification')
        rows.append(row);print('self intersections',p.name,count,flush=True)
        dump(root/'reports/self_intersections_R14.json',dict(passed=all(x['passed'] for x in rows),meshes=rows,
            coplanar_and_numerical_contact_limitations='VCGLib floating-point predicates; not an exact-arithmetic proof'))
        del ms;gc.collect()


def connectivity(root,cache):
    data=json.loads((root/'data/rail_paths_R14.json').read_text());ports=json.loads((root/'data/ports_R14.json').read_text())['ports']
    rows=[];pairs=[];minimum=2.
    for i,line in enumerate(data['lines']):
        lid=line['line_id'];core=load_core(cache/f'{lid}_core.npz')
        components=core.decompose();positive=[q for q in components if q.volume()>1e-8]
        row=dict(line_id=lid,name_ja=line['name_ja'],core_sha256=digest(cache/f'{lid}_core.npz'),
            positive_fluid_components=len(positive),nonpositive_boundary_components=len(components)-len(positive),
            ports=sum(p['line_id']==lid for p in ports),port_membership=[])
        # A signed-distance test on the independent component checks both ends
        # of every nozzle; a port record alone is not connectivity evidence.
        from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray
        if len(positive)==1:
            raw=positive[0].to_mesh64();v=np.asarray(raw.vert_properties)[:,:3];f=np.asarray(raw.tri_verts)
            locator=locator_for(trimesh.Trimesh(v,f,process=False));poly=locator.GetDataSet()
            signed=vtk.vtkImplicitPolyDataDistance();signed.SetInput(poly)
            for p in ports:
                if p['line_id']!=lid:continue
                a=float(signed.EvaluateFunction(p['attachment_enu_mm']));b=float(signed.EvaluateFunction(p['mouth_center_enu_mm']))
                row['port_membership'].append(dict(port_id=p['port_id'],attachment_signed_mm=a,mouth_signed_mm=b,passed=a<0 and b<0))
            del locator,poly,signed,raw,v,f
        row['passed']=len(positive)==1 and row['ports']>=2 and len(row['port_membership'])==row['ports'] and all(p['passed'] for p in row['port_membership'])
        rows.append(row)
        for other in data['lines'][i+1:]:
            c=load_core(cache/f"{other['line_id']}_core.npz")
            gap=float(core.min_gap(c,2.));minimum=min(minimum,gap)
            pairs.append(dict(line_a=lid,line_b=other['line_id'],capped_search_distance_mm=2.,gap_mm=gap,passed=gap>=.5))
            del c;gc.collect()
        print('core',lid,'components',len(positive),'ports',row['ports'],flush=True)
        dump(root/'reports/lumen_connectivity_R14.json',dict(passed=all(r['passed'] for r in rows) and all(p['passed'] for p in pairs),
            lines=rows,pairs=pairs,minimum_interline_lumen_gap_mm=minimum,checked_all_pairs=len(pairs)==91,
            nominal_lumen_diameter_mm=1.2,physical_flow_test_performed=False))
        del core,positive,components;gc.collect()


def passage_and_walls(root):
    path,mesh=design_mesh(root);loc=locator_for(mesh)
    data=json.loads((root/'data/rail_paths_R14.json').read_text());ports=json.loads((root/'data/ports_R14.json').read_text())['ports']
    lines=[];port_rows=[]
    for line in data['lines']:
        fail=[];ray_count=0;wall=[];radii=[];wall_missing=0
        for trace in line['paths']:
            p=np.asarray(trace['points_enu_mm'],float)
            for k,(a,b) in enumerate(zip(p[:-1],p[1:])):
                delta=b-a;length=np.linalg.norm(delta)
                if length<1e-10:continue
                axis=delta/length;n,q=basis(axis)
                for off in [np.zeros(3),.45*n,-.45*n,.45*q,-.45*q]:
                    hit=first_hit(loc,a+off,b+off);ray_count+=1
                    if hit['hit']:fail.append(dict(trace_id=trace['trace_id'],segment=k,offset=off.tolist(),**hit))
                # Sample each original segment midpoint, up to 12 directions.
                # All source bends are therefore sampled, including short ones.
                center=(a+b)/2
                for angle in np.arange(12)*2*np.pi/12:
                    direction=np.cos(angle)*n+np.sin(angle)*q
                    inner=first_hit(loc,center,center+3.5*direction)
                    if not inner['hit']:wall_missing+=1;continue
                    r=float(np.linalg.norm(np.array(inner['point_enu_mm'])-center))
                    radii.append(r)
                    # Retain the first interval of solid after crossing the bore.
                    outer=first_hit(loc,center+(r+.0001)*direction,center+5.*direction)
                    if outer['hit']:
                        outer_r=float(np.linalg.norm(np.array(outer['point_enu_mm'])-center))
                        if outer_r>r:wall.append(outer_r-r)
                    else:wall_missing+=1
        row=dict(line_id=line['line_id'],name_ja=line['name_ja'],route_rays=ray_count,
            offset_radius_mm=.45,obstruction_hits=fail,all_passage_rays_clear=not fail,
            radial_direction_samples=len(radii),minimum_sampled_bore_radius_mm=min(radii) if radii else None,
            minimum_sampled_radial_wall_mm=min(wall) if wall else None,
            wall_samples_with_no_exit_within_5mm=wall_missing,
            sampled_minimum_is_not_a_global_wall_thickness_bound=True)
        lines.append(row);print('rays',line['line_id'],ray_count,'blocked',len(fail),'wall',row['minimum_sampled_radial_wall_mm'],flush=True)
        dump(root/'reports/saved_body_passage_R14.json',dict(body_sha256=digest(path),lines=lines,ports=port_rows,in_progress=True))
    for p in ports:
        a=np.asarray(p['attachment_enu_mm']);mouth=np.asarray(p['mouth_center_enu_mm']);axis=np.asarray(p['axis_enu'])
        internal=probe_rays(loc,mouth,axis,.45,True,a,p['length_mm'],mesh)
        external=probe_rays(loc,mouth,axis,.8,False,body=mesh)
        n,q=basis(axis);checks=[]
        for angle in np.arange(8)*2*np.pi/8:
            checks.append(ray_winding(loc,mesh,mouth+.1*axis+.8*(np.cos(angle)*n+np.sin(angle)*q)))
        port_rows.append(dict(port_id=p['port_id'],line_id=p['line_id'],
            internal_rays=internal,external_rays=external,mouth_empty_checks=checks,
            passed=not any(x['hit'] for x in internal+external) and all(x['empty_by_oriented_ray'] for x in checks)))
    dump(root/'reports/saved_body_passage_R14.json',dict(body_sha256=digest(path),
        body_file=path.name,generator_imported=False,lines=lines,ports=port_rows,
        passed=all(l['all_passage_rays_clear'] for l in lines) and all(p['passed'] for p in port_rows),
        route_rays=sum(l['route_rays'] for l in lines),port_count=len(port_rows),
        minimum_sampled_radial_wall_mm=min(l['minimum_sampled_radial_wall_mm'] for l in lines),
        minimum_sampled_bore_radius_mm=min(l['minimum_sampled_bore_radius_mm'] for l in lines),
        limitation='Finite triangle-intersection rays, 12 radial directions at each original segment midpoint; sampled wall lengths are not a global minimum proof. Photopolymer process and pigment flow are untested.'))


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--cache',type=Path)
    ap.add_argument('--stage',choices=['self','connectivity','passage'],required=True);a=ap.parse_args();root=a.root.resolve()
    if a.stage=='self':self_intersections(root)
    elif a.stage=='connectivity':connectivity(root,a.cache.resolve())
    elif a.stage=='passage':passage_and_walls(root)

if __name__=='__main__':main()
