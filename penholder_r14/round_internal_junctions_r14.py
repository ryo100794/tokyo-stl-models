#!/usr/bin/env python3
"""Bounded local removal of solid cusps between cavities of the SAME R14 line.

Inputs remain untouched. This is a geometric trial, not print qualification.
The old core is retained by union. A local cavity closing rounds thin solid
junction tips. Each added cutter is checked against the unchanged outer
positive and ALL other line cores. Saved STL topology and self intersections
are rechecked. Source centerline coordinates and port coordinates never move.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,hashlib,json,time,sys
from pathlib import Path
import numpy as np,manifold3d as mf,trimesh,pymeshlab as ml
from scipy.spatial import cKDTree
from verify_saved_wall_regions import audit,locator_for,first_hit


def clusters(points,link_mm=2.3):
    p=np.array(points);parent=list(range(len(p)))
    def find(a):
        while a!=parent[a]:parent[a]=parent[parent[a]];a=parent[a]
        return a
    for a,b in cKDTree(p).query_pairs(link_mm):
        a=find(a);b=find(b);parent[a]=b
    groups={}
    for i in range(len(p)):groups.setdefault(find(i),[]).append(i)
    return list(groups.values())


def msh(obj):
    v,f=g.arrays(obj);return trimesh.Trimesh(v,f,process=False)


def critical_rays(path,T,examples):
    mm=trimesh.load_mesh(path,process=True);mm.apply_transform(np.linalg.inv(T));loc=locator_for(mm);rows=[]
    for ex in examples:
        c=np.array(ex['center_mm']);d=np.array(ex['direction']);x=c.copy();hits=[]
        for i in range(10):
            h=first_hit(loc,x,c+20*d)
            if h is None:break
            q,cell=h;hits.append({'distance_from_center_mm':float(np.linalg.norm(q-c)),'point_enu_mm':q.tolist(),'face_id':cell,'face_vertices_enu_mm':mm.triangles[cell].tolist(),'face_normal_enu':mm.face_normals[cell].tolist()});x=q+1e-4*d
        rows.append({'trace_id':ex['trace_id'],'segment':ex['segment'],'old_interval_mm':ex['thickness_mm'],'center_mm':c.tolist(),'direction':d.tolist(),'hits':hits})
    return {'file':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'samples':rows}


def local_sharp_edges(path,T,line,sites):
    """Find geometrically sharp convex material ridges facing same-line fluid.

    Surface kind is labeled by orientation relative to the saved nearest axis.
    This is a diagnostic neighborhood classifier, NOT a manufacturing limit.
    """
    mm=trimesh.load_mesh(path,process=True);mm.apply_transform(np.linalg.inv(T))
    allp=[np.array(x['points_enu_mm']) for x in line['paths']]
    starts=np.concatenate([p[:-1] for p in allp]);delta=np.concatenate([np.diff(p,axis=0) for p in allp]);den=np.maximum((delta*delta).sum(1),1e-24)
    ad=mm.face_adjacency;edge=mm.face_adjacency_edges;mid=mm.vertices[edge].mean(axis=1);angles=mm.face_adjacency_angles;convex=mm.face_adjacency_convex
    candidate=np.zeros(len(ad),dtype=bool)
    for lo,hi in sites:candidate|=np.all((mid>=lo)&(mid<=hi),axis=1)
    ids=np.flatnonzero(candidate & convex & (angles>np.radians(60)))
    entries=[]
    def is_inner(cid):
        q=mm.triangles_center[cid];u=np.clip(((q-starts)*delta).sum(1)/den,0,1);v=q-starts-u[:,None]*delta;j=np.argmin((v*v).sum(1));dist=np.linalg.norm(v[j]);return float(mm.face_normals[cid]@v[j]/max(dist,1e-20))<-.15
    for i in ids:
        a,b=ad[i]
        if not is_inner(a) or not is_inner(b):continue
        entries.append({'edge_vertex_ids':edge[i].tolist(),'edge_vertices_enu_mm':mm.vertices[edge[i]].tolist(),'face_ids':[int(a),int(b)],'normal_angle_degrees':float(np.degrees(angles[i])),'material_included_angle_degrees':float(180-np.degrees(angles[i])),'edge_length_mm':float(np.linalg.norm(np.diff(mm.vertices[edge[i]],axis=0)))})
    return {'method':'Convex material edges in patch bounds, both adjacent facets oriented toward nearest same-line axes, face-normal angle over 60 deg. Threshold is a diagnostic setting, not printer acceptance.','convex_inner_edges_over_60_degrees':len(entries),'convex_inner_edges_over_120_degrees':sum(x['normal_angle_degrees']>120 for x in entries),'maximum_normal_angle_degrees':max([x['normal_angle_degrees'] for x in entries],default=0.),'sharp_edges':entries}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--cache',type=Path,required=True);ap.add_argument('--diagnostic',type=Path,help='defaults to ROOT/data/internal_cusp_seed_R14.json');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--lines',nargs='*');ap.add_argument('--ball-inradius',type=float,default=.25);ap.add_argument('--patch-inset',type=float,default=.0005);ap.add_argument('--outer-guard',type=float,default=.52);a=ap.parse_args()
    if a.out.resolve()==a.cache.resolve():
        raise ValueError('The output directory must differ from the input hull cache')
    if a.diagnostic is None:a.diagnostic=a.root/'data/internal_cusp_seed_R14.json'
    # Runtime helpers come from the selected R14 root, never a fixed workspace.
    sys.path.insert(0,str((a.root/'source').resolve()))
    global g,read_stl_bytes,weld_near_vertices,topology
    import geometry_r14 as g
    from robust_export_r14 import read_stl_bytes,weld_near_vertices
    from build_r14 import topology
    a.out.mkdir(exist_ok=True,parents=True)
    data=json.loads((a.root/'data/rail_paths_R14.json').read_text());ports=json.loads((a.root/'data/ports_R14.json').read_text())['ports'];T=np.array(json.loads((a.root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major']);diag=json.loads(a.diagnostic.read_text())
    for name,expected in diag.get('source_data_sha256',{}).items():
        actual=hashlib.sha256((a.root/'data'/name).read_bytes()).hexdigest()
        if actual!=expected:raise ValueError(f'R14 source data differ from the saved cusp seed: {name}')
    k=mf.Manifold.sphere(1.,16);kt=msh(k);inradius=float(np.min(np.abs(np.einsum('ij,ij->i',kt.face_normals,kt.vertices[kt.faces[:,0]]))));radius=a.ball_inradius/inradius;kernel=mf.Manifold.sphere(radius,16)
    erode_kernel=mf.Manifold.sphere((a.ball_inradius+a.patch_inset)/inradius,16)
    guard_kernel=mf.Manifold.sphere(a.outer_guard/inradius,16)
    result={'method':'same-line local polyhedral morphological closing, then limited to the outer-positive erosion; retain the complete original core by union','ball_kernel_sides':16,'ball_kernel_unit_inradius':inradius,'ball_kernel_circumradius_mm':radius,'ball_kernel_inscribed_radius_mm':a.ball_inradius,'patch_inset_mm':a.patch_inset,'outer_positive_erosion_ball_inradius_mm':a.outer_guard,'local_cleaning_simplify_mm':1e-5,'saved_simplify_mm':.001,'minimum_geometric_guard_mm':.5,'centerlines_moved':False,'port_coordinates_moved':False,'outer_positive_changed':False,'physical_tests_performed':False,'global_internal_minimum_wall_proven':False,'lines':[]}
    result['source_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['seed_sha256']=hashlib.sha256(a.diagnostic.read_bytes()).hexdigest()
    result['source_data_sha256']=diag.get('source_data_sha256',{})
    original={lid:g.load_npz(a.cache/f'{lid}_core.npz') for lid in [f'L{i:02d}' for i in range(1,15)]};final={**original}
    for dr in diag['lines']:
        lid=dr['line_id'];ex=dr['examples']['septum']
        if not ex or a.lines and lid not in a.lines:continue
        begin=time.time();line=next(x for x in data['lines'] if x['line_id']==lid);outer=g.load_npz(a.cache/f'{lid}_outer.npz');core=original[lid];newcore=core;before_solid=outer-core;before_components=before_solid.decompose();ext=mf.Manifold.cube([1000,1000,1000],center=True)-outer
        p=np.array([(np.array(e['entry_mm'])+e['exit_mm'])/2 for e in ex]);sites=[];site_rows=[]
        for k,group in enumerate(clusters(p)):
            lo=p[group].min(axis=0)-2.;hi=p[group].max(axis=0)+2.;roi=mf.Manifold.cube(hi-lo).translate(lo);sites.append((lo,hi));local=core^roi
            print(lid,'site',k,'local triangles',local.num_tri(),flush=True)
            expanded=local.minkowski_sum(kernel)
            # A 0.0005 mm additional erosion puts the patch strictly inside the
            # old core wherever the closing did not add meaningful cavity.
            # This avoids near-coincident duplicate faces in the union. The
            # original core is never eroded or narrowed.
            closed=expanded.minkowski_difference(erode_kernel)
            clean=closed.as_original().simplify(1e-5)
            # A local erosion of the unchanged positive is a conservative
            # allowed domain for new cuts. The kernel contains a true sphere
            # of radius outer_guard, so all accepted patch points remain at
            # least that far from the positive's outside, including near a
            # planned mouth. Existing mouth cutters stay in original core.
            guard_roi=mf.Manifold.cube(hi-lo+2).translate(lo-1)
            safe=(outer^guard_roi).minkowski_difference(guard_kernel).as_original().simplify(1e-5)
            clean=(clean^safe).as_original().simplify(1e-5)
            patch_gap=clean.min_gap(ext,1.)
            newcore=newcore+clean
            g.save_npz(clean,a.out/f'{lid}_site{k:02d}_closed.npz')
            for label,obj in [('before',(outer-core)^roi),('after',(outer-newcore)^roi)]:
                msh(obj).export(a.out/f'{lid}_site{k:02d}_{label}_ENU_DIAGNOSTIC_ONLY.stl')
            site_rows.append({'site':k,'critical_example_indices':group,'roi_enu_mm':[lo.tolist(),hi.tolist()],'source_local_core_volume_mm3':local.volume(),'closed_local_core_volume_mm3':clean.volume(),'closed_local_triangles':clean.num_tri(),'closed_local_component_count':len(clean.decompose()),'whole_patch_to_positive_exterior_gap_mm':patch_gap})
            print(lid,'site',k,'closed clean triangles',clean.num_tri(),'seconds',time.time()-begin,flush=True)
        solid=outer-newcore
        exteriorgap=min(x['whole_patch_to_positive_exterior_gap_mm'] for x in site_rows);foreign=[{'line_id':other,'minimum_gap_mm':newcore.min_gap(c,1.)} for other,c in original.items() if other!=lid]
        row={'line_id':lid,'original_core_sha256':hashlib.sha256((a.cache/f'{lid}_core.npz').read_bytes()).hexdigest(),'original_outer_sha256':hashlib.sha256((a.cache/f'{lid}_outer.npz').read_bytes()).hexdigest(),'sites':site_rows,'extra_void_volume_mm3':newcore.volume()-core.volume(),'whole_patch_to_exterior_gap_mm':exteriorgap,'new_core_gaps_to_other_original_cores':foreign,'core_component_count':len(newcore.decompose()),'solid_component_count_before':len(before_components),'solid_components_before':[{'volume_mm3':x.volume(),'bounds_enu_mm':list(x.bounding_box())} for x in before_components],'solid_component_count':len(solid.decompose()),'source_core_volume_mm3':core.volume(),'new_core_volume_mm3':newcore.volume(),'source_centerline_sha256':hashlib.sha256((a.root/'data/rail_paths_R14.json').read_bytes()).hexdigest(),'ports_sha256':hashlib.sha256((a.root/'data/ports_R14.json').read_bytes()).hexdigest()}
        row['source_geometry_guards_pass']=bool(exteriorgap>=.5 and min(x['minimum_gap_mm'] for x in foreign)>=.5 and row['core_component_count']==1 and row['solid_component_count']==len(before_components))
        print(lid,'guards',row['source_geometry_guards_pass'],'exterior',exteriorgap,'foreign',min(x['minimum_gap_mm'] for x in foreign),'components',row['core_component_count'],row['solid_component_count'],flush=True)
        g.save_npz(newcore,a.out/f'{lid}_core_after.npz');g.save_npz(solid,a.out/f'{lid}_solid_after.npz');final[lid]=newcore
        # Always export the trial; failed guards remain explicit and do not authorize adoption.
        v,f=g.arrays(solid.as_original().simplify(.001));v=v@T[:3,:3].T+T[:3,3];mm,raw=read_stl_bytes(trimesh.Trimesh(v,f,process=False));mm,repair=weld_near_vertices(mm,.00005);mm,raw=read_stl_bytes(mm);path=a.out/f'rail_R14_{lid}_cusp_trial.stl';path.write_bytes(raw)
        row['saved_stl']={'file':path.name,'sha256':hashlib.sha256(raw).hexdigest(),'topology':topology(mm),'numerical_weld':repair}
        ms=ml.MeshSet();ms.add_mesh(ml.Mesh(mm.vertices,mm.faces.astype(np.int32)));ms.apply_filter('compute_selection_by_self_intersections_per_face');row['saved_stl']['self_intersecting_faces']=int(ms.current_mesh().selected_face_number())
        row['saved_stl']['connected_surface_bodies']=int(mm.body_count)
        row['critical_samples_before']=critical_rays(a.root/'models'/f'rail_R14_{lid}.stl',T,ex);row['critical_samples_after']=critical_rays(path,T,ex)
        row['sharp_edges_before']=local_sharp_edges(a.root/'models'/f'rail_R14_{lid}.stl',T,line,sites);row['sharp_edges_after']=local_sharp_edges(path,T,line,sites)
        row['saved_wall_diagnostic']=audit(line,[p for p in ports if p['line_id']==lid],path,T,8,.5)
        row['seconds']=time.time()-begin;result['lines'].append(row);(a.out/'internal_cusp_repair_trials.json').write_text(json.dumps(result,indent=2)+'\n')
        print(lid,'DONE triangles',row['saved_stl']['topology']['triangles'],'selfix',row['saved_stl']['self_intersecting_faces'],'sharp120 before/after',row['sharp_edges_before']['convex_inner_edges_over_120_degrees'],row['sharp_edges_after']['convex_inner_edges_over_120_degrees'],'wall',row['saved_wall_diagnostic']['minimum_containment_normal_mm'],'septum',row['saved_wall_diagnostic']['minimum_septum_radial_mm'],flush=True)
    pairs=[]
    for i,(lid,core) in enumerate(final.items()):
        for other,c in list(final.items())[i+1:]:pairs.append({'line_a':lid,'line_b':other,'minimum_gap_mm':core.min_gap(c,1.)})
    result['all_91_final_core_pairs']=pairs;result['minimum_final_core_pair_gap_mm']=min(x['minimum_gap_mm'] for x in pairs);result['all_final_core_pair_gaps_ge_0_5_mm']=result['minimum_final_core_pair_gap_mm']>=.5
    for row in result['lines']:
        t=row['saved_stl']['topology']
        row['numerical_integration_gates_pass']=bool(row['source_geometry_guards_pass'] and t['watertight'] and t['winding_consistent'] and t['nonmanifold_vertex_links']==0 and t['zero_area_triangles']==0 and t['duplicate_faces']==0 and row['saved_stl']['self_intersecting_faces']==0 and row['saved_stl']['connected_surface_bodies']==row['solid_component_count_before'])
    result['numerical_integration_gates_pass']=bool(result['all_final_core_pair_gaps_ge_0_5_mm'] and all(r['numerical_integration_gates_pass'] for r in result['lines']))
    (a.out/'internal_cusp_repair_trials.json').write_text(json.dumps(result,indent=2)+'\n')
    files=[{'path':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(a.out.iterdir()) if p.is_file() and p.name!='artifact_manifest_R14_internal_junctions.json']
    (a.out/'artifact_manifest_R14_internal_junctions.json').write_text(json.dumps({'numerical_integration_gates_pass':result['numerical_integration_gates_pass'],'physical_tests_performed':False,'files':files},indent=2)+'\n')
    print('ALL PAIRS',result['minimum_final_core_pair_gap_mm'],flush=True)
    if not result['numerical_integration_gates_pass']:
        raise SystemExit('Numerical integration gates failed; inspect the output report before adoption')
if __name__=='__main__':main()
