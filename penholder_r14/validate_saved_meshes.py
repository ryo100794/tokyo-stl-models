#!/usr/bin/env python3
"""Read-only, independent validation of STL bytes. Does not repair or modify inputs."""
from __future__ import annotations
import argparse, hashlib, json, os, sys, time, gc
from pathlib import Path
import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import trimesh
import pymeshlab


def sha256(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def vertex_links(faces, nvertices, edge_ok, orientation_ok):
    # For consistently oriented closed triangular edge manifolds, each directed
    # vertex link is a disjoint union of cycles; exactly one cycle is required.
    # General local undirected link is checked if the edge prerequisite fails.
    starts = faces.reshape(-1)
    nxt = np.roll(faces,-1,axis=1).reshape(-1)
    end = np.roll(faces,-2,axis=1).reshape(-1)
    order = np.argsort(starts,kind='stable')
    counts = np.bincount(starts,minlength=nvertices)
    offsets = np.r_[0,np.cumsum(counts)]
    nxt,end = nxt[order],end[order]
    bad=[]; isolated=0
    for v in range(nvertices):
        a,b = offsets[v:v+2]
        if a==b:
            isolated+=1;continue
        pairs=list(zip(nxt[a:b].tolist(),end[a:b].tolist()))
        if edge_ok and orientation_ok:
            succ=dict(pairs)
            x=pairs[0][0]; seen=set()
            while x not in seen and x in succ:
                seen.add(x);x=succ[x]
            valid=(len(succ)==len(pairs) and len(seen)==len(pairs) and x==pairs[0][0])
        else:
            nbr={}
            for x,y in pairs:
                nbr.setdefault(x,[]).append(y);nbr.setdefault(y,[]).append(x)
            seen=set(); todo=[pairs[0][0]]
            while todo:
                x=todo.pop()
                if x in seen:continue
                seen.add(x);todo.extend(nbr[x])
            valid=all(len(q)==2 for q in nbr.values()) and len(seen)==len(nbr)
        if not valid:bad.append(int(v))
    return {'nonmanifold_vertex_count':len(bad),'nonmanifold_vertex_ids_first_50':bad[:50],
            'nonmanifold_vertex_ids':bad,
            'isolated_vertex_count':isolated,'method':'full local vertex-link cycle traversal after exact-coordinate weld'}


def validate(path, intersection=True, links=True):
    path=Path(path); ts=time.monotonic()
    print('START '+str(path),flush=True)
    if path.suffix.lower()=='.npz':
        data=np.load(path,allow_pickle=False)
        vkey = 'vertices' if 'vertices' in data.files else 'v'
        fkey = 'faces' if 'faces' in data.files else 'f'
        m=trimesh.Trimesh(vertices=data[vkey],faces=data[fkey],process=False)
    else:
        m=trimesh.load_mesh(path,process=False)
    if not isinstance(m,trimesh.Trimesh): raise TypeError('Expected one triangular STL mesh')
    # STL has no shared-index topology. Weld only bit-identical saved coordinates.
    v,remap=np.unique(np.asarray(m.vertices),axis=0,return_inverse=True)
    f=remap[np.asarray(m.faces)]
    del m,remap;gc.collect()
    mesh=trimesh.Trimesh(vertices=v,faces=f,process=False)
    e=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    sorted_e=np.sort(e,axis=1)
    eu,inv,cnt=np.unique(sorted_e,axis=0,return_inverse=True,return_counts=True)
    signs=np.where(e[:,0]<e[:,1],1,-1)
    orientation=np.bincount(inv,weights=signs,minlength=len(cnt))
    boundary=int(np.count_nonzero(cnt==1))
    nonmanifold=int(np.count_nonzero(cnt>2))
    inconsistent=int(np.count_nonzero((cnt==2)&(orientation!=0)))
    graph=coo_matrix((np.ones(len(eu),dtype=np.uint8),(eu[:,0],eu[:,1])),shape=(len(v),len(v))).tocsr()
    nc,lab=connected_components(graph,directed=False,return_labels=True)
    flab=lab[f[:,0]]
    tri=v[f]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    signed_v=np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2]))/6
    volumes=np.bincount(flab,weights=signed_v,minlength=nc)
    areas=np.bincount(flab,weights=area,minlength=nc)
    fn=np.bincount(flab,minlength=nc)
    vn=np.bincount(lab,minlength=nc)
    components=[]
    for i in np.argsort(fn)[::-1]:
        verts=v[lab==i]
        components.append({'index':int(i),'vertices':int(vn[i]),'faces':int(fn[i]),
                           'signed_volume_mm3':float(volumes[i]),'area_mm2':float(areas[i]),
                           'bounds_mm':[verts.min(axis=0).tolist(),verts.max(axis=0).tolist()]})
    result={
        'file':str(path),'basename':path.name,'sha256':sha256(path),'bytes':path.stat().st_size,
        'input_modified':False,'units':'mm as specified by project; STL carries no units',
        'welding':'exact identical saved coordinates only; no tolerance merge, face removal, smoothing or repair',
        'vertices_after_exact_weld':len(v),'faces':len(f),'unique_edges':len(eu),
        'boundary_edges':boundary,'edges_with_more_than_two_faces':nonmanifold,
        'inconsistent_two_face_edge_orientation':inconsistent,
        'zero_area_triangles_exact':int(np.count_nonzero(area==0)),
        'triangles_area_below_1e_minus_12_mm2':int(np.count_nonzero(area<1e-12)),
        'duplicate_faces_ignoring_orientation':int(len(f)-len(np.unique(np.sort(f,axis=1),axis=0))),
        'finite_vertex_coordinates':bool(np.isfinite(v).all()),
        'signed_volume_mm3':float(signed_v.sum()),'surface_area_mm2':float(area.sum()),
        'bounds_mm':[v.min(axis=0).tolist(),v.max(axis=0).tolist()],
        'dimensions_mm':np.ptp(v,axis=0).tolist(),
        'connected_surface_components':int(nc),'components':components,
        'surface_components_caveat':'Surface components are not by themselves solid-component counts: negative oriented shells may bound enclosed cavities. Positive separately bounded shells require nesting/overlap interpretation.',
        'euler_characteristic':int(len(v)-len(eu)+len(f)),
        'trimesh_is_watertight':bool(mesh.is_watertight),'trimesh_is_winding_consistent':bool(mesh.is_winding_consistent),
    }
    del graph,lab,flab,tri,area,signed_v,sorted_e,e,eu,inv,cnt,orientation;gc.collect()
    if links:
        print('VERTEX_LINKS '+path.name,flush=True)
        result['vertex_links']=vertex_links(f,len(v),boundary==0 and nonmanifold==0,inconsistent==0)
    if intersection:
        print('SELF_INTERSECTIONS '+path.name,flush=True)
        ms=pymeshlab.MeshSet();ms.add_mesh(pymeshlab.Mesh(v,f))
        ms.compute_selection_by_self_intersections_per_face()
        selected=ms.current_mesh().face_selection_array()
        result['self_intersection_test']={
            'selected_face_count':int(np.count_nonzero(selected)),
            'selected_face_ids_first_100':np.flatnonzero(selected)[:100].tolist(),
            'method':'PyMeshLab compute_selection_by_self_intersections_per_face on exact-coordinate-welded saved STL',
            'limitation':'Floating-point triangle intersection predicate; zero means none detected, not an exact-arithmetic guarantee. Adjacent/shared-topology contacts and coincident/coplanar degeneracies may need separate assessment.'}
        del ms
    result['elapsed_seconds']=round(time.monotonic()-ts,3)
    print('DONE '+path.name+' '+json.dumps({k:result[k] for k in ['faces','connected_surface_components','boundary_edges','edges_with_more_than_two_faces','elapsed_seconds']})+' '+json.dumps(result.get('self_intersection_test',{})),flush=True)
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('inputs',nargs='+');ap.add_argument('--out',required=True)
    ap.add_argument('--skip-intersections',action='store_true');ap.add_argument('--skip-vertex-links',action='store_true')
    a=ap.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    all_results=[]
    for path in a.inputs:
        try:
            result=validate(path,not a.skip_intersections,not a.skip_vertex_links)
        except Exception as e:
            result={'file':path,'error':repr(e),'validation_status':'incomplete'}
            print('ERROR '+path+' '+repr(e),flush=True)
        (out/(Path(path).stem+'.independent.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
        all_results.append(result);gc.collect()
    report={'schema':'R14-independent-saved-STL-validation-v1','runtime':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'trimesh':trimesh.__version__},'tests':all_results,
      'does_not_establish':['material wall-thickness certification','fluid connectivity or leak-tightness','successful SLA printing','support adequacy','successful physical flow','strength or pressure rating']}
    (out/'independent_saved_meshes_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')

if __name__=='__main__':main()
