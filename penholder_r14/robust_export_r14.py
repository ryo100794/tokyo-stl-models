"""Bounded numerical seam welding for R14 STL serialization.

This does not fill holes or run a surface reconstruction. It only merges
vertices separated by a measured sub-micrometre distance, removes collapsed
faces, and cancels coincident oppositely oriented triangles. All acceptance
checks operate on reloaded STL bytes. Geometry qualification is separate.
"""
from pathlib import Path
import argparse, gc, hashlib, io, json, time
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from geometry_r14 import load_npz, arrays


def weld_near_vertices(mesh, tolerance, eligible_vertices=None):
    v=np.asarray(mesh.vertices).copy(); f=np.asarray(mesh.faces).copy()
    if eligible_vertices is None:
        pairs=cKDTree(v).query_pairs(tolerance, output_type='ndarray')
    else:
        eligible_vertices=np.unique(np.asarray(eligible_vertices,dtype=int))
        local=cKDTree(v[eligible_vertices]).query_pairs(tolerance, output_type='ndarray')
        pairs=eligible_vertices[local]
    parent=np.arange(len(v))
    def find(a):
        while parent[a]!=a:
            parent[a]=parent[parent[a]];a=parent[a]
        return a
    for a,b in pairs:
        a=find(a);b=find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    involved=np.unique(pairs) if len(pairs) else np.empty(0,dtype=int)
    for a in involved:parent[a]=find(a)
    ids,labels=np.unique(parent,return_inverse=True)
    counts=np.bincount(labels)
    w=np.column_stack([np.bincount(labels,weights=v[:,k])/counts for k in range(3)])
    # Quantize now; the measured displacement includes actual STL rounding.
    w=w.astype(np.float32).astype(np.float64)
    displacement=np.linalg.norm(v-w[labels],axis=1)
    max_move=float(displacement.max(initial=0))
    if max_move>tolerance*1.5:
        raise ValueError('Transitive weld exceeds the bounded displacement')
    ff=labels[f]
    deg=(ff[:,0]==ff[:,1])|(ff[:,1]==ff[:,2])|(ff[:,2]==ff[:,0])
    ff=ff[~deg]
    keys=np.sort(ff,axis=1)
    _,first,group,multiplicity=np.unique(keys,axis=0,return_index=True,return_inverse=True,return_counts=True)
    keep=np.ones(len(ff),dtype=bool)
    for k in np.flatnonzero(multiplicity>1):
        positions=np.flatnonzero(group==k)
        faces=ff[positions]
        parity=np.where(((faces[:,0]>faces[:,1]).astype(int)+(faces[:,0]>faces[:,2])+(faces[:,1]>faces[:,2]))%2, -1,1)
        keep[positions]=False
        if parity.sum()!=0:keep[positions[np.flatnonzero(parity==np.sign(parity.sum()))[0]]]=True
    result=trimesh.Trimesh(w,ff[keep],process=False)
    nondeg=result.nondegenerate_faces(height=1e-10)
    result.update_faces(nondeg);result.remove_unreferenced_vertices()
    return result,dict(weld_tolerance_mm=tolerance,near_vertex_pairs=len(pairs),
        merged_vertex_count=len(v)-len(w),maximum_vertex_displacement_mm=max_move,
        collapsed_faces=int(deg.sum()+np.count_nonzero(~nondeg)),
        coincident_faces_removed=int(np.count_nonzero(~keep)),hole_filling_performed=False,
        scope='all vertices' if eligible_vertices is None else 'only vertices of diagnosed self-intersection faces',
        merged_pairs_before_mm=v[pairs].tolist() if eligible_vertices is not None else None)


def read_stl_bytes(mesh):
    raw=trimesh.exchange.stl.export_stl(mesh)
    saved=trimesh.load_mesh(io.BytesIO(raw),file_type='stl',process=True)
    return saved,raw


def basic_topology(m):
    _,count=np.unique(m.edges_sorted,axis=0,return_counts=True)
    return dict(vertices=len(m.vertices),triangles=len(m.faces),
        watertight=bool(m.is_watertight),winding_consistent=bool(m.is_winding_consistent),
        boundary_edges=int((count==1).sum()),nonmanifold_edges=int((count>2).sum()),
        zero_area_triangles=int((m.area_faces<=0).sum()),volume_mm3=float(m.volume))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--transform',type=Path)
    ap.add_argument('--simplify',type=float,default=.005);ap.add_argument('--weld',type=float,default=.00005)
    a=ap.parse_args();start=time.time();obj=load_npz(a.input)
    obj=obj.as_original().simplify(a.simplify)
    v,f=arrays(obj);del obj;gc.collect()
    if a.transform:
        T=np.asarray(json.loads(a.transform.read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
        v=v@T[:3,:3].T+T[:3,3]
    m,raw=read_stl_bytes(trimesh.Trimesh(v,f,process=False));del v,f,raw;gc.collect()
    before=basic_topology(m);print('before',before,flush=True)
    m,repair=weld_near_vertices(m,a.weld)
    m,raw=read_stl_bytes(m);after=basic_topology(m)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_bytes(raw)
    report=dict(input=str(a.input),output=a.output.name,simplification_tolerance_mm=a.simplify,
        before=before,repair=repair,after=after,sha256=hashlib.sha256(raw).hexdigest(),
        seconds=time.time()-start,full_geometry_qualified=False)
    a.output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':main()
