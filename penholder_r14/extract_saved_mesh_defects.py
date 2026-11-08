#!/usr/bin/env python3
"""Export defect locations in saved coordinates and optionally original ENU."""
from pathlib import Path
import argparse,json,time,sys,hashlib,inspect
import numpy as np,trimesh,pymeshlab
from validate_saved_meshes import vertex_links

def main():
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('--out',required=True);ap.add_argument('--topology-report',required=True);ap.add_argument('--print-transform');a=ap.parse_args()
    out=Path(a.out);out.mkdir(exist_ok=True,parents=True);src=Path(a.input);sha=hashlib.sha256(src.read_bytes()).hexdigest()
    old=json.loads(Path(a.topology_report).read_text());assert old['sha256']==sha
    assert old['boundary_edges']==0 and old['edges_with_more_than_two_faces']==0 and old['inconsistent_two_face_edge_orientation']==0
    m=trimesh.load_mesh(src,process=False);v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces]
    T=np.eye(4)
    if a.print_transform:T=np.asarray(json.loads(Path(a.print_transform).read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    Ti=np.linalg.inv(T)
    def enu(points):return np.einsum('...j,ij->...i',points,Ti[:3,:3])+Ti[:3,3]
    options={'return_all_ids':True} if 'return_all_ids' in inspect.signature(vertex_links).parameters else {}
    bad=vertex_links(f,len(v),True,True,**options)['nonmanifold_vertex_ids'];bad=np.asarray(bad,dtype=np.int64)
    np.savez_compressed(out/'nonmanifold_vertices.npz',vertex_ids=bad,coordinates_saved_mm=v[bad],coordinates_enu_mm=enu(v[bad]))
    print('VERTICES',len(bad),flush=True)
    ms=pymeshlab.MeshSet();ms.add_mesh(pymeshlab.Mesh(v,f));ms.compute_selection_by_self_intersections_per_face();sel=np.flatnonzero(ms.current_mesh().face_selection_array())
    tri=v[f[sel]]
    np.savez_compressed(out/'self_intersection_selection.npz',face_ids=sel,vertex_ids=f[sel],triangles_saved_mm=tri,triangles_enu_mm=enu(tri))
    report=dict(source=str(src),source_sha256=sha,saved_coordinates='print XYZ mm' if a.print_transform else 'source ENU mm',print_transform_inverted_only_for_repair_coordinates=bool(a.print_transform),exact_predicate_must_use_saved_coordinates=True,selected_self_intersection_face_count=len(sel),selected_self_intersection_face_ids=sel.tolist(),nonmanifold_vertex_count=len(bad),nonmanifold_vertices=[dict(vertex_id=int(i),coordinates_saved_mm=v[i].tolist(),coordinates_enu_mm=enu(v[i]).tolist()) for i in bad])
    (out/'defects.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print('DONE',len(sel),len(bad),flush=True)

if __name__=='__main__':main()
