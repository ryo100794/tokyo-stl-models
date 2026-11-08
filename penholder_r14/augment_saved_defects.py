"""Add bounded material around diagnosed defects, then require a complete lumen recut.

The diagnostic mesh is never smoothed or hole-filled. Repair seeds are saved in
design coordinates so the operation is explicit and reproducible. This operation
does not itself accept a model; independent saved-file validation is mandatory.
"""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
import trimesh
import pymeshlab as pm
from geometry_r14 import arrays, mf, unite, load_npz, save_npz


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model',type=Path,required=True)
    ap.add_argument('--positive',type=Path,required=True)
    ap.add_argument('--qa',type=Path,required=True)
    ap.add_argument('--coordinates',type=Path)
    ap.add_argument('--triangle-radius',type=float,default=.03)
    ap.add_argument('--vertex-radius',type=float,default=.06)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True);started=time.time()
    raw=trimesh.load_mesh(a.model,process=False)
    v,inv=np.unique(raw.vertices,axis=0,return_inverse=True);f=inv[raw.faces]
    ms=pm.MeshSet();ms.add_mesh(pm.Mesh(v,f));ms.compute_selection_by_self_intersections_per_face()
    ids=np.flatnonzero(ms.current_mesh().face_selection_array())
    triangles=v[f[ids]].copy()
    qa=json.loads(a.qa.read_text())
    expected=hashlib.sha256(a.model.read_bytes()).hexdigest()
    if qa['sha256']!=expected:raise ValueError('QA does not match the saved diagnostic model')
    vertex_ids=qa['vertex_links'].get('nonmanifold_vertex_ids',qa['vertex_links']['nonmanifold_vertex_ids_first_50'])
    if len(vertex_ids)!=qa['vertex_links']['nonmanifold_vertex_count']:raise ValueError('Incomplete defect vertex list')
    points=v[vertex_ids].copy()
    if a.coordinates:
        T=np.asarray(json.loads(a.coordinates.read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major']);invT=np.linalg.inv(T)
        triangles=triangles@invT[:3,:3].T+invT[:3,3]
        points=points@invT[:3,:3].T+invT[:3,3]
    k,_=arrays(mf.Manifold.sphere(a.triangle_radius,8))
    patches=[mf.Manifold.hull_points((t[:,None,:]+k[None,:,:]).reshape(-1,3)) for t in triangles]
    patches += [mf.Manifold.sphere(a.vertex_radius,12).translate(p) for p in points]
    augmentation=unite(patches)
    save_npz(augmentation,a.out/'augmentation.npz')
    np.savez_compressed(a.out/'repair_geometry.npz',triangles_enu_mm=triangles,points_enu_mm=points)
    positive=unite([load_npz(a.positive),augmentation]);save_npz(positive,a.out/'positive.npz')
    report=dict(source_stl=a.model.name,source_stl_sha256=expected,self_intersection_faces=int(len(ids)),
                nonmanifold_vertices=len(points),triangle_augmentation_radius_mm=a.triangle_radius,
                vertex_augmentation_radius_mm=a.vertex_radius,all_core_recuts_required=True,
                augmentation_volume_mm3=augmentation.volume(),seconds=time.time()-started,
                diagnostic_only=True)
    (a.out/'method.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':main()
