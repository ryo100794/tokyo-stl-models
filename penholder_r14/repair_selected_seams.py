"""Bounded float32 seam correction restricted to diagnosed intersecting faces.

No hole filling or global smoothing is performed. The output remains a
candidate and must pass independent topology, intersection and channel tests.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np,trimesh,pymeshlab as pm
from robust_export_r14 import weld_near_vertices,read_stl_bytes,basic_topology

def repair(input_path,output_path,tolerance=.00005,diagnosis=None,pair_seeds=None):
    raw=trimesh.load_mesh(input_path,process=False)
    v,inv=np.unique(raw.vertices,axis=0,return_inverse=True);f=inv[raw.faces]
    mesh=trimesh.Trimesh(v,f,process=False)
    if pair_seeds is not None:
        seeds=json.loads(Path(pair_seeds).read_text())
        positions=np.asarray(seeds['pairs_saved_print_mm']).reshape(-1,3)
        lookup={tuple(p):i for i,p in enumerate(v)}
        if any(tuple(p) not in lookup for p in positions):
            raise ValueError('Prescribed numerical seam positions are absent from the saved mesh')
        eligible=np.asarray([lookup[tuple(p)] for p in positions],dtype=int)
        selected=np.flatnonzero(np.isin(f,eligible).any(axis=1))
    elif diagnosis is None:
        ms=pm.MeshSet();ms.add_mesh(pm.Mesh(v,f));ms.compute_selection_by_self_intersections_per_face()
        selected=np.flatnonzero(ms.current_mesh().face_selection_array())
    else:
        d=json.loads(Path(diagnosis).read_text())
        if d['source_sha256']!=hashlib.sha256(input_path.read_bytes()).hexdigest():
            raise ValueError('Exact diagnosis does not match input STL')
        selected=np.asarray(d['abnormal_face_ids'],dtype=int)
    if pair_seeds is None:eligible=np.unique(f[selected])
    mesh,weld=weld_near_vertices(mesh,tolerance,eligible)
    mesh,data=read_stl_bytes(mesh)
    output_path.parent.mkdir(parents=True,exist_ok=True);output_path.write_bytes(data)
    record=dict(source_file=input_path.name,source_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
                output_file=output_path.name,sha256=hashlib.sha256(data).hexdigest(),
                diagnosed_face_ids=selected.tolist(),diagnosis_source=Path(pair_seeds).name if pair_seeds is not None else 'MeshLab' if diagnosis is None else Path(diagnosis).name,
                weld=weld,topology=basic_topology(mesh),
                independent_validation_required=True,physical_tests_performed=False)
    output_path.with_suffix('.seam_correction.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True);return record

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
    ap.add_argument('--tolerance',type=float,default=.00005);ap.add_argument('--diagnosis',type=Path);ap.add_argument('--pair-seeds',type=Path);a=ap.parse_args()
    repair(a.input,a.output,a.tolerance,a.diagnosis,a.pair_seeds)
