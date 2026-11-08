#!/usr/bin/env python3
"""Reproduce the digitally repaired R14 body from versioned inputs.

The default uses a SHA-pinned positive checkpoint and recalculates the final
subtraction, print-coordinate serialization, boundary selection and diagnosed
seam correction. --rebuild-positive replays all documented local additions
from the recovered base and current rail positives. Boolean face ordering may
vary; acceptance always concerns the actual output bytes and independent QA.
The seed-only positive replay is an experimental diagnostic path: the tested
replay has a nonmanifold saved edge and is rejected by the basic topology gate.
The pinned-checkpoint path is the supported reproduction of the delivered body.
No physical print or pigment-flow qualification is performed by this script.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,gc,hashlib,json,time
from pathlib import Path
import numpy as np,trimesh
from geometry_r14 import arrays,mf,load_npz,save_npz,unite
from robust_export_r14 import read_stl_bytes,basic_topology
from extract_principal_boundary import extract
from repair_selected_seams import repair
from core_signature_r14 import signature as core_signature


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def check_file(base,record):
    p=base/record['file']
    if sha(p)!=record['sha256']:raise ValueError(f'Input checksum mismatch: {p}')
    return p

def apply_recipe(positive,root,recipe,out):
    for i,step in enumerate(recipe['steps']):
        p=check_file(root/'data',step);kind=step['kind'];items=[]
        if kind=='edge_capsules':
            k,_=arrays(mf.Manifold.sphere(step['radius_mm'],step['sphere_segments']))
            segments=np.asarray(json.loads(p.read_text()))
            items=[mf.Manifold.hull_points(np.vstack([k+a,k+b])) for a,b in segments]
            positive=unite([positive]+items)
        elif kind=='triangle_and_point_capsules':
            with np.load(p) as z:
                triangles=z[step['triangle_key']]
                points=z[step['point_key']] if step.get('point_key') else np.empty((0,3))
            k,_=arrays(mf.Manifold.sphere(step['triangle_radius_mm'],step['triangle_sphere_segments']))
            items=[mf.Manifold.hull_points((t[:,None,:]+k[None,:,:]).reshape(-1,3)) for t in triangles]
            if len(points):
                items += [mf.Manifold.sphere(step['point_radius_mm'],step['point_sphere_segments']).translate(x) for x in points]
            augmentation=unite(items);positive=unite([positive,augmentation])
        elif kind=='points_and_bridge':
            with np.load(p) as z:points=z[step['point_key']]
            items=[mf.Manifold.sphere(step['point_radius_mm'],step['point_sphere_segments']).translate(x) for x in points]
            items.append(load_npz(check_file(root/'data',step['bridge'])))
            positive=unite([positive]+items)
        else:raise ValueError(f'Unknown repair step: {kind}')
        if step.get('reload_after'):
            checkpoint=out/f'positive_step_{i:02d}.npz';save_npz(positive,checkpoint)
            positive=load_npz(checkpoint)
        print('REPAIR',i,kind,positive.num_tri(),flush=True)
        del items;gc.collect()
    return positive


def serialize_candidate(obj,T,settings):
    obj=obj.as_original().simplify(settings['simplification_mm']).transform(T[:3])
    for _ in range(settings['float32_roundtrip_count']):
        v,f=arrays(obj)
        obj=mf.Manifold(mf.Mesh64(v.astype(np.float32).astype(float),np.asarray(f,dtype=np.uint64),
                                tolerance=settings['roundtrip_manifold_tolerance_mm']))
        obj=obj.as_original().simplify(settings['roundtrip_simplification_mm'])
    v,f=arrays(obj);m,_=read_stl_bytes(trimesh.Trimesh(v,f,process=False))
    m.update_faces(m.nondegenerate_faces(height=settings['minimum_face_height_mm']))
    m.update_faces(m.unique_faces());m.remove_unreferenced_vertices()
    return m


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--core-cache',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--rebuild-positive',action='store_true');a=p.parse_args();start=time.time()
    root=a.root.resolve();a.out.mkdir(parents=True,exist_ok=True)
    recipe=json.loads((root/'data/body_repair_recipe_R14.json').read_text())
    for name,digest in recipe['source_inputs_sha256'].items():
        if sha(root/name)!=digest:raise ValueError(f'Recovered input changed: {name}')
    if a.rebuild_positive:
        positive=apply_recipe(load_npz(a.cache/'positive.npz'),root,recipe,a.out)
        save_npz(positive,a.out/'replayed_positive.npz')
    else:
        if not recipe.get('positive_checkpoint'):raise ValueError('No pinned positive checkpoint in recipe')
        positive=load_npz(check_file(root/'data',recipe['positive_checkpoint']))
    core_geometry=json.loads((root/'data/body_core_geometry_R14.json').read_text())
    for lid,digest in recipe.get('final_core_sha256',{}).items():
        path=a.core_cache/f'{lid}_core.npz'
        if sha(path)!=digest and core_signature(path)!=core_geometry['oriented_geometry_sha256'][path.name]:
            raise ValueError(f'Final lumen surface differs: {lid}')
    combined=a.core_cache/'all_cores.npz'
    if core_signature(combined)!=core_geometry['oriented_geometry_sha256'][combined.name]:
        raise ValueError('Combined lumen cache does not match the recorded final14 cores')
    cores=load_npz(combined)
    body=positive-cores
    T=np.asarray(json.loads((root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    mesh=serialize_candidate(body,T,recipe['serialization'])
    raw=a.out/'raw_with_internal_shells_print.stl';raw.write_bytes(trimesh.exchange.stl.export_stl(mesh))
    principal=a.out/'principal_before_seam.stl';extract(raw,principal)
    target=a.out/'subway_R14_hollow_candidate_print.stl'
    pair_record=recipe['serialization'].get('prescribed_seam_pairs')
    pair_seeds=check_file(root/'data',pair_record) if pair_record else None
    seam=repair(principal,target,recipe['serialization']['selected_intersection_seam_weld_mm'],pair_seeds=pair_seeds)
    topology=seam['topology']
    basic_pass=bool(topology['watertight'] and topology['winding_consistent']
                    and topology['nonmanifold_edges']==0 and topology['boundary_edges']==0
                    and topology['zero_area_triangles']==0 and topology['volume_mm3']>0)
    record=dict(recipe_sha256=sha(root/'data/body_repair_recipe_R14.json'),
                source_script_sha256=sha(Path(__file__)),
                core_geometry_manifest_sha256=sha(root/'data/body_core_geometry_R14.json'),
                all_core_file_sha256=sha(combined),
                positive_rebuilt_from_seeds=a.rebuild_positive,all_final_cores_recut=True,
                raw_sha256=sha(raw),principal_before_seam_sha256=sha(principal),
                candidate_sha256=sha(target),candidate_file=target.name,
                basic_serialized_topology_passed=basic_pass,
                candidate_status='requires independent acceptance' if basic_pass else 'REJECTED: failed saved-file topology',
                supported_reproduction_mode='SHA-pinned positive checkpoint plus final lumen recut and saved-mesh corrections',
                diagnostic_only=True,independent_acceptance_required=True,
                required_checks=['topology and vertex links','all-triangle exact self intersections',
                                 'removed-surface containment and intended-core separation',
                                 'installed channels and 69 openings','wall regions','stored normals and orientation'],
                physical_tests_performed=False,seconds=time.time()-start)
    (a.out/'body_reproduction_R14.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    if not basic_pass:
        raise SystemExit('Rejected saved candidate: basic topology failed. Retain the diagnostics and use the pinned-checkpoint reproduction path for the validated delivered R14.')

if __name__=='__main__':main()
