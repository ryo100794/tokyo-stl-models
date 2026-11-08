#!/usr/bin/env python3
"""Test every subsidiary enclosed volume against the intended final cutters."""
from pathlib import Path
import argparse,json,time,hashlib,gc
import numpy as np,trimesh,manifold3d as mf
from mesh_validation_utils import subset,load_npz_arrays

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def asmf(v,f):return mf.Manifold(mf.Mesh64(v.astype('float64'),f.astype('uint64')))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('mesh');ap.add_argument('--face-labels',required=True);ap.add_argument('--component-inventory',required=True);ap.add_argument('--core-dir',required=True);ap.add_argument('--print-transform',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
 out=Path(a.out);out.mkdir(exist_ok=True,parents=True);src=Path(a.mesh);start=time.monotonic();invrep=json.loads(Path(a.component_inventory).read_text());assert invrep['sha256']==sha(src)
 m=trimesh.load_mesh(src,process=False);v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];del m,inv;gc.collect()
 T=np.asarray(json.loads(Path(a.print_transform).read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major']);Ti=np.linalg.inv(T);v=v@Ti[:3,:3].T+Ti[:3,3]
 with np.load(a.face_labels) as z:lab=z['labels']
 manifest=[];coreobjs=[]
 for p in sorted(Path(a.core_dir).glob('L[0-9][0-9]_core.npz')):
  cv,cf=load_npz_arrays(p);o=asmf(cv,cf);assert o.status()==mf.Error.NoError;coreobjs.append(o);manifest.append(dict(file=str(p),line_id=p.stem[:3],sha256=sha(p),triangles=o.num_tri(),volume_mm3=o.volume()))
 assert len(coreobjs)==14;cores=mf.Manifold.batch_boolean(coreobjs,mf.OpType.Add);assert cores.status()==mf.Error.NoError
 rows=[]
 for c in invrep['components']:
  i=c['component']
  if i==invrep['primary_component']:continue
  sv,sf=subset(v,f,lab==i);neg=c['signed_volume_mm3']<0
  if neg:sf=sf[:,::-1]
  obj=asmf(sv,sf);row=dict(component=i,original_signed_volume_mm3=c['signed_volume_mm3'],negative_shell_reoriented_to_positive_enclosed_space=neg,query_manifold_status=str(obj.status()))
  if obj.status()!=mf.Error.NoError:row['separation_gate_passed']=False;row['reason']='Cannot interpret this closed surface as a Boolean volume'
  else:
   gap=float(obj.min_gap(cores,5.));inter=obj^cores
   row.update(minimum_gap_to_union_of_final_cores_mm_search5=gap,core_intersection_volume_mm3=inter.volume(),core_intersection_triangles=inter.num_tri(),separation_gate_passed=bool(gap>1e-4 and inter.is_empty()))
   if not row['separation_gate_passed']:
    row['individual_core_checks']=[dict(line_id=manifest[j]['line_id'],minimum_gap_mm_search5=float(obj.min_gap(co,5.)),intersection_volume_mm3=(obj^co).volume()) for j,co in enumerate(coreobjs)]
  rows.append(row)
  if len(rows)%50==0:print('CHECKED',len(rows),time.monotonic()-start,flush=True)
 report=dict(source=str(src),source_sha256=sha(src),coordinates='Saved print XYZ inverted by recorded transform into ENU mm for core comparison',print_transform=T.tolist(),cores=manifest,union_core_triangles=cores.num_tri(),union_core_volume_mm3=cores.volume(),subsidiary_components_tested=len(rows),all_separation_gates_passed=all(r['separation_gate_passed'] for r in rows),minimum_reported_gap_mm=min(r.get('minimum_gap_to_union_of_final_cores_mm_search5',0) for r in rows),components=rows,method='Boolean overlap and minimum gap against the true union of all 14 final double-precision cutters. Negative shells reversed to represent their enclosed cavity volume.',separation_margin_mm=1e-4,limitations=['Geometric separation only; no physical print or flow qualification.','A 5 mm distance value is a capped search result, meaning at least 5 mm.','A failed or unclassifiable separation gate must be reviewed and must not be bypassed solely because a shell is small.','Core hash manifest must match the final cutters used in the delivered generation pipeline.'],elapsed_seconds=time.monotonic()-start)
 (out/'secondary_core_separation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');print('DONE',json.dumps({k:report[k] for k in ['subsidiary_components_tested','all_separation_gates_passed','minimum_reported_gap_mm','elapsed_seconds']}),flush=True)
if __name__=='__main__':main()
