#!/usr/bin/env python3
"""Exhaustive BVH candidate enumeration and exact saved-mesh triangle predicates."""
from pathlib import Path
import argparse,json,time,subprocess,hashlib,gc
import numpy as np,trimesh
from build_exact_engines import ensure_engine
from run_exact_triangle_pairs import CLASSES

def write_indexed(v,f,path):
 with Path(path).open('wb') as s:
  np.asarray([len(v),len(f)],dtype='<u8').tofile(s);np.asarray(v,dtype='<f8').tofile(s);np.asarray(f,dtype='<u4').tofile(s)

def run_indexed(v,f,prefix,threads=3):
 inp=prefix.with_suffix('.indexed_input.bin');out=prefix.with_suffix('.bad_pairs.bin');write_indexed(v,f,inp)
 r=subprocess.run([str(ensure_engine('exact_mesh_bvh')),str(inp),str(out),str(threads)],stdout=subprocess.PIPE,text=True,check=True)
 stats=json.loads(r.stdout.strip().splitlines()[-1]);dtype=np.dtype([('face_a','<u4'),('face_b','<u4'),('code','u1')]);hits=np.fromfile(out,dtype=dtype)
 return stats,hits,inp,out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('mesh');ap.add_argument('--out',required=True);ap.add_argument('--threads',type=int,default=3);ap.add_argument('--keep-input',action='store_true');a=ap.parse_args()
 src=Path(a.mesh);out=Path(a.out);out.parent.mkdir(exist_ok=True,parents=True);start=time.monotonic();sha=hashlib.sha256(src.read_bytes()).hexdigest()
 m=trimesh.load_mesh(src,process=False);v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];assert np.isfinite(v).all();del m,inv;gc.collect()
 inp=out.with_suffix('.indexed_input.bin');hitsfile=out.with_suffix('.bad_pairs.bin');write_indexed(v,f,inp);nf=len(f);nv=len(v);del v,f;gc.collect()
 print('INPUT_READY',nv,nf,flush=True)
 result=subprocess.run([str(ensure_engine('exact_mesh_bvh')),str(inp),str(hitsfile),str(a.threads)],stdout=subprocess.PIPE,text=True,check=True)
 stats=json.loads(result.stdout.strip().splitlines()[-1]);dtype=np.dtype([('face_a','<u4'),('face_b','<u4'),('code','u1')]);hits=np.fromfile(hitsfile,dtype=dtype)
 deg=np.fromfile(str(hitsfile)+'.degenerate.bin',dtype='<u4')
 np.savez_compressed(out.with_suffix('.npz'),face_id_pairs=np.c_[hits['face_a'],hits['face_b']],classification_codes=hits['code'],degenerate_face_ids=deg)
 assert hashlib.sha256(src.read_bytes()).hexdigest()==sha
 report=dict(schema='R14-exhaustive-exact-saved-mesh-intersection-v1',source=str(src),source_sha256=sha,exact_coordinate_weld_only=True,engine=stats,class_counts={CLASSES[k]:n for k,n in enumerate(stats['counts'])},abnormal_pair_count=len(hits),abnormal_face_ids=np.unique(np.r_[hits['face_a'],hits['face_b']]).tolist(),degenerate_face_ids=deg.tolist(),no_abnormal_intersections_and_no_exact_degenerate_faces=bool(len(hits)==0 and len(deg)==0),candidate_completeness='All triangles in the saved mesh; a closed AABB BVH self-traversal enumerates every overlapping triangle-box pair once. No MeshLab selection filter.',arithmetic='All binary64 representations of saved coordinates decoded exactly to integers with a common power of two. GMP arbitrary-precision integer orientations and rational intersection comparisons, without geometric tolerance.',adjacency='Only an intersection wholly contained in the exact shared vertex or shared edge is normal. Other intersections remain failures, including those in a pair that shares one vertex.',does_not_establish=['edge or vertex manifoldness (reported separately)','wall thickness','fluid connectivity','support adequacy','SLA print success','physical flow or drainage'],elapsed_seconds=time.monotonic()-start)
 out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
 if not a.keep_input:inp.unlink()
 print(json.dumps({k:report[k] for k in ['source_sha256','abnormal_pair_count','no_abnormal_intersections_and_no_exact_degenerate_faces','elapsed_seconds']}),flush=True)
if __name__=='__main__':main()
