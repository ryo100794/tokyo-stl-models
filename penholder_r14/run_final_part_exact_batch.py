#!/usr/bin/env python3
"""Run the independent exhaustive exact mesh validator once on R14's 15 parts."""
import argparse,concurrent.futures,hashlib,json,subprocess,sys,time
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--validator-dir',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--jobs',type=int,default=2);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 source_names=['run_exact_mesh_bvh.py','exact_mesh_bvh.cpp','exact_predicate.hpp','build_exact_engines.py','run_exact_triangle_pairs.py'];sources={n:sha(a.validator_dir/n) for n in source_names}
 files=[a.root/'models'/f'rail_R14_L{i:02d}.stl' for i in range(1,15)]+[a.root/'models/flow_coupon_R14.stl'];initial={p.name:sha(p) for p in files};start=time.time();summary={'scope':'Full saved-mesh exact triangle intersection scan of 14 final standalone line STLs and the unchanged coupon; no geometry changes','validator_sources_sha256':sources,'source_script_sha256':sha(__file__),'parts':[],'completed':False}
 def run(p):
  out=a.out/f'{p.stem}_exact_intersections.json';log=a.out/f'{p.stem}.log'
  with log.open('w') as stream:r=subprocess.run([sys.executable,str(a.validator_dir/'run_exact_mesh_bvh.py'),str(p),'--out',str(out),'--threads','1'],stdout=stream,stderr=stream)
  if r.returncode:return {'file':p.name,'returncode':r.returncode,'log':log.name,'input_sha256':initial[p.name],'error':'independent validator process failed'}
  q=json.loads(out.read_text());return {'file':p.name,'report':out.name,'report_sha256':sha(out),'source_sha256':q['source_sha256'],'input_hash_unchanged':sha(p)==initial[p.name]==q['source_sha256'],'abnormal_pair_count':q['abnormal_pair_count'],'degenerate_face_count':len(q['degenerate_face_ids']),'class_counts':q['class_counts'],'no_abnormal_intersections_and_no_exact_degenerate_faces':q['no_abnormal_intersections_and_no_exact_degenerate_faces'],'elapsed_seconds':q['elapsed_seconds']}
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
  pending=[pool.submit(run,p) for p in files]
  for future in concurrent.futures.as_completed(pending):
   row=future.result();summary['parts'].append(row);summary['parts'].sort(key=lambda x:x['file']);(a.out/'final15_exact_intersection_summary_R14.json').write_text(json.dumps(summary,indent=2)+'\n');print(row['file'],'abnormal',row.get('abnormal_pair_count'),'degenerate',row.get('degenerate_face_count'),'error',row.get('error'),flush=True)
 summary['completed']=True;summary['all_source_mesh_bytes_unchanged']=all(x.get('input_hash_unchanged',False) for x in summary['parts']);summary['all_15_pass']=len(summary['parts'])==15 and all(x.get('no_abnormal_intersections_and_no_exact_degenerate_faces',False) for x in summary['parts']);summary['validator_sources_unchanged']=all(sha(a.validator_dir/n)==h for n,h in sources.items());summary['elapsed_seconds']=time.time()-start;(a.out/'final15_exact_intersection_summary_R14.json').write_text(json.dumps(summary,indent=2)+'\n');print('SUMMARY',summary['all_15_pass'],summary['all_source_mesh_bytes_unchanged'],'seconds',summary['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
