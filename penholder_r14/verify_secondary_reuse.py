#!/usr/bin/env python3
"""Permit reuse of core-separation results only for exactly unchanged oriented secondary triangles and identical cutters."""
from pathlib import Path
import argparse,json,hashlib,struct
import numpy as np

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def triangles(path,labels,primary):
 with Path(path).open('rb') as f:f.seek(80);n=struct.unpack('<I',f.read(4))[0]
 dt=np.dtype([('normal','<f4',(3,)),('v','<f4',(3,3)),('attr','<u2')]);assert Path(path).stat().st_size==84+50*n
 mm=np.memmap(path,mode='r',offset=84,dtype=dt,shape=(n,));lab=np.load(labels)['labels'];assert len(lab)==n
 tr=np.asarray(mm['v'][lab!=primary],dtype='<f4').copy();tr[tr==0]=0
 # Canonicalize by cyclic rotation only, retaining triangle orientation.
 first=np.lexsort((tr[:,:,2],tr[:,:,1],tr[:,:,0]),axis=1)[:,0];idx=(first[:,None]+np.arange(3))%3
 tr=np.take_along_axis(tr,idx[:,:,None],axis=1).reshape(-1,9);order=np.lexsort(tuple(tr[:,k] for k in range(8,-1,-1)))
 return np.ascontiguousarray(tr[order])
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--previous-inventory',required=True);ap.add_argument('--previous-labels',required=True);ap.add_argument('--previous-core-report',required=True);ap.add_argument('--new-inventory',required=True);ap.add_argument('--new-labels',required=True);ap.add_argument('--current-core-dir',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
 old=json.loads(Path(a.previous_inventory).read_text());new=json.loads(Path(a.new_inventory).read_text());prev=json.loads(Path(a.previous_core_report).read_text());assert old['sha256']==sha(old['source']) and new['sha256']==sha(new['source']);assert prev['source_sha256']==old['sha256']
 ot=triangles(old['source'],a.previous_labels,old['primary_component']);nt=triangles(new['source'],a.new_labels,new['primary_component'])
 same=ot.shape==nt.shape and np.array_equal(ot,nt);checks=[]
 for c in prev['cores']:
  p=Path(a.current_core_dir)/(c['line_id']+'_core.npz');h=sha(p);checks.append(dict(line_id=c['line_id'],previous_sha256=c['sha256'],current_sha256=h,identical=h==c['sha256']))
 report=dict(schema='R14-exact-secondary-triangle-core-result-reuse-v1',previous_source_sha256=old['sha256'],new_source_sha256=new['sha256'],previous_secondary_triangles=len(ot),new_secondary_triangles=len(nt),oriented_secondary_triangle_multisets_exactly_identical=same,previous_canonical_oriented_triangle_sha256=hashlib.sha256(ot.tobytes()).hexdigest(),new_canonical_oriented_triangle_sha256=hashlib.sha256(nt.tobytes()).hexdigest(),canonicalization='Exact saved float32 coordinates; normalize signed zero; cyclically rotate each triangle only, preserving orientation; lexicographically sort triangle records.',core_checks=checks,previous_core_separation_report=str(Path(a.previous_core_report)),previous_core_report_sha256=sha(a.previous_core_report),core_separation_reuse_gate_passed=bool(same and all(c['identical'] for c in checks) and prev['all_separation_gates_passed']),reused_minimum_gap_mm=prev['minimum_reported_gap_mm'],reused_component_count=prev['subsidiary_components_tested'],limitation='Only the unchanged secondary volumes and unchanged final cutter separation are reused. Containment against the changed principal boundary requires a fresh nesting test; final saved principal topology/intersections/flow remain separate gates.')
 Path(a.out).write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['previous_secondary_triangles','new_secondary_triangles','oriented_secondary_triangle_multisets_exactly_identical','core_separation_reuse_gate_passed','reused_minimum_gap_mm']}),flush=True)
if __name__=='__main__':main()
