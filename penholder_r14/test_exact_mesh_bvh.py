#!/usr/bin/env python3
"""Verify BVH pair completeness against exhaustive exact all-pairs on dyadic fixtures."""
from pathlib import Path
import itertools,json
import numpy as np
from run_exact_triangle_pairs import run
from run_exact_mesh_bvh import run_indexed
R=Path(__file__).parent/'exact_predicate_tests';R.mkdir(exist_ok=True)
rng=np.random.default_rng(20261109);tr=rng.integers(-16,17,size=(100,3,3))/8
# Include shared-edge, shared-vertex, coplanar, duplicate, and exact tiny-gap examples.
a=np.array([[0,0,0],[2,0,0],[0,2,0]],float)
extra=[a,a+[0,0,2.**-80],[[0,0,0],[2,0,0],[0,0,2]],[[0,0,0],[1,1,-1],[1,1,1]],a,[[0,0,0],[2,0,0],[1,1,0]]]
tr=np.concatenate([tr,np.asarray(extra)]);pairs=np.asarray(list(itertools.combinations(range(len(tr)),2)),dtype=np.uint32)
codes,log=run(tr,pairs,R/'bvh_brute_force');expected={(int(p[0]),int(p[1]),int(c)) for p,c in zip(pairs,codes) if c>=2}
v,inv=np.unique(tr.reshape(-1,3),axis=0,return_inverse=True);f=inv.reshape(-1,3);checks=[]
for threads in [1,3]:
 stats,hits,_,_=run_indexed(v,f,R/f'bvh_threads_{threads}',threads)
 actual={(int(x['face_a']),int(x['face_b']),int(x['code'])) for x in hits};checks.append(dict(threads=threads,match=actual==expected,expected_abnormal_pairs=len(expected),actual_abnormal_pairs=len(actual),missing=len(expected-actual),extra=len(actual-expected),engine=stats));assert actual==expected
report=dict(triangles=len(tr),brute_force_pairs=len(pairs),brute_force_engine=log,checks=checks)
(R/'bvh_completeness.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
