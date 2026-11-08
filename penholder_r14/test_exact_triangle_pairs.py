#!/usr/bin/env python3
"""Analytic regression cases for the exact predicate, including shared-vertex crossing."""
from pathlib import Path
import itertools,json
import numpy as np
from run_exact_triangle_pairs import run,CLASSES
O=Path(__file__).parent/'exact_predicate_tests';O.mkdir(exist_ok=True)
a=np.array([[0,0,0],[2,0,0],[0,2,0]],dtype=float)
cases=[
 ('parallel_disjoint',a,a+[0,0,1],0),
 ('expected_noncoplanar_shared_edge',a,[[0,0,0],[2,0,0],[0,0,2]],1),
 ('expected_shared_vertex',a,[[0,0,0],[-2,0,1],[0,-2,1]],1),
 ('proper_crossing',[[-2,-2,0],[2,-2,0],[0,2,0]],[[0,-1,-1],[0,-1,1],[0,1,0]],3),
 ('shared_vertex_and_crossing_elsewhere',a,[[0,0,0],[1,1,-1],[1,1,1]],3),
 ('coplanar_area_overlap',a,[[.5,.5,0],[2,.5,0],[.5,2,0]],4),
 ('expected_coplanar_shared_edge',a,[[0,0,0],[2,0,0],[0,-2,0]],1),
 ('shared_edge_but_area_overlap',a,[[0,0,0],[2,0,0],[1,1,0]],4),
 ('coplanar_shared_vertex_only',a,[[0,0,0],[-1,0,0],[0,-1,0]],1),
 ('unshared_point_touch',a,[[1,1,0],[2,2,0],[1,2,0]],2),
 ('unshared_partial_edge_touch',a,[[.5,0,0],[1.5,0,0],[1,-1,0]],5),
 ('tiny_exact_parallel_gap',a,a+[0,0,2.**-80],0),
 ('degenerate_triangle',[[0,0,0],[1,0,0],[2,0,0]],a,6),
 ('identical_duplicate',a,a,4),
 ('tiny_exact_crossing',a,[[.5,.5,-2.**-80],[.5,.5,2.**-80],[1,1,0]],3),
]
tris=[];pairs=[];expected=[];labels=[]
for name,x,y,code in cases:
 for pa in itertools.permutations(range(3)):
  for pb in itertools.permutations(range(3)):
   i=len(tris);tris.extend([np.asarray(x)[list(pa)],np.asarray(y)[list(pb)]]);pairs.append([i,i+1]);expected.append(code);labels.append(name)
codes,log=run(tris,pairs,O/'analytic_cases')
fail=[dict(case=labels[i],expected=CLASSES[expected[i]],actual=CLASSES[int(codes[i])]) for i in range(len(codes)) if codes[i]!=expected[i]]
r=dict(cases=len(cases),orientation_permutations=36,total_checks=len(codes),passed=int(sum(codes==np.asarray(expected))),failures=fail,engine_stdout=log)
(O/'analytic_cases.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2));assert not fail
