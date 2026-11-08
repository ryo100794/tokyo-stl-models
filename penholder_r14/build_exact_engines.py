#!/usr/bin/env python3
"""Build exact predicates locally; no package installation or network access."""
from pathlib import Path
import argparse,subprocess
R=Path(__file__).resolve().parent

def ensure_engine(name):
 target=R/name;source=R/(name+'.cpp');header=R/'exact_predicate.hpp'
 if not target.exists() or target.stat().st_mtime<max(source.stat().st_mtime,header.stat().st_mtime):
  cmd=['g++','-std=c++17','-O2']
  if name=='exact_mesh_bvh':cmd.append('-fopenmp')
  cmd.extend([str(source),'-lgmpxx','-lgmp','-o',str(target)])
  subprocess.run(cmd,check=True)
 return target
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--engine',choices=['exact_triangle_pairs','exact_mesh_bvh','all'],default='all');a=ap.parse_args()
 for n in (['exact_triangle_pairs','exact_mesh_bvh'] if a.engine=='all' else [a.engine]):print(ensure_engine(n))
