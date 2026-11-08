#!/usr/bin/env python3
"""All-vertex containment plus exact cross-boundary triangle pairs for subsidiary shells."""
from pathlib import Path
import argparse,json,time,hashlib,gc
import numpy as np,trimesh,vtk
from rtree import index
from mesh_validation_utils import subset,locator,closest_vertices
from run_exact_triangle_pairs import run,CLASSES

def main():
 ap=argparse.ArgumentParser();ap.add_argument('mesh');ap.add_argument('--face-labels',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
 out=Path(a.out);out.mkdir(exist_ok=True,parents=True);src=Path(a.mesh);start=time.monotonic()
 m=trimesh.load_mesh(src,process=False);v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];del m,inv;gc.collect()
 with np.load(a.face_labels) as z:lab=z['labels']
 assert len(lab)==len(f)
 counts=np.bincount(lab);mainidx=int(counts.argmax());mainids=np.flatnonzero(lab==mainidx);otherids=np.flatnonzero(lab!=mainidx)
 mv,mf=subset(v,f,lab==mainidx);loc,poly=locator(mv,mf);inside=vtk.vtkSelectEnclosedPoints();inside.SetTolerance(1e-10);inside.Initialize(poly)
 rows=[]
 for i in np.flatnonzero(counts):
  if i==mainidx:continue
  ids=np.unique(f[lab==i].ravel());vv=v[ids]
  cls=np.asarray([int(inside.IsInsideSurface(q)) for q in vv]);ds,cp,c=closest_vertices(vv,loc)
  rows.append(dict(component=int(i),all_boundary_vertices=len(ids),vertices_inside_primary=int(cls.sum()),vertices_outside_primary=int((cls==0).sum()),minimum_vertex_distance_to_primary_mm=float(ds.min()),outside_vertex_ids=ids[cls==0].tolist()))
 print('ALL VERTICES DONE',len(rows),time.monotonic()-start,flush=True)
 tr=v[f];lo=tr.min(axis=1);hi=tr.max(axis=1)
 prop=index.Property();prop.dimension=3
 rt=index.Index(((int(i),tuple(np.r_[lo[i],hi[i]]),None) for i in mainids),properties=prop)
 print('TREE',time.monotonic()-start,flush=True)
 pairs=np.asarray([(int(i),j) for i in otherids for j in rt.intersection(tuple(np.r_[lo[i],hi[i]]))],dtype=np.uint32).reshape(-1,2)
 print('PAIRS',len(pairs),time.monotonic()-start,flush=True)
 allids=np.unique(pairs.ravel());localpairs=np.searchsorted(allids,pairs)
 codes,log=run(tr[allids],localpairs,out/'cross_boundary_pairs')
 np.savez_compressed(out/'cross_boundary_exact_pairs.npz',face_id_pairs=pairs,classification_codes=codes)
 for r in rows:
  mask=lab[pairs[:,0]]==r['component'];bad=mask&(codes>=2)
  r['exact_cross_boundary_candidate_pairs']=int(mask.sum());r['exact_cross_boundary_abnormal_pairs']=int(bad.sum())
  r['containment_gate_passed']=bool(r['vertices_outside_primary']==0 and not bad.any())
 report=dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),primary_component=mainidx,all_subsidiary_components=len(rows),all_subsidiary_boundary_vertices_tested=sum(r['all_boundary_vertices'] for r in rows),exact_cross_boundary_candidate_pairs=len(pairs),exact_counts={CLASSES[k]:int((codes==k).sum()) for k in CLASSES},containment_gate_passed_all=all(r['containment_gate_passed'] for r in rows),components=rows,engine_stdout=log,methods=['Every boundary vertex of every subsidiary component classified by VTK ray parity, tolerance 1e-10 of model diagonal.','Every subsidiary triangle queried against all primary triangle AABBs using closed-coordinate boxes; every candidate tested with GMP exact integer/rational triangle intersection predicate.','Shared contacts count as normal only if fully contained in exact shared vertex/edge.'],limitations=['This report does not establish primary self-intersection freedom; require separate primary topology and exhaustive exact-intersection validation.','No shell is deleted by this script. Noncontained or cross-boundary-intersecting components must not be removed solely by size.','Primary-to-subsidiary checks are exhaustive; primary self-pairs and subsidiary self-pairs are not part of this particular report.'],elapsed_seconds=time.monotonic()-start)
 (out/'nesting_verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
 print('DONE',json.dumps({k:report[k] for k in ['all_subsidiary_components','all_subsidiary_boundary_vertices_tested','exact_cross_boundary_candidate_pairs','exact_counts','containment_gate_passed_all','elapsed_seconds']}),flush=True)
if __name__=='__main__':main()
