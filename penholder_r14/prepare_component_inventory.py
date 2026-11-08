#!/usr/bin/env python3
"""Record all saved STL components and verify optional principal extraction."""
from pathlib import Path
import argparse,json,hashlib,time,gc
import numpy as np,trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('mesh');ap.add_argument('--out',required=True);ap.add_argument('--principal');a=ap.parse_args()
 out=Path(a.out);out.mkdir(exist_ok=True,parents=True);t=time.monotonic();src=Path(a.mesh)
 m=trimesh.load_mesh(src,process=False);v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];del m,inv;gc.collect()
 e=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);g=coo_matrix((np.ones(len(e),dtype=np.uint8),(e[:,0],e[:,1])),shape=(len(v),len(v))).tocsr();nc,lab=connected_components(g,directed=False);flab=lab[f[:,0]];del e,g;gc.collect()
 np.savez_compressed(out/'face_component_labels.npz',labels=flab)
 tri=v[f];vol=np.einsum('ij,ij->i',tri[:,0],np.cross(tri[:,1],tri[:,2]))/6;area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
 counts=np.bincount(flab,minlength=nc);vs=np.bincount(flab,weights=vol,minlength=nc);ars=np.bincount(flab,weights=area,minlength=nc);rows=[]
 for i in np.argsort(counts)[::-1]:
  vv=v[lab==i];rows.append(dict(component=int(i),faces=int(counts[i]),vertices=len(vv),signed_volume_mm3=float(vs[i]),surface_area_mm2=float(ars[i]),bounds_saved_mm=[vv.min(axis=0).tolist(),vv.max(axis=0).tolist()]))
 mainidx=rows[0]['component'];report=dict(source=str(src),sha256=sha(src),coordinate_welding='exact identical saved coordinates only',connected_surface_components=nc,primary_component=mainidx,components=rows,secondary_summary=dict(components=nc-1,positive_count=sum(c['signed_volume_mm3']>0 for c in rows[1:]),negative_count=sum(c['signed_volume_mm3']<0 for c in rows[1:]),positive_volume_mm3=sum(max(0,c['signed_volume_mm3']) for c in rows[1:]),negative_signed_volume_mm3=sum(min(0,c['signed_volume_mm3']) for c in rows[1:]),positive_surface_area_mm2=sum(c['surface_area_mm2'] for c in rows[1:] if c['signed_volume_mm3']>0),negative_surface_area_mm2=sum(c['surface_area_mm2'] for c in rows[1:] if c['signed_volume_mm3']<0)))
 if a.principal:
  p=Path(a.principal);pm=trimesh.load_mesh(p,process=False);pt=pm.triangles;expected=tri[flab==mainidx];same=pt.shape==expected.shape and np.array_equal(pt,expected)
  report['principal_extraction_comparison']=dict(file=str(p),sha256=sha(p),faces=len(pt),saved_triangle_coordinates_and_order_exactly_equal=same)
  if not same:
   # Sort face corner coordinates and complete face rows, preserving the exact coordinate set.
   def canonical(x):
    order=np.lexsort((x[:,:,2],x[:,:,1],x[:,:,0]),axis=1);x=np.take_along_axis(x,order[:,:,None],axis=1).reshape(-1,9);idx=np.lexsort(tuple(x[:,k] for k in range(8,-1,-1)));return x[idx]
   report['principal_extraction_comparison']['unoriented_triangle_coordinate_set_exactly_equal']=bool(pt.shape==expected.shape and np.array_equal(canonical(pt),canonical(expected)))
 report['elapsed_seconds']=time.monotonic()-t;(out/'component_inventory.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['connected_surface_components','secondary_summary','elapsed_seconds']}));print(report.get('principal_extraction_comparison'))
if __name__=='__main__':main()
