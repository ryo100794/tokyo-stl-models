#!/usr/bin/env python3
import argparse,json,hashlib
from pathlib import Path
import numpy as np

def check(path):
 p=Path(path);raw=p.read_bytes();n=int.from_bytes(raw[80:84],'little')
 if len(raw)!=84+50*n:return {'file':str(p),'status':'binary STL layout unsupported by this checker'}
 dt=np.dtype([('normal','<f4',(3,)),('vertices','<f4',(3,3)),('attribute','<u2')])
 rec=np.frombuffer(raw,dtype=dt,offset=84,count=n);v=rec['vertices'].astype(np.float64);normal=rec['normal'].astype(np.float64)
 gn=np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0]);a=np.linalg.norm(gn,axis=1);b=np.linalg.norm(normal,axis=1);good=(a>0)&(b>0)
 cos=np.einsum('ij,ij->i',gn[good],normal[good])/(a[good]*b[good])
 return {'file':str(p),'sha256':hashlib.sha256(raw).hexdigest(),'faces':n,'nonfinite_stored_normals':int((~np.isfinite(normal).all(axis=1)).sum()),'zero_stored_normals':int((b==0).sum()),'zero_geometric_normals':int((a==0).sum()),'stored_normals_opposed_to_vertex_winding':int((cos<0).sum()),'minimum_normal_dot_cosine':float(cos.min()) if len(cos) else None,'cosine_below_0_99999_count':int((cos<.99999).sum()),'input_modified':False,'scope':'binary STL normal records compared to cross product of stored triangle vertices; does not repair winding or certify physical printing'}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('inputs',nargs='+');ap.add_argument('--out',required=True);a=ap.parse_args();r=[check(p) for p in a.inputs]
 Path(a.out).write_text(json.dumps({'checks':r},indent=2)+'\n');print(json.dumps([{k:q.get(k) for k in ['file','stored_normals_opposed_to_vertex_winding','zero_stored_normals','cosine_below_0_99999_count']} for q in r],indent=2))
