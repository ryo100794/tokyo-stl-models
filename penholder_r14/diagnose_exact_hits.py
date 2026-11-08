#!/usr/bin/env python3
"""Locate exact triangle-pair findings against intended core and retained base."""
from pathlib import Path
import argparse,json,hashlib,struct
from fractions import Fraction as F
import numpy as np
from mesh_validation_utils import load_npz_arrays,locator,closest_vertices
from run_exact_triangle_pairs import CLASSES

def coplanar_area(A,B):
 n=np.cross(A[1]-A[0],A[2]-A[0]);n/=np.linalg.norm(n);drop=int(np.argmax(abs(n)));cols=[k for k in range(3) if k!=drop]
 a=[[F(float(q[k])) for k in cols] for q in A];b=[[F(float(q[k])) for k in cols] for q in B]
 def sub(x,y):return [x[0]-y[0],x[1]-y[1]]
 def cross(x,y):return x[0]*y[1]-x[1]*y[0]
 if cross(sub(b[1],b[0]),sub(b[2],b[0]))<0:b=b[::-1]
 poly=a
 for j in range(3):
  s,t=b[j],b[(j+1)%3];edge=sub(t,s)
  def side(x):return cross(edge,sub(x,s))
  out=[]
  if not poly:break
  prev=poly[-1];sp=side(prev)
  for cur in poly:
   sc=side(cur)
   if (sc>=0)!=(sp>=0):
    r=sp/(sp-sc);out.append([prev[k]+r*(cur[k]-prev[k]) for k in range(2)])
   if sc>=0:out.append(cur)
   prev=cur;sp=sc
  poly=out
 ar=abs(sum(poly[i][0]*poly[(i+1)%len(poly)][1]-poly[(i+1)%len(poly)][0]*poly[i][1] for i in range(len(poly))))/2 if len(poly)>2 else F(0)
 return dict(projected_area_exact_mm2_fraction=str(ar),surface_overlap_area_mm2=float(ar)/abs(n[drop]))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('mesh');ap.add_argument('--exact-report',required=True);ap.add_argument('--exact-pairs',required=True);ap.add_argument('--print-transform',required=True);ap.add_argument('--core-dir',required=True);ap.add_argument('--base-dir');ap.add_argument('--base-report');ap.add_argument('--out',required=True);a=ap.parse_args()
 src=Path(a.mesh);out=Path(a.out);out.mkdir(exist_ok=True,parents=True);sha=hashlib.sha256(src.read_bytes()).hexdigest();rep=json.loads(Path(a.exact_report).read_text());assert rep['source_sha256']==sha
 with np.load(a.exact_pairs) as z:pairs=z['face_id_pairs'];codes=z['classification_codes']
 ids=np.unique(pairs);dtype=np.dtype([('normal','<f4',(3,)),('vertices','<f4',(3,3)),('attr','<u2')])
 with src.open('rb') as s:s.seek(80);nf=struct.unpack('<I',s.read(4))[0]
 assert src.stat().st_size==84+50*nf;m=np.memmap(src,mode='r',offset=84,dtype=dtype,shape=(nf,));tr=np.asarray(m['vertices'][ids],dtype=np.float64)
 T=np.asarray(json.loads(Path(a.print_transform).read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major']);Ti=np.linalg.inv(T);enu=np.einsum('...j,ij->...i',tr,Ti[:3,:3])+Ti[:3,3]
 np.savez_compressed(out/'selected_triangles.npz',face_ids=ids,triangles_saved_mm=tr,triangles_enu_mm=enu,exact_face_pairs=pairs,classification_codes=codes)
 centers=enu.mean(axis=1);pts=np.vstack([centers,enu.reshape(-1,3)]);mat=np.full((len(pts),14),np.inf);normals=np.zeros((len(pts),14,3));core_manifest=[]
 for j,p in enumerate(sorted(Path(a.core_dir).glob('L[0-9][0-9]_core.npz'))):
  v,f=load_npz_arrays(p);loc,poly=locator(v,f);ds,cp,ci=closest_vertices(pts,loc);mat[:,j]=ds;nn=np.cross(v[f[ci,1]]-v[f[ci,0]],v[f[ci,2]]-v[f[ci,0]]);nn/=np.linalg.norm(nn,axis=1)[:,None];normals[:,j]=nn;core_manifest.append(dict(line_id=p.stem[:3],sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 base_mat=np.full((len(centers),43),np.inf);base_names={}
 if a.base_dir:
  if a.base_report:base_names={g['index']:g['name'] for g in json.loads(Path(a.base_report).read_text())['groups']}
  for p in sorted(Path(a.base_dir).glob('base_[0-9][0-9].npz')):
   i=int(p.stem[-2:]);v,f=load_npz_arrays(p);lo=np.linalg.norm(np.maximum(np.maximum(v.min(axis=0)-centers,centers-v.max(axis=0)),0),axis=1);which=np.flatnonzero(lo<=base_mat.min(axis=1))
   if len(which):loc,poly=locator(v,f);ds,_,_=closest_vertices(centers[which],loc);base_mat[which,i]=ds
 nn=np.cross(enu[:,1]-enu[:,0],enu[:,2]-enu[:,0]);nn/=np.linalg.norm(nn,axis=1)[:,None];rows=[]
 for i,face in enumerate(ids):
  j=int(mat[i].argmin());k=int(base_mat[i].argmin());pv=len(ids)+3*i+np.arange(3)
  rows.append(dict(face_id=int(face),triangle_saved_mm=tr[i].tolist(),triangle_enu_mm=enu[i].tolist(),centroid_enu_mm=centers[i].tolist(),nearest_core_line=core_manifest[j]['line_id'],centroid_distance_to_core_surface_mm=float(mat[i,j]),vertex_distances_to_same_core_mm=mat[pv,j].tolist(),dot_body_normal_core_outward_normal=float(nn[i]@normals[i,j]),nearest_base_group=k if a.base_dir else None,nearest_base_name=base_names.get(k),centroid_distance_to_nearest_base_surface_mm=float(base_mat[i,k]) if a.base_dir else None))
 prows=[]
 for pair,code in zip(pairs,codes):
  ia,ib=np.searchsorted(ids,pair);A,B=tr[ia],tr[ib];vv=np.linalg.norm(A[:,None,:]-B[None,:,:],axis=2);nz=vv[vv>0];row=dict(face_ids=pair.tolist(),classification=CLASSES[int(code)],exact_common_vertex_count=int((vv==0).sum()),minimum_distinct_vertex_distance_mm=float(nz.min()) if len(nz) else 0,midpoint_enu_mm=((centers[ia]+centers[ib])/2).tolist())
  if code==4:row['coplanar_overlap']=coplanar_area(A,B)
  touching=[]
  for typ,X,Y in [('Avertex_Bedge',A,B),('Bvertex_Aedge',B,A)]:
   for i,x in enumerate(X):
    for j in range(3):
     y=Y[j];delta=Y[(j+1)%3]-y;t=float(np.dot(x-y,delta)/np.dot(delta,delta));tclip=np.clip(t,0,1);dd=float(np.linalg.norm(x-(y+tclip*delta)))
     if dd<1e-10 and 1e-10<t<1-1e-10:touching.append(dict(type=typ,vertex=i,edge_start=j,edge_parameter=t,distance_mm=dd))
  row['vertex_interior_edge_contacts_float_screen']=touching;prows.append(row)
 report=dict(source=str(src),source_sha256=sha,exact_counts=rep.get('class_counts',rep.get('counts')),affected_faces=len(ids),core_manifest=core_manifest,faces=rows,pairs=prows,scope='Coordinate diagnosis only. Exact classifications are inherited from the saved full-mesh integer/rational test. Nearest-surface distances and contact scales are floating-point measurements.')
 (out/'diagnosis.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
 print('FACES',len(rows),'PAIRS',len(prows))
 for r in rows:print(r['face_id'],r['centroid_enu_mm'],'core',r['nearest_core_line'],r['centroid_distance_to_core_surface_mm'],'base',r['nearest_base_name'],r['centroid_distance_to_nearest_base_surface_mm'])
 print('PAIR_SCALES',[(r['face_ids'],r['classification'],r['minimum_distinct_vertex_distance_mm'],r.get('coplanar_overlap')) for r in prows])
if __name__=='__main__':main()
