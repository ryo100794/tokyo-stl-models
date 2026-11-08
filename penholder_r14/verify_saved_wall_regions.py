#!/usr/bin/env python3
"""Diagnostic wall classification on serialized R14 standalone line STLs.

No geometry is changed. Input centerlines/ports only identify intended fluid
regions and openings; all thickness/intersection evidence is from STL triangles.
Finite surface sampling is not a global thickness proof or printing approval.
"""
import argparse, hashlib, json, math, time
from pathlib import Path
import numpy as np
import trimesh
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray


def locator_for(mesh):
    p=vtk.vtkPoints();p.SetData(numpy_to_vtk(np.asarray(mesh.vertices),deep=True))
    f=np.asarray(mesh.faces,dtype=np.int64);c=vtk.vtkCellArray()
    c.SetData(numpy_to_vtkIdTypeArray(np.arange(0,f.size+1,3),deep=True),numpy_to_vtkIdTypeArray(f.ravel(),deep=True))
    poly=vtk.vtkPolyData();poly.SetPoints(p);poly.SetPolys(c)
    loc=vtk.vtkStaticCellLocator();loc.SetDataSet(poly);loc.BuildLocator();return loc


def first_hit(loc,a,b):
    t=vtk.reference(0.);sub=vtk.reference(0);cid=vtk.reference(0);hit=[0.,0.,0.];pc=[0.,0.,0.]
    yes=loc.IntersectWithLine(a,b,1e-7,t,hit,pc,sub,cid)
    return (np.array(hit),int(cid)) if yes else None


def basis(u):
    n=np.cross(u,np.eye(3)[np.argmin(np.abs(u))]);n/=np.linalg.norm(n)
    return n,np.cross(u,n)


def nearest_axis(q,starts,delta,denom):
    t=np.clip(((q-starts)*delta).sum(1)/denom,0,1)
    v=q-starts-t[:,None]*delta
    j=int(np.argmin((v*v).sum(1)))
    return v[j],float(np.linalg.norm(v[j])),j


def escape_via_port(c,d,ports,maxdist):
    for p in ports:
        axis=np.array(p['axis_enu']);mouth=np.array(p['mouth_center_enu_mm']);den=axis@d
        if abs(den)<1e-9:continue
        t=float((mouth-c)@axis/den)
        if t<0 or t>maxdist:continue
        q=c+t*d
        if np.linalg.norm(q-mouth)<=p['mouth_inner_diameter_mm']/2+.03:
            return p['port_id']
    return None


def audit(line,ports,path,T,angles=8,threshold=.5):
    started=time.time();mesh=trimesh.load_mesh(path,process=True);mesh.apply_transform(np.linalg.inv(T));loc=locator_for(mesh)
    normals=np.asarray(mesh.face_normals);rawpaths=[np.array(t['points_enu_mm']) for t in line['paths']]
    starts=np.concatenate([p[:-1] for p in rawpaths]);delta=np.concatenate([np.diff(p,axis=0) for p in rawpaths]);denom=np.maximum((delta*delta).sum(1),1e-24)
    nports=np.array([p['attachment_enu_mm'] for p in ports]);exclude=np.array([p['length_mm']+2 for p in ports]);maxdist=float(np.linalg.norm(mesh.extents))+10
    rows={'line_id':line['line_id'],'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'ordinary_rays':0,'near_port_rays_excluded':0,'beyond_3mm_bore_hits':0,'missing_boundary':0,'escape_through_intended_mouth':0,'same_line_septum_intervals':0,'ambiguous_intervals':0,'containment_radial_intervals':0,'minimum_containment_radial_mm':None,'minimum_containment_normal_mm':None,'minimum_septum_radial_mm':None,'radial_below_threshold_count':0,'normal_below_threshold_count':0,'examples':{'radial_thin':[],'normal_thin':[],'missing':[],'septum':[]},'sampling':{'angles':angles,'source_segment_midpoints':len(starts),'threshold_mm':threshold,'near_port_exclusion_mm':'length + 2.0','finite_sampling_only':True}}
    radial=[];normalwall=[];septa=[]
    def classify(q,cid):
        vec,dist,j=nearest_axis(q,starts,delta,denom)
        cosine=float(normals[cid]@vec/max(dist,1e-20))
        # This is a classification diagnostic, with undecidable samples retained.
        kind='outer' if cosine>.15 else 'inner' if cosine<-.15 else 'ambiguous'
        return kind,dist,cosine
    for trace,p in zip(line['paths'],rawpaths):
        for k,(a,b) in enumerate(zip(p[:-1],p[1:])):
            ab=b-a;ln=np.linalg.norm(ab)
            if ln<1e-10:continue
            center=(a+b)/2
            if np.any(np.linalg.norm(nports-center,axis=1)<exclude):rows['near_port_rays_excluded']+=angles;continue
            u=ab/ln;n,q=basis(u)
            for angle in np.arange(angles)*2*np.pi/angles:
                d=np.cos(angle)*n+np.sin(angle)*q;rows['ordinary_rays']+=1
                first=first_hit(loc,center,center+maxdist*d)
                base={'trace_id':trace['trace_id'],'segment':k,'angle_rad':float(angle),'center_mm':center.tolist(),'direction':d.tolist()}
                if first is None:
                    port=escape_via_port(center,d,ports,maxdist)
                    if port:rows['escape_through_intended_mouth']+=1
                    else:
                        rows['missing_boundary']+=1
                        if len(rows['examples']['missing'])<12:rows['examples']['missing'].append(base)
                    continue
                a1,cell1=first;dist1=float(np.linalg.norm(a1-center));kind1,axisdist1,dot1=classify(a1,cell1)
                if dist1>3:rows['beyond_3mm_bore_hits']+=1
                second=first_hit(loc,a1+1e-4*d,center+maxdist*d)
                if second is None:
                    rows['ambiguous_intervals']+=1;continue
                a2,cell2=second;thick=float(np.linalg.norm(a2-a1));kind2,axisdist2,dot2=classify(a2,cell2)
                detail={**base,'thickness_mm':thick,'entry_mm':a1.tolist(),'exit_mm':a2.tolist(),'entry_face':cell1,'exit_face':cell2,'entry_kind':kind1,'exit_kind':kind2,'entry_axis_distance_mm':axisdist1,'exit_axis_distance_mm':axisdist2,'entry_normal_dot_axis_radial':dot1,'exit_normal_dot_axis_radial':dot2}
                if kind1!='inner':rows['ambiguous_intervals']+=1;continue
                if kind2=='inner':
                    septa.append(thick);rows['same_line_septum_intervals']+=1
                    if thick<threshold and len(rows['examples']['septum'])<12:rows['examples']['septum'].append(detail)
                elif kind2=='outer':
                    radial.append(thick);rows['containment_radial_intervals']+=1
                    if thick<threshold:
                        rows['radial_below_threshold_count']+=1
                        if len(rows['examples']['radial_thin'])<12:rows['examples']['radial_thin'].append(detail)
                else:rows['ambiguous_intervals']+=1
                # Probe from this actual lumen surface in the opposite face-normal
                # direction; radial interval alone is not the physical shortest wall.
                into=-normals[cell1];exitn=first_hit(loc,a1+1e-4*into,a1+maxdist*into)
                if exitn is not None:
                    en,cn=exitn;kn,dn,dotn=classify(en,cn)
                    if kn=='outer':
                        nt=float(np.linalg.norm(en-a1));normalwall.append(nt)
                        if nt<threshold:
                            rows['normal_below_threshold_count']+=1
                            if len(rows['examples']['normal_thin'])<12:rows['examples']['normal_thin'].append({**base,'thickness_mm':nt,'entry_mm':a1.tolist(),'exit_mm':en.tolist(),'entry_face':cell1,'exit_face':cn,'exit_axis_distance_mm':dn})
    rows['minimum_containment_radial_mm']=min(radial) if radial else None
    rows['minimum_containment_normal_mm']=min(normalwall) if normalwall else None
    rows['minimum_septum_radial_mm']=min(septa) if septa else None
    rows['diagnostic_status']='fail' if rows['missing_boundary'] or rows['radial_below_threshold_count'] or rows['normal_below_threshold_count'] else 'sampled_checks_clear_with_limits'
    rows['seconds']=time.time()-started
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--lines',nargs='*');ap.add_argument('--angles',type=int,default=8);a=ap.parse_args()
    root=a.root;data=json.loads((root/'data/rail_paths_R14.json').read_text());ports=json.loads((root/'data/ports_R14.json').read_text())['ports'];T=np.array(json.loads((root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    result={'scope':'diagnostic on serialized standalone line STLs; no whole-body certification','method':'Oriented ray surfaces classified against saved centerline directions; same-line cavity septa separated from outer fluid containment. Missing crossings searched beyond model bounds, with planned mouth escape detected. Normal rays sampled at actual lumen hits.','global_wall_thickness_proven':False,'physical_tests_performed':False,'lines':[]}
    for line in data['lines']:
        if a.lines and line['line_id'] not in a.lines:continue
        row=audit(line,[p for p in ports if p['line_id']==line['line_id']],root/'models'/f"rail_R14_{line['line_id']}.stl",T,a.angles)
        result['lines'].append(row);a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
        print(row['line_id'],'radial',row['minimum_containment_radial_mm'],'normal',row['minimum_containment_normal_mm'],'missing',row['missing_boundary'],'mouth',row['escape_through_intended_mouth'],'septum',row['minimum_septum_radial_mm'],'ambig',row['ambiguous_intervals'],flush=True)
if __name__=='__main__':main()
