#!/usr/bin/env python3
"""Finite checks of the unchanged R14 coupon's actual saved channels/openings.

Uses the coupon's own XYZ millimetres directly; no railway print transform.
No mesh is changed. Numerical samples are not a physical test or a continuous
minimum-wall or leakage certificate.
"""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np,trimesh
from verify_saved_wall_regions import locator_for,first_hit,basis


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def densify(points,step=.4):
    out=[np.array(points[0],float)]
    for a,b in zip(points[:-1],points[1:]):
        a=np.array(a);b=np.array(b);n=max(1,int(np.ceil(np.linalg.norm(b-a)/step)))
        out.extend(a+(b-a)*t/n for t in range(1,n+1))
    return np.array(out)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.parent.mkdir(exist_ok=True,parents=True)
    spec_path=a.root/'data/flow_coupon_R14_specification.json';mesh_path=a.root/'models/flow_coupon_R14.stl';spec=json.loads(spec_path.read_text());initial=sha(mesh_path);mesh=trimesh.load_mesh(mesh_path,process=True);loc=locator_for(mesh);maxdist=float(np.linalg.norm(mesh.extents))+10
    channels=spec['channels'];geoms=[];allports=[]
    for ch in channels:
        p=np.array(ch['centerline_mm']);starts=p[:-1];delta=np.diff(p,axis=0);den=np.maximum(np.sum(delta*delta,axis=1),1e-24);geoms.append((ch['id'],starts,delta,den,ch['inner_diameter_mm']/2))
        for po in ch['ports']:allports.append((ch,po))
    def fluid_face(q,n):
        for cid,starts,delta,den,r in geoms:
            t=np.clip(np.sum((q-starts)*delta,axis=1)/den,0,1);v=q-starts-t[:,None]*delta;j=np.argmin(np.sum(v*v,axis=1));dist=np.linalg.norm(v[j])
            if abs(dist-r)<.025 and n@v[j]/max(dist,1e-20)<-.1:return cid
        for ch,po in allports:
            origin=np.array(po['attachment_mm']);u=np.array(po['axis']);delta=q-origin;s=delta@u;v=delta-s*u;dist=np.linalg.norm(v);r=ch['inner_diameter_mm']/2;length=po['length_mm'];neck=po['straight_neck_length_mm'];mouth=po['inner_mouth_diameter_mm']/2
            if not (-.5<=s<=length+.8):continue
            expected=r if s<=neck else min(mouth,r+(mouth-r)*(s-neck)/(length-neck))
            if abs(dist-expected)<.025 and n@v/max(dist,1e-20)<-.1:return ch['id']
        return None
    def planned_escape(c,d):
        for ch,po in allports:
            u=np.array(po['axis']);m=np.array(po['mouth_center_mm']);den=u@d
            if abs(den)<1e-10:continue
            t=(m-c)@u/den
            if t>0 and np.linalg.norm(c+t*d-m)<=po['inner_mouth_diameter_mm']/2+.015:return ch['id']+'_'+po['end']
        return None
    def mouth_cap(q,n):
        for ch,po in allports:
            u=np.array(po['axis']);v=q-np.array(po['mouth_center_mm']);axial=v@u
            if abs(axial)<.015 and abs(n@u)>.94 and np.linalg.norm(v-axial*u)<=po['outer_mouth_diameter_mm']/2+.02:return True
        return False
    result={'scope':'Actual saved coupon STL in its own XYZ millimetres; railway transform not applied','mesh_file':mesh_path.name,'mesh_sha256':initial,'specification_sha256':sha(spec_path),'specification_recorded_stl_sha256':spec['geometry_qa']['stl_sha256'],'matches_specification_recorded_stl_hash':initial==spec['geometry_qa']['stl_sha256'],'source_script_sha256':sha(__file__),'peer_validator_sha256':sha(Path(__file__).with_name('verify_saved_wall_regions.py')),'physical_tests_performed':False,'global_wall_or_leak_proof':False,'saved_mesh_topology':{'watertight':bool(mesh.is_watertight),'winding_consistent':bool(mesh.is_winding_consistent),'body_count':int(mesh.body_count),'triangles':len(mesh.faces),'bounds_mm':mesh.bounds.tolist()},'path_step_max_mm':.4,'longitudinal_offsets':'center plus four offsets at 0.75 times the channel inner radius','radial_angles':24,'nominal_channel_count':5,'nominal_ports':10,'channels':[]}
    for ch in channels:
        cid=ch['id'];ri=ch['inner_diameter_mm']/2;p=densify(ch['centerline_mm']);row={'id':cid,'inner_diameter_mm':2*ri,'longitudinal_rays':0,'blocked_longitudinal_rays':0,'radial_rays':0,'planned_mouth_escape_rays':0,'unexpected_missing_radial_boundaries':0,'first_surfaces_not_matching_channel':0,'normal_outer_wall_samples':0,'normal_interchannel_wall_samples':0,'normal_internal_same_channel_samples':0,'intentional_mouth_end_face_exits':0,'missing_normal_exits':0,'minimum_outer_normal_wall_mm':None,'minimum_interchannel_normal_wall_mm':None,'outer_wall_samples_below_0_5_mm':0,'interchannel_wall_samples_below_0_6_mm':0,'examples':[],'ports':[]}
        outer_t=[];other_t=[]
        for aa,bb in zip(p[:-1],p[1:]):
            u=(bb-aa)/np.linalg.norm(bb-aa);n,q=basis(u)
            for off in [np.zeros(3),.75*ri*n,-.75*ri*n,.75*ri*q,-.75*ri*q]:
                row['longitudinal_rays']+=1;h=first_hit(loc,aa+off,bb+off)
                if h is not None:
                    row['blocked_longitudinal_rays']+=1
                    if len(row['examples'])<8:row['examples'].append({'issue':'blocked longitudinal sample','start':(aa+off).tolist(),'end':(bb+off).tolist(),'hit':h[0].tolist(),'face_id':h[1]})
            c=(aa+bb)/2
            for theta in np.arange(24)*2*np.pi/24:
                d=np.cos(theta)*n+np.sin(theta)*q;row['radial_rays']+=1;h=first_hit(loc,c,c+maxdist*d)
                if h is None:
                    if planned_escape(c,d):row['planned_mouth_escape_rays']+=1
                    else:row['unexpected_missing_radial_boundaries']+=1
                    continue
                hit,face=h;normal=mesh.face_normals[face]
                if fluid_face(hit,normal)!=cid:
                    row['first_surfaces_not_matching_channel']+=1
                    if len(row['examples'])<8:row['examples'].append({'issue':'unexpected first radial surface','center':c.tolist(),'hit':hit.tolist(),'face_id':face})
                    continue
                h2=first_hit(loc,hit-1e-4*normal,hit-maxdist*normal)
                if h2 is None:row['missing_normal_exits']+=1;continue
                end,face2=h2;nn=mesh.face_normals[face2];thickness=float(np.linalg.norm(end-hit));exit_kind=fluid_face(end,nn)
                if mouth_cap(end,nn):row['intentional_mouth_end_face_exits']+=1;continue
                if exit_kind==cid:row['normal_internal_same_channel_samples']+=1
                elif exit_kind:
                    row['normal_interchannel_wall_samples']+=1;other_t.append(thickness)
                    if thickness<.6:row['interchannel_wall_samples_below_0_6_mm']+=1
                else:
                    row['normal_outer_wall_samples']+=1;outer_t.append(thickness)
                    if thickness<.5:row['outer_wall_samples_below_0_5_mm']+=1
        for po in ch['ports']:
            attach=np.array(po['attachment_mm']);mouth=np.array(po['mouth_center_mm']);u=np.array(po['axis']);n,q=basis(u);pr={'end':po['end'],'through_nozzle_rays':0,'blocked_through_nozzle_rays':0,'mouth_disk_rays':0,'blocked_mouth_disk_rays':0,'examples':[]}
            for r,kind in [(.75*ri,'through_nozzle'),(.90*po['inner_mouth_diameter_mm']/2,'mouth_disk')]:
                offsets=[np.zeros(3)]+[r*(np.cos(t)*n+np.sin(t)*q) for t in np.arange(8)*2*np.pi/8]
                for off in offsets:
                    aa,bb=(attach-.1*u+off,mouth+.3*u+off) if kind=='through_nozzle' else (mouth-.05*u+off,mouth+.05*u+off)
                    pr[kind+'_rays']+=1;h=first_hit(loc,aa,bb)
                    if h is not None:
                        pr['blocked_'+kind+'_rays']+=1
                        pr['examples'].append({'start':aa.tolist(),'end':bb.tolist(),'hit':h[0].tolist(),'face_id':h[1]})
            pr['sampled_open']=pr['blocked_through_nozzle_rays']==0 and pr['blocked_mouth_disk_rays']==0;row['ports'].append(pr)
        row['minimum_outer_normal_wall_mm']=min(outer_t) if outer_t else None;row['minimum_interchannel_normal_wall_mm']=min(other_t) if other_t else None;row['finite_checks_clear']=not any(row[k] for k in ['blocked_longitudinal_rays','unexpected_missing_radial_boundaries','first_surfaces_not_matching_channel','missing_normal_exits','outer_wall_samples_below_0_5_mm','interchannel_wall_samples_below_0_6_mm']) and all(p['sampled_open'] for p in row['ports']);result['channels'].append(row)
        print(cid,'path rays',row['longitudinal_rays'],'blocked',row['blocked_longitudinal_rays'],'normalwall',row['minimum_outer_normal_wall_mm'],'unexpected',row['first_surfaces_not_matching_channel'],'missing',row['unexpected_missing_radial_boundaries'],'ports',all(p['sampled_open'] for p in row['ports']),flush=True)
    crossing=spec['crossing'];x,y=crossing['crossing_xy_mm'];low=next(c for c in channels if c['id']==crossing['channel_lower']);high=next(c for c in channels if c['id']==crossing['channel_upper']);za=low['centerline_mm'][0][2];zb=high['centerline_mm'][0][2];barriers=[]
    for dx in [-.25,0.,.25]:
        for dy in [-.25,0.,.25]:
            aa=np.array([x+dx,y+dy,za]);bb=np.array([x+dx,y+dy,zb]);h=first_hit(loc,aa,bb);entry=None;exit_=None
            if h:
                entry,fi=h;h2=first_hit(loc,entry+np.array([0,0,1e-4]),bb)
                if h2:exit_,fj=h2
            barriers.append({'xy_mm':[x+dx,y+dy],'material_barrier_present':entry is not None and exit_ is not None,'measured_vertical_material_gap_mm':float(exit_[2]-entry[2]) if entry is not None and exit_ is not None else None,'entry_mm':entry.tolist() if entry is not None else None,'exit_mm':exit_.tolist() if exit_ is not None else None,'entry_face':int(fi) if entry is not None else None,'exit_face':int(fj) if exit_ is not None else None})
    result['crossing_material_barrier_samples']=barriers;result['crossing_minimum_sampled_gap_mm']=min([b['measured_vertical_material_gap_mm'] for b in barriers if b['measured_vertical_material_gap_mm'] is not None],default=None)
    result['mesh_bytes_unchanged']=sha(mesh_path)==initial;result['summary']={'channels_tested':len(result['channels']),'ports_tested':sum(len(c['ports']) for c in result['channels']),'longitudinal_rays':sum(c['longitudinal_rays'] for c in result['channels']),'blocked_longitudinal_rays':sum(c['blocked_longitudinal_rays'] for c in result['channels']),'radial_directions':sum(c['radial_rays'] for c in result['channels']),'minimum_outer_normal_wall_sample_mm':min(c['minimum_outer_normal_wall_mm'] for c in result['channels'] if c['minimum_outer_normal_wall_mm'] is not None),'mouth_disk_rays':sum(p['mouth_disk_rays'] for c in result['channels'] for p in c['ports']),'through_nozzle_rays':sum(p['through_nozzle_rays'] for c in result['channels'] for p in c['ports']),'all_10_sampled_open':all(p['sampled_open'] for c in result['channels'] for p in c['ports']),'all_channel_finite_checks_clear':all(c['finite_checks_clear'] for c in result['channels']),'all_crossing_barriers_meet_0_6_mm':all(b['material_barrier_present'] and b['measured_vertical_material_gap_mm']>=.6 for b in barriers),'physical_test_performed':False}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['summary']),flush=True)
if __name__=='__main__':main()
