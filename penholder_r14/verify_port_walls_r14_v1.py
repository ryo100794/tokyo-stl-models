#!/usr/bin/env python3
"""Finite wall checks at R14 nozzle attachments, necks and tapers.

Read serialized STLs. Use intended cutter surfaces only to label fluid-facing
hits; measure distances on saved triangles. This is not a continuous global
wall proof or a physical-print qualification. Planned mouth annuli are
separately identified and excluded from sidewall minimum statistics.
"""
import argparse,hashlib,json,math,time
from pathlib import Path
import numpy as np,trimesh,vtk
from verify_saved_wall_regions import locator_for,first_hit,escape_via_port

SECTIONS=[-.25,-.125,0.,.125,.25,.5,1.,1.5,1.9,1.99,2.,2.01,2.1,2.4,2.8,3.2,3.4,3.55,3.59]
ANGLES=48
CORE_SURFACE_LABEL_TOL=.008


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def from_npz(path):
    with np.load(path) as z:return trimesh.Trimesh(z['vertices'],z['faces'],process=False)

def nearest(loc,q):
    p=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);dd=vtk.reference(0.);loc.FindClosestPoint(q,p,cell,sub,dd)
    return float(dd)**.5,np.array(p),int(cell)

def mouth_face(q,normal,ports):
    for p in ports:
        u=np.array(p['axis_enu']);m=np.array(p['mouth_center_enu_mm']);delta=q-m
        if abs(delta@u)<.015 and abs(normal@u)>.94 and np.linalg.norm(delta-(delta@u)*u)<=p['mouth_outer_diameter_mm']/2+.03:return p['port_id']
    return None

def category(s):return 'attachment' if s<=.25 else 'neck' if s<1.99 else 'neck_taper_transition' if s<=2.1 else 'taper'


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--outer-cache',type=Path,required=True);ap.add_argument('--core-cache',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--whole',type=Path,help='optional whole-body print-coordinate STL, checked separately from parts');a=ap.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    ports=json.loads((a.root/'data/ports_R14.json').read_text())['ports'];T=np.array(json.loads((a.root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major']);invT=np.linalg.inv(T)
    c=math.cos(math.pi/24);neck=.6*c;taper=.6*c/math.sqrt(1+(.25*c)**2)
    result={'method':'19 axial sections per nozzle, 48 angular directions (both polygon vertices and facet mid-directions), then face-normal material probes from actual saved inner-surface intersections. Faces matching intended core within 0.008 mm are labeled fluid-facing; all measured distances come from saved STL triangles.','source_script_sha256':sha(Path(__file__)),'peer_validator_sha256':sha(Path(__file__).with_name('verify_saved_wall_regions.py')),'source_paths_sha256':sha(a.root/'data/rail_paths_R14.json'),'source_ports_sha256':sha(a.root/'data/ports_R14.json'),'source_coordinates_sha256':sha(a.root/'data/coordinate_systems_R14.json'),'sections_from_attachment_mm':SECTIONS,'angles':ANGLES,'threshold_mm':.5,'ray_restart_offset_mm':.0001,'core_surface_label_tolerance_mm':CORE_SURFACE_LABEL_TOL,'global_wall_proven':False,'physical_tests_performed':False,'analytic_sidewall_bounds':{'neck_parallel_facet_gap_mm':neck,'taper_parallel_facet_gap_mm':taper,'formula_neck':'(1.2-0.6) cos(pi/24)','formula_taper':'(1.2-0.6) cos(pi/24) / sqrt(1 + (0.25 cos(pi/24))^2)','scope':'Unmodified source 24-sided neck and taper side facets. Does not bound the intentional mouth edge, arbitrary material chords, other cutters or export perturbation. Attachment overlap is bounded separately from the actual outer-positive surface distance at the attachment.'},'parts':[],'whole':[]}
    cores={};core_locs={};outers={};outer_locs={}
    for lid in sorted(set(p['line_id'] for p in ports)):
        cores[lid]=from_npz(a.core_cache/f'{lid}_core.npz');core_locs[lid]=locator_for(cores[lid]);outers[lid]=from_npz(a.outer_cache/f'{lid}_outer.npz');outer_locs[lid]=locator_for(outers[lid])
    whole_mesh=None;whole_locator=None
    if a.whole:
        whole_mesh=trimesh.load_mesh(a.whole,process=True);whole_mesh.apply_transform(invT);whole_locator=locator_for(whole_mesh);result['whole_stl_sha256']=sha(a.whole);result['whole_stl_file']=a.whole.name
    for lid in sorted(cores):
        line_ports=[p for p in ports if p['line_id']==lid];path=a.root/'models'/f'rail_R14_{lid}.stl';mesh=trimesh.load_mesh(path,process=True);mesh.apply_transform(invT);loc=locator_for(mesh)
        for variant,mm,locator in [('parts',mesh,loc)]+([('whole',whole_mesh,whole_locator)] if a.whole else []):
            maxdist=float(np.linalg.norm(mm.extents))+10
            for p in line_ports:
                start=time.time();attach=np.array(p['attachment_enu_mm']);u=np.array(p['axis_enu']);u/=np.linalg.norm(u);helper=np.eye(3)[np.argmin(abs(u))];x=np.cross(helper,u);x/=np.linalg.norm(x);y=np.cross(u,x)
                anchor_radius=nearest(outer_locs[lid],attach)[0];inside=bool(outers[lid].contains([attach])[0]);overlap=float(p.get('attachment_overlap_mm',.25));back_bound=anchor_radius-math.sqrt(overlap**2+.6**2)
                row={'line_id':lid,'port_id':p['port_id'],'part_sha256':sha(path),'outer_positive_sha256':sha(a.outer_cache/f'{lid}_outer.npz'),'intended_core_sha256':sha(a.core_cache/f'{lid}_core.npz'),'rays':0,'normal_outer_wall_samples':0,'normal_interline_wall_samples':0,'same_line_internal_exit_samples':0,'planned_mouth_face_exits':0,'planned_mouth_escape_rays':0,'missing_boundary_rays':0,'first_hits_not_matching_own_core':0,'wrong_entry_orientation':0,'missing_normal_exit':0,'outer_wall_below_0_5_mm':0,'interline_wall_below_0_5_mm':0,'minimum_outer_normal_wall_mm':None,'minimum_interline_normal_wall_mm':None,'by_section_region':{k:[] for k in ['attachment','neck','neck_taper_transition','taper']},'outer_positive_anchor_ball':{'attachment_inside_positive':inside,'positive_watertight':bool(outers[lid].is_watertight),'positive_winding_consistent':bool(outers[lid].is_winding_consistent),'distance_attachment_to_positive_surface_mm':anchor_radius,'maximum_extent_of_backward_core_from_attachment_mm':math.sqrt(overlap**2+.6**2),'backward_overlap_source_containment_lower_bound_mm':back_bound if inside and outers[lid].is_watertight else None,'scope':'Source-positive containing ball minus maximum overlap-core radius, valid for s in [-attachment_overlap,0]. Final saved sidewall samples and mouth exclusions are reported separately.'},'worst_outer_sample':None,'examples':{'thin':[],'unexpected':[],'same_line_internal_exit':[]}}
                outer_values=[];other_values=[]
                for s in SECTIONS:
                    if s>=p['length_mm']:continue
                    center=attach+s*u
                    for ai,theta in enumerate(np.arange(ANGLES)*2*np.pi/ANGLES):
                        direction=np.cos(theta)*x+np.sin(theta)*y;row['rays']+=1;meta={'section_mm':s,'angle_index':ai,'angle_rad':float(theta),'center_enu_mm':center.tolist(),'radial_direction_enu':direction.tolist()}
                        first=first_hit(locator,center,center+maxdist*direction)
                        if first is None:
                            mouth=escape_via_port(center,direction,line_ports,maxdist)
                            if mouth:row['planned_mouth_escape_rays']+=1
                            else:
                                row['missing_boundary_rays']+=1
                                if len(row['examples']['unexpected'])<8:row['examples']['unexpected'].append({**meta,'issue':'no material boundary beyond model bounds'})
                            continue
                        q,face=first;n=mm.face_normals[face];core_distance=nearest(core_locs[lid],q)[0]
                        if core_distance>CORE_SURFACE_LABEL_TOL:
                            row['first_hits_not_matching_own_core']+=1
                            if len(row['examples']['unexpected'])<8:row['examples']['unexpected'].append({**meta,'issue':'first hit does not match expected same-line core','core_distance_mm':core_distance,'point_mm':q.tolist(),'face_id':face})
                            continue
                        if n@direction>1e-4:
                            row['wrong_entry_orientation']+=1;continue
                        second=first_hit(locator,q-.0001*n,q-maxdist*n)
                        if second is None:row['missing_normal_exit']+=1;continue
                        end,exitface=second;thickness=float(np.linalg.norm(end-q));en=mm.face_normals[exitface];mouth=mouth_face(end,en,ports if variant=='whole' else line_ports)
                        detail={**meta,'normal_thickness_mm':thickness,'entry_mm':q.tolist(),'exit_mm':end.tolist(),'entry_face':face,'exit_face':exitface,'entry_face_vertices_mm':mm.triangles[face].tolist(),'exit_face_vertices_mm':mm.triangles[exitface].tolist(),'entry_normal_enu':n.tolist(),'exit_normal_enu':en.tolist()}
                        if mouth:row['planned_mouth_face_exits']+=1;continue
                        own_distance=nearest(core_locs[lid],end)[0]
                        if own_distance<=CORE_SURFACE_LABEL_TOL:
                            row['same_line_internal_exit_samples']+=1
                            if thickness<.5 and len(row['examples']['same_line_internal_exit'])<8:row['examples']['same_line_internal_exit'].append(detail)
                            continue
                        other_line=None
                        if variant=='whole':
                            for other,cl in core_locs.items():
                                if other!=lid and nearest(cl,end)[0]<=CORE_SURFACE_LABEL_TOL:other_line=other;break
                        if other_line:
                            row['normal_interline_wall_samples']+=1;other_values.append(thickness)
                            if thickness<.5:
                                row['interline_wall_below_0_5_mm']+=1
                                if len(row['examples']['thin'])<8:row['examples']['thin'].append({**detail,'exit_fluid_line':other_line})
                        else:
                            row['normal_outer_wall_samples']+=1;outer_values.append(thickness);row['by_section_region'][category(s)].append(thickness)
                            if row['worst_outer_sample'] is None or thickness<row['worst_outer_sample']['normal_thickness_mm']:row['worst_outer_sample']=detail
                            if thickness<.5:
                                row['outer_wall_below_0_5_mm']+=1
                                if len(row['examples']['thin'])<8:row['examples']['thin'].append(detail)
                row['minimum_outer_normal_wall_mm']=min(outer_values) if outer_values else None;row['minimum_interline_normal_wall_mm']=min(other_values) if other_values else None
                row['by_section_region']={k:{'sample_count':len(v),'minimum_normal_wall_mm':min(v) if v else None} for k,v in row['by_section_region'].items()}
                row['sampled_checks_clear']=not any(row[k] for k in ['missing_boundary_rays','first_hits_not_matching_own_core','wrong_entry_orientation','missing_normal_exit','outer_wall_below_0_5_mm','interline_wall_below_0_5_mm'])
                row['seconds']=time.time()-start;result[variant].append(row)
                print(variant,p['port_id'],'min',row['minimum_outer_normal_wall_mm'],'interline',row['minimum_interline_normal_wall_mm'],'thin',row['outer_wall_below_0_5_mm'],'unexpected',row['first_hits_not_matching_own_core'],'missing',row['missing_boundary_rays'],'overlap_bound',back_bound,flush=True)
            (a.out/'port_wall_samples_R14.json').write_text(json.dumps(result,indent=2)+'\n')
    for variant in ['parts','whole']:
        rows=result[variant]
        if not rows:continue
        result[f'{variant}_summary']={'ports':len(rows),'rays':sum(x['rays'] for x in rows),'normal_outer_wall_samples':sum(x['normal_outer_wall_samples'] for x in rows),'normal_interline_wall_samples':sum(x['normal_interline_wall_samples'] for x in rows),'minimum_outer_normal_wall_mm':min(x['minimum_outer_normal_wall_mm'] for x in rows if x['minimum_outer_normal_wall_mm'] is not None),'minimum_interline_normal_wall_mm':min([x['minimum_interline_normal_wall_mm'] for x in rows if x['minimum_interline_normal_wall_mm'] is not None],default=None),'minimum_backward_overlap_source_bound_mm':min(x['outer_positive_anchor_ball']['backward_overlap_source_containment_lower_bound_mm'] for x in rows if x['outer_positive_anchor_ball']['backward_overlap_source_containment_lower_bound_mm'] is not None),'all_ports_sampled_checks_clear':all(x['sampled_checks_clear'] for x in rows),'total_outer_wall_samples_below_0_5_mm':sum(x['outer_wall_below_0_5_mm'] for x in rows),'total_interline_wall_samples_below_0_5_mm':sum(x['interline_wall_below_0_5_mm'] for x in rows),'total_missing_boundaries':sum(x['missing_boundary_rays'] for x in rows),'total_unexpected_first_surfaces':sum(x['first_hits_not_matching_own_core'] for x in rows),'total_same_line_internal_exits':sum(x['same_line_internal_exit_samples'] for x in rows),'planned_mouth_face_exits':sum(x['planned_mouth_face_exits'] for x in rows)}
    (a.out/'port_wall_samples_R14.json').write_text(json.dumps(result,indent=2)+'\n');print('SUMMARIES',json.dumps({k:v for k,v in result.items() if k.endswith('_summary')}),flush=True)
if __name__=='__main__':main()
