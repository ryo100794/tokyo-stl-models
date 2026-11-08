#!/usr/bin/env python3
"""Orientation and gravity drainage screening from saved R14 data, not print qualification."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import trimesh

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 root=Path(a.root);out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 co=json.loads((root/'data/coordinate_systems_R14.json').read_text());po=json.loads((root/'data/ports_R14.json').read_text());ro=json.loads((root/'data/rail_paths_R14.json').read_text())
 T=np.array(co['geographic_ENU_mm_to_print_XYZ_mm_row_major']);rot=T[:3,:3];shift=T[:3,3]
 ports=[]
 for q in po['ports']:
  axis=np.array(q['axis_enu'])@rot.T;axis/=np.linalg.norm(axis)
  mouth=np.array(q['mouth_center_enu_mm'])@rot.T+shift;attach=np.array(q['attachment_enu_mm'])@rot.T+shift
  inner_d=q['mouth_inner_diameter_mm'];
  bottom=mouth[2]-inner_d/2*np.sqrt(max(0,1-axis[2]**2))
  ports.append({'port_id':q['port_id'],'line_id':q['line_id'],'trace_id':q['trace_id'],'kind':q['kind'],'print_mouth_center_mm':mouth.tolist(),'print_attachment_mm':attach.tolist(),'print_axis_unit':axis.tolist(),'axis_angle_from_positive_build_Z_degrees':float(np.degrees(np.arccos(np.clip(axis[2],-1.,1.)))),'mouth_minus_attachment_z_mm':float(mouth[2]-attach[2]),'nominal_mouth_lowest_edge_z_mm':float(bottom),'nominal_mouth_inner_diameter_mm':inner_d})
 drains=[]
 for ln in ro['lines']:
  for tr in ln['paths']:
   pts=np.array(tr['points_enu_mm'])@rot.T+shift
   ps=[q for q in ports if q['trace_id']==tr['trace_id'] and q['line_id']==ln['line_id']]
   radius=float(co['tube_nominal_lumen_diameter_mm'])/2
   # Plane perpendicular to local centerline: vertical radial extent is r*sqrt(1-t_z^2)
   tangent=np.gradient(pts,axis=0);tangent/=np.maximum(np.linalg.norm(tangent,axis=1)[:,None],1e-20)
   nominal_bottom=pts[:,2]-radius*np.sqrt(np.maximum(0,1-tangent[:,2]**2))
   lowest_route=float(nominal_bottom.min());ix=int(nominal_bottom.argmin())
   lowest_mouth=min((q['nominal_mouth_lowest_edge_z_mm'] for q in ps),default=None)
   drains.append({'line_id':ln['line_id'],'trace_id':tr['trace_id'],'port_count':len(ps),'nominal_route_lowest_inner_surface_z_mm':lowest_route,'lowest_route_sample_index':ix,'lowest_route_center_xyz_mm':pts[ix].tolist(),'lowest_nominal_mouth_edge_z_mm':lowest_mouth,'nominal_low_point_below_all_same_trace_mouths_mm':max(0,lowest_mouth-lowest_route) if lowest_mouth is not None else None,'screen':'potential retained liquid in this fixed orientation under gravity alone' if lowest_mouth is not None and lowest_mouth>lowest_route+.05 else 'not excluded by the simple elevation screen; full drainage unproven'})
 line_drains=[]
 for ln in ro['lines']:
  ds=[d for d in drains if d['line_id']==ln['line_id']]
  ps=[d for d in ports if d['line_id']==ln['line_id']]
  low=min(d['nominal_route_lowest_inner_surface_z_mm'] for d in ds)
  mouth=min(d['nominal_mouth_lowest_edge_z_mm'] for d in ps)
  line_drains.append({'line_id':ln['line_id'],'lowest_nominal_route_inner_surface_z_mm':low,'lowest_nominal_mouth_edge_z_mm':mouth,'nominal_low_point_below_all_same_line_mouths_mm':max(0,mouth-low),'interpretation':'positive value identifies a fixed-orientation gravity-drainage concern even if same-line traces interconnect; no residual volume or physical-flow claim'})
 coupon=json.loads((root/'data/flow_coupon_R14_specification.json').read_text())
 coupon_ports=[q for ch in coupon['channels'] for q in ch['ports']]
 coupon_screen={'sha256':digest(root/'data/flow_coupon_R14_specification.json'),'own_print_axes':coupon['coordinate_system'],'rail_transform_applied':False,'total_ports':len(coupon_ports),'all_ports_face_positive_z':all(q['axis'][2]>0.99999 for q in coupon_ports),'qualification':'Upright vertical mouths do not demonstrate gravity-only complete drainage; inversion and flushing procedure remain unverified.'}
 orientation=[]
 for f in sorted((root/'models').glob('*.stl')):
  m=trimesh.load_mesh(f,process=False);norm=m.face_normals;area=m.area_faces
  # Report outward normals only. Inner lumen faces are included: this is not a slicer support classification.
  orientation.append({'file':f.name,'sha256':digest(f),'dimensions_xyz_mm':np.ptp(m.vertices,axis=0).tolist(),'min_saved_z_mm':float(m.vertices[:,2].min()),'max_saved_z_mm':float(m.vertices[:,2].max()),'area_with_outward_normal_negative_z_mm2':float(area[norm[:,2]<-1e-8].sum()),'area_with_outward_normal_within_45_degrees_of_negative_z_mm2':float(area[norm[:,2]<-np.sqrt(.5)].sum()),'above_fraction_of_total_surface':float(area[norm[:,2]<-np.sqrt(.5)].sum()/area.sum()),'support_qualification':'NOT EVALUATED: surface-normal screen includes inner lumen and has no layer, resin, support, printer or peel-force model'})
 z=np.array([q['print_axis_unit'][2] for q in ports]);angles=np.array([q['axis_angle_from_positive_build_Z_degrees'] for q in ports])
 report={'schema':'R14-print-orientation-drainage-screen-v2','inputs':{f:digest(root/'data'/f) for f in ['coordinate_systems_R14.json','ports_R14.json','rail_paths_R14.json']},'print_transform':T.tolist(),'orthogonal_rotation_check':bool(np.allclose(rot.T@rot,np.eye(3))),'rotation_determinant':float(np.linalg.det(rot)),'already_applied_to_saved_STL':co['print_transform_baked_into_vertex_coordinates'],'ports_summary':{'total':len(ports),'upward_axis_z_greater_1e_minus_6':int((z>1e-6).sum()),'horizontal_axis_abs_z_at_most_1e_minus_6':int((abs(z)<=1e-6).sum()),'downward_axis_z_less_negative_1e_minus_6':int((z<-1e-6).sum()),'minimum_angle_from_positive_build_Z_degrees':float(angles.min()),'maximum_angle_from_positive_build_Z_degrees':float(angles.max())},'ports':ports,'nominal_trace_drainage_screens':drains,'nominal_line_drainage_screens':line_drains,'coupon_orientation_screen':coupon_screen,'mesh_orientation_screens':orientation,'scope_limitations':['Port direction and nominal elevation use supplied port/route data, not measured physical openings.','Different traces may connect geometrically. The same-trace elevation screen is conservative guidance, not a proof of the residual liquid volume.','Upright gravity screening does not represent inversion, pressure flushing, air evacuation, resin viscosity or capillary effects.','Route elevation does not establish saved-mesh lumen connectivity; use independent mesh and installed-channel evidence tied to the same file SHA256.','No machine profile, material data, slicer supports or experimentally validated minimum wall/channel dimensions were supplied.','Actual SLA printing, washing, curing, leakage and pigment flow remain untested.']}
 (out/'print_orientation_and_drainage_screen_R14.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
 print(json.dumps({'ports':report['ports_summary'],'trace_count':len(drains),'traces_with_elevation_trap':sum((q['nominal_low_point_below_all_same_trace_mouths_mm'] or 0)>.05 for q in drains),'max_potential_head_mm':max((q['nominal_low_point_below_all_same_trace_mouths_mm'] or 0) for q in drains),'by_line':[{ 'line':q['line_id'],'trace':q['trace_id'],'head':q['nominal_low_point_below_all_same_trace_mouths_mm']} for q in drains]},ensure_ascii=False))

if __name__=='__main__':main()
