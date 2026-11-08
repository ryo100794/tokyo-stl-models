#!/usr/bin/env python3
"""Verify channels in every serialized, independent railway STL part."""
import json,hashlib,argparse
from pathlib import Path
import numpy as np
import trimesh
from mesh_ray_checks_r13 import locator_for,first_hit,basis,probe_rays

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);root=ap.parse_args().root
    matrix=np.asarray(json.loads((root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    inv=np.linalg.inv(matrix);data=json.loads((root/'data/rail_paths_R14.json').read_text());ports=json.loads((root/'data/ports_R14.json').read_text())['ports'];rows=[]
    for line in data['lines']:
        path=root/'models'/f"rail_R14_{line['line_id']}.stl";mesh=trimesh.load_mesh(path,process=True);mesh.apply_transform(inv);loc=locator_for(mesh)
        blocked=[];count=0;wall=[];bore=[];ordinary_wall=[];ordinary_bore=[];thin_examples=[];missing_examples=[];missing=0
        line_ports=[p for p in ports if p['line_id']==line['line_id']]
        for trace in line['paths']:
            p=np.asarray(trace['points_enu_mm'])
            for i,(a,b) in enumerate(zip(p[:-1],p[1:])):
                delta=b-a;length=np.linalg.norm(delta)
                if length<1e-10:continue
                n,q=basis(delta/length)
                for off in [np.zeros(3),.45*n,-.45*n,.45*q,-.45*q]:
                    hit=first_hit(loc,a+off,b+off);count+=1
                    if hit['hit']:blocked.append(dict(trace_id=trace['trace_id'],segment=i,**hit))
                for angle in np.arange(8)*np.pi/4:
                    d=np.cos(angle)*n+np.sin(angle)*q;c=(a+b)/2
                    near_port=any(np.linalg.norm(c-np.asarray(p['attachment_enu_mm']))<p['length_mm']+2. for p in line_ports)
                    h=first_hit(loc,c,c+3*d)
                    if not h['hit']:
                        if not near_port:
                            missing+=1
                            if len(missing_examples)<20:missing_examples.append(dict(trace_id=trace['trace_id'],segment=i,point_mm=c.tolist(),direction=d.tolist()))
                        continue
                    r=float(np.linalg.norm(np.asarray(h['point_enu_mm'])-c));bore.append(r)
                    if not near_port:ordinary_bore.append(r)
                    h2=first_hit(loc,c+(r+.0001)*d,c+5*d)
                    if h2['hit']:
                        thickness=float(np.linalg.norm(np.asarray(h2['point_enu_mm'])-c))-r
                        wall.append(thickness)
                        if not near_port and .5<r<.7:
                            ordinary_wall.append(thickness)
                            if thickness<.5:
                                thin_examples.append(dict(radial_interval_mm=thickness,trace_id=trace['trace_id'],segment=i,point_mm=c.tolist(),direction=d.tolist(),inner_triangle=h['triangle_cell_id'],outer_triangle=h2['triangle_cell_id']))
        opening=[]
        for p in ports:
            if p['line_id']!=line['line_id']:continue
            a=np.asarray(p['attachment_enu_mm']);mouth=np.asarray(p['mouth_center_enu_mm']);axis=np.asarray(p['axis_enu'])
            internal=probe_rays(loc,mouth,axis,.45,True,a,p['length_mm'],mesh)
            external=probe_rays(loc,mouth,axis,.8,False,body=mesh)
            opening.append(dict(port_id=p['port_id'],internal_rays_clear=not any(x['hit'] for x in internal),
                approach_rays_clear=not any(x['hit'] for x in external)))
        row=dict(file=path.name,line_id=line['line_id'],name_ja=line['name_ja'],sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            passage_rays=count,parallel_offset_mm=.45,obstruction_hits=blocked,ports=opening,
            minimum_sampled_bore_radius_mm=min(bore),minimum_sampled_radial_wall_mm=min(wall),
            minimum_sampled_radial_wall_away_from_ports_mm=min(ordinary_wall) if ordinary_wall else None,
            ordinary_wall_sample_count=len(ordinary_wall),
            radial_samples_with_no_boundary_within_3mm=missing,
            missing_boundary_examples=missing_examples,
            thinnest_radial_interval_examples=sorted(thin_examples,key=lambda x:x['radial_interval_mm'])[:20],
            provisional_wall_review_threshold_mm=.5,
            wall_review_passed=bool(ordinary_wall and min(ordinary_wall)>=.5 and missing==0),
            wall_scope='Raw radial intervals can approach zero at open lips; ordinary samples exclude distance < nozzle length + 2 mm from any attachment and require a first bore radius between 0.5 and 0.7 mm.',
            passed=not blocked and all(p['internal_rays_clear'] and p['approach_rays_clear'] for p in opening))
        rows.append(row);print(line['line_id'],row['passed'],'blocked',len(blocked),'wall',row['minimum_sampled_radial_wall_mm'],flush=True)
        (root/'reports/individual_passages_R14.json').write_text(json.dumps(dict(
            passed=all(x['passed'] for x in rows),all_14_completed=len(rows)==14,
            scope='Saved individual railway parts; does not certify the installed whole-body lumen.',
            lines=rows,physical_tests_performed=False,wall_sampling_is_not_global_minimum_proof=True),ensure_ascii=False,indent=2)+'\n')

if __name__=='__main__':main()
