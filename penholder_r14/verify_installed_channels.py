"""Read-only channel/port checks on one exact saved R14 STL.

This validator does not import the generator, repair triangles, or fill holes.
Rays include the original polyline bends and a 1.1 mm diameter passage probe.
Finite rays do not prove complete wall thickness or physical fluid performance.
"""
from pathlib import Path
import argparse,json,hashlib,time
import numpy as np,trimesh
from mesh_ray_checks_r13 import locator_for,first_hit,basis,probe_rays,ray_winding


def validate(root,path,enu=False):
    start=time.time();m=trimesh.load_mesh(path,process=True)
    if not enu:
        T=np.asarray(json.loads((root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
        m.apply_transform(np.linalg.inv(T))
    loc=locator_for(m);data=json.loads((root/'data/rail_paths_R14.json').read_text())
    ports=json.loads((root/'data/ports_R14.json').read_text())['ports'];lines=[];openings=[]
    for l in data['lines']:
        count=0;bad=[]
        for tr in l['paths']:
            p=np.asarray(tr['points_enu_mm'])
            for k,(a,b) in enumerate(zip(p[:-1],p[1:])):
                d=b-a
                if np.linalg.norm(d)<1e-10:continue
                n,q=basis(d/np.linalg.norm(d))
                offsets=[np.zeros(3)]+[.55*(np.cos(x)*n+np.sin(x)*q) for x in np.arange(8)*np.pi/4]
                for off in offsets:
                    hit=first_hit(loc,a+off,b+off);count+=1
                    if hit['hit']:bad.append(dict(trace_id=tr['trace_id'],segment=k,offset_mm=off.tolist(),**hit))
        lines.append(dict(line_id=l['line_id'],name_ja=l['name_ja'],rays=count,
                          offset_radius_mm=.55,obstruction_count=len(bad),obstructions=bad,passed=len(bad)==0))
        print('PASSAGE',l['line_id'],count,len(bad),flush=True)
    for p in ports:
        a=np.asarray(p['attachment_enu_mm']);mouth=np.asarray(p['mouth_center_enu_mm']);axis=np.asarray(p['axis_enu'])
        internal=probe_rays(loc,mouth,axis,.55,True,a,p['length_mm'],m)
        approach=probe_rays(loc,mouth,axis,.90,False,body=m)
        n,q=basis(axis)
        empty=[ray_winding(loc,m,mouth+.1*axis+.9*(np.cos(x)*n+np.sin(x)*q))
               for x in np.arange(8)*np.pi/4]
        passed=not any(x['hit'] for x in internal+approach) and all(x['empty_by_oriented_ray'] for x in empty)
        openings.append(dict(port_id=p['port_id'],line_id=p['line_id'],internal_rays=internal,
                             approach_rays=approach,aperture_empty_probes=empty,passed=passed))
    return dict(schema='R14-saved-installed-channel-check-v1',file=path.name,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),input_modified=False,
                generator_imported=False,lines=lines,ports=openings,
                route_rays=sum(x['rays'] for x in lines),port_count=len(openings),
                passed=all(x['passed'] for x in lines+openings),seconds=time.time()-start,
                scope='Saved whole-body triangle surface; original polyline segment and 8 offset rays; nozzle and outside approach rays; oriented aperture emptiness probes.',
                limitations=['Finite rays do not establish continuous global wall thickness.','Only the encoded straight outside approach for each nozzle is tested.','SLA printing, washing, drainage, pigment flow, strength and leakage are untested.'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--model',type=Path);p.add_argument('--enu',action='store_true');p.add_argument('--out',type=Path)
    a=p.parse_args();model=a.model or a.root/'models/subway_R14_hollow_print.stl'
    out=a.out or a.root/'reports/installed_channels_R14.json';out.parent.mkdir(parents=True,exist_ok=True)
    result=validate(a.root,model,a.enu);out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('RESULT',result['passed'],'ports',sum(x['passed'] for x in result['ports']),'/69','seconds',result['seconds'],flush=True)


if __name__=='__main__':main()
