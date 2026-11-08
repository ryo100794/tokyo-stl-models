#!/usr/bin/env python3
"""Reproduce the R14 REVIEW-ONLY individual rail meshes from saved inputs.

This does not claim to finish the one-piece pen holder. It keeps all R13
inputs unchanged, writes a fresh output directory, and records failed wall
checks. Floating-point Boolean/Poisson outputs are not byte reproducible
across dependency versions. Use requirements_rebuild_R14.txt.
"""
import argparse,json,shutil,subprocess,sys,gc
from pathlib import Path
import numpy as np
import trimesh
import pymeshlab as pm
import pymeshfix

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    original=Path(__file__).resolve().parents[1];root=a.out.resolve()
    if root.exists() and any(root.iterdir()):raise ValueError('Output must be new or empty')
    root.mkdir(parents=True,exist_ok=True)
    for name in ['data','reference','source']:shutil.copytree(original/name,root/name,dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','cpu_mesh_raster'))
    for name in ['models','reports','previews']:(root/name).mkdir()
    source=root/'source';cache=root/'_work';cache.mkdir()
    def run(script,*args):
        with (root/'reports'/f'{script}.log').open('w') as out:
            subprocess.run([sys.executable,str(source/script),'--root',str(root),*map(str,args)],stdout=out,stderr=subprocess.STDOUT,check=True)
    for stage in ['rails','base','assemble','parts']:run('build_r14.py','--cache',cache,'--stage',stage)
    sys.path.insert(0,str(source));import geometry_r14 as g,rail_components_r14 as rails
    from build_r14 import export_stl
    from restore_lumens_r14 import path_cutter
    data=json.loads((root/'data/rail_paths_R14.json').read_text());ports=json.loads((root/'data/ports_R14.json').read_text())['ports']
    T=np.asarray(json.loads((root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    def overlap_path(points,radius,sides=24,tolerance=.01):
        points=g.rdp(points,tolerance);sphere=g.mf.Manifold.sphere(radius+.01,sides)
        return g.unite([g.cylinder_between(a,b,radius,sides=sides) for a,b in zip(points[:-1],points[1:])]+[sphere.translate(q) for q in points])
    rails.capsule_path=overlap_path
    special=root/'_special';special.mkdir()
    for lid in ['L06','L12']:
        line=next(l for l in data['lines'] if l['line_id']==lid);outer,core=rails.build_line(line,ports,special);solid=outer-core
        target=root/'models'/f'rail_R14_{lid}.stl'
        if lid=='L06':
            export_stl(solid,target,T,tolerances=(.005,.01,.02))
            m=trimesh.load_mesh(target,process=True);fix=pymeshfix.MeshFix(m.vertices,m.faces,verbose=False);fix.repair(joincomp=False,remove_smallest_components=False)
            trimesh.Trimesh(fix.points,fix.faces,process=False).export(target)
        else:
            v,f=g.arrays(solid.simplify(.01));m=trimesh.Trimesh(v,f,process=False);m.update_faces(m.nondegenerate_faces(height=1e-8));m.remove_unreferenced_vertices()
            ms=pm.MeshSet();ms.add_mesh(pm.Mesh(vertex_matrix=m.vertices,v_normals_matrix=m.vertex_normals))
            ms.generate_surface_reconstruction_screened_poisson(depth=12,fulldepth=5,pointweight=8,preclean=True,threads=2)
            ms.meshing_decimation_quadric_edge_collapse(targetfacenum=240000,preservetopology=True,preservenormal=True,qualitythr=.3,autoclean=True)
            out=ms.current_mesh();m=trimesh.Trimesh(out.vertex_matrix(),out.face_matrix(),process=False)
            positive=g.as_manifold(m.vertices,m.faces);cutters=[path_cutter(t['points_enu_mm'],.62) for t in line['paths']]
            for port in ports:
                if port['line_id']!=lid:continue
                length=port['length_mm'];cutters.append(g.ring_profile_solid(port['attachment_enu_mm'],port['axis_enu'],[(-.25,.62),(2.,.62),(length,1.02),(length+.8,1.02)]))
            export_stl(positive-g.unite(cutters),target,T,tolerances=(.005,.01,.015,.02))
        del solid,outer,core;gc.collect()
    coupon=root/'_coupon';subprocess.run([sys.executable,str(source/'build_flow_coupon_R14.py'),'--out',str(coupon)],check=True)
    shutil.copy2(coupon/'flow_coupon_R14.stl',root/'models/flow_coupon_R14.stl')
    shutil.copy2(coupon/'flow_coupon_R14_specification.json',root/'data/flow_coupon_R14_specification.json')
    run('verify_parts_r14.py');run('verify_r14.py','--stage','self');run('render_r14.py')
    review=json.loads((root/'reports/individual_passages_R14.json').read_text())
    print(json.dumps({'generation_finished':True,'whole_penholder_complete':False,
        'all_route_rays_clear':all(x['passed'] for x in review['lines']),
        'all_provisional_wall_checks_passed':all(x['wall_review_passed'] for x in review['lines']),
        'physical_print_and_flow_tests_performed':False},indent=2))

if __name__=='__main__':main()
