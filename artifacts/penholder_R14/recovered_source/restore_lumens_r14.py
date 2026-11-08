#!/usr/bin/env python3
"""Re-cut preserved R13 centre lines through the reconstructed closed body.

The reconstruction is never accepted as evidence of channel clearance. A
1.24 mm nominal cutter is explicitly subtracted afterwards; this allows a
0.02 mm radial numerical clearance while retaining the 2.4 mm outer envelope.
"""
import argparse,json,gc
from pathlib import Path
import numpy as np
import trimesh
from geometry_r14 import *
from build_r14 import export_stl,dump

def path_cutter(points,radius):
    p=rdp(points,.01)
    sphere=mf.Manifold.sphere(radius+.01,24)
    return unite([cylinder_between(a,b,radius,sides=24) for a,b in zip(p[:-1],p[1:])]+[sphere.translate(q) for q in p])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--cache',type=Path,required=True);a=ap.parse_args();root=a.root.resolve();cache=a.cache.resolve();cache.mkdir(exist_ok=True)
    data=json.loads((root/'data/rail_paths_R14.json').read_text());ports=json.loads((root/'data/ports_R14.json').read_text())['ports']
    T=np.asarray(json.loads((root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    mesh=trimesh.load_mesh(root/'models/subway_R14_hollow_print.stl',process=True);mesh.apply_transform(np.linalg.inv(T));body=as_manifold(mesh.vertices,mesh.faces);del mesh;gc.collect()
    for line in data['lines']:
        lid=line['line_id'];path=cache/f'{lid}_core.npz'
        if path.exists():core=load_npz(path)
        else:
            items=[path_cutter(t['points_enu_mm'],.62) for t in line['paths']]
            for p in ports:
                if p['line_id']!=lid:continue
                length=p['length_mm'];items.append(ring_profile_solid(p['attachment_enu_mm'],p['axis_enu'],[(-.25,.62),(2.,.62),(length,1.02),(length+.8,1.02)]))
            core=unite(items);save_npz(core,path);del items
        print('cut preserved path',lid,flush=True);body=body-core;body=body.simplify(.001);del core;gc.collect()
        save_npz(body,cache/'body_recut.npz')
    print('EXPORT RECUT BODY',body.num_tri(),flush=True)
    info=export_stl(body,root/'models/subway_R14_hollow_print.stl',T,tolerances=(.01,.015,.02))
    info.update(nominal_final_rail_bore_mm=1.24,nominal_straight_wall_mm=.58,nominal_mouth_bore_mm=2.04,
                reconstructed_body_followed_by_explicit_lumen_subtraction=True)
    dump(root/'reports/body_recut_R14.json',info)
    mesh=trimesh.load_mesh(root/'models/subway_R14_hollow_print.stl',process=True);mesh.units='mm'
    (root/'models/subway_R14_hollow_print.3mf').write_bytes(trimesh.exchange.threemf.export_3MF(mesh,compresslevel=6))

if __name__=='__main__':main()
