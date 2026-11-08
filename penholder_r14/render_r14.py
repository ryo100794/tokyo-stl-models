#!/usr/bin/env python3
"""CPU previews from the actual serialized R14 meshes, in design ENU mm."""
import argparse,json
from pathlib import Path
import numpy as np
import trimesh
from PIL import Image,ImageDraw,ImageFont
from cpu_mesh_render import render_cpu

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args();root=a.root
    matrix=np.asarray(json.loads((root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    inverse=np.linalg.inv(matrix)
    data=json.loads((root/'data/rail_paths_R14.json').read_text());parts=[]
    for line in data['lines']:
        mesh=trimesh.load_mesh(root/'models'/f"rail_R14_{line['line_id']}.stl",process=True);mesh.apply_transform(inverse)
        color=line['source_material']['pbrMetallicRoughness']['baseColorFactor']
        parts.append(dict(vertices=np.asarray(mesh.vertices),faces=np.asarray(mesh.faces),color=color))
    camera=dict(center=[0,0,7],eye=[175,-240,205],up=[0,0,1],scale=130)
    output=root/'previews';output.mkdir(exist_ok=True)
    im=render_cpu(parts,camera,1500,1100,[246,248,250],supersample=2)
    draw=ImageDraw.Draw(im);font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',24)
    draw.text((35,25),'R14 | 14 rebuilt routes | nominal bore 1.2 mm | physical tests pending',fill=(32,43,55),font=font)
    im.save(output/'R14_railways.png')
    body=root/'models/subway_R14_hollow_print.stl'
    if body.exists():
        m=trimesh.load_mesh(body,process=True);m.apply_transform(inverse)
        im=render_cpu([dict(vertices=m.vertices,faces=m.faces,color=[.78,.80,.82,1])],camera,1500,1100,[246,248,250],supersample=2)
        ImageDraw.Draw(im).text((35,25),'R14 | hollow pen holder | physical tests pending',fill=(32,43,55),font=font);im.save(output/'R14_overview.png')
    m=trimesh.load_mesh(root/'models/flow_coupon_R14.stl',process=True)
    im=render_cpu([dict(vertices=m.vertices,faces=m.faces,color=[.53,.73,.81,1])],dict(center=[27,34,3],eye=[112,-65,125],up=[0,0,1],scale=48),1200,950,[246,248,250],supersample=2)
    ImageDraw.Draw(im).text((30,24),'R14 | flow coupon | bores 1.0 / 1.2 / 1.4 mm',fill=(32,43,55),font=font);im.save(output/'R14_flow_coupon.png')

if __name__=='__main__':main()
