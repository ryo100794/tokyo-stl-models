"""Render exact opaque triangles in software with an orthographic depth buffer."""
from pathlib import Path
import struct
import subprocess
import tempfile
import numpy as np
from PIL import Image


def render_cpu(parts,camera,width,height,background,supersample=2):
    root=Path(__file__).resolve().parent
    source=root/'cpu_mesh_raster.cpp';binary=root/'cpu_mesh_raster'
    if not binary.exists() or binary.stat().st_mtime<source.stat().st_mtime:
        subprocess.run(['g++','-O3','-std=c++17',str(source),'-o',str(binary)],check=True)
    w,h=width*supersample,height*supersample
    center=np.asarray(camera['center'],float);eye=np.asarray(camera['eye'],float);up=np.asarray(camera['up'],float)
    look=center-eye;look/=np.linalg.norm(look)
    right=np.cross(look,up);right/=np.linalg.norm(right)
    up=np.cross(right,look)
    basis=np.column_stack([right,up,-look])
    px=h/(2*camera['scale'])
    with tempfile.TemporaryDirectory(prefix='raster_',dir=root) as temp:
        input_path=Path(temp)/'scene.bin';output_path=Path(temp)/'pixels.rgb'
        with input_path.open('wb') as stream:
            stream.write(struct.pack('<III3B',w,h,len(parts),*background))
            for part in parts:
                vertices=np.asarray(part['vertices'],float);faces=np.asarray(part['faces'],np.uint32)
                normals=part.get('normals')
                if normals is None or len(normals)!=len(vertices):
                    normals=np.zeros_like(vertices)
                    tri=vertices[faces];cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
                    for corner in range(3):np.add.at(normals,faces[:,corner],cross)
                    normals/=np.maximum(np.linalg.norm(normals,axis=1,keepdims=True),1e-30)
                screen=(vertices-center)@basis
                screen[:,:2]*=px;screen[:,0]+=w/2;screen[:,1]=h/2-screen[:,1]
                packed=np.column_stack([screen,np.asarray(normals)@basis]).astype('<f4')
                color=part.get('color',[.58,.66,.71,1])[:3]
                stream.write(struct.pack('<II3f',len(vertices),len(faces),*color))
                stream.write(packed.tobytes());stream.write(faces.astype('<u4').tobytes())
        subprocess.run([str(binary),str(input_path),str(output_path)],check=True)
        pixels=np.fromfile(output_path,dtype=np.uint8).reshape(h,w,3)
        result=Image.fromarray(pixels)
        if supersample!=1:result=result.resize((width,height),Image.Resampling.LANCZOS)
    return result
