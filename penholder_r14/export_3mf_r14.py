"""Export millimetre 3MF from the delivered STL and verify the loaded surface."""
import argparse,hashlib,json,zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np,trimesh

def canonical(mesh):
    v,remap=np.unique(np.asarray(mesh.vertices),axis=0,return_inverse=True)
    f=remap[np.asarray(mesh.faces)]
    # Only cyclic rotations are allowed; reversing a face remains different.
    start=np.argmin(f,axis=1);f=np.take_along_axis(f,(np.arange(3)[None,:]+start[:,None])%3,axis=1)
    f=f[np.lexsort(f[:,::-1].T)]
    return v,f

def export(source,target,report):
    mesh=trimesh.load_mesh(source,process=True);mesh.units='mm'
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_bytes(trimesh.exchange.threemf.export_3MF(mesh,compresslevel=6))
    restored=trimesh.load(target,force='mesh',process=False)
    va,fa=canonical(mesh);vb,fb=canonical(restored)
    same=bool(np.array_equal(va,vb) and np.array_equal(fa,fb))
    with zipfile.ZipFile(target) as z:
        models=[n for n in z.namelist() if n.endswith('.model')]
        units=[ET.fromstring(z.read(n)).attrib.get('unit','millimeter') for n in models]
    record=dict(source_stl=source.name,source_stl_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                file=target.name,sha256=hashlib.sha256(target.read_bytes()).hexdigest(),bytes=target.stat().st_size,
                units=units,triangle_count=len(restored.faces),vertex_coordinates_exactly_preserved=bool(np.array_equal(va,vb)),
                oriented_triangle_set_exactly_preserved=same,source_transform_already_baked=True,
                file_conversion_only=True,physical_tests_performed=False)
    report.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)
    if not same or units!=['millimeter']:raise ValueError('3MF surface or unit equivalence failed')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('target',type=Path);ap.add_argument('--report',type=Path,required=True);a=ap.parse_args()
    export(a.source,a.target,a.report)
