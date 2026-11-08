"""Extract a connected saved boundary without changing any retained triangle.

This is an intermediate geometric operation, not acceptance by itself.
Other components must be audited for containment and intended-fluid exclusion;
the extracted surface must pass all independent channel and topology checks.
"""
from pathlib import Path
import argparse,json,hashlib,time
import numpy as np,trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def extract(path,output):
    started=time.time();m=trimesh.load_mesh(path,process=False)
    v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces]
    e=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    graph=coo_matrix((np.ones(len(e),dtype=np.uint8),(e[:,0],e[:,1])),shape=(len(v),len(v))).tocsr()
    count,labels=connected_components(graph,directed=False);face_labels=labels[f[:,0]]
    counts=np.bincount(face_labels);principal=int(counts.argmax())
    q=trimesh.Trimesh(v,f[face_labels==principal],process=False);q.remove_unreferenced_vertices()
    output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes(trimesh.exchange.stl.export_stl(q))
    report=dict(source_file=path.name,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                output_file=output.name,sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                method='Principal exact-coordinate connected boundary; retained triangles and coordinates are unchanged.',
                source_surface_components=int(count),retained_faces=len(q.faces),removed_faces=len(f)-len(q.faces),
                input_modified=False,requires_subsidiary_component_audit=True,
                print_or_flow_qualification=False,seconds=time.time()-started)
    output.with_suffix('.extraction.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();extract(a.input,a.output)
