#!/usr/bin/env python3
"""Run exact integer/rational pair predicate on a saved selection of mesh faces."""
from pathlib import Path
import argparse,json,time,subprocess,hashlib
import numpy as np
from rtree import index
from build_exact_engines import ensure_engine
ROOT=Path(__file__).parent
CLASSES={0:'disjoint',1:'expected_shared_vertex_or_edge_contact',2:'other_point_contact',3:'noncoplanar_proper_crossing_segment',4:'coplanar_area_overlap',5:'other_segment_contact',6:'degenerate_triangle'}

def run(tri,pairs,outprefix):
    tri=np.asarray(tri,dtype='<f8');pairs=np.asarray(pairs,dtype='<u4')
    inp=outprefix.with_suffix('.exact_input.bin');rawout=outprefix.with_suffix('.exact_classes.bin')
    with inp.open('wb') as f:
        np.asarray([len(tri),len(pairs)],dtype='<u8').tofile(f);tri.tofile(f);pairs.tofile(f)
    r=subprocess.run([str(ensure_engine('exact_triangle_pairs')),str(inp),str(rawout)],capture_output=True,text=True,check=True)
    codes=np.fromfile(rawout,dtype=np.uint8)
    if len(codes)!=len(pairs):raise RuntimeError('Incomplete exact output')
    return codes,r.stdout.strip()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('selected_npz');ap.add_argument('--out',required=True);ap.add_argument('--pair-source');a=ap.parse_args()
    start=time.monotonic();src=Path(a.selected_npz);out=Path(a.out);out.parent.mkdir(exist_ok=True,parents=True)
    with np.load(src) as d:
        tri=d['triangles_saved_mm' if 'triangles_saved_mm' in d.files else 'triangles_enu_mm'].copy();fid=d['face_ids'].copy()
    if a.pair_source:
        with np.load(a.pair_source) as d:ids=d['face_id_pairs']
        pairs=np.searchsorted(fid,ids)
        assert np.array_equal(fid[pairs],ids)
    else:
        lo=tri.min(axis=1);hi=tri.max(axis=1);pr=index.Property();pr.dimension=3
        rt=index.Index(((i,tuple(np.r_[lo[i],hi[i]]),None) for i in range(len(tri))),properties=pr)
        pairs=np.asarray([(i,j) for i in range(len(tri)) for j in rt.intersection(tuple(np.r_[lo[i],hi[i]])) if j>i],dtype=np.uint32)
    codes,log=run(tri,pairs,out)
    np.savez_compressed(out.with_suffix('.npz'),face_id_pairs=fid[pairs],classification_codes=codes)
    count={CLASSES[i]:int((codes==i).sum()) for i in CLASSES}
    bad=(codes>=2)&(codes<=6)
    examples=[]
    for code in [3,4,5,2,6]:
        for i in np.flatnonzero(codes==code)[:10]:examples.append(dict(classification=CLASSES[code],face_ids=fid[pairs[i]].tolist(),triangles_saved_coordinates_mm=tri[pairs[i]].tolist()))
    report=dict(schema='R14-exact-triangle-pair-predicate-v1',source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),selected_input_faces=len(tri),AABB_candidate_pairs=len(pairs),counts=count,abnormal_intersection_face_count=len(np.unique(fid[pairs[bad]])),abnormal_intersection_face_ids=np.unique(fid[pairs[bad]]).tolist(),proper_crossing_face_ids=np.unique(fid[pairs[codes==3]]).tolist(),examples=examples,arithmetic='Every finite binary64 coordinate decoded exactly; all positions multiplied by a common power of two to integers; GMP arbitrary-precision orient3D/orient2D signs and rational interval/segment comparisons. No geometric epsilon.',normal_adjacency='A pair is excluded only when its full intersection is contained in its exact shared vertex or shared edge. A shared vertex does not exempt an intersection elsewhere.',engine_stdout=log,limitations=['Pair domain is the saved selected-face set; this does not prove completeness for all faces of the source mesh.','Abnormal point/segment contact is distinguished from proper crossing and coplanar area overlap.','No wall-thickness, print or physical-flow qualification.'],elapsed_seconds=time.monotonic()-start)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['selected_input_faces','AABB_candidate_pairs','counts','abnormal_intersection_face_count','engine_stdout','elapsed_seconds']}),flush=True)

if __name__=='__main__':main()
