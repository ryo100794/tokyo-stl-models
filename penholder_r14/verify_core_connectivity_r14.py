#!/usr/bin/env python3
"""Check final intended fluid cores, all69 port memberships and all91 core gaps.

These float64 cutter checks are paired with independent saved-STL passage
tests. They do not establish physical pigment flow or printer suitability.
"""
import argparse, gc, hashlib, json, os, time
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import trimesh
import manifold3d as mf
import vtk
from mesh_ray_checks_r13 import locator_for,first_hit,basis,ray_winding,probe_rays


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
def load_core(p):
    with np.load(p) as d:
        return mf.Manifold(mf.Mesh64(d['vertices'].astype(np.float64),d['faces'].astype(np.uint64)))


def design_mesh(root):
    path=root/'models/subway_R14_hollow_print.stl'
    m=trimesh.load_mesh(path,process=True)
    matrix=np.asarray(json.loads((root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    m.apply_transform(np.linalg.inv(matrix))
    return path,m


def connectivity(root,cache):
    data=json.loads((root/'data/rail_paths_R14.json').read_text());ports=json.loads((root/'data/ports_R14.json').read_text())['ports']
    rows=[];pairs=[];minimum=2.
    for i,line in enumerate(data['lines']):
        lid=line['line_id'];core=load_core(cache/f'{lid}_core.npz')
        components=core.decompose();positive=[q for q in components if q.volume()>1e-8]
        row=dict(line_id=lid,name_ja=line['name_ja'],core_sha256=digest(cache/f'{lid}_core.npz'),
            positive_fluid_components=len(positive),nonpositive_boundary_components=len(components)-len(positive),
            ports=sum(p['line_id']==lid for p in ports),port_membership=[])
        # A signed-distance test on the independent component checks both ends
        # of every nozzle; a port record alone is not connectivity evidence.
        from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray
        if len(positive)==1:
            raw=positive[0].to_mesh64();v=np.asarray(raw.vert_properties)[:,:3];f=np.asarray(raw.tri_verts)
            locator=locator_for(trimesh.Trimesh(v,f,process=False));poly=locator.GetDataSet()
            signed=vtk.vtkImplicitPolyDataDistance();signed.SetInput(poly)
            for p in ports:
                if p['line_id']!=lid:continue
                a=float(signed.EvaluateFunction(p['attachment_enu_mm']));b=float(signed.EvaluateFunction(p['mouth_center_enu_mm']))
                row['port_membership'].append(dict(port_id=p['port_id'],attachment_signed_mm=a,mouth_signed_mm=b,passed=a<0 and b<0))
            del locator,poly,signed,raw,v,f
        row['passed']=len(positive)==1 and row['ports']>=2 and len(row['port_membership'])==row['ports'] and all(p['passed'] for p in row['port_membership'])
        rows.append(row)
        for other in data['lines'][i+1:]:
            c=load_core(cache/f"{other['line_id']}_core.npz")
            gap=float(core.min_gap(c,2.));minimum=min(minimum,gap)
            pairs.append(dict(line_a=lid,line_b=other['line_id'],capped_search_distance_mm=2.,gap_mm=gap,passed=gap>=.5))
            del c;gc.collect()
        print('core',lid,'components',len(positive),'ports',row['ports'],flush=True)
        dump(root/'reports/lumen_connectivity_R14.json',dict(passed=all(r['passed'] for r in rows) and all(p['passed'] for p in pairs),
            lines=rows,pairs=pairs,minimum_interline_lumen_gap_mm=minimum,checked_all_pairs=len(pairs)==91,
            nominal_lumen_diameter_mm=1.2,physical_flow_test_performed=False))
        del core,positive,components;gc.collect()



if __name__=='__main__':
    p=argparse.ArgumentParser(description='Check final 14 intended cores and their 69 openings; saved-body passage is separate.')
    p.add_argument('--root',type=Path,required=True);p.add_argument('--cache',type=Path,required=True);a=p.parse_args()
    connectivity(a.root,a.cache)
