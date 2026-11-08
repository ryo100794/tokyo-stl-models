"""Export hull-generated individual R14 tubes without global mesh repair."""
from pathlib import Path
import argparse,hashlib,json,gc,time
import numpy as np
import trimesh
from geometry_r14 import load_npz, arrays
from robust_export_r14 import read_stl_bytes,weld_near_vertices,basic_topology
from build_r14 import topology


def export_one(obj,target,transform,simplifications=(.001,.002,.005,.01)):
    attempts=[]
    for tolerance in simplifications:
        source=obj.as_original().simplify(tolerance)
        v,f=arrays(source);v=v@transform[:3,:3].T+transform[:3,3]
        mesh,_=read_stl_bytes(trimesh.Trimesh(v,f,process=False))
        before=basic_topology(mesh)
        mesh,weld=weld_near_vertices(mesh,.00005)
        mesh,raw=read_stl_bytes(mesh)
        qa=basic_topology(mesh)
        attempt=dict(simplification_tolerance_mm=tolerance,before=before,weld=weld,after=qa)
        attempts.append(attempt)
        if qa['watertight'] and qa['winding_consistent'] and qa['zero_area_triangles']==0 and qa['volume_mm3']>0:
            complete=topology(mesh)
            if complete['nonmanifold_vertex_links']==0:
                target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
                return dict(file=target.name,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),
                            topology=complete,export_attempts=attempts,
                            geometry_method='polyhedral capsule Minkowski sweep; Boolean outer minus core',
                            surface_reconstruction_used=False,hole_filling_used=False,
                            physical_tests_performed=False)
        del source,v,f,mesh,raw;gc.collect()
    raise ValueError(f'No valid STL serialization for {target.name}: {attempts}')


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--cache',type=Path,required=True)
    a=p.parse_args();T=np.asarray(json.loads((a.root/'data/coordinate_systems_R14.json').read_text())['geographic_ENU_mm_to_print_XYZ_mm_row_major'])
    lines=json.loads((a.root/'data/rail_paths_R14.json').read_text())['lines'];records=[]
    for line in lines:
        lid=line['line_id'];start=time.time()
        part=load_npz(a.cache/f'{lid}_outer.npz')-load_npz(a.cache/f'{lid}_core.npz')
        row=export_one(part,a.root/'models'/f'rail_R14_{lid}.stl',T)
        row.update(line_id=lid,name_ja=line['name_ja'],seconds=time.time()-start)
        records.append(row);del part;gc.collect()
        (a.root/'reports/repaired_part_exports_R14.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
        print(lid,row['topology']['triangles'],row['seconds'],row['sha256'],flush=True)


if __name__=='__main__':main()
