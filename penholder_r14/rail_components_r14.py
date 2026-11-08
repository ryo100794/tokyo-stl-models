"""R14 line solids derived from the recovered R13 generator.

Saved R13 centre lines, marker adjustments and 69 nozzle placements are used
without rerouting. Tube segments use the shared R14 sweep implementation.
"""
from __future__ import annotations
import time,hashlib,json
import numpy as np
from geometry_r14 import *
RO=1.2
RI=.6
STATION_RADIUS=1.42
PORT_LENGTH=3.6
PORT_MOUTH_RADIUS=1.0
PORT_OUTER_RADIUS=1.6
PORT_NECK_LENGTH=2.0
PORT_ATTACHMENT_OVERLAP=.25

def distance_point_segments(point,points):
    p=np.asarray(points,float);a=p[:-1];d=np.diff(p,axis=0)
    t=np.clip(np.sum((point-a)*d,axis=1)/np.maximum(np.sum(d*d,axis=1),1e-15),0,1)
    dist=np.linalg.norm(a+t[:,None]*d-point,axis=1)
    k=int(np.argmin(dist))
    return float(dist[k]),k,float(t[k])

def recover_displacement(line,point):
    """Interpolate the original-to-adjusted path displacement at a station."""
    best=None
    for path in line['paths']:
        original=np.asarray(path.get('original_points_enu_mm',path['points_enu_mm']))
        updated=np.asarray(path['points_enu_mm'])
        d,i,t=distance_point_segments(np.asarray(point),original)
        shift=(updated[i]-original[i])*(1-t)+(updated[i+1]-original[i+1])*t
        if best is None or d<best[0]:best=(d,shift,path['trace_id'],i,t)
    return best

def color_for(line):
    m=line['source_material']
    return m.get('pbrMetallicRoughness',{}).get('baseColorFactor',[.6,.6,.6,1])

def nozzle(spec):
    p=np.asarray(spec['attachment_enu_mm']);u=np.asarray(spec['axis_enu']);length=spec['length_mm']
    outer=ring_profile_solid(p,u,[(0,RO),(PORT_NECK_LENGTH,RO),
        (length,PORT_OUTER_RADIUS)])
    overlap=float(spec.get('attachment_overlap_mm',PORT_ATTACHMENT_OVERLAP))
    core_body=ring_profile_solid(p,u,[(-overlap,RI),(PORT_NECK_LENGTH,RI),
        (length,PORT_MOUTH_RADIUS),(length+.8,PORT_MOUTH_RADIUS)])
    # The rail sweep already surrounds the attachment point.  The nozzle
    # intersects it volumetrically; a short overlap also avoids a coplanar
    # cutter cap at the outer inlet face.  No duplicate tangent sphere is used.
    return outer,core_body

def build_line(line,ports,cache):
    t=time.time();lid=line['line_id']
    outer_path=cache/f'{lid}_outer.npz';core_path=cache/f'{lid}_core.npz'
    receipt=cache/f'{lid}_cache_key.json'
    config=dict(algorithm='polyhedral_capsule_hull_v1',outer_radius_mm=RO,
                inner_radius_mm=RI,line=line,ports=[p for p in ports if p['line_id']==lid],
                geometry_source_sha256=hashlib.sha256(Path(__file__).with_name('geometry_r14.py').read_bytes()).hexdigest(),
                component_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    key=hashlib.sha256(json.dumps(config,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    if outer_path.exists() and core_path.exists() and receipt.exists():
        saved=json.loads(receipt.read_text())
        if saved.get('input_and_code_sha256')==key:
            return load_npz(outer_path),load_npz(core_path)
    positive=[];negative=[]
    for path in line['paths']:
        p=path['points_enu_mm']
        print('path',path['trace_id'],'points',len(p),'reduced',len(rdp(p)),flush=True)
        positive.append(capsule_path(p,RO));negative.append(capsule_path(p,RI))
    for station in line['stations']:
        c=np.asarray(station['display_center_enu_mm']);distance,shift,*_=recover_displacement(line,c)
        marker=np.asarray(station['marker_centerline_enu_mm'])
        if distance>1.0:
            # Preserve the two inherited outlying decorative markers; never
            # invent a pigment branch or named station coordinate for them.
            positive.append(capsule_path(marker,.9))
            continue
        marker=marker+shift
        positive.append(capsule_path(marker,STATION_RADIUS,tolerance=0.))
        if np.linalg.norm(shift)>.03:
            positive.append(capsule_path([c,c+shift],.68,tolerance=0.))
    for port in ports:
        if port['line_id']!=lid:continue
        outer,inner=nozzle(port);positive.append(outer);negative.append(inner)
    outer,core=unite(positive),unite(negative)
    assert outer.status()==mf.Error.NoError and core.status()==mf.Error.NoError
    save_npz(outer,outer_path);save_npz(core,core_path)
    receipt.write_text(json.dumps({'input_and_code_sha256':key,'algorithm':config['algorithm']},indent=2)+'\n')
    print('line done',lid,'triangles',outer.num_tri(),core.num_tri(),'seconds',time.time()-t,flush=True)
    return outer,core
