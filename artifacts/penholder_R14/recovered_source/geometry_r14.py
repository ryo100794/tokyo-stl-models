"""Geometry utilities for the pigment-channel edition of the Tokyo rail holder.

All calculation coordinates use millimetres, X east, Y north, Z up.  The
functions build positive solids and cutters separately so that every retained
station support and landscape is cut by the same lumen before export.
"""
from __future__ import annotations

import gc
import hashlib
import json
from pathlib import Path
import struct

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import manifold3d as mf
import trimesh


def read_glb(path):
    """Read triangle groups from the verified identity-node R12 reference."""
    raw = Path(path).read_bytes()
    assert struct.unpack_from('<III', raw) == (0x46546C67, 2, len(raw))
    chunks, offset = {}, 12
    while offset < len(raw):
        count, kind = struct.unpack_from('<II', raw, offset)
        chunks[kind] = raw[offset + 8:offset + 8 + count]
        offset += 8 + count
    doc, blob = json.loads(chunks[0x4E4F534A]), chunks[0x004E4942]
    def get(index):
        acc = doc['accessors'][index]
        view = doc['bufferViews'][acc['bufferView']]
        dtype = np.dtype({5121:'u1',5123:'<u2',5125:'<u4',5126:'<f4'}[acc['componentType']])
        columns = {'SCALAR':1,'VEC2':2,'VEC3':3,'VEC4':4}[acc['type']]
        return np.ndarray((acc['count'],columns), dtype=dtype, buffer=blob,
            offset=view.get('byteOffset',0)+acc.get('byteOffset',0),
            strides=(view.get('byteStride',dtype.itemsize*columns),dtype.itemsize)).copy()
    for n in doc['nodes']:
        assert 'matrix' not in n or np.allclose(np.asarray(n['matrix']).reshape(4,4),np.eye(4))
        assert np.allclose(n.get('translation',[0,0,0]),0)
        assert np.allclose(n.get('scale',[1,1,1]),1)
        assert np.allclose(n.get('rotation',[0,0,0,1]),[0,0,0,1])
    result=[]
    for mesh in doc['meshes']:
        assert len(mesh['primitives'])==1
        p=mesh['primitives'][0]
        v=get(p['attributes']['POSITION']).astype(float)[:,[0,2,1]]*[1000,-1000,1000]
        f=get(p['indices']).reshape(-1,3)
        mat=doc['materials'][p['material']]
        result.append({'name':mesh['name'],'vertices':v,'faces':f,
            'color':mat.get('pbrMetallicRoughness',{}).get('baseColorFactor',[.7,.7,.7,1])})
    return result,doc,hashlib.sha256(raw).hexdigest()


def as_manifold(vertices, faces, merge=False):
    """Convert indexed triangles, optionally weld the normal-split relief seam."""
    mesh=mf.Mesh64(np.asarray(vertices,dtype=np.float64),np.asarray(faces,dtype=np.uint64))
    if merge:
        mesh.merge()
    obj=mf.Manifold(mesh)
    if obj.status()!=mf.Error.NoError:
        raise ValueError(f'Manifold input status: {obj.status()}')
    return obj


def unite(parts):
    """Perform a true Boolean union, not an intersecting triangle concatenation."""
    parts=[p for p in parts if not p.is_empty()]
    if not parts:
        return mf.Manifold()
    return mf.Manifold.batch_boolean(parts,mf.OpType.Add)


def closed_group(group, contact_margin_mm=0.):
    """Weld display seams then union the individual closed input components."""
    try:
        obj=as_manifold(group['vertices'],group['faces'])
    except ValueError:
        obj=as_manifold(group['vertices'],group['faces'],merge=True)
    # Explicit indexed-component splitting avoids quadratic allocation in the
    # generic decompose path for the many tiny access-lattice members.
    v,f=arrays(obj)
    edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    graph=coo_matrix((np.ones(len(edges),dtype=np.uint8),(edges[:,0],edges[:,1])),
        shape=(len(v),len(v))).tocsr()
    count,labels=connected_components(graph,directed=False)
    if count==1 and contact_margin_mm == 0:return obj
    face_labels=labels[f[:,0]]
    order=np.argsort(face_labels,kind='stable')
    blocks=np.split(order,np.flatnonzero(np.diff(face_labels[order]))+1)
    parts=[]
    for block in blocks:
        ids,inv=np.unique(f[block].ravel(),return_inverse=True)
        vv=v[ids];ff=inv.reshape(-1,3)
        if contact_margin_mm:
            # The inherited lattice has tangent primitive contacts. Turn them
            # into positive overlaps before Boolean union, instead of patching
            # holes in a final STL. This explicitly expands each primitive by
            # a bounded vertex-normal displacement, recorded in the build.
            local=trimesh.Trimesh(vv,ff,process=False)
            if local.volume <= 0:
                raise ValueError('Contact expansion requires positive primitive components')
            vv=vv+contact_margin_mm*np.asarray(local.vertex_normals)
        parts.append(as_manifold(vv,ff))
    return unite(parts)


def arrays(obj):
    mesh=obj.to_mesh64()
    return np.asarray(mesh.vert_properties)[:,:3].copy(),np.asarray(mesh.tri_verts).copy()


def save_npz(obj,path):
    # Keep Boolean intersections in double precision until a delivery format
    # explicitly requires float32.  Repeated float32 checkpoint round trips
    # can collapse tiny gaps into coincident sheets.
    mesh=obj.to_mesh64()
    v=np.asarray(mesh.vert_properties)[:,:3].copy()
    f=np.asarray(mesh.tri_verts,dtype=np.uint32).copy()
    path = Path(path)
    temporary = path.with_name(path.name + '.pending')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream,vertices=v,faces=f)
    import zipfile
    with zipfile.ZipFile(temporary) as archive:
        if archive.testzip() is not None:
            raise ValueError('Incomplete mesh checkpoint')
    temporary.replace(path)


def load_npz(path):
    with np.load(path) as d:
        v=d['vertices']; f=d['faces']
    if v.dtype.itemsize==8:
        obj=mf.Manifold(mf.Mesh64(np.asarray(v,dtype=np.float64),
            np.asarray(f,dtype=np.uint64)))
        if obj.status()!=mf.Error.NoError:
            raise ValueError(f'Manifold checkpoint status: {obj.status()}')
        return obj
    return as_manifold(v,f)


def ring_profile_solid(origin,axis,profile,sides=24):
    """One closed solid with shared rings and no internal butt-joint faces.

    ``profile`` contains increasing (distance, radius) pairs along ``axis``.
    Rings at a neck, taper or mouth share indices with both adjacent spans;
    only the first and last rings are capped.
    """
    origin=np.asarray(origin,float);axis=np.asarray(axis,float)
    axis=axis/np.linalg.norm(axis)
    profile=np.asarray(profile,float)
    assert len(profile)>=2 and np.all(np.diff(profile[:,0])>0)
    assert np.all(profile[:,1]>0)
    helper=np.eye(3)[np.argmin(np.abs(axis))]
    x=np.cross(helper,axis);x/=np.linalg.norm(x)
    y=np.cross(axis,x)
    angles=np.arange(sides)*2*np.pi/sides
    radial=np.cos(angles)[:,None]*x+np.sin(angles)[:,None]*y
    vertices=np.concatenate([origin+d*axis+r*radial for d,r in profile])
    start_center=len(vertices);end_center=start_center+1
    vertices=np.vstack([vertices,origin+profile[0,0]*axis,origin+profile[-1,0]*axis])
    faces=[]
    for ring in range(len(profile)-1):
        a=ring*sides;b=(ring+1)*sides
        for i in range(sides):
            j=(i+1)%sides
            faces.extend([[a+i,a+j,b+j],[a+i,b+j,b+i]])
    last=(len(profile)-1)*sides
    for i in range(sides):
        j=(i+1)%sides
        faces.extend([[start_center,j,i],[end_center,last+i,last+j]])
    obj=mf.Manifold(mf.Mesh64(np.asarray(vertices,dtype=np.float64),
        np.asarray(faces,dtype=np.uint64)))
    assert obj.status()==mf.Error.NoError and obj.volume()>0
    return obj


def cylinder_between(a,b,r0,r1=None,sides=24):
    a,b=np.asarray(a,float),np.asarray(b,float)
    d=b-a
    length=float(np.linalg.norm(d))
    if length<1e-8:
        return mf.Manifold()
    z=d/length
    axis=np.eye(3)[np.argmin(np.abs(z))]
    x=np.cross(axis,z); x/=np.linalg.norm(x)
    y=np.cross(z,x)
    transform=np.column_stack([x,y,z,a])
    if r1 is None or r1 == r0:
        # Use the R14 indexed, consistently oriented shared tube API. Long or
        # sharply bent routes are a union of capped segments plus round joints;
        # a single ring sweep at a cusp can overlap itself.
        from subway_n02_welded_tube_model import make_welded_tube
        m = make_welded_tube([a.tolist(), b.tolist()], r0, sides, closed=False)
        return as_manifold(m.vertices, m.faces)
    return mf.Manifold.cylinder(length,r0,r1,sides).transform(transform)


def rdp(points,tolerance=.003):
    """Reduce nearly collinear points with a bounded geometric chord error."""
    p=np.asarray(points,float)
    if len(p)<=2:
        return p
    keep={0,len(p)-1}; stack=[(0,len(p)-1)]
    while stack:
        start,end=stack.pop()
        if end<=start+1: continue
        d=p[end]-p[start]; dd=np.dot(d,d)
        t=np.clip((p[start+1:end]-p[start])@d/max(dd,1e-30),0,1)
        distance=np.linalg.norm(p[start+1:end]-(p[start]+t[:,None]*d),axis=1)
        k=int(np.argmax(distance))
        if distance[k]>tolerance:
            mid=start+1+k;keep.add(mid);stack.extend([(start,mid),(mid,end)])
    return p[sorted(keep)]


def capsule_path(points,radius,sides=24,tolerance=.01):
    """Union cylinders and round joints for a continuous, closed sweep cutter."""
    p=rdp(points,tolerance)
    sphere=mf.Manifold.sphere(radius,sides)
    items=[cylinder_between(a,b,radius,sides=sides) for a,b in zip(p[:-1],p[1:])]
    items += [sphere.translate(q) for q in p]
    return unite(items)


def sample_path(points,step=.30):
    p=np.asarray(points,float)
    lens=np.linalg.norm(np.diff(p,axis=0),axis=1)
    s=np.r_[0,np.cumsum(lens)]
    samples=np.unique(np.r_[np.arange(0,s[-1],step),s[-1]])
    return np.column_stack([np.interp(samples,s,p[:,k]) for k in range(3)])


def path_point(points,svalue):
    p=np.asarray(points,float)
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(p,axis=0),axis=1))]
    i=int(np.clip(np.searchsorted(s,svalue,side='right')-1,0,len(p)-2))
    t=np.clip((svalue-s[i])/max(s[i+1]-s[i],1e-12),0,1)
    tangent=p[i+1]-p[i]; tangent/=np.linalg.norm(tangent)
    return p[i]*(1-t)+p[i+1]*t,tangent,i,float(t)


def write_stl(obj,path):
    """Export and re-read a binary STL in mm, retaining true lumen surfaces.

    A 0.00001 mm simplification removes numerical Boolean slivers.  Any
    zero-area or duplicate triangles created by float32 STL welding are
    removed explicitly; no holes are patched or filled.
    """
    export_object=obj.simplify(1e-5)
    v,f=arrays(export_object)
    mesh=trimesh.Trimesh(v,f,process=False)
    Path(path).write_bytes(trimesh.exchange.stl.export_stl(mesh))
    del v,f,mesh,export_object
    gc.collect()
    saved=trimesh.load_mesh(path,process=True)
    initial_faces=len(saved.faces)
    saved.update_faces(saved.nondegenerate_faces())
    saved.update_faces(saved.unique_faces())
    saved.remove_unreferenced_vertices()
    assert saved.is_watertight and saved.is_winding_consistent, 'Saved STL is not a closed consistently oriented mesh'
    cleaned_faces=len(saved.faces)
    saved._cache.clear()
    Path(path).write_bytes(trimesh.exchange.stl.export_stl(saved))
    del saved
    gc.collect()
    saved=trimesh.load_mesh(path,process=True)
    assert saved.is_watertight and saved.is_winding_consistent, 'Re-exported STL failed topology check'
    return {'vertices':len(saved.vertices),'triangles':len(saved.faces),'volume_mm3':float(saved.volume),
        'watertight':bool(saved.is_watertight),'winding_consistent':bool(saved.is_winding_consistent),
        'manifold_status':str(obj.status()),'saved_stl_reloaded':True,
        'export_simplification_tolerance_mm':1e-5,
        'post_stl_removed_degenerate_or_duplicate_triangles':initial_faces-len(saved.faces),
        'bounds_mm':saved.bounds.tolist(),'dimensions_mm':saved.extents.tolist(),
        'stl_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest()}


def write_glb(groups,path,extras):
    """Write compact indexed GLB, m units and glTF Y-up, from ENU mm groups."""
    from crease_normals import split_normals
    doc={'asset':{'version':'2.0','generator':'Tokyo Rail Holder R13'},
        'scene':0,'scenes':[{'nodes':[]}],'nodes':[],'meshes':[],'materials':[],
        'buffers':[{'byteLength':0}],'bufferViews':[],'accessors':[],'extras':extras}
    binary=bytearray()
    def add(values,component,kind,target):
        nonlocal binary
        while len(binary)%4:binary+=b'\0'
        a=np.ascontiguousarray(values)
        view=len(doc['bufferViews']); start=len(binary);binary+=a.tobytes()
        doc['bufferViews'].append({'buffer':0,'byteOffset':start,'byteLength':a.nbytes,'target':target})
        acc={'bufferView':view,'componentType':component,'count':len(a),'type':kind}
        if kind=='VEC3':acc.update(min=a.min(0).tolist(),max=a.max(0).tolist())
        index=len(doc['accessors']);doc['accessors'].append(acc);return index
    for g in groups:
        if len(g['faces'])==0:continue
        v=(np.asarray(g['vertices'])[:,[0,2,1]]*[.001,.001,-.001]).astype('<f4')
        f=np.asarray(g['faces'],dtype='<u4')
        v,f,normals=split_normals(v,f,38.)
        v=np.asarray(v,dtype='<f4')
        f=np.asarray(f,dtype='<u4')
        normals=np.asarray(normals,dtype='<f4')
        pi=add(v,5126,'VEC3',34962);ni=add(normals,5126,'VEC3',34962)
        fi=add(f.ravel(),5125,'SCALAR',34963)
        i=len(doc['meshes']); color=list(g.get('color',[.75,.75,.75,1]))
        if len(color)==3:color.append(1)
        doc['materials'].append({'name':g['name'],'pbrMetallicRoughness':{
            'baseColorFactor':color,'metallicFactor':0,'roughnessFactor':.52},'doubleSided':False})
        doc['meshes'].append({'name':g['name'],'primitives':[{
            'attributes':{'POSITION':pi,'NORMAL':ni},'indices':fi,'material':i}]})
        doc['nodes'].append({'name':g['name'],'mesh':i})
        doc['scenes'][0]['nodes'].append(i)
    doc['buffers'][0]['byteLength']=len(binary)
    j=json.dumps(doc,ensure_ascii=False,separators=(',',':')).encode()
    j+=b' '*((-len(j))%4);binary+=b'\0'*((-len(binary))%4)
    result=struct.pack('<III',0x46546C67,2,12+8+len(j)+8+len(binary))
    result+=struct.pack('<II',len(j),0x4E4F534A)+j
    result+=struct.pack('<II',len(binary),0x004E4942)+binary
    Path(path).write_bytes(result)


def object_group(name,obj,color):
    v,f=arrays(obj)
    return {'name':name,'vertices':v,'faces':f,'color':color}
