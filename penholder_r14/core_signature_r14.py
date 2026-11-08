"""Hash an oriented NPZ triangle surface independently of vertex/face ordering."""
import hashlib,struct
import numpy as np

def signature(path):
    with np.load(path) as z:
        vertices=np.asarray(z['vertices'],dtype=np.float64)
        faces=np.asarray(z['faces'],dtype=np.int64)
    used,face_remap=np.unique(faces,return_inverse=True)
    v,remap=np.unique(vertices[used],axis=0,return_inverse=True)
    f=remap[face_remap.reshape(-1,3)]
    start=np.argmin(f,axis=1)
    f=np.take_along_axis(f,(np.arange(3)[None,:]+start[:,None])%3,axis=1)
    f=f[np.lexsort((f[:,2],f[:,1],f[:,0]))]
    h=hashlib.sha256();h.update(b'R14 oriented triangle geometry v1\0')
    h.update(struct.pack('<QQ',len(v),len(f)))
    h.update(v.astype('<f8').tobytes());h.update(f.astype('<u4').tobytes())
    return h.hexdigest()
