"""Independent ray queries over saved triangle coordinates.

These helpers do not import the model generator or change the input geometry.
A VTK static triangle locator supports passage rays and oriented far-ray probes.
"""
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray


def locator_for(mesh):
    points=vtk.vtkPoints();points.SetData(numpy_to_vtk(np.asarray(mesh.vertices),deep=True))
    faces=np.asarray(mesh.faces,dtype=np.int64)
    cells=vtk.vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0,len(faces)*3+1,3,dtype=np.int64),deep=True),
                  numpy_to_vtkIdTypeArray(faces.ravel(),deep=True))
    poly=vtk.vtkPolyData();poly.SetPoints(points);poly.SetPolys(cells)
    loc=vtk.vtkStaticCellLocator();loc.SetDataSet(poly);loc.BuildLocator()
    return loc


def first_hit(locator,a,b):
    t=vtk.reference(0.);sub=vtk.reference(0);cell=vtk.reference(0)
    hit=[0.,0.,0.];pc=[0.,0.,0.]
    yes=locator.IntersectWithLine(a,b,1e-6,t,hit,pc,sub,cell)
    return {'hit':bool(yes),'fraction':float(t),'point_enu_mm':hit if yes else None,
            'triangle_cell_id':int(cell) if yes else None}



def basis(axis):
    n=np.cross(axis,[0.,0.,1.])
    if np.linalg.norm(n)<.1:n=np.cross(axis,[0.,1.,0.])
    n/=np.linalg.norm(n)
    return n,np.cross(axis,n)


def probe_rays(locator,mouth,axis,radius=.8,inner=False,attachment=None,length=None,body=None):
    n,b=basis(axis);rows=[]
    for r in [0.,radius]:
        for angle in ([0.] if r==0 else np.linspace(0,2*np.pi,8,endpoint=False)):
            offset=r*(n*np.cos(angle)+b*np.sin(angle))
            # Begin inside the rail/port overlap, so the attachment plane is
            # tested as well as the neck, taper, and mouth. The final cutter
            # extends 0.25 mm behind the attachment.
            start=(attachment-.10*axis if inner else mouth+.1*axis)+offset
            end=(mouth+.15*axis if inner else mouth+2.*axis)+offset
            hit=first_hit(locator,start,end)
            if hit['hit'] and body is not None:
                cell=hit['triangle_cell_id'];tri=body.vertices[body.faces[cell]].astype(float)
                area=float(np.linalg.norm(np.cross(tri[1]-tri[0],tri[2]-tri[0]))/2)
                hit['triangle_area_mm2']=area
            if inner and hit['hit']:
                distance=float((np.asarray(hit['point_enu_mm'])-attachment)@axis)
                hit['distance_along_nozzle_mm']=distance
                hit['near_neck_joint_plane']=abs(distance-2.)<2e-4
                hit['near_mouth_plane']=abs(distance-length)<2e-4
            rows.append(dict(hit,radius_mm=r,angle_rad=float(angle)))
    return rows


def ray_winding(locator,body,q):
    """Oriented intersection winding with paired zero-thickness faces cancelled."""
    direction=np.array([.324717957,.567143290,.756828461]);direction/=np.linalg.norm(direction)
    points=vtk.vtkPoints();ids=vtk.vtkIdList()
    locator.IntersectWithLine(q,q+1000*direction,1e-7,points,ids)
    events=[]
    for i in range(ids.GetNumberOfIds()):
        tri=body.vertices[body.faces[ids.GetId(i)]].astype(float)
        normal=np.cross(tri[1]-tri[0],tri[2]-tri[0]);dot=float(normal@direction)
        if abs(dot)<1e-14:continue
        t=float((np.asarray(points.GetPoint(i))-q)@direction)
        events.append((t,1 if dot>0 else -1))
    events.sort();groups=[]
    for t,sign in events:
        if groups and abs(t-groups[-1][0])<1e-5:groups[-1][1].add(sign)
        else:groups.append([t,{sign}])
    winding=sum(sum(signs) for _,signs in groups)
    return {'oriented_far_ray_winding':int(winding),'triangle_hits':len(events),'surface_crossing_clusters':len(groups),
            'empty_by_oriented_ray':winding==0}

