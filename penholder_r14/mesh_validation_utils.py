"""Portable array and read-only mesh-query helpers for independent validation."""
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray

def load_npz_arrays(path):
 with np.load(path,allow_pickle=False) as a:return a['vertices'].copy(),a['faces'].copy()
def subset(v,f,ids):
 sel=f[ids];vi,inv=np.unique(sel.reshape(-1),return_inverse=True)
 return v[vi],inv.reshape(-1,3)
def locator(v,f):
 pts=vtk.vtkPoints();pts.SetData(numpy_to_vtk(v,deep=True))
 cells=vtk.vtkCellArray();cells.SetData(numpy_to_vtkIdTypeArray(np.arange(len(f)+1,dtype=np.int64)*3,deep=True),numpy_to_vtkIdTypeArray(f.astype(np.int64).ravel(),deep=True))
 poly=vtk.vtkPolyData();poly.SetPoints(pts);poly.SetPolys(cells)
 loc=vtk.vtkStaticCellLocator();loc.SetDataSet(poly);loc.BuildLocator();return loc,poly
def closest_vertices(points,loc):
 ds=[];cps=[];cs=[]
 for p in points:
  cp=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);dd=vtk.reference(0.)
  loc.FindClosestPoint(p,cp,cell,sub,dd);ds.append(float(dd)**.5);cps.append(cp);cs.append(int(cell))
 return np.asarray(ds),np.asarray(cps),np.asarray(cs)
