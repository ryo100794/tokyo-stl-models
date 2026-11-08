"""Split display vertices at sharp edges without altering any triangle positions."""
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray,vtk_to_numpy


def split_normals(vertices,faces,angle=38.):
    vertices=np.asarray(vertices);faces=np.asarray(faces,dtype=np.int64)
    points=vtk.vtkPoints();points.SetData(numpy_to_vtk(vertices,deep=True))
    cells=vtk.vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0,len(faces)*3+1,3,dtype=np.int64),deep=True),
                  numpy_to_vtkIdTypeArray(faces.ravel(),deep=True))
    data=vtk.vtkPolyData();data.SetPoints(points);data.SetPolys(cells)
    normal=vtk.vtkPolyDataNormals();normal.SetInputData(data);normal.SetFeatureAngle(angle)
    normal.SplittingOn();normal.ConsistencyOff();normal.AutoOrientNormalsOff()
    normal.ComputePointNormalsOn();normal.ComputeCellNormalsOff();normal.Update()
    out=normal.GetOutput()
    result_vertices=vtk_to_numpy(out.GetPoints().GetData()).copy()
    result_faces=vtk_to_numpy(out.GetPolys().GetData()).reshape(-1,4)[:,1:].astype(np.uint32)
    result_normals=vtk_to_numpy(out.GetPointData().GetNormals()).copy()
    if result_faces.shape!=faces.shape:raise AssertionError('Normal filter changed triangle count')
    for start in range(0,len(faces),50000):
        if not np.array_equal(vertices[faces[start:start+50000]],result_vertices[result_faces[start:start+50000]]):
            raise AssertionError('Normal filter changed triangle positions or winding')
    return result_vertices,result_faces,result_normals
