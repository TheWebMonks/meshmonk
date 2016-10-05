import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa

inputfilename = "/src/data/paris/558/skin95.stl"
outputfilename = "/out/skin95out.stl"

########
# READ #
########

reader = vtk.vtkSTLReader()
reader.SetFileName(inputfilename)
reader.Update()

polydata = vtk.vtkPolyData()
polydata = reader.GetOutput() #type vtkPolyData
points = polydata.GetPoints()

"""
NOTE: normals, vectors or scalars cannot be retrieved from the polydata object if this is not available in the file
that has been loaded by the reader. Typically, this will require applying some filters to calculate that data (especially
normals). So for now, we only know the locations (x,y,z) and the connectivity.
"""
#pointdata = polydata.GetPointData() #type vtkPointData
#vectors = pointdata.GetVectors() #type vtkDataArray
#normals = pointdata.GetNormals() #type vtkDataArray


################
# VTK TO NUMPY #
################
coordinates = points.GetData()
coordinatesPy = vtk_to_numpy(coordinates)

print(coordinatesPy.ndim)
print(coordinatesPy.shape)
print(coordinatesPy.size)
print(coordinatesPy.dtype)
print(coordinatesPy.itemsize)
