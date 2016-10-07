import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

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
nPoints = coordinatesPy.shape[0];

polygons = polydata.GetPolys().GetData()
polygonsPy = vtk_to_numpy(polygons)
"""
The very first number of the polygonsPy-array gives us the number of vertices per face.
We're assuming all the faces have the same number of vertices! So for a triangular mesh,
this would look something like:
polygonsPy = [ 3 12 11 24 3 15 27 47 3 24 65 54 ...]

Where each fourth number is a '3' which means each polygon has three vertices. Behind each 3 then are
three numbers that give the indices of the vertices which belong to that triangle.

So let's reshape the polygonsPy array so we would get:
polygonsPy = [
12 11 24
15 27 47
24 65 54
...
]
"""
nCorners = polygonsPy[0]
nPolygons = polygonsPy.size/(nCorners+1) #the size of polygonsPy is nPolygons * (nCorners +1)
polygonsPy = polygonsPy.reshape(nPolygons,nCorners+1)
polygonsPy = np.delete(polygonsPy,0,1) #delete the first column that just contains nCorners in each element



print(coordinatesPy.ndim)
print(coordinatesPy.shape)
print(coordinatesPy.size)
print(coordinatesPy.dtype)
print(coordinatesPy.itemsize)
