import vtk
import numpy as np
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

"""
Let's create some test data: a pyramid with only 4 nodes.
"""
coordinatesPy = np.array([(0,0,1),(1,-1,0),(-1,-1,0),(0,1,0)], dtype = float)
connectivityPy = np.array([(0,1,2),(0,2,3),(0,3,1),(1,2,3)], dtype = int)


"""
Now we'll attempt to put this data into a vtkPolyData object so that we'll be able to perform VTK filters
on it later, when we need to.
"""

#Deduce some parameters
nPoints = coordinatesPy.shape[0]
nFaces = connectivityPy.shape[0]
nCorners = connectivityPy.shape[1] #the number of corners for the mesh's polygon faces

###
#Set the coordinates of each point
###
coordinatesVtk = numpy_to_vtk(coordinatesPy)
pointsVtk = vtk.vtkPoints()
pointsVtk.SetData(coordinatesVtk)

###
#Set the connectivity. We have to rework our connectiviy matrix a bit, though.
###
#Reworking the connectivity matrix to the format that VTK requires
nCornersArray = np.ones((nFaces,1),dtype = int)*nCorners
connectivityPy = np.hstack((nCornersArray, connectivityPy))
connectivityPy = np.ravel(connectivityPy,order='C') #flatten the matrix to C-style array
connectivityVtk = numpy_to_vtkIdTypeArray(connectivityPy)

#Building cells object that the vtkPolyData object needs, where we set the connectivity.
cellsVtk = vtk.vtkCellArray()
cellsVtk.SetCells(nFaces,connectivityVtk)

###
#Insert the Coordinates and Faces into VTK's PolyData object
polydata = vtk.vtkPolyData()
polydata.SetPoints(pointsVtk)
polydata.SetPolys(cellsVtk)


##########
#CHECK: Output what we have to control what we did.
##########
writer = vtk.vtkXMLPolyDataWriter();
writer.SetFileName("/out/polydata_test.vtp");
writer.SetInputData(polydata)
writer.Write()

#=====Let's see if we can obtain the normals now=====#
# normalsFilter = vtk.vtkPolyDataNormals()
# normalsFilter.ComputePointNormalsOn()
# normalsFilter.ComputeCellNormalsOff()
# normalsFilter.SetInputData(polydata)
# polydata2
# normalsFilter.SetOutput(polydata)
