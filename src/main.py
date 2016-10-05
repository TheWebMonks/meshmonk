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

##########
# RENDER #
##########

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(reader.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Create a rendering window and renderer
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)

# Create a renderwindowinteractor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# Assign actor to the renderer
ren.AddActor(actor)

# Enable user interface interactor
iren.Initialize()
renWin.Render()
iren.Start()

################
# VTK TO NUMPY #
################
#nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
#nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)


#########
# WRITE #
#########
stlWriter = vtk.vtkSTLWriter()
stlWriter.SetFileName(outputfilename)
stlWriter.SetInputConnection(reader.GetOutputPort())
stlWriter.Write()
