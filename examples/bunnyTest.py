#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
from openmesh import *
#Importing the rest of utilities
import numpy as np


#=====Play around with the Stanford bunny=====#
#Load
bunny = TriMesh()
read_mesh(bunny, "/home/jonatan/kuleuven-algorithms/src/data/bunny90.obj")

#Modify
#let's loop over all the vertices and move them one unit vector along their normal direction
oldPos = TriMesh.Point(0,0,0)
newPos = TriMesh.Point(0,0,0)
normal = TriMesh.Normal(0,0,0)
bunny.request_vertex_normals()
bunny.request_face_normals() #we need the face normals to calculate the vertex normals
bunny.update_normals()
bunny.release_face_normals(); #after the calculation we release the face normals since we don't need them anymore
for vh in bunny.vertices():
    oldPos = bunny.point(vh)
    #print "oldPos" + ' %(x).2f %(y).2f %(z).2f' % {"x": oldPos[0], "y": oldPos[1], "z": oldPos[2]}
    normal = bunny.normal(vh)
    #print "normal" + ' %(x).2f %(y).2f %(z).2f' % {"x": normal[0], "y": normal[1], "z": normal[2]}
    newPos = oldPos + 0.01*normal
    #print "newPos" + ' %(x).2f %(y).2f %(z).2f' % {"x": newPos[0], "y": newPos[1], "z": newPos[2]}
    bunny.set_point(vh,newPos)


#Save
write_mesh(bunny, "/home/jonatan/kuleuven-algorithms/src/data/fucked_up_bunny.obj")
