#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
from openmesh import *
#Importing the rest of utilities
import numpy as np



#=====Create a new triangle mesh and populate it with test data=====#
#For reference, the OpenMesh python tutorial can be found at http://www.openmesh.org/media/Documentations/OpenMesh-Doc-Latest/a00020.html
pyramid = TriMesh()

#Add the corners of a simple pyramid. The function call returns the handle to the vertex, hence we name it 'vh' (vertex handle)
vh0 = pyramid.add_vertex(TriMesh.Point(0, 0, 1))
vh1 = pyramid.add_vertex(TriMesh.Point(1, -1, 0))
vh2 = pyramid.add_vertex(TriMesh.Point(-1, -1, 0))
vh3 = pyramid.add_vertex(TriMesh.Point(0, 1, 0))

#Define the faces of the pyramid
fh0 = pyramid.add_face(vh0, vh1, vh2)
fh1 = pyramid.add_face(vh0, vh2, vh3)
fh2 = pyramid.add_face(vh0, vh3, vh1)
fh2 = pyramid.add_face(vh1, vh2, vh3)

#Save
write_mesh(pyramid, "/home/jonatan/kuleuven-algorithms/src/data/pyramid.obj")

#Modify
#let's loop over all the vertices and move them one unit vector along their normal direction
oldPos = TriMesh.Point(0,0,0)
newPos = TriMesh.Point(0,0,0)
normal = TriMesh.Normal(0,0,0)
pyramid.request_vertex_normals()
pyramid.request_face_normals()
pyramid.update_normals()
pyramid.release_face_normals()
for vh in pyramid.vertices():
    oldPos = pyramid.point(vh)
    #print "oldPos" + ' %(x).2f %(y).2f %(z).2f' % {"x": oldPos[0], "y": oldPos[1], "z": oldPos[2]}
    normal = pyramid.normal(vh)
    #print "normal" + ' %(x).2f %(y).2f %(z).2f' % {"x": normal[0], "y": normal[1], "z": normal[2]}
    newPos = oldPos + normal*0.01
    #print "newPos" + ' %(x).2f %(y).2f %(z).2f' % {"x": newPos[0], "y": newPos[1], "z": newPos[2]}
    pyramid.set_point(vh,newPos)


#Save
write_mesh(pyramid, "/home/jonatan/kuleuven-algorithms/src/data/fucked_up_pyramid.obj")

##########
#DATA EXTRACTION
##########
#=====Let's see if we can extract all the data we need from a mesh=====#
#parameters needed for initialization
nVertices = pyramid.n_vertices()
nFaces = pyramid.n_faces()
maxValence = 7 #assumption that no vertex will be connected that more than 7 other vertices

#Extract positions, normals and connectivity
positions = np.zeros((nVertices,3), dtype = float)
normals = np.zeros((nVertices,3), dtype = float)
neighbourIndices = -1 * np.ones((nVertices, maxValence), dtype = int)

i = 0
for i, vh in enumerate(pyramid.vertices()):
    positions[i,0] = pyramid.point(vh)[0]
    positions[i,1] = pyramid.point(vh)[1]
    positions[i,2] = pyramid.point(vh)[2]
    normals[i,0] = pyramid.normal(vh)[0]
    normals[i,1] = pyramid.normal(vh)[1]
    normals[i,2] = pyramid.normal(vh)[2]

    for j, vvh in enumerate(pyramid.vv_iter(vh)):
        neighbourIndices[i,j] = pyramid.point(vvh).idx()

#TODO: verder het stukje hier net boven: zodat we indices van neighbours kunnen opslaan. Als te lang duurt, gewoon skippen want lijkt makkelijker in c++. Dan andersom de keten doen: van numpy naar OpenMesh. En dan wta algoritmes tussenin proberen!

print(positions)
print(normals)
print(neighbourIndices)
