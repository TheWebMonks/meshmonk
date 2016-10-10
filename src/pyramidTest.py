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
fh3 = pyramid.add_face(vh1, vh3, vh2)

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
maxValenceAssumed = 7 #assumption that no vertex will be connected that more than 7 other vertices

#Extract positions, normals and connectivity
positions = np.zeros((nVertices,3), dtype = float)
normals = np.zeros((nVertices,3), dtype = float)
faces = np.zeros((nFaces,3), dtype = int)
neighbourIndices = -1 * np.ones((nVertices, maxValenceAssumed), dtype = int)

face = TriMesh.Face(0,0,0)
for i, fh in enumerate(pyramid.faces()): #loop over faces
    for j,fvh in enumerate(pyramid.fv(fh)): #for current face, loop over all vertices (fv = face-vertex iterator)
        faces[i,j] = fvh.idx()
        #TODO: DEZE LOOP DEBUGGEN

maxValence = 0 #we save the maximum valence so that we can release array memory later
for i, vh in enumerate(pyramid.vertices()):
    positions[i,0] = pyramid.point(vh)[0]
    positions[i,1] = pyramid.point(vh)[1]
    positions[i,2] = pyramid.point(vh)[2]
    normals[i,0] = pyramid.normal(vh)[0]
    normals[i,1] = pyramid.normal(vh)[1]
    normals[i,2] = pyramid.normal(vh)[2]
    for j, vvh in enumerate(pyramid.vv(vh)):
        neighbourIndices[i,j] = vvh.idx()
        if (j+1) > maxValence:
            maxValence = j+1

#release memory for neighbourIndices matrix
neighbourIndices = np.delete(neighbourIndices,range(maxValence,maxValenceAssumed),1) #delete the first column that just contains nCorners in each element

##########
#DATA INSERTION
##########
#=====Let's see if we can insert all the data we need into both an new mesh=====#
"""
The data formats we assume to be working with to represent a triangular mesh is the following:
-positions: (nVertices x 3) numpy matrix of floats/doubles with each row the x,y,z coordinates for each vertex
-normals: (nVertices x 3) numpy matrix of floats/doubles with each row the x,y,z values for each vertex (!) normal (not the face normal)
-neighbours: (nVertices x maxNumNeighbours) numpy matrix of integers with each row the indices of a vertex's neighbouring vertices (elements equal to -1 are used to fill up the matrix)
-faces: (nFaces x 3) numpy of integers with each row the three indices of the vertices that belong to that face.

NOTE: neighbours and faces contains similar information. However, we provision the possibility that the neighbouring vertices
can later be determined not only through connecting edges, but can also be found through a real neighbourhood search.
"""
pyramidNew = TriMesh()

pos = TriMesh.Point(0,0,0)
vh_list = []
for i in range(nVertices):
    pos[0] = positions[i,0]
    pos[1] = positions[i,1]
    pos[2] = positions[i,2]
    vh = pyramidNew.add_vertex(pos)
    vh_list.append(vh)

id = [ 0 for i in range(3)]
for i in range(nFaces):
    id[0] = faces[i,0]
    id[1] = faces[i,1]
    id[2] = faces[i,2]
    fh = pyramidNew.add_face(vh_list[id[0]], vh_list[id[1]], vh_list[id[2]])
    print(id)

#Save
write_mesh(pyramidNew, "/home/jonatan/kuleuven-algorithms/src/data/fucked_up_pyramid.obj")
