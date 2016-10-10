#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
from openmesh import *
#Importing the rest of utilities
import numpy as np

"""
In this example, we will perform rigid registration of the original
(downsampled) bunny to the fucked up bunny.

#####
#Some definitions and info:
#####

MESHES
-floatingMesh: the mesh that we will register, and hence will be transformed.
This mesh is both an input and an output object.
-targetMesh: the mesh we're trying to register to, hence 'target'. This mesh is
just an input object.
-correspondingMesh: for each vertex of the floating mesh, we're gonna look for
corresponding points in the target mesh. So what we're going to do is, make a
copy of the floating mesh, and store the corresponding points by inserting the
corresponding positions into this mesh copy. This mesh is an internal object.

WEIGHTS
In the computation of the correspondences and the transformation, the
contribution for each vertex of both the floating mesh and the target mesh can
be weighed. E.g. when there are abnormal regions in one of the meshes, this can
be detected and flagged by assigning a weight equal to zero for the vertices
that belong to that region.
-floatingWeights: weights for the vertices of the floating mesh. This array is
both an input and an output object.
-targetWeights: weights for the vertices of the target mesh. This array is both
an input and an output object.

PROCESS PARAMETERS
-maxNumIterations: a hard limit for the maximum number of iterations in the ICP
iterative loop.
-maxNumNeighbours: the maximum number of vertices that can be considered one
vertex's neighbours. (Neighbours can be determined in various ways. They are
important for some operations where neighbouring vertices play a role in the
computation.)
"""

##########
# SET PARAMETERS
##########
maxNumIterations = 10
maxNumNeighbours = 10 #

##########
# PREPARE DATA
##########
floatingMeshPath = "/home/jonatan/kuleuven-algorithms/src/data/bunny90.obj"
targetMeshPath = "/home/jonatan/kuleuven-algorithms/src/data/fucked_up_bunny.obj"

##Load meshes
floatingMesh = TriMesh()
read_mesh(floatingMesh, floatingMeshPath)

targetMesh = TriMesh()
read_mesh(targetMesh, targetMeshPath)

###construct the corresponding mesh by deep copying the floating mesh
correspondingMesh = TriMesh()
read_mesh(correspondingMesh, floatingMeshPath)

##Convert to the matrices we need
###Obtain info for initialization
numFloatingVertices = floatingMesh.n_vertices()
numFloatingFaces = floatingMesh.n_faces()
numTargetVertices = targetMesh.n_vertices()
numTargetFaces = targetMesh.n_faces()

###Initialize positions
floatingPositions = np.zeros((numFloatingVertices,3), dtype = float)
correspondingPositions = np.zeros((numFloatingVertices,3), dtype = float)
targetPositions = np.zeros((numTargetVertices,3), dtype = float)
### Initialize normals
floatingNormals = np.zeros((numFloatingVertices,3), dtype = float)
correspondingNormals = np.zeros((numFloatingVertices,3), dtype = float)
targetNormals = np.zeros((numTargetVertices,3), dtype = float)
### Initialize neighbours (these matrices are filled with '-1' elements)
floatingNeighbours = -1 * np.ones((numFloatingVertices, maxNumNeighbours), dtype = int)
correspondingNeighbours = floatingNeighbours #a shallow copy suffices
targetNeighbours = -1 * np.ones((numTargetVertices, maxNumNeighbours), dtype = int)
