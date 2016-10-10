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
-kappa: mahalanobis distance for outlier determination (e.g. kappa = 2 means
that 2 standard deviations of the distribution is the middle of the transition
zone from inlier to outlier classification)
"""

##########
# SET PARAMETERS
##########
maxNumIterations = 1
maxNumNeighbours = 20 #
kappa = 3

##########
# PREPARE DATA
##########


##Load meshes
floatingMeshPath = "/home/jonatan/kuleuven-algorithms/src/data/bunny90.obj"
targetMeshPath = "/home/jonatan/kuleuven-algorithms/src/data/fucked_up_bunny.obj"
floatingMesh = TriMesh()
read_mesh(floatingMesh, floatingMeshPath)
targetMesh = TriMesh()
read_mesh(targetMesh, targetMeshPath)

###construct the corresponding mesh by deep copying the floating mesh
correspondingMesh = TriMesh()
read_mesh(correspondingMesh, floatingMeshPath)

###make sure the vertex normals are there
floatingMesh.request_vertex_normals()
floatingMesh.request_face_normals()
floatingMesh.update_normals()
floatingMesh.release_face_normals()
targetMesh.request_vertex_normals()
targetMesh.request_face_normals()
targetMesh.update_normals()
targetMesh.release_face_normals()
correspondingMesh.request_vertex_normals()
correspondingMesh.request_face_normals()
correspondingMesh.update_normals()
correspondingMesh.release_face_normals()

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
### Initialize weights
floatingWeights = np.ones((numFloatingVertices,1), dtype = float)
targetWeights = np.ones((numFloatingVertices,1), dtype = float)

###Extract data from meshes
for i, vh in enumerate(floatingMesh.vertices()):
    floatingPositions[i,0] = floatingMesh.point(vh)[0]
    floatingPositions[i,1] = floatingMesh.point(vh)[1]
    floatingPositions[i,2] = floatingMesh.point(vh)[2]
    floatingNormals[i,0] = floatingMesh.normal(vh)[0]
    floatingNormals[i,1] = floatingMesh.normal(vh)[1]
    floatingNormals[i,2] = floatingMesh.normal(vh)[2]
    for j,vvh in enumerate(floatingMesh.vv(vh)):
        floatingNeighbours[i,j] = vvh.idx()
        #warning: make sure you limit this loop to number of maxNumNeighbours

for i, vh in enumerate(targetMesh.vertices()):
    targetPositions[i,0] = targetMesh.point(vh)[0]
    targetPositions[i,1] = targetMesh.point(vh)[1]
    targetPositions[i,2] = targetMesh.point(vh)[2]
    targetNormals[i,0] = targetMesh.normal(vh)[0]
    targetNormals[i,1] = targetMesh.normal(vh)[1]
    targetNormals[i,2] = targetMesh.normal(vh)[2]
    for j,vvh in enumerate(targetMesh.vv(vh)):
        targetNeighbours[i,j] = vvh.idx()
        #warning: make sure you limit this loop to number of maxNumNeighbours

##########
# ICP
##########
"""

"""
for iteration in range(0,maxNumIterations):
    ##1) Determine Nearest neighbours. We'll simply use index correspondences for the bunny
    correspondingPositions = targetPositions
    correspondingNormals = targetNormals
    ##2) Determine weights. A weight related to the gaussian distance distribution suffices.
    ### Update the distribution parameters
    sigmaNumerator = 0
    sigmaDenominator = 0
    for i in range(numFloatingVertices):
        distance = np.linalg.norm(correspondingPositions[i,] - floatingPositions[i,])
        sigmaNumerator = sigmaNumerator + floatingWeights[i] * distance * distance
        sigmaDenominator = sigmaDenominator + floatingWeights[i]
        print(distance)
        print(sigmaNumerator)
        print(sigmaDenominator)
    sigma = np.sqrt(sigmaNumerator/sigmaDenominator)
    lambdaa = 1.0/(np.sqrt(2 * 3.14159) * sigma) * np.exp(-0.5 * kappa * kappa)
    print(sigma)
    print(lambdaa)
    ### Recalculate the weights
    for i, vh in enumerate(floatingMesh.vertices()):
        distance = np.linalg.norm(correspondingPositions[i,] - floatingPositions[i,])
        inlierProbability = 1.0/(np.sqrt(2 * 3.14159) * sigma) * np.exp(-0.5 * np.square(distance/sigma))
        floatingWeights[i] = inlierProbability / (inlierProbability + lambdaa)
        print(distance)
        print(inlierProbability)
        print(floatingWeights[i])
    ##3) Determine and update transformation.
    ###3.1 Compute the centroids of the Floating and Corresponding mesh
    floatingCentroid = np.array([0.0,0.0,0.0], dtype = float)
    correspondingCentroid = np.array([0.0,0.0,0.0], dtype = float)
    weightSum = 0.0
    for i in range(numFloatingVertices):
        floatingCentroid = floatingCentroid + floatingWeights[i] * floatingPositions[i,]
        correspondingCentroid = correspondingCentroid + floatingWeights[i] * correspondingPositions[i,]
        weightSum = weightSum + floatingWeights[i]
    floatingCentroid = floatingCentroid / weightSum
    correspondingCentroid = correspondingCentroid / weightSum
    ###3.2 Compute the Cross Variance matrix
    crossVarianceMatrix = numpy.zeros((3,3), dtype = float)
    for i in range(numFloatingVertices):
        #CrossVarMat = sum(weight[i] * floatingPosition[i] * correspondingPosition[i]_Transposed)
        crossVarianceMatrix = crossVarianceMatrix + floatingWeights[i] * np.outer(floatingPositions[i,:], correspondingPositions[i,:])
    crossVarianceMatrix = crossVarianceMatrix / weightSum - np.outer(floatingCentroid, correspondingCentroid)
    ###3.3 Compute the Anti-Symmetric matrix
    antiSymmetricMatrix = crossVarianceMatrix - crossVarianceMatrix.transpose()
    ###3.4 Use the cyclic elements of the Anti-Symmetric matrix to construct delta
    delta = ones(())
    ##4) Apply transformation
