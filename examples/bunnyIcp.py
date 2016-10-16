#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
from openmesh import *
#Importing the rest of utilities
import numpy as np
from scipy import linalg

floatingMeshPath = "/home/jonatan/kuleuven-algorithms/src/data/bunny90.obj"
targetMeshPath = "/home/jonatan/kuleuven-algorithms/src/data/fucked_up_bunny2.obj"
outputMeshName = "/home/jonatan/kuleuven-algorithms/src/data/registered_bunny.obj"


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
-kappaa: mahalanobis distance for outlier determination (e.g. kappaa = 2 means
that 2 standard deviations of the distribution is the middle of the transition
zone from inlier to outlier classification)
-adjustScale: whether or not to change the size of the floating mesh during ICP
"""

##########
# SET PARAMETERS
##########
maxNumIterations = 1
maxNumNeighbours = 20 #
kappaa = 3
adjustScale = False

##########
# PREPARE DATA
##########


##Load meshes
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
floatingPositions = np.zeros((3,numFloatingVertices), dtype = float)
correspondingPositions = np.zeros((3,numFloatingVertices), dtype = float)
targetPositions = np.zeros((3,numTargetVertices), dtype = float)
### Initialize normals
floatingNormals = np.zeros((3,numFloatingVertices), dtype = float)
correspondingNormals = np.zeros((3,numFloatingVertices), dtype = float)
targetNormals = np.zeros((3,numTargetVertices), dtype = float)
### Initialize neighbours (these matrices are filled with '-1' elements)
floatingNeighbours = -1.0 * np.ones((maxNumNeighbours, numFloatingVertices), dtype = int)
correspondingNeighbours = floatingNeighbours #a shallow copy suffices
targetNeighbours = -1.0 * np.ones((maxNumNeighbours, numTargetVertices), dtype = int)
### Initialize weights
floatingWeights = np.ones((numFloatingVertices), dtype = float)
targetWeights = np.ones((numFloatingVertices), dtype = float)

###Extract data from meshes
for i, vh in enumerate(floatingMesh.vertices()):
    floatingPositions[0,i] = floatingMesh.point(vh)[0]
    floatingPositions[1,i] = floatingMesh.point(vh)[1]
    floatingPositions[2,i] = floatingMesh.point(vh)[2]
    floatingNormals[0,i] = floatingMesh.normal(vh)[0]
    floatingNormals[1,i] = floatingMesh.normal(vh)[1]
    floatingNormals[2,i] = floatingMesh.normal(vh)[2]
    for j,vvh in enumerate(floatingMesh.vv(vh)):
        floatingNeighbours[j,i] = vvh.idx()
        #warning: make sure you limit this loop to number of maxNumNeighbours

for i, vh in enumerate(targetMesh.vertices()):
    targetPositions[0,i] = targetMesh.point(vh)[0]
    targetPositions[1,i] = targetMesh.point(vh)[1]
    targetPositions[2,i] = targetMesh.point(vh)[2]
    targetNormals[0,i] = targetMesh.normal(vh)[0]
    targetNormals[1,i] = targetMesh.normal(vh)[1]
    targetNormals[2,i] = targetMesh.normal(vh)[2]
    for j,vvh in enumerate(targetMesh.vv(vh)):
        targetNeighbours[j,i] = vvh.idx()
        #warning: make sure you limit this loop to number of maxNumNeighbours

##########
# ICP
##########
"""

"""
#for iteration in range(0,maxNumIterations):
##1) Determine Nearest neighbours. We'll simply use index correspondences for the bunny
correspondingPositions = targetPositions
##2) Determine weights. A weight related to the gaussian distance distribution suffices.
### Update the distribution parameters
sigmaa = 0.0
lambdaa = 0.0
sigmaNumerator = 0.0
sigmaDenominator = 0.0
for i in range(numFloatingVertices):
    distance = np.linalg.norm(correspondingPositions[:,i] - floatingPositions[:,i])
    sigmaNumerator = sigmaNumerator + floatingWeights[i] * distance * distance
    sigmaDenominator = sigmaDenominator + floatingWeights[i]

sigmaa = np.sqrt(sigmaNumerator/sigmaDenominator)
lambdaa = 1.0/(np.sqrt(2 * 3.14159) * sigmaa) * np.exp(-0.5 * kappaa * kappaa)
### Recalculate the weights
for i in range(numFloatingVertices):
    distance = np.linalg.norm(correspondingPositions[:,i] - floatingPositions[:,i])
    inlierProbability = 1.0/(np.sqrt(2 * 3.14159) * sigmaa) * np.exp(-0.5 * np.square(distance/sigmaa))
    #floatingWeights[i] = inlierProbability / (inlierProbability + lambdaa) #TODO: I left the weight calculation out for now

##3) Determine and update transformation.
###3.1 Compute the centroids of the Floating and Corresponding mesh
floatingCentroid = np.array([0.0,0.0,0.0], dtype = float)
correspondingCentroid = np.array([0.0,0.0,0.0], dtype = float)
weightSum = 0.0
for i in range(numFloatingVertices):
    floatingCentroid = floatingCentroid + floatingWeights[i] * floatingPositions[:,i]
    correspondingCentroid = correspondingCentroid + floatingWeights[i] * correspondingPositions[:,i]
    weightSum = weightSum + floatingWeights[i]

floatingCentroid = floatingCentroid / weightSum
correspondingCentroid = correspondingCentroid / weightSum
###3.2 Compute the Cross Variance matrix
crossVarianceMatrix = np.zeros((3,3), dtype = float)
for i in range(numFloatingVertices):
    #CrossVarMat = sum(weight[i] * floatingPosition[i] * correspondingPosition[i]_Transposed)
    crossVarianceMatrix = crossVarianceMatrix + floatingWeights[i] * np.outer(floatingPositions[:,i], correspondingPositions[:,i])

crossVarianceMatrix = crossVarianceMatrix / weightSum - np.outer(floatingCentroid, correspondingCentroid)
###3.3 Compute the Anti-Symmetric matrix
antiSymmetricMatrix = crossVarianceMatrix - crossVarianceMatrix.transpose()
###3.4 Use the cyclic elements of the Anti-Symmetric matrix to construct delta
delta = np.zeros(3)
delta[0] = antiSymmetricMatrix[1,2];
delta[1] = antiSymmetricMatrix[2,0];
delta[2] = antiSymmetricMatrix[0,1];
###3.5 Compute Q
Q = np.zeros((4,4), dtype = float)
Q[0,0] = np.trace(crossVarianceMatrix)
####take care with transposition of delta here. Depending on your library, it might be the opposite than this
Q[1:4,0] = delta
Q[0,1:4] = delta #in some libraries, you might have to transpose delta here
Q[1:4,1:4] = crossVarianceMatrix + crossVarianceMatrix.transpose() - np.identity(3, dtype = float)*np.trace(crossVarianceMatrix)
###3.6 we now compute the rotation quaternion by finding the eigenvector of
#Q of its largest eigenvalue
eigenValues, eigenVectors = np.linalg.eig(Q)
indexMaxVal = 0
maxVal = 0
for i, eigenValue in enumerate(eigenValues):
    if abs(eigenValue) > maxVal:
        maxVal = abs(eigenValue)
        indexMaxVal = i

rotQ = eigenVectors[:,indexMaxVal]
###3.7 Now we can construct the rotation matrix
rotMatTemp = np.zeros((3,3), dtype = float)
####diagonal elements
rotMatTemp[0,0] = np.square(rotQ[0]) + np.square(rotQ[1]) - np.square(rotQ[2]) - np.square(rotQ[3])
rotMatTemp[1,1] = np.square(rotQ[0]) + np.square(rotQ[2]) - np.square(rotQ[1]) - np.square(rotQ[3])
rotMatTemp[2,2] = np.square(rotQ[0]) + np.square(rotQ[3]) - np.square(rotQ[1]) - np.square(rotQ[2])
####remaining elements
rotMatTemp[1,0] = 2.0 * (rotQ[1] * rotQ[2] + rotQ[0] * rotQ[3])
rotMatTemp[0,1] = 2.0 * (rotQ[1] * rotQ[2] - rotQ[0] * rotQ[3])
rotMatTemp[2,0] = 2.0 * (rotQ[1] * rotQ[3] - rotQ[0] * rotQ[2])
rotMatTemp[0,2] = 2.0 * (rotQ[1] * rotQ[3] + rotQ[0] * rotQ[2])
rotMatTemp[2,1] = 2.0 * (rotQ[2] * rotQ[3] + rotQ[0] * rotQ[1])
rotMatTemp[1,2] = 2.0 * (rotQ[2] * rotQ[3] - rotQ[0] * rotQ[1])
### 3.8 Adjust the scale (optional)
scaleFactor = 1.0 #>1 to grow ; <1 to shrink
if adjustScale:
    numerator = 0.0
    denominator = 0.0
    for i in range(numFloatingVertices):
        centeredFloatingPosition = rotMatTemp.dot(floatingPositions[:,i] - floatingCentroid)
        centeredCorrespondingPosition = correspondingPositions[:,i] - correspondingCentroid
        numerator = numerator + floatingWeights[i] * np.dot(centeredCorrespondingPosition, centeredFloatingPosition)
        denominator = denominator + floatingWeights[i] * np.dot(centeredFloatingPosition, centeredFloatingPosition)
    scaleFactor = numerator / denominator

### 3.9 Compute the remaining translation necessary between the centroids
translation = correspondingCentroid - scaleFactor * rotMatTemp.dot(floatingCentroid)
### 3.10 Let's (finally) compute the entire transformation matrix!
translationMatrix = np.identity(4, dtype = float)
rotationMatrix = np.identity(4, dtype = float)
translationMatrix[0:3,3] = translation
rotationMatrix[0:3,0:3] = scaleFactor * rotMatTemp
#First rotate, then translate (transformations using homogeneous transformation
#matrices are executed from right to left).
transformationMatrix = translationMatrix.dot(rotationMatrix)

##4) Apply transformation
oldFloatingPositions = floatingPositions.copy()

floatingPosition = np.ones((4), dtype = float)
for i in range(numFloatingVertices):
    floatingPosition[0:3] = floatingPositions[:,i].copy()
    print(floatingPosition)
    floatingPositions[:,i] = transformationMatrix.dot(floatingPosition)[0:3]
    print(floatingPositions[:,i])

#print(oldFloatingPositions[100:110,:])
#print(floatingPositions[100:110,:])
#print(oldFloatingPositions[100:110] - floatingPositions[100:110])
##########
#EXPORT DATA
##########
##New vertex positions
newPosition = TriMesh.Point(0.0,0.0,0.0)
for i, vh in enumerate(floatingMesh.vertices()):
    newPosition[0] = floatingPositions[0,i]
    newPosition[1] = floatingPositions[1,i]
    newPosition[2] = floatingPositions[2,i]
    floatingMesh.set_point(vh,newPosition)

##New vertex normals
floatingMesh.request_vertex_normals()
floatingMesh.request_face_normals()
floatingMesh.update_normals()
floatingMesh.release_face_normals()

##Save the mesh
write_mesh(floatingMesh, outputMeshName)
