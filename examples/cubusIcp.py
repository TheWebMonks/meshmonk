#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
from openmesh import *
#Importing the rest of utilities
import numpy as np
from scipy import linalg

from context import registration

"""
In this example, we will perform rigid registration of a self-made object,
which are basically one pyramid tied to another pyramid upside down.

#####
#Some definitions and info:
#####

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

# Obtain info for initialization
numFloatingVertices = 8
numFloatingFaces = 8

# Initialize matrices
floatingPositions = np.zeros((numFloatingVertices,3), dtype = float)
correspondingPositions = np.zeros((numFloatingVertices,3), dtype = float)
targetPositions = np.zeros((numFloatingVertices,3), dtype = float)
floatingNormals = np.zeros((numFloatingVertices,3), dtype = float)
correspondingNormals = np.zeros((numFloatingVertices,3), dtype = float)
targetNormals = np.zeros((numFloatingVertices,3), dtype = float)
floatingFeatures = np.zeros((numFloatingVertices,6), dtype = float)
correspondingFeatures = np.zeros((numFloatingVertices,6), dtype = float)
targetFeatures = np.zeros((numFloatingVertices,6), dtype = float)
# Initialize weights
floatingWeights = np.ones((numFloatingVertices), dtype = float)

# Put in the right positions and normals
floatingPositions[0,:] = [1.0,1.0,1.0]
floatingPositions[1,:] = [1.0,-1.0,1.0]
floatingPositions[2,:] = [-1.0,-1.0,1.0]
floatingPositions[3,:] = [-1.0,1.0,1.0]
floatingPositions[4,:] = [1.0,1.0,-1.0]
floatingPositions[5,:] = [1.0,-1.0,-1.0]
floatingPositions[6,:] = [-1.0,-1.0,-1.0]
floatingPositions[7,:] = [-1.0,1.0,-1.0]
## The normals are pointing in the same direction as the position vectors, but
## they have to be normalized to be unit vectors.
for i in range(numFloatingVertices):
    length = linalg.norm(floatingPositions[i,:])
    floatingNormals[i,:] = floatingPositions[i,:] / length

# The features are the concatenation of positions and normals
floatingFeatures = np.hstack((floatingPositions,floatingNormals))

"""
The original transformation which we will apply to the floating data to obtain
our target data. We will try to recover this matrix through ICP.
"""
originalTransformation = np.array([[-1.0/np.sqrt(6.0), -1.0/np.sqrt(6.0), 2.0/np.sqrt(6.0), 0.0],[1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0), 0.0, 0.0],[1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0), 0.0],[0.0, 0.0, 0.0, 1.0]])

position = np.array([0.0,0.0,0.0,1.0])
normal = np.array([0.0,0.0,0.0,1.0])
for i in range(numFloatingVertices):
    position[0:3] = floatingFeatures[i,0:3]
    normal[0:3] = floatingFeatures[i,3:6]
    targetPositions[i,:] = originalTransformation.dot(position)[0:3]
    targetNormals[i,:] = originalTransformation.dot(normal)[0:3]


##########
# ICP
##########
"""

"""
#for iteration in range(0,maxNumIterations):
##1) Determine Nearest neighbours. We'll simply use index correspondences for the bunny
registration.core.wknn_affinity(floatingPositions, features2, affinity12, k = 3)
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
transformationMatrix = rotationMatrix.dot(translationMatrix)

##4) Apply transformation
oldFloatingPositions = floatingPositions.copy()

floatingPosition = np.ones((4), dtype = float)
for i in range(numFloatingVertices):
    floatingPosition[0:3] = floatingPositions[:,i].copy()
    print(floatingPosition)
    floatingPositions[:,i] = transformationMatrix.dot(floatingPosition)[0:3]
    print(floatingPositions[:,i])
