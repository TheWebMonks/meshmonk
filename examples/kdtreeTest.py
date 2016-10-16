#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
from openmesh import *
#Importing the rest of utilities
import numpy as np
from numpy import linalg
from scipy import spatial

"""
In this file, we're going to test the use of a kd-tree in python.
"""

#PARAMETERS
k = 3
leafsize = 1 #should be between 5 and 50 or so
eps = 0.01 #margin of error on distance allowed
p = 2 #Minkowski-norm: 2 for Euclidean
maxDist = 10
scaleNormals = 1.0 #by scaling the normals, their relative importance in the k-nn search can be de-/increased



#Build a simple 3D cubus
##Obtain info for initialization
numTargetVertices = 8
numFloatingVertices = 8
targetPositions = np.zeros((numTargetVertices,3), dtype = float)
floatingPositions = np.zeros((numFloatVertices,3), dtype = float)
targetNormals = np.zeros((numTargetVertices,3), dtype = float)
floatingNormals = np.zeros((numFloatVertices,3), dtype = float)
targetFeatures = np.zeros((numTargetVertices,6), dtype = float)
floatingFeatures = np.zeros((numFloatVertices,6), dtype = float)
##Put in the right positions
targetPositions[0,:] = [1.0,1.0,1.0]
targetPositions[1,:] = [1.0,-1.0,1.0]
targetPositions[2,:] = [-1.0,-1.0,1.0]
targetPositions[3,:] = [-1.0,1.0,1.0]
targetPositions[4,:] = [1.0,1.0,-1.0]
targetPositions[5,:] = [1.0,-1.0,-1.0]
targetPositions[6,:] = [-1.0,-1.0,-1.0]
targetPositions[7,:] = [-1.0,1.0,-1.0]
targetNormals[0,:] = [1.0,1.0,1.0]
targetNormals[1,:] = [1.0,-1.0,1.0]
targetNormals[2,:] = [-1.0,-1.0,1.0]
targetNormals[3,:] = [-1.0,1.0,1.0]
targetNormals[4,:] = [1.0,1.0,-1.0]
targetNormals[5,:] = [1.0,-1.0,-1.0]
targetNormals[6,:] = [-1.0,-1.0,-1.0]
targetNormals[7,:] = [-1.0,1.0,-1.0]
norms = np.linalg.norm(targetNormals, axis=-1)[:, np.newaxis] #normalize normals
targetNormals = targetNormals / norms
targetFeatures = np.hstack((targetPositions,scaleNormals*targetNormals))
##Set up inlier probabilities and vertex statuses.
"""
Difference between inlier probability and vertex status:
-inlier prob is a measure for how likely a vertex is an inlier. The values are
in a continuous range between 0.0 and 1.0. The value can be determined
internally during the registration process, in various ways.
-vertex status: A vertex can be flagged externally as 1 or 0 (binary). This is
not adjusted during registration, but is used to prevent these vertices from
influencing computations of correspondences, outliers and the transformation.

The status should be used for example to pre-screen for erroneous data. This
flagging does not change in the registration process.
The inlier probabilities can take the status info into account, but they are not
the same thing.

E.g. when a target vertex is marked 'status 0', every floating vertex that is
attracted to this target vertex is automatically marked an outlier itself.
However, being attracted to nodes with 'inlier weight 0.0' does not make a
vertex an outlier.
"""
targetInlierProb = np.ones((numTargetVertices), dtype = float)
targetStatus = np.ones((numTargetVertices), dtype = int)
floatingInlierProb = np.ones((numTargetVertices), dtype = float)
floatingStatus = np.ones((numTargetVertices), dtype = int)

##We'll determine our floating object by rotating the cubus around the z-axis
rotationMatrix = np.array([[0.955336489125606, 0.29552020666134, 0.0],[-0.29552020666134, 0.955336489125606, 0.0], [0.0, 0.0, 1.0]])
floatingPositions = targetPositions.dot(rotationMatrix.transpose())
floatingNormals = targetNormals.dot(rotationMatrix.transpose())
floatingFeatures = np.hstack((floatingPositions,scaleNormals*floatingNormals))


#construct the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.cKDTree.html#scipy.spatial.cKDTree)
kdTreeTarget = spatial.cKDTree(targetFeatures, leafsize)
kdTreeFloat = spatial.cKDTree(floatingFeatures, leafsize)

#query the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query)
distFloatToTarget, idxFloatToTarget = kdTreeTarget.query(floatingFeatures, k, eps, p, maxDist)
distTargetToFloat, idxTargetToFloat = kdTreeFloat.query(targetFeatures, k, eps, p, maxDist)

#Construct the affinity matrices
affinityFloatToTarget = np.zeros((numFloatingVertices, numTargetVertices), dtype = float)
affinityTargetToFloat = np.zeros((numTargetVertices, numFloatingVertices), dtype = float)

#loop over the floating vertices to determine their affinity with the target ones
for i in range(numFloatingVertices):
    ##Loop over k nearest neighbours found in the target vertices
    for j, idx in enumerate(idxFloatToTarget[i,:]):
        #idx now indicates the index of the target vertex that was determined a
        #neighbour of the current floating vertex
        distance = distFloatToTarget[i,j]
        if distance < 0.001:
            distance = 0.001 #numeric stability
        weight = 1.0/(distance * distance)
        #this weight should be inserted in the affinity matrix element that links
        #floating point at index i with target point at index 'idx'
        affinityFloatToTarget[i,idx] = weight
        #finally, check the status of the neighbour target vertex and set the
        #inlier probability of this floating vertex to zero if status was zero.
        if targetStatus[idx] == 0:
            floatingInlierProb[i] = 0

#loop over the target vertices to determine their affinity with the floating ones
for i in range(numTargetVertices):
    ##Loop over k nearest neighbours found in the floating vertices
    for j, idx in enumerate(idxTargetToFloat[i,:]):
        #idx now indicates the index of the target vertex that was determined a
        #neighbour of the current floating vertex
        distance = distTargetToFloat[i,j]
        if distance < 0.001:
            distance = 0.001 #numeric stability
        weight = 1.0/(distance * distance)
        #this weight should be inserted in the affinity matrix element that links
        #floating point at index i with target point at index 'idx'
        affinityTargetToFloat[i,idx] = weight
        #finally, check the status of the neighbour floating vertex and set the
        #inlier probability of this target vertex to zero if status was zero.
        if floatingStatus[idx] == 0:
            targetInlierProb[i] = 0

#normalize both affinity matrices (the sum of each row has to equal 1.0)
rowSums = affinityFloatToTarget.sum(axis=1, keepdims=True)
affinityFloatToTarget = affinityFloatToTarget / rowSums
rowSums = affinityTargetToFloat.sum(axis=1, keepdims=True)
affinityTargetToFloat = affinityTargetToFloat / rowSums

#Sum the affinity matrices and normalize again
affinityFloatToTarget = affinityFloatToTarget + affinityTargetToFloat
rowSums = affinityFloatToTarget.sum(axis=1, keepdims=True)
affinityFloatToTarget = affinityFloatToTarget / rowSums

"""
Now we have constructed the affinity matrix using push-pull forces, let's
construct the corresponding positions.
"""
correspondingFeatures = affinityFloatToTarget.dot(targetFeatures)
#extract positions and normals (and normalize the normals)
(correspondingPositions,correspondingNormals) = np.hsplit(correspondingFeatures,2)
norms = np.linalg.norm(correspondingNormals, axis=-1)[:, np.newaxis] #normalize normals
correspondingNormals = correspondingNormals / norms
