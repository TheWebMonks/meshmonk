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



#Build a simple 3D cubus
##Obtain info for initialization
targetPositions = np.zeros((numTargetVertices,3), dtype = float)
floatPositions = np.zeros((numFloatVertices,3), dtype = float)
targetWeights = np.ones((numTargetVertices,3), dtype = float)
floatWeights = np.ones((numTargetVertices,3), dtype = float)
##Put in the right positions
targetPositions[0,:] = [1.0,1.0,1.0]
targetPositions[1,:] = [1.0,-1.0,1.0]
targetPositions[2,:] = [-1.0,-1.0,1.0]
targetPositions[3,:] = [-1.0,1.0,1.0]
targetPositions[4,:] = [1.0,1.0,-1.0]
targetPositions[5,:] = [1.0,-1.0,-1.0]
targetPositions[6,:] = [-1.0,-1.0,-1.0]
targetPositions[7,:] = [-1.0,1.0,-1.0]
numTargetVertices = targetPositions.shape[0]
numFloatingVertices = numTargetVertices

##We'll determine our floating object by rotating the cubus around the z-axis
rotationMatrix = np.array([[0.955336489125606, 0.29552020666134, 0.0],[-0.29552020666134, 0.955336489125606, 0.0], [0.0, 0.0, 1.0]])
floatPositions = targetPositions.dot(rotationMatrix.transpose())


#construct the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.cKDTree.html#scipy.spatial.cKDTree)
kdTreeTarget = spatial.cKDTree(targetPositions, leafsize)
kdTreeFloat = spatial.cKDTree(floatPositions, leafsize)

#query the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query)
distFloatToTarget, idxFloatToTarget = kdTreeTarget.query(floatPositions, k, eps, p, maxDist)
distTargetToFloat, idxTargetToFloat = kdTreeFloat.query(targetPositions, k, eps, p, maxDist)

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

#normalize both affinity matrices (the sum of each row has to equal 1.0)
rowSums = affinityFloatToTarget.sum(axis=1, keepdims=True)
affinityFloatToTarget = affinityFloatToTarget / rowSums
rowSums = affinityTargetToFloat.sum(axis=1, keepdims=True)
affinityTargetToFloat = affinityTargetToFloat / rowSums

#Sum the affinity matrices and normalize again
affinityFloatToTarget = affinityFloatToTarget + affinityTargetToFloat
rowSums = affinityFloatToTarget.sum(axis=1, keepdims=True)
affinityFloatToTarget = affinityFloatToTarget / rowSums
