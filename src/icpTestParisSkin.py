#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
#if problems with libstdc++ : http://askubuntu.com/questions/575505/glibcxx-3-4-20-not-found-how-to-fix-this-error
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
from openmesh import *
#Importing the rest of utilities
import numpy as np
from numpy import linalg
from scipy import spatial

"""
In this file, we're going to test the combination of our first ICP version
together with the use of weighted k-NN (wknn) correspondences on (downsampled)
real life data.
"""
##########
#PARAMETERS
##########
# DATA
floatingMeshPath = "/home/jonatan/kuleuven-algorithms/src/data/paris/558/skin95.obj"
targetMeshPath = "/home/jonatan/kuleuven-algorithms/src/data/paris/559/skin95.obj"
outputMeshName = "/home/jonatan/kuleuven-algorithms/src/data/paris/results/icp_558_to_559.obj"
# CORRESPONDENCES
k = 3
leafsize = 15 #should be between 5 and 50 or so
eps = 0.0001 #margin of error on distance allowed
p = 2 #Minkowski-norm: 2 for Euclidean
maxDist = 1000
scaleNormals = 1.0 #by scaling the normals, their relative importance in the k-nn search can be de-/increased
# OUTLIER DETECTION
kappaa = 3
# ICP PROCESS
maxNumIterations = 10
adjustScale = False
"""
################################################################################
##########
#DATA PREPARATION
##########
################################################################################
"""
#INITIALIZATION
##Load meshes
floatingMesh = TriMesh()
targetMesh = TriMesh()
correspondingMesh = TriMesh()
read_mesh(floatingMesh, floatingMeshPath)
read_mesh(targetMesh, targetMeshPath)
read_mesh(correspondingMesh, floatingMeshPath) #copy of floating mesh

##make sure the vertex normals are there
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

##Obtain some parameters for initialization
numFloatingVertices = floatingMesh.n_vertices()
numFloatingFaces = floatingMesh.n_faces()
numTargetVertices = targetMesh.n_vertices()
numTargetFaces = targetMesh.n_faces()

##Initialize positions, normals and feature matrices
targetPositions = np.zeros((numTargetVertices,3), dtype = float)
floatingPositions = np.zeros((numFloatingVertices,3), dtype = float)
correspondingPositions = np.zeros((numFloatingVertices,3), dtype = float)
targetNormals = np.zeros((numTargetVertices,3), dtype = float)
floatingNormals = np.zeros((numFloatingVertices,3), dtype = float)
correspondingNormals = np.zeros((numFloatingVertices,3), dtype = float)
targetFeatures = np.zeros((numTargetVertices,6), dtype = float)
floatingFeatures = np.zeros((numFloatingVertices,6), dtype = float)
correspondingFeatures = np.zeros((numFloatingVertices,6), dtype = float)
## Initialize inlier probabilities and statuses
targetInlierProb = np.ones((numTargetVertices), dtype = float)
floatingInlierProb = np.ones((numFloatingVertices), dtype = float)
targetStatus = np.ones((numTargetVertices), dtype = int)
floatingStatus = np.ones((numFloatingVertices), dtype = int)

##Extract data from meshes
###floating mesh positions and normals
for i, vh in enumerate(floatingMesh.vertices()):
    floatingPositions[i,0] = floatingMesh.point(vh)[0]
    floatingPositions[i,1] = floatingMesh.point(vh)[1]
    floatingPositions[i,2] = floatingMesh.point(vh)[2]
    floatingNormals[i,0] = floatingMesh.normal(vh)[0]
    floatingNormals[i,1] = floatingMesh.normal(vh)[1]
    floatingNormals[i,2] = floatingMesh.normal(vh)[2]

###concatenate the two to build the floating mesh feature matrix
floatingFeatures = np.hstack((floatingPositions,scaleNormals * floatingNormals))

###target mesh positions and normals
for i, vh in enumerate(targetMesh.vertices()):
    targetPositions[i,0] = targetMesh.point(vh)[0]
    targetPositions[i,1] = targetMesh.point(vh)[1]
    targetPositions[i,2] = targetMesh.point(vh)[2]
    targetNormals[i,0] = targetMesh.normal(vh)[0]
    targetNormals[i,1] = targetMesh.normal(vh)[1]
    targetNormals[i,2] = targetMesh.normal(vh)[2]

###concatenate the two to build the target mesh feature matrix
targetFeatures = np.hstack((targetPositions,scaleNormals * targetNormals))


"""
################################################################################
##########
# ICP
##########
################################################################################
"""

for iteration in range(0,maxNumIterations): #TODO: un-comment later
    """
    ##########
    #CORRESPONDENCES
    ##########
    """
    #1) Determine correspondences as weighted nearest neighbours
    ##construct the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.cKDTree.html#scipy.spatial.cKDTree)
    kdTreeTarget = spatial.cKDTree(targetFeatures, leafsize)
    kdTreeFloating = spatial.cKDTree(floatingFeatures, leafsize)
    ##query the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query)
    distFloatingToTarget, idxFloatingToTarget = kdTreeTarget.query(floatingFeatures, k, eps, p, maxDist)
    distTargetToFloating, idxTargetToFloating = kdTreeFloating.query(targetFeatures, k, eps, p, maxDist)
    ##Construct the affinity matrices
    ###Initialize the matrices
    affinityFloatingToTarget = np.zeros((numFloatingVertices, numTargetVertices), dtype = float)
    affinityTargetToFloating = np.zeros((numTargetVertices, numFloatingVertices), dtype = float)
    ###loop over the floating vertices to determine their affinity with the target ones
    for i in range(numFloatingVertices):
        ##Loop over k nearest neighbours found in the target vertices
        for j, idx in enumerate(idxFloatingToTarget[i,:]):
            #idx now indicates the index of the target vertex that was determined a
            #neighbour of the current floating vertex
            distance = distFloatingToTarget[i,j]
            if distance < 0.001:
                distance = 0.001 #numeric stability
            weight = 1.0/(distance * distance)
            #this weight should be inserted in the affinity matrix element that links
            #floating point at index i with target point at index 'idx'
            affinityFloatingToTarget[i,idx] = weight
            #finally, check the status of the neighbour target vertex and set the
            #inlier probability of this floating vertex to zero if status was zero.
            if targetStatus[idx] == 0:
                floatingInlierProb[i] = 0

    ##loop over the target vertices to determine their affinity with the floating ones
    for i in range(numTargetVertices):
        ##Loop over k nearest neighbours found in the floating vertices
        for j, idx in enumerate(idxTargetToFloating[i,:]):
            #idx now indicates the index of the target vertex that was determined a
            #neighbour of the current floating vertex
            distance = distTargetToFloating[i,j]
            if distance < 0.001:
                distance = 0.001 #numeric stability
            weight = 1.0/(distance * distance)
            #this weight should be inserted in the affinity matrix element that links
            #floating point at index i with target point at index 'idx'
            affinityTargetToFloating[i,idx] = weight
            #finally, check the status of the neighbour floating vertex and set the
            #inlier probability of this target vertex to zero if status was zero.
            if floatingStatus[idx] == 0:
                targetInlierProb[i] = 0

    ##normalize both affinity matrices (the sum of each row has to equal 1.0)
    rowSums = affinityFloatingToTarget.sum(axis=1, keepdims=True)
    affinityFloatingToTarget = affinityFloatingToTarget / rowSums
    rowSums = affinityTargetToFloating.sum(axis=1, keepdims=True)
    affinityTargetToFloating = affinityTargetToFloating / rowSums
    ##Sum the affinity matrices and normalize again
    affinityFloatingToTarget = affinityFloatingToTarget + affinityTargetToFloating.transpose()
    rowSums = affinityFloatingToTarget.sum(axis=1, keepdims=True)
    affinityFloatingToTarget = affinityFloatingToTarget / rowSums
    ##Now we have constructed the affinity matrix using push-pull forces, let's
    ##construct the corresponding positions.
    correspondingFeatures = affinityFloatingToTarget.dot(targetFeatures)
    ##extract positions and normals (and normalize the normals)
    (correspondingPositions,correspondingNormals) = np.hsplit(correspondingFeatures,2)
    norms = np.linalg.norm(correspondingNormals, axis=-1)[:, np.newaxis] #normalize normals
    correspondingNormals = correspondingNormals / norms
    """
    IMPORTANT NOTE: For the ease of implementing the outlier detection and
    transformation in the ICP algorithms, we transpose our our matrices here so
    the each column corresponds to a single vertex and the rows to a different
    co-ordinate or feature.
    """
    floatingPositions = floatingPositions.transpose()
    targetPositions = targetPositions.transpose()
    correspondingPositions = correspondingPositions.transpose()
    floatingNormals = floatingNormals.transpose()
    targetNormals = targetNormals.transpose()
    correspondingNormals = correspondingNormals.transpose()
    floatingFeatures = floatingFeatures.transpose()
    targetFeatures = targetFeatures.transpose()
    correspondingFeatures = correspondingFeatures.transpose()
    """
    ##########
    #OUTLIER DETECTION
    ##########
    """
    #2) Determine weights. A weight related to the gaussian distance distribution suffices.
    ## Update the distribution parameters
    sigmaa = 0.0
    lambdaa = 0.0
    sigmaNumerator = 0.0
    sigmaDenominator = 0.0
    for i in range(numFloatingVertices):
        distance = np.linalg.norm(correspondingFeatures[:,i] - floatingFeatures[:,i])
        sigmaNumerator = sigmaNumerator + floatingInlierProb[i] * distance * distance
        sigmaDenominator = sigmaDenominator + floatingInlierProb[i]

    sigmaa = np.sqrt(sigmaNumerator/sigmaDenominator)
    lambdaa = 1.0/(np.sqrt(2 * 3.14159) * sigmaa) * np.exp(-0.5 * kappaa * kappaa)

    ### Recalculate the weights
    for i in range(numFloatingVertices):
        distance = np.linalg.norm(correspondingFeatures[:,i] - floatingFeatures[:,i])
        inlierProbability = 1.0/(np.sqrt(2 * 3.14159) * sigmaa) * np.exp(-0.5 * np.square(distance/sigmaa))
        floatingInlierProb[i] = inlierProbability / (inlierProbability + lambdaa) #TODO: I left the weight calculation out for now
    """
    ##########
    #TRANSFORMATION
    ##########
    """
    ##3) Determine and update transformation.
    ###3.1 Compute the centroids of the Floating and Corresponding mesh
    floatingCentroid = np.array([0.0,0.0,0.0], dtype = float)
    correspondingCentroid = np.array([0.0,0.0,0.0], dtype = float)
    weightSum = 0.0
    for i in range(numFloatingVertices):
        floatingCentroid = floatingCentroid + floatingInlierProb[i] * floatingPositions[:,i]
        correspondingCentroid = correspondingCentroid + floatingInlierProb[i] * correspondingPositions[:,i]
        weightSum = weightSum + floatingInlierProb[i]

    floatingCentroid = floatingCentroid / weightSum
    correspondingCentroid = correspondingCentroid / weightSum
    ###3.2 Compute the Cross Variance matrix
    crossVarianceMatrix = np.zeros((3,3), dtype = float)

    for i in range(numFloatingVertices):
        #CrossVarMat = sum(weight[i] * floatingPosition[i] * correspondingPosition[i]_Transposed)
        crossVarianceMatrix = crossVarianceMatrix + floatingInlierProb[i] * np.outer(floatingPositions[:,i], correspondingPositions[:,i])
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
    rotationMatrix = np.zeros((3,3), dtype = float)
    ####diagonal elements
    rotationMatrix[0,0] = np.square(rotQ[0]) + np.square(rotQ[1]) - np.square(rotQ[2]) - np.square(rotQ[3])
    rotationMatrix[1,1] = np.square(rotQ[0]) + np.square(rotQ[2]) - np.square(rotQ[1]) - np.square(rotQ[3])
    rotationMatrix[2,2] = np.square(rotQ[0]) + np.square(rotQ[3]) - np.square(rotQ[1]) - np.square(rotQ[2])
    ####remaining elements
    rotationMatrix[1,0] = 2.0 * (rotQ[1] * rotQ[2] + rotQ[0] * rotQ[3])
    rotationMatrix[0,1] = 2.0 * (rotQ[1] * rotQ[2] - rotQ[0] * rotQ[3])
    rotationMatrix[2,0] = 2.0 * (rotQ[1] * rotQ[3] - rotQ[0] * rotQ[2])
    rotationMatrix[0,2] = 2.0 * (rotQ[1] * rotQ[3] + rotQ[0] * rotQ[2])
    rotationMatrix[2,1] = 2.0 * (rotQ[2] * rotQ[3] + rotQ[0] * rotQ[1])
    rotationMatrix[1,2] = 2.0 * (rotQ[2] * rotQ[3] - rotQ[0] * rotQ[1])
    ### 3.8 Adjust the scale (optional)
    scaleFactor = 1.0 #>1 to grow ; <1 to shrink
    if adjustScale:
        numerator = 0.0
        denominator = 0.0
        for i in range(numFloatingVertices):
            centeredFloatingPosition = rotationMatrix.dot(floatingPositions[:,i] - floatingCentroid)
            centeredCorrespondingPosition = correspondingPositions[:,i] - correspondingCentroid
            numerator = numerator + floatingInlierProb[i] * np.dot(centeredCorrespondingPosition, centeredFloatingPosition)
            denominator = denominator + floatingInlierProb[i] * np.dot(centeredFloatingPosition, centeredFloatingPosition)
        scaleFactor = numerator / denominator
    ### 3.9 Compute the remaining translation necessary between the centroids
    translation = correspondingCentroid - scaleFactor * rotationMatrix.dot(floatingCentroid)
    ### 3.10 Let's (finally) compute the entire transformation matrix!
    translationMatrix = np.identity(4, dtype = float)
    scaledRotationMatrix = np.identity(4, dtype = float)
    translationMatrix[0:3,3] = translation
    scaledRotationMatrix[0:3,0:3] = scaleFactor * rotationMatrix
    ####First rotate, then translate (transformations using homogeneous
    ####transformation matrices are executed from right to left).
    transformationMatrix = translationMatrix.dot(scaledRotationMatrix)
    ##4) Apply transformation to update the floating surface
    ###4.1 Apply total transformation to the positions
    position = np.ones((4), dtype = float)
    for i in range(numFloatingVertices):
        position[0:3] = floatingPositions[:,i].copy()
        floatingPositions[:,i] = transformationMatrix.dot(position)[0:3]
    ###4.2 Apply the rotation to the normals
    normal = np.ones((3), dtype = float)
    for i in range(numFloatingVertices):
        normal = floatingNormals[:,i].copy()
        floatingNormals[:,i] = rotationMatrix.dot(normal)
    ####Normalize the normals again (avoiding accumulation of rounding errors)
    norms = np.linalg.norm(floatingNormals, axis=-1)[:, np.newaxis] #normalize normals
    floatingNormals = floatingNormals / norms
    ###4.3 Update floating feature matrix
    floatingFeatures[0:3,:] = floatingPositions
    floatingFeatures[3:6,:] = scaleNormals * floatingNormals
    """
    IMPORTANT NOTE: We have to reverse the transposition that we did earlier. We did
    that to make the implementation of the outlier detection and transformation in
    the ICP algorithm easier, but now we have to undo that.
    """
    floatingPositions = floatingPositions.transpose()
    targetPositions = targetPositions.transpose()
    correspondingPositions = correspondingPositions.transpose()
    floatingNormals = floatingNormals.transpose()
    targetNormals = targetNormals.transpose()
    correspondingNormals = correspondingNormals.transpose()
    floatingFeatures = floatingFeatures.transpose()
    targetFeatures = targetFeatures.transpose()
    correspondingFeatures = correspondingFeatures.transpose()

"""
################################################################################
##########
#EXPORT DATA
##########
################################################################################
"""
##New vertex positions
newPosition = TriMesh.Point(0.0,0.0,0.0)
for i, vh in enumerate(floatingMesh.vertices()):
    newPosition[0] = floatingPositions[i,0]
    newPosition[1] = floatingPositions[i,1]
    newPosition[2] = floatingPositions[i,2]
    floatingMesh.set_point(vh,newPosition)

##New vertex normals
floatingMesh.request_vertex_normals()
floatingMesh.request_face_normals()
floatingMesh.update_normals()
floatingMesh.release_face_normals()

##Save the mesh
write_mesh(floatingMesh, outputMeshName)
