"""
   Copyright 2016 Jonatan Snyders

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
import openmesh
#Importing the rest of utilities
import numpy
from numpy import linalg
from scipy import spatial







def nearest_neighbours(features1, features2, k = 3, leafsize = 15, eps = 0.0001, p = 2, maxDistance = 1000):
    """
    GOAL
    This function searches for the k nearest neighbours in the features2-set for
    each element in the features1 set. It outputs the indices of each neighbour
    and the distances between each element of features1 and its neighbours.

    INPUT
    -features1:
    -features2:

    PARAMETERS
    -k(= 3): number of nearest neighbours
    -leafsize(= 15): should be between 5 and 50 or so
    -eps(=0.0001): margin of error on distance allowed
    -p(=2): Minkowski-norm: 2 for Euclidean
    -maxDistance(=1000)

    RETURNS
    -distances
    -neighbourIndices.
    """
    # Obtain info and initialize outputs
    numElements1 = features1.shape[0]
    kdTree = spatial.cKDTree(features2, leafsize)
    # Query the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query)
    distances, neighbourIndices = kdTree.query(features1, k, eps, p, maxDistance)

    return distances, neighbourIndices

def gaussian_vector_interpolation(position, fieldPositions, fieldVectors, fieldWeights, sigmaa):
    """
    GOAL
    This function does 'Gaussian interpolation' of vectors in an irregular
    vector field to estimate the vector of a certain position.

    INPUT
    -position:
    The postion where you want to calculate the vector value at
    -fieldPositions:
    Positions of the vectors of the vector field. The dimensions of this matrix
    are numVectors x 3
    -fieldVectors:
    Values of the vectors of the vector field. The dimensions of this matrix
    are numVectors x numFeatures
    -fieldWeights:
    The contribution of each vector field vector can be weighed.

    PARAMETERS
    -sigmaa:
    The sigma used in the Gaussian averaging.

    RETURNS
    -interpolatedVector:
    Computed through Gaussian averaging, taking additional weights (see
    fieldWeights) into account.
    """
    # Obtain info and initialize outputs
    numVectors = fieldVectors.shape[0]
    numFeatures = fieldVectors.shape[1]
    interpolatedVector = numpy.zeros((numFeatures), dtype = float)

    # Loop over the field vectors
    sumWeights = 0.0
    for i, vector in enumerate(fieldVectors):
        ## Compute distance
        distanceVector = fieldPositions[i] - position
        distance = linalg.norm(distanceVector)
        ## From the distance, we can compute the gaussian weight
        gaussianWeight = numpy.exp(-0.5 * numpy.square(distance / sigmaa))
        ## Let's take the user-defined weight into account
        combinedWeight = gaussianWeight * fieldWeights[i]
        ## Weigh the current vector field vector and sum it
        interpolatedVector = interpolatedVector + combinedWeight * vector
        ## Add the weight to the total sum of weights for normalization after.
        sumWeights = sumWeights + combinedWeight

    # Normalize the vector
    interpolatedVector = interpolatedVector / sumWeights

    return interpolatedVector


def gaussian_scalar_interpolation(position, fieldPositions, fieldScalars, fieldWeights, sigmaa):
    """
    GOAL
    This function does 'Gaussian interpolation' of scalars in an irregular
    scalar field to estimate the scalar of a certain position.

    INPUT
    -position:
    The postion where you want to calculate the scalar value at
    -fieldPositions:
    Positions of the nodes of the scalar field. The dimensions of this matrix
    are numNodes x 3
    -fieldScalars:
    Scalar values of the nodes of the scalar field. The dimensions of this array
    are numNodes x 1
    -fieldWeights:
    The contribution of each scalar field node can be weighed.

    PARAMETERS
    -sigmaa:
    The sigma used in the Gaussian averaging.

    RETURNS
    -interpolatedScalar:
    Computed through Gaussian averaging, taking additional weights (see
    fieldWeights) into account.
    """
    # Obtain info and initialize outputs
    numNodes = fieldScalars.shape[0]
    interpolatedScalar = 0.0

    # Loop over the field nodes
    sumWeights = 0.0
    for i, scalar in enumerate(fieldScalars):
        ## Compute distance
        distanceVector = fieldPositions[i] - position
        distance = linalg.norm(distanceVector)
        ## From the distance, we can compute the gaussian weight
        gaussianWeight = numpy.exp(-0.5 * numpy.square(distance / sigmaa))
        ## Let's take the user-defined weight into account
        combinedWeight = gaussianWeight * fieldWeights[i]
        ## Weigh the current node's scalar and sum it
        interpolatedScalar = interpolatedScalar + combinedWeight * scalar
        ## Add the weight to the total sum of weights for normalization after.
        sumWeights = sumWeights + combinedWeight

    # Normalize the vector
    interpolatedScalar = interpolatedScalar / sumWeights

    return interpolatedScalar

def gaussian_smoothing_displacement_field(floatingPositions, unregulatedDisplacementField, regulatedDisplacementField, floatingWeights, numNeighbours, gaussianSigma):
    """
    GOAL
    This function performs gaussian smoothing on a displacement field.

    INPUT
    -floatingPositions
    -unregulatedDisplacementField:
    The displacement field that should be regulated/smoothed.
    -floatingWeights

    OUTPUT
    -regulatedDisplacementField:
    The resulting regulated displacement field.


    PARAMETERS
    -numNeighbours:
    For the smoothing, the nearest neighbours for each floating positions have
    to be found. The number should be high enough so that every significant
    contribution (up to a distance of e.g. 3*gaussianSigma) is included. But low
    enough to keep computational speed high.
    -gaussianSigma:
    The value for sigma of the gaussian used for the smoothing.

    RETURNS
    """
    # Info & Initialization
    numFloatingVertices = floatingPositions.shape[0]
    # Determine for each displacement the (closely) neighbouring displacements
    distances, neighbourIndices = nearest_neighbours(floatingPositions, floatingPositions, numNeighbours, 15, 0.0001, 2, 1000.0)
    # Use the neighbouring displacements to smooth each individual displacement
    for i in range(numFloatingVertices):
        position = floatingPositions[i,:]
        neighbourPositions = numpy.zeros((numNeighbours,3), dtype = float)
        neighbourDisplacements = numpy.zeros((numNeighbours,3), dtype = float)
        neighbourWeights = numpy.zeros((numNeighbours), dtype = float)
        ## For the current displacement, get all the neighbouring displacements,
        ## positions, and weights (needed for Gaussian smoothing!).
        for j in range(numNeighbours):
            neighbourIndex = neighbourIndices[i,j]
            neighbourPositions[j,:] = floatingPositions[neighbourIndex,:]
            neighbourDisplacements[j,:] = unregulatedDisplacementField[neighbourIndex,:]
            neighbourWeights[j] = floatingWeights[neighbourIndex]

        ## Gaussian averaging of neighbouring displacements
        regulatedDisplacement = gaussian_vector_interpolation(position, neighbourPositions, neighbourDisplacements, neighbourWeights, gaussianSigma)
        regulatedDisplacementField[i,:] = regulatedDisplacement


def openmesh_normals_from_positions(mesh, newPositions):
    """
    GOAL
    This function recalculates the normals for a given mesh, when there are new
    positions for the nodes.

    INPUT
    -mesh:
    this has to be a mesh of openmesh's TriMesh type
    -newPositions:
    the new positions that have to be inserted into the mesh to recompute the
    vertex normals.

    PARAMETERS

    RETURNS
    -newNormals:
    The new vertex normals that were recomputed given the new positions and mesh
    structure.
    """
    # Info and Initialization
    numVertices = newPositions.shape[0]
    newNormals = numpy.zeros((numVertices,3), dtype = float)
    # Put positions back into the floating mesh
    newPosition = openmesh.TriMesh.Point(0.0,0.0,0.0)
    for i, vh in enumerate(mesh.vertices()):
        newPosition[0] = newPositions[i,0]
        newPosition[1] = newPositions[i,1]
        newPosition[2] = newPositions[i,2]
        mesh.set_point(vh,newPosition)

    # Let openmesh recalculate the vertex normals
    mesh.request_vertex_normals()
    mesh.request_face_normals()
    mesh.update_normals()
    mesh.release_face_normals()

    # Get the normals back out from the openmesh TriMesh
    for i, vh in enumerate(mesh.vertices()):
        newNormals[i,0] = mesh.normal(vh)[0]
        newNormals[i,1] = mesh.normal(vh)[1]
        newNormals[i,2] = mesh.normal(vh)[2]

    return newNormals
