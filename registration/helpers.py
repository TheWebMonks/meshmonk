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
from openmesh import *
#Importing the rest of utilities
import numpy as np
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
    interpolatedVector = np.zeros((numFeatures), dtype = float)

    # Loop over the field vectors
    sumWeights = 0.0
    for i, vector in enumerate(fieldVectors):
        ## Compute distance
        distanceVector = fieldPositions[i] - position
        distance = linalg.norm(distanceVector)
        ## From the distance, we can compute the gaussian weight
        gaussianWeight = np.exp(-0.5 * np.square(distance / sigmaa))
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
        gaussianWeight = np.exp(-0.5 * np.square(distance / sigmaa))
        ## Let's take the user-defined weight into account
        combinedWeight = gaussianWeight * fieldWeights[i]
        ## Weigh the current node's scalar and sum it
        interpolatedScalar = interpolatedScalar + combinedWeight * scalar
        ## Add the weight to the total sum of weights for normalization after.
        sumWeights = sumWeights + combinedWeight

    # Normalize the vector
    interpolatedScalar = interpolatedScalar / sumWeights

    return interpolatedScalar
