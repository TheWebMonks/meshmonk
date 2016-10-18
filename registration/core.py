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
#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
#if problems with libstdc++ : http://askubuntu.com/questions/575505/glibcxx-3-4-20-not-found-how-to-fix-this-error
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
from openmesh import *
#Importing the rest of utilities
import numpy as np
from numpy import linalg
from scipy import spatial
import helpers

def affinity_to_correspondences(features2, weights2, flags2, affinity, correspondingFeatures, correspondingWeights, correspondingFlags, flagRoundingLimit = 0.8):
    """
    # GOAL
    This function computes all the corresponding features, weights and flags,
    when the affinity for a set of features, weights and flags is given.

    # INPUTS
    -features2
    -weights2
    -flags2
    -affinity

    # PARAMETERS
    -flagRoundingLimit:
    Flags are binary. Anything over this double is rounded up, whereas anything
    under it is rounded down. A suggested cut-off is 0.8, which means that if
    an element flagged zero contributes 20 percent or to the affinity, the
    corresponding element should be flagged a zero as well.

    # OUTPUT
    -correspondingFeatures
    -correspondingWeights
    -correspondingFlags

    # RETURNS
    """
    correspondingFeatures = affinity.dot(features2)
    correspondingWeights = affinity.dot(weights2)
    correspondingFlags = affinity.dot(flags2)
    # Flags are binary. We will round them down if lower than the flag rounding
    # limit (see explanation in parameter description).
    for i,flag in enumerate(correspondingFlags):
        if flag > flagRoundingLimit:
            correspondingFlags[i] = 1.0
        else:
            correspondingFlags[i] = 0.0

def wknn_affinity(features1, features2, weights2, affinity12, k = 3):
    """
    # GOAL
    For each element in features1, you're going to determine affinity weights
    which link it to elements in features2. This is based on the Euclidean
    distance of k nearest neighbours found for each element of features1.

    # INPUTS
    -features1
    -features2
    -weights2


    # PARAMETERS
    -k(=3): number of nearest neighbours

    # OUTPUT
    -indices12
    -affinity12

    # RETURNS
    """
    # Obtain some required info
    numElements1 = features1.shape[0]
    numElements2 = features2.shape[0]
    # Determine the nearest neighbours
    helpers.nearest_neighbours(features1, features2, distances, neighbourIndices, k):
    # Compute the affinity matrix
    ## Initialize the matrix
    affinity12 = np.zeros((numElements1, numElements2), dtype = float)
    ## Loop over the first feature set to determine their affinity with the
    ## second set.
    for i in range(numVertices1):
        ### Loop over k nearest neighbours
        for j, idx in enumerate(indices12[i,:]):
            ### idx now indicates the index of the second set element that was
            ## determined a neighbour of the current first set element
            distance = distances12[i,j]
            if distance < 0.001:
                distance = 0.001 #numeric stability
            ## The affinity is 1 over the squared distance, weighed by the
            ## corresponding weight of the nearest neighbour (weights2)
            affinity = weights2[idx] * 1.0/(distance * distance)
            if affinity < 0.0001:
                affinity = 0.0001 #numeric stability for the normalization
            ### this should be inserted in the affinity matrix element that
            ### links the first set element index i with second element set at
            ### index 'idx'
            affinity12[i,idx] = affinity

    ## Normalize the affinity matrix (the sum of each row has to equal 1.0)
    rowSums = affinity12.sum(axis=1, keepdims=True)
    affinity12 = affinity12 / rowSums
    return True

def fuse_affinities(affinity12, affinity21, weights1, weights2, affinityTotal):
    """
    # GOAL
    Fuse the two affinity matrices together, taking weights into account.

    # INPUTS
    -affinity12
    -affinity21: dimensions should be transposed of affinity12
    -weights1
    -weights2

    # PARAMETERS
    -k(=3): number of nearest neighbours

    # OUTPUT
    -affinityTotal

    # RETURNS
    """
    # Obtain Info
    numElements1 = affinity12.shape[0]
    numElements2 = affinity12.shape[1]
    # Initialize total affinity matrix
    for i in range(numElements1):
        w1 = weights1[i]
        w2 = weights2[i]
        sumWeights = w1 + w2
        if sumWeights > 0.0001:
            #TODO: should tranpose affinity21?
            affinityTotal[i,:] = (w2 * affinity12[i,:] + w1 * affinity21.tranpose()[i,:]) / sumWeights
        else:
            affinityTotal[i,:] = np.zeros((numElements2), dtype = float)

    # Normaly, the affinity matrix should be normalized now. TODO: verify this
    #rowSums = affinityFloatingToTarget.sum(axis=1, keepdims=True)
    #affinityFloatingToTarget = affinityFloatingToTarget / rowSums
    return True

def inlier_detection(features, correspondingFeatures, neighbourIndices):
    """
    # GOAL
    Determine which elements are inliers ('normal') and which are outliers
    ('abnormal'). There are different ways to to that, but they should all
    result in a scalar assigmnent that represents the probability of an element
    being an inlier.

    # INPUTS
    -features
    -correspondingFeatures
    -neighbourIndices


    # PARAMETERS


    # OUTPUT
    -correspondingFeatures
    -indices12
    -indices21

    # RETURNS
    """
    # Flagged element handling
    ## Get indices of flagged elements
    ## Find which elements

    numFeatures = features.shape[0]
    sigmaa = 0.0
    lambdaa = 0.0
    sigmaNumerator = 0.0
    sigmaDenominator = 0.0
    for i in range(numFeatures):
        distance = np.linalg.norm(correspondingFeatures[i,:] - features[i,:])
        sigmaNumerator = sigmaNumerator + floatingInlierProb[i] * distance * distance
        sigmaDenominator = sigmaDenominator + floatingInlierProb[i]

        sigmaa = np.sqrt(sigmaNumerator/sigmaDenominator)
        lambdaa = 1.0/(np.sqrt(2 * 3.14159) * sigmaa) * np.exp(-0.5 * kappaa * kappaa)

        ### Recalculate the weights
        for i in range(numFloatingVertices):
            distance = np.linalg.norm(correspondingFeatures[:,i] - floatingFeatures[:,i])
            inlierProbability = 1.0/(np.sqrt(2 * 3.14159) * sigmaa) * np.exp(-0.5 * np.square(distance/sigmaa))
            floatingInlierProb[i] = inlierProbability / (inlierProbability + lambdaa) #TODO: I left the weight calculation out for now

    return True
