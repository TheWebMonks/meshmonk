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



def wknn_affinity(features1, features2, affinity12, k = 3):
    """
    # GOAL
    For each element in features1, you're going to determine affinity weights
    which link it to elements in features2. This is based on the Euclidean
    distance of k nearest neighbours found for each element of features1.

    # INPUTS
    -features1
    -features2

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
    helpers.nearest_neighbours(features1, features2, distances, neighbourIndices, k)
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
            ## The affinity is 1 over the squared distance
            affinity = 1.0 / (distance * distance)
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

def fuse_affinities(affinity12, affinity21, affinityFused):
    """
    # GOAL
    Fuse the two affinity matrices together.

    # INPUTS
    -affinity12
    -affinity21: dimensions should be transposed of affinity12

    # PARAMETERS

    # OUTPUT
    -affinityFused

    # RETURNS
    """
    # Fusing is done by simple averaging
    affinityFused = affinity12 + affinity21.transpose()
    rowSums = affinityFused.sum(axis=1, keepdims=True)
    affinityFused = affinityFused / rowSums

    return True

def affinity_to_correspondences(features2, flags2, affinity, correspondingFeatures, correspondingFlags, flagRoundingLimit = 0.9):
    """
    # GOAL
    This function computes all the corresponding features and flags,
    when the affinity for a set of features and flags is given.

    # INPUTS
    -features2
    -flags2
    -affinity

    # PARAMETERS
    -flagRoundingLimit:
    Flags are binary. Anything over this double is rounded up, whereas anything
    under it is rounded down. A suggested cut-off is 0.9, which means that if
    an element flagged zero contributes 10 percent or to the affinity, the
    corresponding element should be flagged a zero as well.

    # OUTPUT
    -correspondingFeatures
    -correspondingFlags

    # RETURNS
    """
    correspondingFeatures = affinity.dot(features2)
    correspondingFlags = affinity.dot(flags2)
    # Flags are binary. We will round them down if lower than the flag rounding
    # limit (see explanation in parameter description).
    for i,flag in enumerate(correspondingFlags):
        if flag > flagRoundingLimit:
            correspondingFlags[i] = 1.0
        else:
            correspondingFlags[i] = 0.0

def inlier_detection(features, correspondingFeatures, correspondingFlags, inlierProbability):
    """
    # GOAL
    Determine which elements are inliers ('normal') and which are outliers
    ('abnormal'). There are different ways to to that, but they should all
    result in a scalar assigmnent that represents the probability of an element
    being an inlier.

    # INPUTS
    -features
    -correspondingFeatures
    -correspondingFlags

    # PARAMETERS

    # OUTPUT
    -inlierProbability

    # RETURNS
    """
    # Info
    numElements = features.shape[0]

    # Flag based inlier/outlier classification
    for i,flag in enumerate(correspondingFlags):
        if flag < 0.5:
            inlierProbability[i] = 0.0

    # Distance based inlier/outlier classification
    ## Re-calculate the parameters sigma and lambda
    sigmaa = 0.0
    lambdaa = 0.0
    sigmaNumerator = 0.0
    sigmaDenominator = 0.0
    for i in range(numElements):
        distance = np.linalg.norm(correspondingFeatures[i,:] - features[i,:])
        sigmaNumerator = sigmaNumerator + inlierProbability[i] * distance * distance
        sigmaDenominator = sigmaDenominator + inlierProbability[i]

    sigmaa = np.sqrt(sigmaNumerator/sigmaDenominator)
    lambdaa = 1.0/(np.sqrt(2 * 3.14159) * sigmaa) * np.exp(-0.5 * kappaa * kappaa)
    ## Recalculate the distance-based probabilities
    for i in range(numElements):
        distance = np.linalg.norm(correspondingFeatures[:,i] - features[:,i])
        probability = 1.0/(np.sqrt(2 * 3.14159) * sigmaa) * np.exp(-0.5 * np.square(distance/sigmaa))
        inlierProbability[i] = probability / (probability + lambdaa) #TODO: I left the weight calculation out for now

    #Gradient Based inlier/outlier classification
    for i in range(numElements):
        normal = features[i,3:6]
        correspondingNormal = features[i,3:6]
        ## Dot product gives an idea of how well they point in the same
        ## direction. This gives a weight between -1.0 and +1.0
        dotProduct = normal.dot(correspondingNormal)
        ## Rescale this result so that it's continuous between 0.0 and +1.0
        probability = dotProduct / 2.0 + 0.5
        inlierProbability[i] = inlierProbability[i] * probability

    return True
