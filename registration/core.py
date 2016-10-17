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

def symmetric_wknn_correspondences(features1, features2, correspondingFeatures, indices12, indices21, k = 3):
    """
    # GOAL
    For each element in features1, you're going to determine corresponding
    features by looking at the elements present in features2. This is done with
    'symmetric weighted k nearest neighbours'. Concretely:
    ## K Nearest Neighbours
    The k nearest neighbours for each element in features1 are sought in
    features2 using the Euclidean distance metric
    ## Weighted:
    The k nearest neighbours are combined through weighted averaging
    to obtain a correspondence
    ## Symmetric:
    The final correspondences in features2 for features1 are determined by
    combining weighted k-nn correspondences in features2 for features1 and
    weighted k-nn correspondences in features1 for features2. This can be though
    of as combining push and pull forces.

    # INPUTS
    -features1
    -features2


    # PARAMETERS


    # OUTPUT
    -correspondingFeatures
    -indices12
    -indices21

    # RETURNS
    """
    helpers.nearest_neighbours(features1, features2, distances, neighbourIndices, k):
    helpers.nearest_neighbours(features1, features2, distances, neighbourIndices, k):
    # Construct the affinity matrices
    ## Initialize the matrices
    affinity12 = np.zeros((numVertices1, numVertices2), dtype = float)
    affinity21 = np.zeros((numVertices2, numVertices1), dtype = float)
    ## Loop over the first feature set to determine their affinity with the
    ## second set.
    for i in range(numVertices1):
        ### Loop over k nearest neighbours found in the second set
        for j, idx in enumerate(indices12[i,:]):
            ### idx now indicates the index of the second set element that was
            ## determined a neighbour of the current first set element
            distance = distances12[i,j]
            if distance < 0.001:
                distance = 0.001 #numeric stability
            weight = 1.0/(distance * distance)
            ### this weight should be inserted in the affinity matrix element that
            ### links the first set element index i with second element set at
            ### index 'idx'
            affinity12[i,idx] = weight

    ## Loop over the second set elements to determine their affinity with the
    ## first set elements
    for i in range(numVertices2):
        ### Loop over k nearest neighbours found in the first set elements
        for j, idx in enumerate(indices21[i,:]):
            ### idx now indicates the index of the first set element that was
            ### determined a neighbour of the current second set element
            distance = distances21[i,j]
            if distance < 0.001:
                distance = 0.001 #numeric stability
            weight = 1.0/(distance * distance)
            ### This weight should be inserted in the affinity matrix element
            ### that links the second set element at index i with the first set
            ### element at index 'idx'
            affinity21[i,idx] = weight

    ## Normalize both affinity matrices (the sum of each row has to equal 1.0)
    rowSums = affinity12.sum(axis=1, keepdims=True)
    affinity12 = affinity12 / rowSums
    rowSums = affinity21.sum(axis=1, keepdims=True)
    affinity21 = affinity21 / rowSums
    ## Sum the affinity matrices and normalize again
    affinity12 = affinity12 + affinity21.transpose()
    rowSums = affinity12.sum(axis=1, keepdims=True)
    affinity12 = affinity12 / rowSums
    ## Now we have constructed the affinity matrix using push-pull forces, let's
    ## construct the corresponding features.
    correspondingFeatures = affinity12.dot(targetFeatures)

def inlier_detection(features, correspondingFeatures, indices12, indices21, k = 3):
    """
    # GOAL
    For each element in features1, you're going to determine corresponding
    features by looking at the elements present in features2. This is done with
    'symmetric weighted k nearest neighbours'. Concretely:
    ## K Nearest Neighbours
    The k nearest neighbours for each element in features1 are sought in
    features2 using the Euclidean distance metric
    ## Weighted:
    The k nearest neighbours are combined through weighted averaging
    to obtain a correspondence
    ## Symmetric:
    The final correspondences in features2 for features1 are determined by
    combining weighted k-nn correspondences in features2 for features1 and
    weighted k-nn correspondences in features1 for features2. This can be though
    of as combining push and pull forces.

    # INPUTS
    -features1
    -features2


    # PARAMETERS


    # OUTPUT
    -correspondingFeatures
    -indices12
    -indices21

    # RETURNS
    """
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
