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







def nearest_neighbours(features1, features2, distances, neighbourIndices, k = 3, leafsize = 15, eps = 0.0001, p = 2, maxDistance = 1000):
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

    OUTPUT
    -distances
    -neighbourIndices
    -kd-tree (see 'RETURN')

    RETURNS
    -kd-tree: This function returns the kd-tree that is built internally.
    """
    kdTree = spatial.cKDTree(features2, leafsize)
    ##query the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query)
    distances, neighbourIndices = kdTree.query(features1, k, eps, p, maxDistance)

    return kdTree
