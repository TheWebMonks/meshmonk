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

from context import registration

"""
In this example, we will perform nonrigid registration of the stanford bunny.
"""

##########
# SET PARAMETERS
##########
# Data I/O
floatingMeshPath = "/home/jonatan/kuleuven-algorithms/examples/data/fucked_up_bunny.obj"
targetMeshPath = "/home/jonatan/kuleuven-algorithms/examples/data/bunny90.obj"
resultingMeshPath = "/home/jonatan/kuleuven-algorithms/examples/data/bunnyNonRigid.obj"
# Correspondences
wknnNumNeighbours = 3
# Inlier Detection
kappaa = 3
adjustScale = False
# Transformation
numNeighbourDisplacements = 6
sigmaSmoothing = 0.01
numViscousSmoothingIterationsList = [10, 5, 2, 1]
numElasticSmoothingIterationsList = [10, 5, 2, 1]

##########
# PREPARE DATA
##########
# Load from file
floatingMesh = openmesh.TriMesh()
targetMesh = openmesh.TriMesh()
openmesh.read_mesh(floatingMesh, floatingMeshPath)
openmesh.read_mesh(targetMesh, targetMeshPath)

# Decimate the mesh
# create decimater and module handle
decimator = openmesh.PolyMeshDecimater(floatingMesh)
modHandle = openmesh.PolyMeshModQuadricHandle()
# add modules
decimator.add(modHandle)
decimator.module(modHandle).set_max_err(0.001)

# decimate
decimator.initialize()
decimator.decimate_to(100)

floatingMesh.garbage_collection()
openmesh.write_mesh(m, "/home/jonatan/kuleuven-algorithms/examples/data/openmesh_decimated_bunny.obj")


# Obtain info and initialize matrices
numFloatingVertices = floatingMesh.n_vertices()
numTargetVertices = targetMesh.n_vertices()
## Initialize weights and flags
floatingWeights = numpy.ones((numFloatingVertices), dtype = float)
targetWeights = numpy.ones((numTargetVertices), dtype = float)
floatingFlags = numpy.ones((numFloatingVertices), dtype = float)
targetFlags = numpy.ones((numTargetVertices), dtype = float)

# Obtain the floating and target mesh features (= positions and normals)
floatingFeatures = registration.helpers.openmesh_to_numpy_features(floatingMesh)
targetFeatures = registration.helpers.openmesh_to_numpy_features(targetMesh)
correspondingFeatures = numpy.zeros((floatingFeatures.shape), dtype = float)


##########
# NONRIGID REGISTRATION
##########
"""

"""
## Initialize
originalFloatingPositions = numpy.copy(floatingFeatures[:,0:3])
regulatedDisplacementField = numpy.zeros((numFloatingVertices,3), dtype = float)

##TODO: HERE WE SHOULD START THE ANNEALING SCHEME
for numViscousSmoothingIterations, numElasticSmoothingIterations in zip(numViscousSmoothingIterationsList, numElasticSmoothingIterationsList):
    ## 1) Determine Nearest neighbours.
    affinity = registration.core.wknn_affinity(floatingFeatures, targetFeatures, wknnNumNeighbours)
    ###TODO: REMOVE THIS HACK WHERE WE SET AFFINITY TO UNITY (BECAUSE BUNNIES CORRESPOND BY INDEX)
    affinity = numpy.identity(numFloatingVertices, dtype = float)
    correspondingFeatures, correspondingFlags = registration.core.affinity_to_correspondences(targetFeatures, targetFlags, affinity, 0.5)
    ## 2) Determine inlier weights.
    floatingWeights = registration.core.inlier_detection(floatingFeatures, correspondingFeatures, correspondingFlags, floatingWeights, kappaa)
    ## 3) Determine and update transformation.
    ### Compute a viscoelastic transformation
    registration.core.compute_viscoelastic_transformation(floatingFeatures[:,0:3], correspondingFeatures[:,0:3], floatingWeights, regulatedDisplacementField, numNeighbourDisplacements, sigmaSmoothing, numViscousSmoothingIterations, numElasticSmoothingIterations)

    ### Apply the computed displacement field
    floatingFeatures[:,0:3] = originalFloatingPositions + regulatedDisplacementField

    ## 4) Re-calculate the mesh's properties (like normals e.g.)
    floatingFeatures[:,3:6] = registration.helpers.openmesh_normals_from_positions(floatingMesh, floatingFeatures[:,0:3])

##########
# EXPORT DATA
##########
# Save the mesh
openmesh.write_mesh(floatingMesh, resultingMeshPath)
