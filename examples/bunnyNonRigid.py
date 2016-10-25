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
# Correspondences
wknnNumNeighbours = 3
# Inlier Detection
kappaa = 3
adjustScale = False
floatingMeshPath = "/home/jonatan/kuleuven-algorithms/examples/data/fucked_up_bunny.obj"
targetMeshPath = "/home/jonatan/kuleuven-algorithms/examples/data/bunny90.obj"
resultingMeshPath = "/home/jonatan/kuleuven-algorithms/examples/data/bunnyNonRigid.obj"
# Transformation
numNeighbourForces = 6
numNeighbourDisplacements = 6
sigmaSmoothing = 1.0 # no idea on right value

##########
# PREPARE DATA
##########
# Load from file
floatingMesh = openmesh.TriMesh()
targetMesh = openmesh.TriMesh()
openmesh.read_mesh(floatingMesh, floatingMeshPath)
openmesh.read_mesh(targetMesh, targetMeshPath)

# Obtain info and initialize matrices
numFloatingVertices = floatingMesh.n_vertices()
numTargetVertices = targetMesh.n_vertices()
floatingPositions = numpy.zeros((numFloatingVertices,3), dtype = float)
targetPositions = numpy.zeros((numTargetVertices,3), dtype = float)
floatingNormals = numpy.zeros((numFloatingVertices,3), dtype = float)
targetNormals = numpy.zeros((numTargetVertices,3), dtype = float)
# Initialize weights and flags
floatingWeights = numpy.ones((numFloatingVertices), dtype = float)
targetWeights = numpy.ones((numTargetVertices), dtype = float)
floatingFlags = numpy.ones((numFloatingVertices), dtype = float)
targetFlags = numpy.ones((numTargetVertices), dtype = float)

# Put in the right positions and normals
## Make sure normals are present
floatingMesh.request_vertex_normals()
floatingMesh.request_face_normals()
floatingMesh.update_normals()
floatingMesh.release_face_normals()
for i, vh in enumerate(floatingMesh.vertices()):
    floatingPositions[i,0] = floatingMesh.point(vh)[0]
    floatingPositions[i,1] = floatingMesh.point(vh)[1]
    floatingPositions[i,2] = floatingMesh.point(vh)[2]
    floatingNormals[i,0] = floatingMesh.normal(vh)[0]
    floatingNormals[i,1] = floatingMesh.normal(vh)[1]
    floatingNormals[i,2] = floatingMesh.normal(vh)[2]

targetMesh.request_vertex_normals()
targetMesh.request_face_normals()
targetMesh.update_normals()
targetMesh.release_face_normals()
for i, vh in enumerate(targetMesh.vertices()):
    targetPositions[i,0] = targetMesh.point(vh)[0]
    targetPositions[i,1] = targetMesh.point(vh)[1]
    targetPositions[i,2] = targetMesh.point(vh)[2]
    targetNormals[i,0] = targetMesh.normal(vh)[0]
    targetNormals[i,1] = targetMesh.normal(vh)[1]
    targetNormals[i,2] = targetMesh.normal(vh)[2]

# The features are the concatenation of positions and normals
floatingFeatures = numpy.hstack((floatingPositions,floatingNormals))
targetFeatures = numpy.hstack((targetPositions,targetNormals))


##########
# NONRIGID REGISTRATION
##########
"""

"""
## Initialize
initialFloatingPositions = numpy.copy(floatingPositions)
regulatedDisplacementField = numpy.zeros((numFloatingVertices,3), dtype = float)

##TODO: HERE WE SHOULD START THE ITERATIONS

## 1) Determine Nearest neighbours.
affinity = registration.core.wknn_affinity(floatingFeatures, targetFeatures, wknnNumNeighbours)
###TODO: REMOVE THIS CHEAT WHERE WE SET AFFINITY TO UNITY (BECAUSE BUNNIES CORRESPOND BY INDEX)
affinity = numpy.identity(numFloatingVertices, dtype = float)
correspondingFeatures, correspondingFlags = registration.core.affinity_to_correspondences(targetFeatures, targetFlags, affinity, 0.5)
## 2) Determine inlier weights.
floatingWeights = registration.core.inlier_detection(floatingFeatures, correspondingFeatures, correspondingFlags, floatingWeights, kappaa)
## 3) Determine and update transformation.
### We only use the positions here, not the normals. So let's separate them.
floatingPositions = floatingFeatures[:,0:3]
targetPositions = targetFeatures[:,0:3]
correspondingPositions = correspondingFeatures[:,0:3]
floatingNormals = floatingFeatures[:,3:6]
targetNormals = targetFeatures[:,3:6]
correspondingNormals = correspondingFeatures[:,3:6]
### The 'Force Field' is what drives the deformation: the difference between
### the floating vertices and their correspondences
forceField = correspondingPositions - floatingPositions
regulatedForceField = numpy.zeros((numFloatingVertices,3), dtype = float)
### Let's regulate the Force Field
registration.helpers.gaussian_smoothing_displacement_field(floatingPositions, forceField, regulatedForceField, floatingWeights, numNeighbourForces, sigmaSmoothing)

### Add the regulated Force Field to the Displacement Field
displacementField = regulatedDisplacementField + regulatedForceField
### Regulate the Displacement Field (same to what we did with the Force Field)
registration.helpers.gaussian_smoothing_displacement_field(floatingPositions, displacementField, regulatedDisplacementField, floatingWeights, numNeighbourDisplacements, sigmaSmoothing)

### Apply the regulated Displacement Field
floatingPositions = initialFloatingPositions + regulatedDisplacementField
### Re-calculate the normals
floatingNormals = registration.helpers.openmesh_normals_from_positions(floatingMesh, floatingPositions)
### Re-combine the new positions and normals into the floating feature set
floatingFeatures = numpy.hstack((floatingPositions,floatingNormals)) #TODO: UITPUT BUNNY BEKIJKEN?

##########
# EXPORT DATA
##########
# Save the mesh
openmesh.write_mesh(floatingMesh, resultingMeshPath)
