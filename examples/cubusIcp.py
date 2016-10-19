#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
import openmesh as om
#Importing the rest of utilities
import numpy as np
from scipy import linalg

from context import registration

"""
In this example, we will perform rigid registration of a self-made object,
which are basically one pyramid tied to another pyramid upside down.

#####
#Some definitions and info:
#####

PROCESS PARAMETERS
-maxNumIterations: a hard limit for the maximum number of iterations in the ICP
iterative loop.
-maxNumNeighbours: the maximum number of vertices that can be considered one
vertex's neighbours. (Neighbours can be determined in various ways. They are
important for some operations where neighbouring vertices play a role in the
computation.)
-kappaa: mahalanobis distance for outlier determination (e.g. kappaa = 2 means
that 2 standard deviations of the distribution is the middle of the transition
zone from inlier to outlier classification)
-adjustScale: whether or not to change the size of the floating mesh during ICP
"""

##########
# SET PARAMETERS
##########
# ICP management
maxNumIterations = 1
maxNumNeighbours = 20 #
# Correspondences
wknnNumNeighbours = 1
# Inlier Detection
kappaa = 3
adjustScale = False
floatingMeshPath = "/home/jonatan/kuleuven-algorithms/examples/data/cubus.obj"
targetMeshPath = "/home/jonatan/kuleuven-algorithms/examples/data/cubus_transformed.obj"
resultingMeshPath = "/home/jonatan/kuleuven-algorithms/examples/data/cubus_registered.obj"

##########
# PREPARE DATA
##########
# Load from file
floatingMesh = om.TriMesh()
targetMesh = om.TriMesh()
om.read_mesh(floatingMesh, floatingMeshPath)
om.read_mesh(targetMesh, targetMeshPath)

# Obtain info and initialize matrices
numFloatingVertices = floatingMesh.n_vertices()
numTargetVertices = targetMesh.n_vertices()
floatingPositions = np.zeros((numFloatingVertices,3), dtype = float)
targetPositions = np.zeros((numTargetVertices,3), dtype = float)
floatingNormals = np.zeros((numFloatingVertices,3), dtype = float)
targetNormals = np.zeros((numTargetVertices,3), dtype = float)
# Initialize weights and flags
floatingWeights = np.ones((numFloatingVertices), dtype = float)
targetWeights = np.ones((numTargetVertices), dtype = float)
floatingFlags = np.ones((numFloatingVertices), dtype = float)
targetFlags = np.ones((numTargetVertices), dtype = float)
#targetFlags[3] = 0.0 #randomly attributed

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
floatingFeatures = np.hstack((floatingPositions,floatingNormals))
targetFeatures = np.hstack((targetPositions,targetNormals))


##########
# ICP
##########
"""

"""
#for iteration in range(0,maxNumIterations):
##1) Determine Nearest neighbours. We'll simply use index correspondences for the bunny
affinity = registration.core.wknn_affinity(floatingFeatures, targetFeatures, wknnNumNeighbours)
correspondingFeatures, correspondingFlags = registration.core.affinity_to_correspondences(targetFeatures, targetFlags, affinity, 0.5)
##2) Determine weights. A weight related to the gaussian distance distribution suffices.
### Update the distribution parameters
floatingWeights = registration.core.inlier_detection(floatingFeatures, correspondingFeatures, correspondingFlags, floatingWeights, 3.0)
##3) Determine and update transformation.
transformationMatrix = registration.core.rigid_transformation(floatingFeatures, correspondingFeatures, floatingWeights, False)



##########
# EXPORT DATA
##########
##New vertex positions
newPosition = om.TriMesh.Point(0.0,0.0,0.0)
for i, vh in enumerate(floatingMesh.vertices()):
    newPosition[0] = floatingFeatures[i,0]
    newPosition[1] = floatingFeatures[i,1]
    newPosition[2] = floatingFeatures[i,2]
    floatingMesh.set_point(vh,newPosition)

##New vertex normals
floatingMesh.request_vertex_normals()
floatingMesh.request_face_normals()
floatingMesh.update_normals()
floatingMesh.release_face_normals()

##Save the mesh
om.write_mesh(floatingMesh, resultingMeshPath)
