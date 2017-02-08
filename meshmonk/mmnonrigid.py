#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 14:12:26 2017

@author: Jonatan Snyders
"""

import openmesh
#Importing the rest of utilities
import numpy
import time
import argparse
import helpers
import core
#from context import meshmonk

"""
In this example, we will perform nonrigid registration of the stanford bunny.
"""

##########
# ARGUMENT PARSING
##########

parser = argparse.ArgumentParser()
# Add arguments for input and output data
parser.add_argument('float', help = 'The absolute path for the floating mesh.')
parser.add_argument('target', help = 'The absolute path for the target mesh.')
parser.add_argument('result', help = 'The absolute path for the nonrigidly registered mesh result.')

# Add parameter arguments
parser.add_argument("--params", help = 'Optional - the absolute path for a parameter file.')

#Parse the arguments
args = parser.parse_args()
## Parameters
if args.params:
    print ("Parameter file used - overwriting standard parameters.\n")
else:
    print ("No parameter file given - using standard parameters.\n")

#Print the arguments for the user
print("Floating mesh path: \n\"" + args.float + "\"")
print("Target mesh path: \n\"" + args.target + "\"")
print("Resulting mesh path: \n\"" + args.result + "\"")


#floatingMeshPath = "/home/jonatan/projects/kuleuven-algorithms/examples/data/fucked_up_bunny.obj"
#targetMeshPath = "/home/jonatan/projects/kuleuven-algorithms/examples/data/bunny90.obj"
#resultingMeshPath = "/home/jonatan/projects/kuleuven-algorithms/examples/data/bunnyNonRigid.obj"



# Correspondences
wknnNumNeighbours = 3
# Inlier Detection
kappaa = 3.0
adjustScale = False
# Transformation



##########
# PREPARE DATA
##########
# Load from file
floatingMesh = openmesh.TriMesh()
targetMesh = openmesh.TriMesh()
openmesh.read_mesh(floatingMesh, args.float)
openmesh.read_mesh(targetMesh, args.target)


## Decimate the mesh
## create decimater and module handle
#decimator = openmesh.PolyMeshDecimater(floatingMesh)
#modHandle = openmesh.PolyMeshModQuadricHandle()
## add modules
#decimator.add(modHandle)
#decimator.module(modHandle).set_max_err(0.001)
#
## decimate
#decimator.initialize()
#decimator.decimate_to(100)
#
#floatingMesh.garbage_collection()
#openmesh.write_mesh(floatingMesh, "/home/jonatan/kuleuven-algorithms/examples/data/openmesh_decimated_bunny.obj")


# Obtain info and initialize matrices
numFloatingVertices = floatingMesh.n_vertices()
numTargetVertices = targetMesh.n_vertices()
## Initialize weights and flags
floatingWeights = numpy.ones((numFloatingVertices), dtype = float)
targetWeights = numpy.ones((numTargetVertices), dtype = float)
floatingFlags = numpy.ones((numFloatingVertices), dtype = float)
targetFlags = numpy.ones((numTargetVertices), dtype = float)

# Obtain the floating and target mesh features (= positions and normals)
floatingFeatures = helpers.openmesh_to_numpy_features(floatingMesh)
targetFeatures = helpers.openmesh_to_numpy_features(targetMesh)
correspondingFeatures = numpy.zeros((floatingFeatures.shape), dtype = float)
correspondingFlags = numpy.ones((floatingFlags.shape), dtype = float)


##########
# NONRIGID REGISTRATION
##########
"""

"""
## Initialize
originalFloatingPositions = numpy.copy(floatingFeatures[:,0:3])
regulatedDisplacementField = numpy.zeros((numFloatingVertices,3), dtype = float)
symCorrespondenceFilter = core.SymCorrespondenceFilter(floatingFeatures,
                                                              floatingFlags,
                                                              targetFeatures,
                                                              targetFlags,
                                                              correspondingFeatures,
                                                              correspondingFlags,
                                                              wknnNumNeighbours)
## Set up inlier filter
inlierFilter = core.InlierFilter(floatingFeatures, correspondingFeatures,
                                              correspondingFlags, floatingWeights,
                                              kappaa)


##TODO: HERE WE SHOULD START THE ANNEALING SCHEME
numViscousSmoothingIterationsList = [55, 34, 21, 13, 8, 5, 3, 2, 1, 1]
numElasticSmoothingIterationsList = [55, 34, 21, 13, 8, 5, 3, 2, 1, 1]
## Set up transformation filter
numNeighbourDisplacements = 10
sigmaSmoothing = 10.0
transformationFilter = core.ViscoElasticFilter(floatingFeatures,
                                                            correspondingFeatures,
                                                            floatingWeights,
                                                            10,
                                                            sigmaSmoothing,
                                                            numViscousIterations = 1,
                                                            numElasticIterations = 1)
iteration = 0
for numViscousSmoothingIterations, numElasticSmoothingIterations in zip(numViscousSmoothingIterationsList, numElasticSmoothingIterationsList):
    timeStart = time.time()
    ## 1) Determine Nearest neighbours.
    symCorrespondenceFilter.set_floating_features(floatingFeatures, floatingFlags)
    symCorrespondenceFilter.update()
    ## 2) Determine inlier weights.
    inlierFilter.update()
    ## 3) Determine and update transformation.
    transformationFilter.set_parameters(numNeighbourDisplacements,
                                            sigmaSmoothing,
                                            numViscousSmoothingIterations,
                                            numElasticSmoothingIterations)
    transformationFilter.update()

    ## 4) Re-calculate the mesh's properties (like normals e.g.)
    floatingFeatures[:,3:6] = helpers.openmesh_normals_from_positions(floatingMesh, floatingFeatures[:,0:3])
    
    timeEnd = time.time()
    print "Iteration " + str(iteration) + " took " + str(timeEnd-timeStart) 
    iteration = iteration + 1
##########
# EXPORT DATA
##########
# Save the mesh
openmesh.write_mesh(floatingMesh, args.result)
print "Exported result."
