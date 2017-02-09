#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 14:12:26 2017

@author: Jonatan Snyders
"""

#Importing the rest of utilities
import argparse
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



nonrigidTransformer = core.RegistrationManager(args.float, args.target, args.result, 'nonrigid')
nonrigidTransformer.update()
