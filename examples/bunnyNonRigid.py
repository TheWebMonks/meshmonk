#import sys
#sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
#Importing the rest of utilities

#from context import meshmonk
from meshmonk import core

"""
In this example, we will perform nonrigid registration of the stanford bunny.
"""

##########
# SET PARAMETERS
##########
# Data I/O
floatingMeshPath = "/home/jonatan/projects/meshmonk/examples/data/fucked_up_bunny.obj"
targetMeshPath = "/home/jonatan/projects/meshmonk/examples/data/bunny90.obj"
resultingMeshPath = "/home/jonatan/projects/meshmonk/examples/data/bunnyNonRigid.obj"

nonrigidTransformer = core.RegistrationManager(floatingMeshPath, targetMeshPath, resultingMeshPath, 'full')
nonrigidTransformer.update()
