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
import openmesh as om
import numpy as np
from scipy import linalg

from context import registration

##########
# TEST VECTOR INTERPOLATION
##########
# Set up test data
position = np.array([0.0,0.0,0.0], dtype = float)
fieldPositions = np.array([[1.0,1.0,0.0],[1.0,-1.0,0.0],[-1.0,-1.0,0.0],[-1.0,1.0,0.0]], dtype = float)
fieldVectors = np.array([[1.0,1.0,1.0],[1.0,-1.0,1.0],[-1.0,-1.0,1.0],[-1.0,1.0,1.0]], dtype = float)
fieldWeights = np.ones((4), dtype = float)
sigmaa = 1.0
expectedResult = np.array([0.0,0.0,1.0], dtype = float)
interpolatedVector = registration.helpers.gaussian_vector_interpolation(position, fieldPositions, fieldVectors, fieldWeights, sigmaa)
# Check if the result equals the expected result.
if linalg.norm(interpolatedVector - expectedResult) < 0.00001:
    print('gaussian_vector_interpolation is successful!')
else:
    print('A problem has occured in gaussian_vector_interpolation!')

##########
# TEST SCALAR INTERPOLATION
##########
# Set up test data
position = np.array([0.0,0.0,0.0], dtype = float)
fieldPositions = np.array([[1.0,1.0,0.0],[1.0,-1.0,0.0],[-1.0,-1.0,0.0],[-1.0,1.0,0.0]], dtype = float)
fieldScalars = np.array([1.0,2.0,3.0,4.0], dtype = float)
fieldWeights = np.ones((4), dtype = float)
sigmaa = 1.0
expectedResult = 2.5
interpolatedScalar = registration.helpers.gaussian_scalar_interpolation(position, fieldPositions, fieldScalars, fieldWeights, sigmaa)
# Check if the result equals the expected result.
if linalg.norm(interpolatedScalar - expectedResult) < 0.00001:
    print('gaussian_scalar_interpolation is successful!')
else:
    print('A problem has occured in gaussian_scalar_interpolation!')
