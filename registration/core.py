#Linking to the OpenMesh bindings. (Check out http://stackoverflow.com/questions/33637373/openmesh-with-python-3-4 for troubleshooting)
#if problems with libstdc++ : http://askubuntu.com/questions/575505/glibcxx-3-4-20-not-found-how-to-fix-this-error
import sys
sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
from openmesh import *
#Importing the rest of utilities
import numpy
from numpy import linalg
from scipy import spatial
from scipy.interpolate import interp1d
from abc import ABCMeta, abstractmethod
import helpers



def wknn_affinity(features1, features2, neighbourDistances, neighbourIndices, k = 3):
    """
    # GOAL
    For each element in features1, you're going to determine affinity weights
    with its neighbouring elements in features2. This is based on the Euclidean
    distance of k nearest neighbours found for each element of features1.

    # INPUTS
    -features1
    -features2
    -neighbourDistances
    -neighbourIndices

    # PARAMETERS
    -k(=3): number of nearest neighbours

    # RETURNS
    -affinity12
    """
    # Obtain some required info and initialize data structures
    numElements1 = features1.shape[0]
    numElements2 = features2.shape[0]
    affinity12 = numpy.zeros((numElements1, numElements2), dtype = float)
    # Compute the affinity matrix
    ## Loop over the first feature set to determine their affinity with the
    ## second set.
    for i in range(numElements1):
        if k == 1:
            idx = neighbourIndices[i]
            dist = neighbourDistances[i]
            if dist < 0.001:
                dist = 0.001 #numeric stability
            ## The affinity is 1 over the squared distance
            affinityElement = 1.0 / (dist * dist)
            if affinityElement < 0.0001:
                affinityElement = 0.0001 #numeric stability for the normalization
            affinity12[i,idx] = affinityElement
        else: ### Loop over k nearest neighbours
            for j, idx in enumerate(neighbourIndices[i,:]):
                ### idx now indicates the index of the second set element that was
                ## determined a neighbour of the current first set element
                distance = neighbourDistances[i,j]
                if distance < 0.001:
                    distance = 0.001 #numeric stability
                ## The affinity is 1 over the squared distance
                affinityElement = 1.0 / (distance * distance)
                if affinityElement < 0.0001:
                    affinityElement = 0.0001 #numeric stability for the normalization
                ### this should be inserted in the affinity matrix element that
                ### links the first set element index i with second element set at
                ### index 'idx'
                affinity12[i,idx] = affinityElement

    ## Normalize the affinity matrix (the sum of each row has to equal 1.0)
    rowSums = affinity12.sum(axis=1, keepdims=True)
    affinity12 = affinity12 / rowSums
    return affinity12

    
def fuse_affinities(affinity12, affinity21):
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

    return affinityFused

def affinity_to_correspondences(features2, flags2, affinity, flagRoundingLimit = 0.9):
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

    # RETURNS
    -correspondingFeatures
    -correspondingFlags
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

    return correspondingFeatures, correspondingFlags

class GenericCorrespondenceFilter(object):
    # A generic template class for correspondence filters
    
    __metaclass__ = ABCMeta
    
    def __init__(self, inFloatingFeatures, inFloatingFlags,
                 inTargetFeatures, inTargetFlags,
                 outCorrespondingFeatures, outCorrespondingFlags):
        self.inFloatingFeatures = inFloatingFeatures
        self.inFloatingFlags = inFloatingFlags
        self.inTargetFeatures = inTargetFeatures
        self.inTargetFlags = inTargetFlags
        self.outCorrespondingFeatures = outCorrespondingFeatures
        self.outCorrespondingFlags = outCorrespondingFlags
        self._numFloatingElements = self.inFloatingFeatures.shape[0]
        self._numTargetElements = self.inTargetFeatures.shape[0]
        self._affinity = numpy.zeros([self._numFloatingElements,
                                      self._numTargetElements],
                                    dtype = float)
        self._outputUpToDate = False
        
    def set_floating_features(self, inFloatingFeatures, inFloatingFlags):
        self.inFloatingFeatures = inFloatingFeatures
        self.inFloatingFlags = inFloatingFlags
        self._outputUpToDate = False
        
    def set_target_features(self, inTargetFeatures, inTargetFlags):
        self.inTargetFeatures = inTargetFeatures
        self.inTargetFlags = inTargetFlags
        self._outputUpToDate = False
        
    def set_corresponding_features(self, outCorrespondingFeatures, outCorrespondingFlags):
        self.outCorrespondingFeatures = outCorrespondingFeatures
        self.outCorrespondingFlags = outCorrespondingFlags
        self._outputUpToDate = False
    
    def _get_affinity(self):
        # Use with care!
        return self._affinity
        
    @abstractmethod
    def _construct_affinity(self):
        # This is the main method each correspondence filter should implement
        # by itself!
        pass
    
    def _affinity_to_correspondences(self):
        paramFlagRoundingLimit = 0.9
        
        correspondingFeatures, correspondingFlags = \
            affinity_to_correspondences(self.inTargetFeatures,
                                        self.inTargetFlags,
                                        self._affinity, 
                                        paramFlagRoundingLimit)
        # Copy result into the user requested output
        for i in range(self._numFloatingElements):
            self.outCorrespondingFeatures[i,:] = correspondingFeatures[i,:]
            self.outCorrespondingFlags[i] = correspondingFlags[i]
            
        self._outputUpToDate = True
    
    def update(self):
        # Construct the affinity
        self._construct_affinity()
        # Use the affinity matrix to determine corresponding features (& flags)
        self._affinity_to_correspondences()
    
class CorrespondenceFilter(GenericCorrespondenceFilter):
    #NOTE: we don't want this object to hold any data, just references to 
    # input and output data, and perform processing on inputs to put those in 
    # those outputs.
    
    def __init__(self, inFloatingFeatures, inFloatingFlags,
                 inTargetFeatures, inTargetFlags,
                 outCorrespondingFeatures, outCorrespondingFlags,
                 paramNumNeighbours = 3):
        GenericCorrespondenceFilter.__init__(self, inFloatingFeatures, inFloatingFlags,
                                             inTargetFeatures, inTargetFlags,
                                             outCorrespondingFeatures, outCorrespondingFlags)
        self.paramNumNeighbours = paramNumNeighbours
        self._kdTree = None
        self._kdTreeUpToDate = False
        self._neighbourDistances = None
        self._neighbourIndices = None
        
        
    def set_target_features(self, inTargetFeatures, inTargetFlags):
        # When the target features are changed, this filter needs to 
        # reconstruct its kd-tree.
        GenericCorrespondenceFilter.set_target_features(self, inTargetFeatures,
                                                        inTargetFlags)
        self._kdTreeUpToDate = False
        

    def _nearest_neighbours(self):
        # init kd-tree
        if (self._kdTreeUpToDate == False):
            self._kdTree = spatial.cKDTree(self.inTargetFeatures, 15)
            self._kdTreeUpToDate = True
        # Query the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query)
        ## set some required parameters
        paramEps = 0.0001
        paramP = 2
        paramMaxDistance = 1000
        ## actual querying for each floating position
        self._neighbourDistances, self._neighbourIndices = \
                    self._kdTree.query(self.inFloatingFeatures,
                                       self.paramNumNeighbours,
                                       paramEps,
                                       paramP,
                                       paramMaxDistance)
        
    
    def _construct_affinity(self):
        # Compute nearest neighbours
        self._nearest_neighbours()
        # Compute affinity from floating to neighbouring target features
        self._affinity = wknn_affinity(self.inFloatingFeatures,
                                       self.inTargetFeatures,
                                       self._neighbourDistances,
                                       self._neighbourIndices,
                                       self.paramNumNeighbours)
    
        
        
        
class SymCorrespondenceFilter(GenericCorrespondenceFilter):
    
    def __init__(self, inFloatingFeatures, inFloatingFlags,
                 inTargetFeatures, inTargetFlags,
                 outCorrespondingFeatures, outCorrespondingFlags,
                 paramNumNeighbours = 3):
        self.paramNumNeighbours = paramNumNeighbours
        GenericCorrespondenceFilter.__init__(self, inFloatingFeatures, inFloatingFlags,
                                             inTargetFeatures, inTargetFlags,
                                             outCorrespondingFeatures, outCorrespondingFlags)
        self.pushFilter = CorrespondenceFilter(inFloatingFeatures, inFloatingFlags,
                                               inTargetFeatures, inTargetFlags,
                                               [], [],
                                               self.paramNumNeighbours)
        # The pull filter does the reverse of the push filter. It looks for
        # correspondences for the target features in the floating features set
        self.pullFilter = CorrespondenceFilter(inTargetFeatures, inTargetFlags,
                                               inFloatingFeatures, inFloatingFlags,
                                               [], [],
                                               self.paramNumNeighbours)
    
    def set_floating_features(self, inFloatingFeatures, inFloatingFlags):
        self.pushFilter.set_floating_features(inFloatingFeatures, inFloatingFlags)
        self.pullFilter.set_target_features(inFloatingFeatures, inFloatingFlags)
        self._outputUpToDate = False
        
    def set_target_features(self, inTargetFeatures, inTargetFlags):
        self.pushFilter.set_target_features(inTargetFeatures, inTargetFlags)
        self.pullFilter.set_floating_features(inTargetFeatures, inTargetFlags)
        self._outputUpToDate = False
        
    def _construct_affinity(self):
        # Compute affinities in both directions
        self.pushFilter._construct_affinity()
        self.pullFilter._construct_affinity()
        
        # Retrieve the affinities so that we can merge them
        pushAffinity = self.pushFilter._get_affinity()
        pullAffinity = self.pullFilter._get_affinity()
        
        # Merge the affinities
        self._affinity = fuse_affinities(pushAffinity, pullAffinity)
        
    
    
def inlier_detection(features, correspondingFeatures, correspondingFlags, oldProbability, kappaa = 3.0):
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
    -inlierProbability:
    Is processed here and returned when the function ends.

    # PARAMETERS
    -kappaa(=3): Mahalanobis distance that determines cut-off in- vs outliers

    # OUTPUTS

    # RETURNS
    -inlierProbability
    """
    # Info & Initialization
    numElements = features.shape[0]
    # Flag based inlier/outlier classification
    ## Do an element-wise multiplication
    newProbability = correspondingFlags.astype(dtype = float, copy = True)

    # Distance based inlier/outlier classification
    ## Re-calculate the parameters sigma and lambda
    sigmaa = 0.0
    lambdaa = 0.0
    sigmaNumerator = 0.0
    sigmaDenominator = 0.0
    for i in range(numElements):
        distance = numpy.linalg.norm(correspondingFeatures[i,:] - features[i,:])
        sigmaNumerator = sigmaNumerator + newProbability[i] * distance * distance
        sigmaDenominator = sigmaDenominator + newProbability[i]

    sigmaa = numpy.sqrt(sigmaNumerator/sigmaDenominator)
    lambdaa = 1.0/(numpy.sqrt(2 * 3.14159) * sigmaa) * numpy.exp(-0.5 * kappaa * kappaa)
    ## Recalculate the distance-based probabilities
    for i in range(numElements):
        distance = numpy.linalg.norm(correspondingFeatures[i,:] - features[i,:])
        probability = 1.0/(numpy.sqrt(2 * 3.14159) * sigmaa) * numpy.exp(-0.5 * numpy.square(distance/sigmaa))
        probability = probability / (probability + lambdaa)
        newProbability[i] = newProbability[i] * probability

    #Gradient Based inlier/outlier classification
    for i in range(numElements):
        normal = features[i,3:6]
        correspondingNormal = correspondingFeatures[i,3:6]
        ## Dot product gives an idea of how well they point in the same
        ## direction. This gives a weight between -1.0 and +1.0
        dotProduct = normal.dot(correspondingNormal)
        ## Rescale this result so that it's continuous between 0.0 and +1.0
        probability = dotProduct / 2.0 + 0.5
        newProbability[i] = newProbability[i] * probability

    return newProbability

class InlierFilter(object):
    """
    This filter object determines inlier weights
    """
    def __init__(self, inFeatures, inCorrespondingFeatures, inCorrespondingFlags,
                 outWeights, paramKappaa = 3.0):
        self.inFeatures = inFeatures
        self.inCorrespondingFeatures = inCorrespondingFeatures
        self.inCorrespondingFlags = inCorrespondingFlags
        self.outWeights = outWeights
        self.paramKappaa = paramKappaa
        
    def update(self):
        # Compute new weights
        newWeights = inlier_detection(self.inFeatures, self.inCorrespondingFeatures,
                                      self.inCorrespondingFlags, self.outWeights,
                                      self.paramKappaa)
        # Copy new weights into the output
        for i in range(len(self.outWeights)):
            self.outWeights[i] = newWeights[i]
    
def rigid_transformation(floatingFeatures, correspondingFeatures, floatingWeights, allowScaling = False):
    """
    # GOAL
    This function computes the rigid transformation between a set a features and
    a set of corresponding features. Each correspondence can be weighed between
    0.0 and 1.0.
    The features are automatically transformed, and the function returns the
    transformation matrix that was used.

    # INPUTS
    -floatingFeatures
    -correspondingFeatures
    -floatingWeights

    # PARAMETERS
    -allowScaling:
    Whether or not to allow scaling.

    # RETURNS
    -transformationMatrix
    """
    # Info and initialization
    numElements = floatingFeatures.shape[0]
    numFeatures = floatingFeatures.shape[1]
    floatingFeaturesT = [] #numFeatures x numElements matrix
    correspondingFeaturesT = [] #numFeatures x numElements matrix
    # Tranpose the data if necessary
    if numElements > numFeatures: #this should normally be the case
        if numFeatures < 6:
            print('Warning: small set of features? Transposition during computation of the transformation might have gone wrong.')
        floatingFeaturesT = floatingFeatures.transpose()
        correspondingFeaturesT = correspondingFeatures.transpose()
    else:
        print('Warning: input of rigid transformation expects rows to correspond with elements, not features, and to have more elements than features per element.')

    floatingCentroid = numpy.zeros((3), dtype = float)
    correspondingCentroid = numpy.zeros((3), dtype = float)
    weightSum = 0.0
    for i in range(numElements):
        floatingCentroid = floatingCentroid + floatingWeights[i] * floatingFeaturesT[0:3,i]
        correspondingCentroid = correspondingCentroid + floatingWeights[i] * correspondingFeaturesT[0:3,i]
        weightSum = weightSum + floatingWeights[i]

    floatingCentroid = floatingCentroid / weightSum
    correspondingCentroid = correspondingCentroid / weightSum
    ###3.2 Compute the Cross Variance matrix
    crossVarianceMatrix = numpy.zeros((3,3), dtype = float)
    for i in range(numElements):
        #CrossVarMat = sum(weight[i] * floatingPosition[i] * correspondingPosition[i]_Transposed)
        crossVarianceMatrix = crossVarianceMatrix + floatingWeights[i] * numpy.outer(floatingFeaturesT[0:3,i], correspondingFeaturesT[0:3,i])

    crossVarianceMatrix = crossVarianceMatrix / weightSum - numpy.outer(floatingCentroid, correspondingCentroid)
    ###3.3 Compute the Anti-Symmetric matrix
    antiSymmetricMatrix = crossVarianceMatrix - crossVarianceMatrix.transpose()
    ###3.4 Use the cyclic elements of the Anti-Symmetric matrix to construct delta
    delta = numpy.zeros(3)
    delta[0] = antiSymmetricMatrix[1,2];
    delta[1] = antiSymmetricMatrix[2,0];
    delta[2] = antiSymmetricMatrix[0,1];
    ###3.5 Compute Q
    Q = numpy.zeros((4,4), dtype = float)
    Q[0,0] = numpy.trace(crossVarianceMatrix)
    ####take care with transposition of delta here. Depending on your library, it might be the opposite than this
    Q[1:4,0] = delta
    Q[0,1:4] = delta #in some libraries, you might have to transpose delta here
    Q[1:4,1:4] = crossVarianceMatrix + crossVarianceMatrix.transpose() - numpy.identity(3, dtype = float)*numpy.trace(crossVarianceMatrix)
    ###3.6 we now compute the rotation quaternion by finding the eigenvector of
    #Q of its largest eigenvalue
    eigenValues, eigenVectors = numpy.linalg.eig(Q)
    indexMaxVal = 0
    maxVal = 0
    for i, eigenValue in enumerate(eigenValues):
        if abs(eigenValue) > maxVal:
            maxVal = abs(eigenValue)
            indexMaxVal = i

    rotQ = eigenVectors[:,indexMaxVal]
    ###3.7 Now we can construct the rotation matrix
    rotMatTemp = numpy.zeros((3,3), dtype = float)
    ####diagonal elements
    rotMatTemp[0,0] = numpy.square(rotQ[0]) + numpy.square(rotQ[1]) - numpy.square(rotQ[2]) - numpy.square(rotQ[3])
    rotMatTemp[1,1] = numpy.square(rotQ[0]) + numpy.square(rotQ[2]) - numpy.square(rotQ[1]) - numpy.square(rotQ[3])
    rotMatTemp[2,2] = numpy.square(rotQ[0]) + numpy.square(rotQ[3]) - numpy.square(rotQ[1]) - numpy.square(rotQ[2])
    ####remaining elements
    rotMatTemp[1,0] = 2.0 * (rotQ[1] * rotQ[2] + rotQ[0] * rotQ[3])
    rotMatTemp[0,1] = 2.0 * (rotQ[1] * rotQ[2] - rotQ[0] * rotQ[3])
    rotMatTemp[2,0] = 2.0 * (rotQ[1] * rotQ[3] - rotQ[0] * rotQ[2])
    rotMatTemp[0,2] = 2.0 * (rotQ[1] * rotQ[3] + rotQ[0] * rotQ[2])
    rotMatTemp[2,1] = 2.0 * (rotQ[2] * rotQ[3] + rotQ[0] * rotQ[1])
    rotMatTemp[1,2] = 2.0 * (rotQ[2] * rotQ[3] - rotQ[0] * rotQ[1])
    ### 3.8 Adjust the scale (optional)
    scaleFactor = 1.0 #>1 to grow ; <1 to shrink
    if allowScaling:
        numerator = 0.0
        denominator = 0.0
        for i in range(numElements):
            centeredFloatingPosition = rotMatTemp.dot(floatingFeaturesT[0:3,i] - floatingCentroid)
            centeredCorrespondingPosition = correspondingFeaturesT[0:3,i] - correspondingCentroid
            numerator = numerator + floatingWeights[i] * numpy.dot(centeredCorrespondingPosition, centeredFloatingPosition)
            denominator = denominator + floatingWeights[i] * numpy.dot(centeredFloatingPosition, centeredFloatingPosition)
        scaleFactor = numerator / denominator

    ### 3.9 Compute the remaining translation necessary between the centroids
    translation = correspondingCentroid - scaleFactor * rotMatTemp.dot(floatingCentroid)
    ### 3.10 Let's (finally) compute the entire transformation matrix!
    translationMatrix = numpy.identity(4, dtype = float)
    rotationMatrix = numpy.identity(4, dtype = float)
    translationMatrix[0:3,3] = translation
    rotationMatrix[0:3,0:3] = scaleFactor * rotMatTemp
    transformationMatrix = rotationMatrix.dot(translationMatrix)

    ##4) Apply transformation
    vector = numpy.ones((4), dtype = float)
    for i in range(numElements):
        vector[0:3] = floatingFeaturesT[0:3,i].copy()
        floatingFeatures[i,0:3] = transformationMatrix.dot(vector)[0:3]

    return transformationMatrix

def compute_viscoelastic_transformation(currentFloatingPositions, correspondingPositions, floatingWeights, toBeUpdatedDisplacementField, numNeighbourDisplacements, sigmaSmoothing = 1.0, numViscousSmoothingIterations = 1, numElasticSmoothingIterations = 1):
    """
    # GOAL
    This function computes the rigid transformation between a set a features and
    a set of corresponding features. Each correspondence can be weighed between
    0.0 and 1.0.
    The features are automatically transformed, and the function returns the
    transformation matrix that was used.

    # INPUTS
    -floatingPositions
    -correspondingPositions
    -floatingWeights

    # OUTPUTS
    -toBeUpdatedDisplacementField:
    The current displacement field that should be updated.

    # PARAMETERS
    -numNeighbourDisplacements:
    For the regularization, the nearest neighbours for each floating positions
    have to be found. The number should be high enough so that every significant
    contribution (up to a distance of e.g. 3*sigmaSmoothing) is included. But low
    enough to keep computational speed high.
    -sigmaSmoothing:
    The value for sigma of the gaussian used for the regularization.
    -numViscousSmoothingIterations:
    Number of times the viscous deformation is smoothed.
    -numElasticSmoothingIterations:
    Number of times the elastic deformation is smoothed

    # RETURNS
    """
    # Info and Initialization
    numFloatingVertices = currentFloatingPositions.shape[0]

    # Viscous Part
    ## The 'Force Field' is what drives the deformation: the difference between
    ## the floating vertices and their correspondences. By regulating it,
    ## viscous behaviour is obtained.
    ### Compute the Force Field
    unregulatedForceField = correspondingPositions - currentFloatingPositions
    ### Regulate the Force Field (Gaussian smoothing, iteratively)
    regulatedForceField = numpy.zeros((numFloatingVertices,3), dtype = float)
    for i in range(numViscousSmoothingIterations):
        helpers.gaussian_smoothing_vector_field(currentFloatingPositions, unregulatedForceField, regulatedForceField, floatingWeights, numNeighbourDisplacements, sigmaSmoothing)
        unregulatedForceField = regulatedForceField

    # Elastic Part
    ## Add the regulated Force Field to the current Displacement Field that has
    ## to be updated.
    unregulatedDisplacementField = toBeUpdatedDisplacementField + regulatedForceField
    ## Regulate the new Displacement Field (Gaussian smoothing, iteratively)
    for i in range(numElasticSmoothingIterations):
        helpers.gaussian_smoothing_vector_field(currentFloatingPositions, unregulatedDisplacementField, toBeUpdatedDisplacementField, floatingWeights, numNeighbourDisplacements, sigmaSmoothing)
        unregulatedDisplacementField = toBeUpdatedDisplacementField

class TransformationFilter(object):
    
    __metaclass__ = ABCMeta
    
    def __init__(self, ioFloatingFeatures, inCorrespondingFeatures,
                 inFloatingWeights):
        self.ioFloatingFeatures = ioFloatingFeatures
        self.inCorrespondingFeatures = inCorrespondingFeatures 
        self.inFloatingWeights = inFloatingWeights
    
    @abstractmethod
    def update(self):
        pass
    
class RigidTransformationFilter(TransformationFilter):
    
    def update(self):
        rigid_transformation(self.ioFloatingFeatures, self.inCorrespondingFeatures,
                             self.inFloatingWeights, allowScaling = False)
        
class ViscoElasticFilter(TransformationFilter):
    
    def __init__(self, ioFloatingFeatures, inCorrespondingFeatures,
                 inFloatingWeights, paramNumNeighbours, paramSigma = 1.0,
                 numViscousIterations = 1, numElasticIterations = 1):
        self.ioFloatingFeatures = ioFloatingFeatures
        self.inCorrespondingFeatures = inCorrespondingFeatures 
        self.inFloatingWeights = inFloatingWeights
        self.paramNumNeighbours = paramNumNeighbours
        self.paramSigma = paramSigma
        self.numViscousIterations = numViscousIterations
        self.numElasticIterations = numElasticIterations
        self._numFloatingElements = self.ioFloatingFeatures.shape[0]
        self.ioDisplacementField = numpy.zeros([self._numFloatingElements,3],
                                               dtype = float)
        self._oldDisplacementField = None
        self._vectorFieldSmoother = helpers.VectorFieldSmoother(self.ioFloatingFeatures[:,0:3],
                                                                self.ioDisplacementField,
                                                                self.inFloatingWeights,
                                                                self.ioDisplacementField,
                                                                self.paramNumNeighbours,
                                                                self.paramSigma)
        
    
    def set_parameters(self, numNeighbours = 10, sigma = 1.0,
                       numViscousIterations = 1, numElasticIterations = 1):
        self.paramNumNeighbours = numNeighbours
        self.paramSigma = sigma
        self.numViscousIterations = numViscousIterations
        self.numElasticIterations = numElasticIterations
        
    def _update_transformation(self):
        # Viscous Part
        ## The 'Force Field' is what drives the deformation: the difference between
        ## the floating vertices and their correspondences. By regulating it,
        ## viscous behaviour is obtained.
        ## 1. Compute the Force Field
        unregulatedForceField = self.inCorrespondingFeatures[:,0:3] - self.ioFloatingFeatures[:,0:3]
        
        ## 2. Set up the vector field smoother
        regulatedForceField = numpy.zeros(unregulatedForceField.shape, dtype = float)
        
        self._vectorFieldSmoother.set_output(regulatedForceField)
        
        ## 3. Use the Vector Field Smoother iteratively to regulate the Force Field
        for i in range(self.numViscousIterations):
            ### Smooth unregulatedForceField (input) and save result in regulatedForceField (output)
            self._vectorFieldSmoother.set_input(self.ioFloatingFeatures[:,0:3],
                                            unregulatedForceField,
                                            self.inFloatingWeights)
            self._vectorFieldSmoother.update()
            ### Copy output back into the input for next iteration
            unregulatedForceField = numpy.copy(regulatedForceField)
    
            
        # Elastic Part
        ## Before computing the new displacement field, save the old one (which
        ## we will need later to apply the difference between the new and the
        ## old displacement field to the floating points)
        self._oldDisplacementField = numpy.copy(self.ioDisplacementField)
        ## 1. Add the regulated Force Field to the current Displacement Field
        ## that has to be updated.
        unregulatedDisplacementField = self.ioDisplacementField + regulatedForceField
        
        ## 2. Set up the Vector Field Smoother
        self._vectorFieldSmoother.set_input(self.ioFloatingFeatures[:,0:3],
                                            unregulatedDisplacementField,
                                            self.inFloatingWeights)
        self._vectorFieldSmoother.set_output(self.ioDisplacementField)
        
        ## 3. Use the Smoother iteratively to regulate the new Displacement Field
        for i in range(self.numElasticIterations):
            ### Smooth unregulatedDisplacementField (input) and save result in self.ioDisplacementField (output)
            self._vectorFieldSmoother.update()
            ### Copy output back into the input for next iteration
            unregulatedDisplacementField = numpy.copy(self.ioDisplacementField)
        
    def _apply_transformation(self):        
        # Apply the differential displacement field to the current floating
        # positions.
        for i in range(self._numFloatingElements):
            self.ioFloatingFeatures[i,0:3] = self.ioFloatingFeatures[i,0:3] - \
                                            self._oldDisplacementField[i,:] + \
                                            self.ioDisplacementField[i,:]

        # Release the memory for the old displacement field
        self._oldDisplacementField = []

    def update(self):
        self._update_transformation()
        self._apply_transformation()
        
        
        
