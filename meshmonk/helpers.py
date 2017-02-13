import sys
# sys.path.append('/home/jonatan/projects/OpenMesh/build/Build/python')
import openmesh
#Importing the rest of utilities
import numpy
from numpy import linalg
from scipy import spatial
from scipy.interpolate import interp1d, LinearNDInterpolator







def nearest_neighbours(features1, features2, k = 3, leafsize = 15, eps = 0.0001, p = 2, maxDistance = 1000):
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

    RETURNS
    -distances
    -neighbourIndices.
    """
    # Obtain info and initialize outputs
    kdTree = spatial.cKDTree(features2, leafsize)
    # Query the kd-tree (http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query)
    distances, neighbourIndices = kdTree.query(features1, k, eps, p, maxDistance)

    return distances, neighbourIndices

def weighted_average_vectors(self, vectors, weights):
    """
    GOAL
    This function does does a weighted average of input vectors and weights

    INPUT
    -vectors:
    -weights:
    Contains a list of weights that correspond to each vector.

    PARAMETERS

    RETURNS
    -averagedVector
    """
    # Compute weighted sum of vectors
    ## Loop over vectors and weights, weigh each vector, sum them, and sum the individual weights
    sumVector = numpy.zeros((1,vectors.shape[1]), dtype = float)
    sumWeights = 0.0
    for i in range(vectors.shape[0]):
        weight = weights[i]
        vector = vectors[i]
        sumVector += weight * vector
        sumWeights += weight
    
    ## Divide the sum of weighted vectors by the sum of weights
    averagedVector = sumVector / sumWeights
    return averagedVector
    
def gaussian_vector_interpolation(position, fieldPositions, fieldVectors, fieldWeights, sigmaa):
    """
    GOAL
    This function does 'Gaussian interpolation' of vectors in an irregular
    vector field to estimate the vector of a certain position.

    INPUT
    -position:
    The postion where you want to calculate the vector value at
    -fieldPositions:
    Positions of the vectors of the vector field. The dimensions of this matrix
    are numVectors x 3
    -fieldVectors:
    Values of the vectors of the vector field. The dimensions of this matrix
    are numVectors x numFeatures
    -fieldWeights:
    The contribution of each vector field vector can be weighed.

    PARAMETERS
    -sigmaa:
    The sigma used in the Gaussian averaging.

    RETURNS
    -interpolatedVector:
    Computed through Gaussian averaging, taking additional weights (see
    fieldWeights) into account.
    """
    # Obtain info and initialize outputs
    numVectors = fieldVectors.shape[0]
    numFeatures = fieldVectors.shape[1]
    interpolatedVector = numpy.zeros((numFeatures), dtype = float)

    # Loop over the field vectors
    sumWeights = 0.0
    for i, vector in enumerate(fieldVectors):
        ## Compute distance
        distanceVector = fieldPositions[i] - position
        distance = linalg.norm(distanceVector)
        ## From the distance, we can compute the gaussian weight
        gaussianWeight = numpy.exp(-0.5 * numpy.square(distance / sigmaa))
        ## Let's take the user-defined weight into account
        combinedWeight = gaussianWeight * fieldWeights[i]
        ## Weigh the current vector field vector and sum it
        interpolatedVector = interpolatedVector + combinedWeight * vector
        ## Add the weight to the total sum of weights for normalization after.
        sumWeights = sumWeights + combinedWeight

    # Normalize the vector
    interpolatedVector = interpolatedVector / sumWeights

    return interpolatedVector


def gaussian_scalar_interpolation(position, fieldPositions, fieldScalars, fieldWeights, sigmaa):
    """
    GOAL
    This function does 'Gaussian interpolation' of scalars in an irregular
    scalar field to estimate the scalar of a certain position.

    INPUT
    -position:
    The postion where you want to calculate the scalar value at
    -fieldPositions:
    Positions of the nodes of the scalar field. The dimensions of this matrix
    are numNodes x 3
    -fieldScalars:
    Scalar values of the nodes of the scalar field. The dimensions of this array
    are numNodes x 1
    -fieldWeights:
    The contribution of each scalar field node can be weighed.

    PARAMETERS
    -sigmaa:
    The sigma used in the Gaussian averaging.

    RETURNS
    -interpolatedScalar:
    Computed through Gaussian averaging, taking additional weights (see
    fieldWeights) into account.
    """
    # Obtain info and initialize outputs
    numNodes = fieldScalars.shape[0]
    interpolatedScalar = 0.0

    # Loop over the field nodes
    sumWeights = 0.0
    for i, scalar in enumerate(fieldScalars):
        ## Compute distance
        distanceVector = fieldPositions[i] - position
        distance = linalg.norm(distanceVector)
        ## From the distance, we can compute the gaussian weight
        gaussianWeight = numpy.exp(-0.5 * numpy.square(distance / sigmaa))
        ## Let's take the user-defined weight into account
        combinedWeight = gaussianWeight * fieldWeights[i]
        ## Weigh the current node's scalar and sum it
        interpolatedScalar = interpolatedScalar + combinedWeight * scalar
        ## Add the weight to the total sum of weights for normalization after.
        sumWeights = sumWeights + combinedWeight

    # Normalize the vector
    interpolatedScalar = interpolatedScalar / sumWeights

    return interpolatedScalar

def gaussian_smoothing_vector_field(fieldPositions, fieldVectors, regulatedFieldVectors, fieldWeights, numNeighbours, gaussianSigma):
    """
    GOAL
    This function performs gaussian smoothing on an entire vector field.

    INPUT
    -fieldPositions
    -fieldVectors:
    The vector field that should be smoothed.
    -fieldWeights

    OUTPUT
    -regulatedFieldVectors:
    The resulting regulated displacement field.


    PARAMETERS
    -numNeighbours:
    For the smoothing, the nearest neighbours for each floating positions have
    to be found. The number should be high enough so that every significant
    contribution (up to a distance of e.g. 3*gaussianSigma) is included. But low
    enough to keep computational speed high.
    -gaussianSigma:
    The value for sigma of the gaussian used for the smoothing.

    RETURNS
    """
    # Info & Initialization
    numNodes = fieldPositions.shape[0]
    # Determine for each field node the (closely) neighbouring nodes
    distances, neighbourIndices = nearest_neighbours(fieldPositions, fieldPositions, numNeighbours, 15, 0.0001, 2, 1000.0)
    # Use the neighbouring field vectors to smooth each individual field vector
    for i in range(numNodes):
        position = fieldPositions[i,:]
        neighbourPositions = numpy.zeros((numNeighbours,3), dtype = float)
        neighbourVectors = numpy.zeros((numNeighbours,3), dtype = float)
        neighbourWeights = numpy.zeros((numNeighbours), dtype = float)
        ## For the current displacement, get all the neighbouring positions,
        ## vectors, and weights (needed for Gaussian smoothing!).
        for j in range(numNeighbours):
            neighbourIndex = neighbourIndices[i,j]
            neighbourPositions[j,:] = fieldPositions[neighbourIndex,:]
            neighbourVectors[j,:] = fieldVectors[neighbourIndex,:]
            neighbourWeights[j] = fieldWeights[neighbourIndex]

        ## Gaussian averaging of neighbouring displacements
        smoothedVector = gaussian_vector_interpolation(position, neighbourPositions, neighbourVectors, neighbourWeights, gaussianSigma)
        regulatedFieldVectors[i,:] = smoothedVector

class GaussianEvaluator(object):
    """
    This class evaluates a gaussian through 1D interpolation of a few gaussian
    node-value pairs. This is meant to evaluate a Gaussian cheaper than 
    computing 1/(2Pi*sigma) * exp(...).
    """
    
    def __init__(self, sigma):
        self.sigma = sigma
        self.xValues = sigma * numpy.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1000.0], dtype = float)
        self.numValues = self.xValues.size
        self.yValues = numpy.zeros((self.numValues), dtype = float)
        self.interpolator = None
        self._construct_interpolator()
        
    def _construct_interpolator(self):  
        # Determine the peak value of the gaussian for the given sigma:
        peak = 1.0/(numpy.sqrt(2.0 * 3.14159) * self.sigma)
        # Now determine the gaussian values for each x-value
        for i in range(self.numValues):
            xVal = self.xValues[i]
            self.yValues[i] = peak * numpy.exp(-(xVal**2/(2*self.sigma**2)))
            
        # Build interpolator using these discrete value pairs
        ## Use scipy's "interp1d" function
        self.interpolator = interp1d(self.xValues, self.yValues, kind = 'linear')
    
    def set_params(self, sigma):
        self.sigma = sigma
        self._construct_interpolator()
        
    def compute(self, xValue):
        # We interpolate the gaussian value for the given distance
        interpolatedValue = self.interpolator(xValue)
        
        return interpolatedValue
        


# TODO: WRITE FUNCTION THAT TAKES INPUT VECTOR FIELD AND OUTPUTS INTERPOLATED
# VECTOR FIELD.
# THEN: REWRITE IT AS A CLASS

def interpolate_vector_field(inSourcePositions, inSourceVectors,
                             inTargetPositions, ioTargetVectors):
    # Info & Initialization
    numTargetPoints = inTargetPositions.shape[0]
    numDimensions = inSourceVectors.shape[1]
    
    # We have to do the interpolation for each dimension separately!
    for dim in range(numDimensions):
        ## Build interpolator (from scipy package)
        interpolator = LinearNDInterpolator(inSourcePositions, inSourceVectors[:,dim],
                                            fill_value=0.0)
        
        ## For each target positions, use the interpolator
        for i in range(numTargetPoints):
            position = inTargetPositions[i,:]
            ioTargetVectors[i,dim] = interpolator.__call__(position)

def interpolate_scalar_field(inSourcePositions, inSourceScalars,
                             inTargetPositions, ioTargetScalars):
    # Info & Initialization
    numTargetPoints = inTargetPositions.shape[0]
    
    ## Build interpolator (from scipy package)
    interpolator = LinearNDInterpolator(inSourcePositions, inSourceScalars,
                                        fill_value=0.0)
    
    ## For each target positions, use the interpolator
    for i in range(numTargetPoints):
        position = inTargetPositions[i,:]
        ioTargetScalars[i] = interpolator.__call__(position)
        
class VectorFieldSmoother(object):
    
    def __init__(self, inPositions, inVectors, inWeights, outVectors, numNeighbours, sigma):
        self.inPositions = inPositions
        self.inVectors = inVectors
        self._numVectors = inVectors.shape[0]
        self.inWeights = inWeights # additional weight assigned to each vector by user (e.g. inlier weight)
        self.outVectors = outVectors # The output containing the smoothed vectors!
        self.numNeighbours = numNeighbours
        self.sigma = sigma
        self._neighbourIndices = None
        self._neighbourIndicesUpToDate = False
        self._neighbourDistances = None
        self._neighbourWeights = numpy.zeros((self._numVectors, self.numNeighbours), dtype = float)
        self._neighbourWeightsUpToDate = False
        self._gaussianEvaluator = GaussianEvaluator(self.sigma)
        
    def set_input(self, inPositions, inVectors, inWeights):
        self.inPositions = inPositions
        self.inVectors = inVectors
        self._numVectors = inVectors.shape[0]
        self.inWeights = inWeights # additional weight assigned to each vector by user (e.g. inlier weight)
        #self._neighbourWeightsUpToDate = False #TODO: add this line or not?
        
    def set_output(self, outVectors):
        self.outVectors = outVectors
        
    def _update_neighbour_indices(self):
        if (self._neighbourIndicesUpToDate == False):
            self._neighbourDistances, self._neighbourIndices =  \
                        nearest_neighbours(self.inPositions,
                                                    self.inPositions,
                                                    self.numNeighbours,
                                                    15, 0.0001, 2, 1000)
        self._neighbourIndicesUpToDate = True
    
    def _update_neighbour_weights(self):
        if (self._neighbourWeightsUpToDate == False):
            # First update the indices and distances of neighbours.
            self._update_neighbour_indices()
            
            # Given the distance to each neighbour, we can compute the gaussian weight
            for i in range(self._numVectors):
                for j in range(self.numNeighbours):
                    distance = self._neighbourDistances[i,j]
                    gaussianWeight = self._gaussianEvaluator.compute(distance)
                    self._neighbourWeights[i,j] = gaussianWeight
        
        self._neighbourWeightsUpToDate = True
        
        # We can free the memory belonging to the neighbour distances. If they
        # don't change, we wouldn't update the weight neither. If they change,
        # they have to be recomputed anyway.
        self._neighbourDistances = []

    def update(self):
        # Update the gaussian weights
        self._update_neighbour_weights()
        
        # Loop over the vectors
        for i in range(self._numVectors):
            ## Initialize output as zero
            self.outVectors[i,:] = numpy.zeros((1,self.inVectors.shape[1]), dtype = float)
            sumWeights = 0.0
            ## For this vector, loop over all its neighbours (to make a
            ## weighted average of the neighbouring vectors)
            for j in range(self.numNeighbours):
                ### Compute the weight of the current neighbour as the combo
                ### of a user defined weight with its gaussian weight.
                neighbourIndex = self._neighbourIndices[i,j]
                combinedWeight = self.inWeights[neighbourIndex] * self._neighbourWeights[i,j]
                ### Weigh the neighbour's vector and add it to the sum
                neighbourVector = self.inVectors[neighbourIndex]
                self.outVectors[i,:] += combinedWeight * neighbourVector
                sumWeights += combinedWeight
            
            ## Divide the sum of vectors by the sum of weights
            self.outVectors[i,:] /= sumWeights

        return self.outVectors
    
    
def openmesh_normals_from_positions(mesh, newPositions):
    """
    GOAL
    This function recalculates the normals for a given mesh, when there are new
    positions for the nodes.

    INPUT
    -mesh:
    this has to be a mesh of openmesh's TriMesh type
    -newPositions:
    the new positions that have to be inserted into the mesh to recompute the
    vertex normals.

    PARAMETERS

    RETURNS
    -newNormals:
    The new vertex normals that were recomputed given the new positions and mesh
    structure.
    """
    # Info and Initialization
    numVertices = newPositions.shape[0]
    newNormals = numpy.zeros((numVertices,3), dtype = float)
    # Put positions back into the floating mesh
    newPosition = openmesh.TriMesh.Point(0.0,0.0,0.0)
    for i, vh in enumerate(mesh.vertices()):
        newPosition[0] = newPositions[i,0]
        newPosition[1] = newPositions[i,1]
        newPosition[2] = newPositions[i,2]
        mesh.set_point(vh,newPosition)

    # Let openmesh recalculate the vertex normals
    mesh.request_vertex_normals()
    mesh.request_face_normals()
    mesh.update_normals()
    mesh.release_face_normals()

    # Get the normals back out from the openmesh TriMesh
    for i, vh in enumerate(mesh.vertices()):
        newNormals[i,0] = mesh.normal(vh)[0]
        newNormals[i,1] = mesh.normal(vh)[1]
        newNormals[i,2] = mesh.normal(vh)[2]

    return newNormals

def openmesh_to_numpy_features(mesh):
    """
    GOAL
    This function converts an openmesh data structure to the feature
    representation in numpy arrays required by the registration framework.

    INPUT
    -mesh:
    this has to be a mesh of openmesh's TriMesh type

    PARAMETERS

    RETURNS
    -features:
    a numVertices x 6 numpy array where the first three columns are made up of
    the positions of the vertices, and the last three columns the normals of
    those vertices.
    """
    # Info and Initialization
    numVertices = mesh.n_vertices()
    features = numpy.zeros((numVertices,6), dtype = float)

    # Let openmesh recalculate the vertex normals (typically, only the face
    # normals are present).
    mesh.request_vertex_normals()
    mesh.request_face_normals()
    mesh.update_normals()
    mesh.release_face_normals()

    # Extract the vertex positions and normals
    for i, vertexHandle in enumerate(mesh.vertices()):
        features[i,0] = mesh.point(vertexHandle)[0]
        features[i,1] = mesh.point(vertexHandle)[1]
        features[i,2] = mesh.point(vertexHandle)[2]
        features[i,3] = mesh.normal(vertexHandle)[0]
        features[i,4] = mesh.normal(vertexHandle)[1]
        features[i,5] = mesh.normal(vertexHandle)[2]

    return features

def openmesh_to_matrices(inMesh):
    """
    GOAL
    This function converts an openmesh trimesh structure to the matrix
    representations used in the registration framework (matrices of vertex
    features and face indices)

    INPUT
    -inMesh:
    this has to be a mesh of openmesh's TriMesh type

    PARAMETERS

    RETURNS
    -outFeatures:
    a numVertices x 6 numpy array where the first three columns are made up of
    the positions of the vertices, and the last three columns the normals of
    those vertices.
    -outIndices:
    a numFaces x 3 numpy array with each element containing the index of a
    vertex belonging to that face.
    """
    # Info and Initialization
    numVertices = inMesh.n_vertices()
    numFaces = inMesh.n_faces()
    outFeatures = numpy.zeros((numVertices,6), dtype = float)
    outIndices = numpy.zeros((numFaces,3), dtype = int)

    # Let openmesh recalculate the vertex normals (typically, only the face
    # normals are present).
    inMesh.request_vertex_normals()
    inMesh.request_face_normals()
    inMesh.update_normals()
    inMesh.release_face_normals()

    # Extract the vertex positions and normals
    i = 0
    for vertexHandle in inMesh.vertices():
        outFeatures[i,0] = inMesh.point(vertexHandle)[0]
        outFeatures[i,1] = inMesh.point(vertexHandle)[1]
        outFeatures[i,2] = inMesh.point(vertexHandle)[2]
        outFeatures[i,3] = inMesh.normal(vertexHandle)[0]
        outFeatures[i,4] = inMesh.normal(vertexHandle)[1]
        outFeatures[i,5] = inMesh.normal(vertexHandle)[2]
        i += 1

    # Extract the face indices
    i = 0
    for faceHandle in inMesh.faces():
        ## Loop over the vertices of the current face
        j = 0
        for vertexHandle in inMesh.fv(faceHandle):
            ## Save the vertex index
            outIndices[i,j] = vertexHandle.idx()
            j += 1

    return outFeatures, outIndices

def matrices_to_openmesh(inVertices, inFaces):
    """
    GOAL
    This function converts matrices into an OpenMesh TriMesh structure.

    INPUT
    -inVertices:
    A matrix containing the x,y,z coordinates for each vertex
    -inFaces:
    A matrix containing the indices for all the nodes of each face

    PARAMETERS

    RETURNS
    -an OpenMesh Trimesh
    """
    # Info and Initialization
    numVertices = inVertices.shape[0]
    numFaces = inFaces.shape[0]
    outMesh = openmesh.TriMesh()
    
    # Add each vertex to the mesh
    vertexHandles = [None] * numVertices
    for i in range(numVertices):
        ## Add current vertex
        vertexHandle = outMesh.add_vertex(openmesh.TriMesh.Point(inVertices[i,0], inVertices[i,1], inVertices[i,2]))
        ## Save the vertex handle that was returned, since we will need it
        ## to construct the faces of the mesh
        vertexHandles[i] = vertexHandle
    
    # Add each face to the mesh
    for i in range(numFaces):
        ## Get vertex indices of current face
        vertexIndex0 = inFaces[i,0]
        vertexIndex1 = inFaces[i,1]
        vertexIndex2 = inFaces[i,2]
        ## Get the corresponding vertex handles (generated when adding the 
        ## vertices to the mesh)
        vertexHandle0 = vertexHandles[vertexIndex0]
        vertexHandle1 = vertexHandles[vertexIndex1]
        vertexHandle2 = vertexHandles[vertexIndex2]
        ## Add a face to the mesh using these vertex handles
        faceHandle = outMesh.addface(vertexHandle0,vertexHandle1,vertexHandle2)
    
    # Return the result
    return outMesh

class DataImporter(object):
    """
    This class takes filenames as input and prepares data that is required
    by the registration framework
    """
    def __init__(self, floatingMeshPath, targetMeshPath):
        self.floatingMeshPath = floatingMeshPath
        self.targetMeshPath = targetMeshPath
        self.floatingMesh = openmesh.TriMesh()
        self.targetMesh = openmesh.TriMesh()
        self.floatingFeatures = []
        self.targetFeatures = []
        self.update()
        
    def update(self):
        print("Importing data....")
        # Load from file
        openmesh.read_mesh(self.floatingMesh, self.floatingMeshPath)
        openmesh.read_mesh(self.targetMesh, self.targetMeshPath)
        # Obtain the floating and target mesh features (= positions and normals)
        self.floatingFeatures = openmesh_to_numpy_features(self.floatingMesh)
        self.targetFeatures = openmesh_to_numpy_features(self.targetMesh)
        print("Data imported.")
        
        
class DataExporter(object):
    """
    This class takes filenames as input and prepares data that is required
    by the registration framework
    """
    def __init__(self,resultingFeatures, resultingMesh, resultingMeshPath):
        self.resultingFeatures = resultingFeatures
        self.resultingMesh = resultingMesh
        self.resultingMeshPath = resultingMeshPath
        self.update()
        
    def update(self):
        print('Exporting data...')
        # Insert positions into the mesh structure
        openmesh_normals_from_positions(self.resultingMesh, self.resultingFeatures[:,0:3])
        # Write the mesh to file
        openmesh.write_mesh(self.resultingMesh, self.resultingMeshPath)
        print('Data exported.')