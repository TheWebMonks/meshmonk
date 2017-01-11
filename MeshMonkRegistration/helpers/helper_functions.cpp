
#include <helper_functions.hpp>

namespace registration {

void normalize_sparse_matrix(SparseMat &ioMat) {
    /*
    # GOAL
    Normalize each row of the inputted sparse matrix.

    # INPUTS
    -ioMat

    # PARAMETERS

    # OUTPUT
    -ioMat
    */

    //# Normalize the rows of the affinity matrix
    //## Initialize a vector that will contain the sum of each row
    const size_t numRows = ioMat.rows();
    std::vector<float> sumRows(numRows, 0.0f);

    //## Loop over the rows and compute the total sum of its elements
    for (size_t i = 0 ; i < ioMat.outerSize() ; i++) {
        for (SparseMat::InnerIterator innerIt(ioMat,i) ; innerIt ; ++innerIt) {
            //### Get index of current row we're in
            const unsigned int currentRowIndex = innerIt.row();
            //### Add current element to the sum of this row.
            const float currentElement = innerIt.value();
            sumRows[currentRowIndex] += currentElement;
        }
    }

    //## Loop over the rows and divide each row by the sum of its elements
    for (size_t i = 0 ; i < ioMat.outerSize() ; i++) {
        for (SparseMat::InnerIterator innerIt(ioMat,i) ; innerIt ; ++innerIt) {
            //### Get index and sum of the current row we're in
            const unsigned int currentRowIndex = innerIt.row();
            const float currentRowSum = sumRows[currentRowIndex];

            //### Get current element and divide it by the sum of its row.
            const float currentElement = innerIt.value();
            innerIt.valueRef() = currentElement / currentRowSum;
        }
    }
}//end normalize_sparse_matrix()


template <typename VecMatType>
void radius_nearest_neighbours(const VecMatType &inQueriedPoints,
                                const VecMatType &inSourcePoints,
                                MatDynInt &outNeighbourIndices,
                                MatDynFloat &outNeighbourSquaredDistances,
                                const float paramRadius,
                                const size_t paramLeafsize){
    /*
    GOAL
    This function searches in a radius for the nearest neighbours in
    'inSourcePoints' for each element in the 'inQueriedPoints' set. It outputs
    the indices of each neighbour and the squared (!) distances between each
    element of 'inQueriedPoints' and its neighbours.

    INPUT
    -inQueriedPoints:
    -inSourcePoints:

    PARAMETERS
    -paramRadius(= 3.0): the radius in which neighbours are queried
    -paramLeafsize(= 15): should be between 5 and 50 or so

    OUTPUT
    -outNeighbourIndices
    -outNeighbourSquaredDistances.
    NOTE: Not for every queried elements, the same number of neighbours will be
    found. So, 'empty' elements in the output are represented by '-1'.

    */

    //# Info and Initialization
    const size_t dimension = inSourcePoints.cols();
    typedef nanoflann::KDTreeEigenMatrixAdaptor<VecMatType>  NanoKDTree;
    const size_t numSourceElements = inSourcePoints.rows();
    const size_t numQueriedElements = inQueriedPoints.rows();


    //# Construct kd-tree
    NanoKDTree kdTree(dimension, inSourcePoints, paramLeafsize);
    kdTree.index->buildIndex();

    //# Query the kd-tree
    //## Loop over the queried features
    //### Initialize variables we'll need during the loop
    std::vector<float> queriedFeature(dimension);
    std::vector< std::vector<size_t> > neighbourIndices;
    std::vector< std::vector<float> > neighbourSquaredDistances;
    size_t maxNumNeighboursFound = 0;

    //### Execute loop
    for (size_t i = 0 ; i < numQueriedElements ; ++i ) {
        //### convert input features to 'queriedFeature' std::vector structure
        //### (required by nanoflann's kd-tree).
        for (size_t j = 0 ; j < dimension ; ++j) {
            queriedFeature[j] = inQueriedPoints(i,j);
        }

        //### Query the kd-tree
        nanoflann::SearchParams searchParams;
        std::vector<std::pair<size_t, float> > resultIndexAndDistancePairs;
        nanoflann::RadiusResultSet<float,size_t> resultSet(paramRadius,resultIndexAndDistancePairs);
        const size_t nMatches = kdTree.index->radiusSearchCustomCallback(&queriedFeature[0],resultSet,searchParams);
        std::cout << "Number of matches: " << nMatches << std::endl;

        //### Copy the result into neighbourIndices and neighbourSquaredDistances
        //### by first copying it into single std::vector<float/int> and then
        //### copying that std::vector into those std::vector<std::vector<>>.
        std::vector<size_t> indices;
        std::vector<float> squaredDistances;
        for (size_t j = 0 ; j < nMatches ; ++j) {
            indices.push_back(resultIndexAndDistancePairs[j].first);
            squaredDistances.push_back(resultIndexAndDistancePairs[j].second);
        }
        neighbourIndices.push_back(indices);
        neighbourSquaredDistances.push_back(squaredDistances);


        //### Check if we found more neighbours than before
        if (nMatches > maxNumNeighboursFound) {
            maxNumNeighboursFound = nMatches;
        }
    }

    //# Copy results into output
    //## Now that we have found all the neighbours and know what the maximum
    //## number of neighbours found was, we can finally initialize the output
    //## matrix size.
    outNeighbourIndices = -1 * MatDynInt::Ones(numQueriedElements, maxNumNeighboursFound);
    outNeighbourSquaredDistances = -1.0 * MatDynFloat::Ones(numQueriedElements, maxNumNeighboursFound);

    //## Loop over results
    for (size_t i = 0 ; i < numQueriedElements ; i++ ) {
        int numNeighbours = neighbourIndices[i].size();

        //### Loop over neighbours
        for (size_t j = 0 ; j < numNeighbours ; j++ ) {
            outNeighbourIndices(i,j) = neighbourIndices[i][j];
            outNeighbourSquaredDistances(i,j) = neighbourSquaredDistances[i][j];
        }
    }

}//end k_nearest_neighbours()








void openmesh_to_eigen_features(const TriMesh &inputMesh,
                                FeatureMat &outputFeatures){
    /*
    GOAL
    This function converts an openmesh data structure to the feature
    representation in eigen matrices arrays required by the registration
    framework.

    INPUT
    -inputMesh:
    this has to be a mesh of openmesh's TriMesh type

    PARAMETERS

    OUTPUT
    -outputFeatures:
    a numVertices x 6 Eigen dense matrixwhere the first three columns are made
    up of the positions of the vertices, and the last three columns the normals
    of those vertices.
    */

    //# Info and Initialization
    const int numVertices = inputMesh.n_vertices();
    outputFeatures = FeatureMat::Zero(numVertices,NUM_FEATURES);

    //# Extract the vertex positions and normals
    unsigned int i = 0;
    TriMesh::Point position(0.0f,0.0f,0.0f);
    TriMesh::Normal normal(0.0f,0.0f,0.0f);
    TriMesh::VertexIter vertexIt(inputMesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(inputMesh.vertices_end());
    for ( ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
        //## Get position and normal
        position = inputMesh.point(vertexIt);
        normal = inputMesh.normal(vertexIt);
        //## Insert position and normal into 'outputFeatures'
        outputFeatures(i,0) = position[0];
        outputFeatures(i,1) = position[1];
        outputFeatures(i,2) = position[2];
        outputFeatures(i,3) = normal[0];
        outputFeatures(i,4) = normal[1];
        outputFeatures(i,5) = normal[2];
    }
}


void load_obj_to_eigen_features(const std::string inObjFilename,
                                TriMesh &outMesh,
                                FeatureMat &outFeatureMatrix){
    /*
    GOAL
    This function loads an OBJ file given by a filename and converts it to
    OpenMesh's TriMesh type and feature matrices we use internally.

    INPUT
    -inObjFilename:

    PARAMETERS

    OUTPUT
    -outMesh
    -outFeatureMatrix
    */

    OpenMesh::IO::_OBJReader_(); //NOTE: this issue should only appear when STATICALLY linking OpenMesh
    if (!OpenMesh::IO::read_mesh(outMesh,inObjFilename)){
        std::cerr << "Read error \n";
        exit(1);
    };

    // Update Normals
    outMesh.request_face_normals();
    outMesh.request_vertex_normals();
    outMesh.update_face_normals();
    outMesh.update_vertex_normals();

    // Convert to Eigen Matrices
    openmesh_to_eigen_features(outMesh, outFeatureMatrix);

}


void eigen_features_to_openmesh(const FeatureMat &inFeatures,
                                TriMesh &outMesh){
    /*
    GOAL
    This function converts a feature representation of a mesh using Eigen
    matrices into a mesh representation using OpenMesh.

    INPUT
    -inFeatures:
    a numVertices x 6 Eigen dense matrix where the first three columns are made
    up of the positions of the vertices, and the last three columns the normals
    of those vertices.

    PARAMETERS

    OUTPUT
    -outMesh:
    this has to be a mesh of openmesh's TriMesh type
    */

    //# Info and Initialization
    const int numRows = inFeatures.rows();
    const int numVertices = outMesh.n_vertices();
     if (numRows != numVertices) {
    std::cerr << "Number of rows does not correspond with number of vertices when"
    << " calling eigen_features_to_openmesh()" << std::endl;
    }

    //# Put positions and normals back into mesh
    TriMesh::Point updatedPosition(0.0f,0.0f,0.0f);
    TriMesh::Normal updatedNormal(0.0f,0.0f,0.0f);
    unsigned int i = 0;
    TriMesh::VertexIter vertexIt(outMesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(outMesh.vertices_end());
    for ( ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
        //## Get positions and normals from the feature matrix
        updatedPosition[0] = inFeatures(i,0);
        updatedPosition[1] = inFeatures(i,1);
        updatedPosition[2] = inFeatures(i,2);
        updatedNormal[0] = inFeatures(i,3);
        updatedNormal[1] = inFeatures(i,4);
        updatedNormal[2] = inFeatures(i,5);
        //## Insert position and normal into mesh
        outMesh.set_point(vertexIt,updatedPosition);
        outMesh.set_normal(vertexIt, updatedNormal);
    }
}


void write_eigen_features_to_obj(const FeatureMat &inFeatures,
                                TriMesh &inMesh,
                                const std::string inObjFilename){
    /*
    GOAL
    This function converts a feature representation of a mesh using Eigen
    matrices into a mesh representation using OpenMesh and then writes the
    result to an obj file.

    INPUT
    -inFeatures:
    a numVertices x 6 Eigen dense matrix where the first three columns are made
    up of the positions of the vertices, and the last three columns the normals
    of those vertices.
    -inMesh
    -inObjFilename

    PARAMETERS

    OUTPUT

    */

    //# Info and Initialization
    const int numRows = inFeatures.rows();
    const int numVertices = inMesh.n_vertices();
     if (numRows != numVertices) {
    std::cerr << "Number of rows does not correspond with number of vertices when"
    << " calling write_eigen_features_to_obj()" << std::endl;
    }

    //# Convert feature matrix to OpenMesh structure
    eigen_features_to_openmesh(inFeatures, inMesh);

    //# Write resulting mesh to OBJ file
    OpenMesh::IO::_OBJWriter_();
    if (!OpenMesh::IO::write_mesh(inMesh, inObjFilename))
    {
        std::cerr << "write error\n";
        exit(1);
    }
}//end write_eigen_features_to_obj()


void update_normals_for_altered_positions(TriMesh &ioMesh,
                                        FeatureMat &ioFeatures){
    /*
    GOAL
    This function takes the positions that are in the first three columns of
    'ioFeatures', inserts them into the 'ioMesh', and uses OpenMesh's
    functionality to recompute the vertex normals given those new vertex
    positions.
    Note(!): the function arguments 'ioMesh' and 'ioFeatures' are both input and
    output here! That's because part of each is required as necessary input, yet
    also another part of each is updated.

    INPUT
    -ioMesh:
    A TriMesh-type mesh. This is a required input because the ioMesh's structure
    (the way vertices are connected to each other) is needed to recompute
    the normals for new vertex positions given by the matrix 'ioFeatures'.
    -ioFeatures:
    a numVertices x 6 Eigen dense matrix. The first three columns give the
    required input: updated vertex positions.

    PARAMETERS

    OUTPUT
    -ioMesh
    A TriMesh-type mesh. Both its positions and normals are updated from the
    positions given by the 'ioFeatures' matrix.
    -ioFeatures:
    a numVertices x 6 Eigen dense matrix. The last three columns contain the
    vertex normals and are updated in this function call.
    */

    //# Info and Initialization
    const int numVertices = ioMesh.n_vertices();
    const int numRows = ioFeatures.rows();
    if (numRows != numVertices) {
    std::cerr << "Number of rows does not correspond with number of vertices when"
    << " calling eigen_features_to_openmesh()" << std::endl;
    }

    //# Insert positions from 'ioFeatures' into 'ioMesh'
    eigen_features_to_openmesh(ioFeatures, ioMesh);

    //# Given the new positions, let ioMesh update its normals
    ioMesh.request_face_normals();
    ioMesh.request_vertex_normals();
    ioMesh.update_normals();
    ioMesh.release_face_normals();

    //# Copy ioMesh's new normals into 'ioFeatures' last three columns.
    unsigned int i = 0;
    TriMesh::Normal updatedNormal(0.0f,0.0f,0.0f);
    TriMesh::VertexIter vertexIt(ioMesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(ioMesh.vertices_end());
    for ( ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
        //## Get normal
        updatedNormal = ioMesh.normal(vertexIt);
        //## Insert normal into ioFeatures
        ioFeatures(i,3) = updatedNormal[0];
        ioFeatures(i,4) = updatedNormal[1];
        ioFeatures(i,5) = updatedNormal[2];
    }

}

}//namespace registration
