
#include "helper_functions.hpp"

namespace registration {


void fuse_affinities(SparseMat &ioAffinity1,
                    const SparseMat &inAffinity2){
    /*
    # GOAL
    Fuse the two affinity matrices together. The result is normalized
    and inserted into ioAffinity1.

    # INPUTS
    -ioAffinity1
    -inAffinity2: dimensions should be transposed of ioAffinity1

    # PARAMETERS

    # OUTPUT
    -ioAffinity1

    # RETURNS
    """
    */
    //# Info and Initialization
    const size_t numRows1 = ioAffinity1.rows();
    const size_t numCols1 = ioAffinity1.cols();
    const size_t numRows2 = inAffinity2.rows();
    const size_t numCols2 = inAffinity2.cols();
    //## Safety check for input sizes
    if((numRows1 != numCols2) || (numCols1 != numRows2)) {
        std::cerr << "The sizes of the inputted matrices in fuse_affinities are wrong. "
        << "Their sizes should be the transpose of each other!" << std::endl;
        std::cout << " Affinity 1 : num rows - " << numRows1 << " | num cols - " << numCols1 << std::endl;
        std::cout << " Affinity 2 : num rows - " << numRows2 << " | num cols - " << numCols2 << std::endl;
    }

    //# Fusing is done by simple averaging
    //## Sum matrices. Note: because Eigen's Sparse matrices require storage
    //## orders to match (row- or column-major) we need to construct a new
    //## temporary matrix of the tranpose of inAffinity2.
    ioAffinity1 += SparseMat(inAffinity2.transpose());
    //## Normalize result
    normalize_sparse_matrix(ioAffinity1);

}


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





void convert_mesh_to_matrices(const TriMesh &inMesh,
                                FeatureMat &outFeatures,
                                FacesMat &outFaces){
    /*
    GOAL
    This function converts an openmesh data structure to a feature matrix
    and a matrix containing the indices of the vertices belonging to each face

    INPUT
    -inMesh:
    this has to be a mesh of openmesh's TriMesh type

    PARAMETERS

    OUTPUT
    -outFeatures:
    a numVertices x 6 Eigen dense matrix where the first three columns are made
    up of the positions of the vertices, and the last three columns the normals
    of those vertices.
    -outFaces:
    a numFaces x 3 Eigen dense matrix where each row contains the corresponding
    indices of the vertices belonging to that face.
    */

    //# Info and Initialization
    const int numVertices = inMesh.n_vertices();
    const int numFaces = inMesh.n_faces();
    outFeatures = FeatureMat::Zero(numVertices,NUM_FEATURES);
    outFaces = FacesMat::Zero(numFaces,3);

    //# Convert the vertex normals and positions to eigen feature matrices
    convert_mesh_to_matrices(inMesh, outFeatures);

    //# Build the matrix containing the faces
    //## Initialize the face iterator
    TriMesh::FaceIter faceIt(inMesh.faces_begin());
    TriMesh::FaceIter faceEnd(inMesh.faces_end());
    //## Loop over every face
    for (size_t i = 0 ; faceIt != faceEnd ; i++, faceIt++) {

        //## For this face, loop over all its vertices
        //### Initialize the face-vertex iterator
        TriMesh::FaceVertexIter fvIt(inMesh, faceIt);
        //### Loop
        for (size_t j = 0 ; fvIt ; j++, fvIt++) {
            outFaces(i,j) = fvIt->idx();
        }
    }

}//end convert_mesh_to_matrices()


void convert_mesh_to_matrices(const TriMesh &inputMesh,
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
}//end convert_mesh_to_matrices()



void convert_mesh_to_matrices(const TriMesh &inMesh,
                                FeatureMat &outFeatures,
                                FacesMat &outFaces,
                                VecDynFloat &outFlags) {
    /*
    GOAL
    This function converts an openmesh data structure to a feature matrix,
    a faces matrix and a matrix containing the flags (values between 0.0 and
    1.0 that were assigned to each vertex).

    INPUT
    -inMesh:
    this has to be a mesh of openmesh's TriMesh type

    PARAMETERS

    OUTPUT
    -outFeatures:
    a numVertices x 6 Eigen dense matrix where the first three columns are made
    up of the positions of the vertices, and the last three columns the normals
    of those vertices.
    -outFaces:
    a numFaces x 3 Eigen dense matrix where each row contains the corresponding
    indices of the vertices belonging to that face.
    -outFlags
    */
    //# Info & Initialization
    size_t numVertices = inMesh.n_vertices();
    outFlags = VecDynFloat::Zero(numVertices);
    OpenMesh::VPropHandleT<float> flags;
    bool flagsExist = inMesh.get_property_handle(flags, "flags");
    if (!flagsExist)
    {
        std::cerr << "Tried to access the 'flags' property of input mesh - couldn't find handle\n";
        exit(1);
    }

    //# Convert to features and faces
    convert_mesh_to_matrices(inMesh, outFeatures, outFaces);

    //# Build the matrix containing the flags
    //## Initialize the vertex iterator
    TriMesh::VertexIter vertexIt(inMesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(inMesh.vertices_end());
    //## Loop over every face
    for (size_t i = 0 ; vertexIt != vertexEnd ; i++, vertexIt++) {
            outFlags[i] = inMesh.property(flags, vertexIt);
    }
}


void convert_matrices_to_mesh(const Vec3Mat &inPositions,
                                const FacesMat &inFaces,
                                TriMesh &outMesh){
    /*
    GOAL
    This function converts a mesh representation using eigen matrices
    into a mesh representation using OpenMesh's TriMesh.

    INPUT
    -inPositions:
    a numVertices x 3 Eigen dense matrix where the three columns are made
    up of the positions of the vertices.
    -inFaces:
    a numFaces x 3 Eigen dense matrix where each row contains the corresponding
    indices of the vertices belonging to that face.

    PARAMETERS

    OUTPUT
    -outMesh:
    this has to be a mesh of openmesh's TriMesh type. The function expects this
    to be an empty mesh !!
    */

    //# Info & Initialization
    const size_t numElements = inPositions.rows();
    FeatureMat features = FeatureMat::Zero(numElements,registration::NUM_FEATURES);

    //# Convert position matrix to feature matrix.
    features.leftCols(3) = inPositions;

    //# Call conversion function.
    convert_matrices_to_mesh(features, inFaces, outMesh);

}

void convert_matrices_to_mesh(const FeatureMat &inFeatures,
                            const FacesMat &inFaces,
                            TriMesh &outMesh) {
    /*
    GOAL
    This function converts a mesh representation using eigen matrices
    into a mesh representation using OpenMesh's TriMesh.

    INPUT
    -inFeatures:
    a numVertices x 6 Eigen dense matrix where the first three columns are made
    up of the positions of the vertices, and the last three columns the normals
    of those vertices.
    -inFaces:
    a numFaces x 3 Eigen dense matrix where each row contains the corresponding
    indices of the vertices belonging to that face.

    PARAMETERS

    OUTPUT
    -outMesh:
    this has to be a mesh of openmesh's TriMesh type. The function expects this
    to be an empty mesh !!
    */

    //# Info and Initialization
    const size_t numVertices = inFeatures.rows();
    const size_t numFaces = inFaces.rows();
    if (outMesh.n_vertices() > 0) { std::cerr<< "convert_eigen_to_openmesh expects an empty mesh as input!" << std::endl; }

    //# Add each vertex and save the vertex handles for adding the faces later
    TriMesh::VertexHandle vertexHandle;
    std::vector<TriMesh::VertexHandle>  vertexHandles;
    for (size_t i = 0 ; i < numVertices ; i++) {
        vertexHandle = outMesh.add_vertex(TriMesh::Point(inFeatures(i,0),inFeatures(i,1),inFeatures(i,2)));
        vertexHandles.push_back(vertexHandle);
    }

    //# Loop over the faces and add them all to the mesh
    for (size_t i = 0 ; i < numFaces ; i++) {
        //## Make a list of the vertex handles belonging to this face
        std::vector<TriMesh::VertexHandle>  faceVertexHandles;
        //### Vertex 0
        size_t vertexIndex = inFaces(i,0);
        vertexHandle = vertexHandles[vertexIndex];
        faceVertexHandles.push_back(vertexHandle);
        //### Vertex 1
        vertexIndex = inFaces(i,1);
        vertexHandle = vertexHandles[vertexIndex];
        faceVertexHandles.push_back(vertexHandle);
        //### Vertex 2
        vertexIndex = inFaces(i,2);
        vertexHandle = vertexHandles[vertexIndex];
        faceVertexHandles.push_back(vertexHandle);

        //## Add the list of vertex handles as a face to the mesh
        outMesh.add_face(faceVertexHandles);
    }

    //# Update the vertex normals of the mesh (while taking care not to flip them)
    update_normals_safely(inFeatures, outMesh);

}//end convert_matrices_to_mesh()




void update_normals_safely(const FeatureMat &features, TriMesh &mesh){
    /*
    GOAL
    Before updating the normals we're gonna check the current average normal.
    After recomputing the normals, we'll check the average normal again. If the
    dot product between them is negative, we assume the normals are flipped and
    we'll return flip them again.

    INPUT
    -features:
    the features should contain the normals of the mesh BEFORE using the mesh to
    update those normals.
    -mesh:
    The mesh should should already have the positions updated, but not yet the normals.

    PARAMETERS

    OUTPUT
    -mesh:
    the output is the mesh with safely updated normals
    */

    //# Info & Initialization
    const size_t numVertices = features.rows();

    //# Compute normal before updating.
    Eigen::Vector3f avgNormalBefore(0.0f,0.0f,0.0f);
    for(size_t i = 0 ; i < numVertices ; i++){
        avgNormalBefore += features.row(i).tail(3);
    }
    avgNormalBefore /= numVertices;

    //# Make sure the mesh has computed/updated its vertex normals
    mesh.request_face_normals();
    mesh.request_vertex_normals();
    mesh.update_normals();
    mesh.release_face_normals();

    //# Compute average normal after updating
    Eigen::Vector3f avgNormalAfter(0.0f,0.0f,0.0f);
    TriMesh::Normal normal(0.0f,0.0f,0.0f);
    TriMesh::VertexIter vertexIt(mesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(mesh.vertices_end());
    for (size_t i = 0 ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
        //## Get normal
        normal = mesh.normal(vertexIt);
        //## Insert position and normal into 'outputFeatures'
        avgNormalAfter(0) += normal[0];
        avgNormalAfter(1) += normal[1];
        avgNormalAfter(2) += normal[2];
    }
    avgNormalAfter /= numVertices;

    //# See if the normals are flipped
    bool normalsAreFlipped = false;
    //## compute the dot product of the average normals
    float dotProduct = avgNormalBefore.dot(avgNormalAfter);
    if (dotProduct < 0.0f){
        normalsAreFlipped = true;
    }

    //# If the normals are flipped, flip them again
    vertexIt = mesh.vertices_begin();
    if (normalsAreFlipped) {
        for (size_t i = 0 ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
            //## Set the vertex normal to itself but flipped.
            mesh.set_normal(vertexIt, -1.0 * mesh.normal(vertexIt));
        }
    }
}

void convert_matrices_to_mesh(const FeatureMat &inFeatures,
                            const FacesMat &inFaces,
                            const VecDynFloat &inFlags,
                            TriMesh &outMesh) {
    /*
    GOAL
    This function converts a mesh representation using eigen matrices
    into a mesh representation using OpenMesh's TriMesh. Additionally,
    flags (scalar between 0.0 and 1.0) are assigned to each vertex.

    INPUT
    -inFeatures:
    a numVertices x 6 Eigen dense matrix where the first three columns are made
    up of the positions of the vertices, and the last three columns the normals
    of those vertices.
    -inFaces:
    a numFaces x 3 Eigen dense matrix where each row contains the corresponding
    indices of the vertices belonging to that face.
    -inFlags

    PARAMETERS

    OUTPUT
    -outMesh:
    this has to be a mesh of openmesh's TriMesh type. The function expects this
    to be an empty mesh !!
    */

    //# Use previously implemented function to build a mesh using features and faces matrices
    convert_matrices_to_mesh(inFeatures, inFaces, outMesh);

    //# Add the flags to the mesh vertices
    OpenMesh::VPropHandleT<float> flags;
    outMesh.add_property(flags, "flags");
    outMesh.property(flags).set_persistent(true);
    TriMesh::VertexIter vertexIt(outMesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(outMesh.vertices_end());
    for (size_t i = 0 ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
        outMesh.property(flags, vertexIt) = inFlags[i];
    }

}//end convert_matrices_to_mesh()


void convert_eigen_to_openmesh(const FeatureMat &inFeatures,
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

void load_obj_to_eigen(const std::string inObjFilename,
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
    convert_mesh_to_matrices(outMesh, outFeatureMatrix);

}





void write_eigen_to_obj(const FeatureMat &inFeatures,
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
    convert_eigen_to_openmesh(inFeatures, inMesh);

    //# Write resulting mesh to OBJ file
    OpenMesh::IO::_OBJWriter_();
    if (!OpenMesh::IO::write_mesh(inMesh, inObjFilename))
    {
        std::cerr << "write error\n";
        exit(1);
    }
}//end write_eigen_features_to_obj()

bool import_data(const std::string inFloatingMeshPath,
                 const std::string inTargetMeshPath,
                 FeatureMat &outFloatingFeatures,
                 FeatureMat &outTargetFeatures,
                 FacesMat &outFloatingFaces,
                 FacesMat &outTargetFaces){
    std::cout << "Importing Data..." << std::endl;

    //# Safety check (do the input files exist?)
    bool inputfilesExist = true;
    //## Floating Mesh
    std::ifstream infile(inFloatingMeshPath);
    if (infile.good() != true){
        inputfilesExist = false;
        std::cerr << "Floating mesh file does not exist" << std::endl;
    }
    infile.close();
    //## Target Mesh
    infile.open(inTargetMeshPath);
    if (infile.good() != true){
        inputfilesExist = false;
        std::cerr << "Target mesh file does not exist" << std::endl;
    }
    infile.close();
    if (!inputfilesExist){
        std::cout<< "DataImporter can't update - input files not found!" << std::endl;
        return inputfilesExist;
    }

    //# Load OpenMesh meshes
    //## Floating Mesh
    TriMesh floatingMesh;
    OpenMesh::IO::_OBJReader_();
    if (!OpenMesh::IO::read_mesh(floatingMesh,inFloatingMeshPath)){
        std::cerr << "Read error \n";
        exit(1);
    };
    //## Target Mesh
    TriMesh targetMesh;
    if (!OpenMesh::IO::read_mesh(targetMesh,inTargetMeshPath)){
        std::cerr << "Read error \n";
        exit(1);
    };

    //# Update the mesh vertex normals
    floatingMesh.request_face_normals();
    floatingMesh.request_vertex_normals();
    floatingMesh.update_face_normals();
    floatingMesh.update_vertex_normals();
    targetMesh.request_face_normals();
    targetMesh.request_vertex_normals();
    targetMesh.update_face_normals();
    targetMesh.update_vertex_normals();

    //# Convert OpenMesh meshes to Eigen matrices
    convert_mesh_to_matrices(floatingMesh, outFloatingFeatures, outFloatingFaces);
    convert_mesh_to_matrices(targetMesh, outTargetFeatures, outTargetFaces);


    std::cout << "Imported Data" << std::endl;

    return true;
}//end import_data()


bool export_data(FeatureMat &inResultFeatures,
                 FacesMat &inResultFaces,
                 const std::string inResultMeshPath) {
    std::cout << "Exporting Data..." << std::endl;

    //# Convert the matrices to a mesh
    TriMesh resultMesh;
    convert_matrices_to_mesh(inResultFeatures, inResultFaces, resultMesh);
    //# Write the mesh to file
    OpenMesh::IO::_OBJWriter_();
    if (!OpenMesh::IO::write_mesh(resultMesh, inResultMeshPath))
    {
        std::cerr << "write error\n";
        exit(1);
    }
    else {std::cout << "Data Exported." << std::endl;}

    return true;
}

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
    TriMesh::Point updatedPosition(0.0f,0.0f,0.0f);
    unsigned int i = 0;
    TriMesh::VertexIter vertexIt(ioMesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(ioMesh.vertices_end());
    for ( ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
        //## Get positions and normals from the feature matrix
        updatedPosition[0] = ioFeatures(i,0);
        updatedPosition[1] = ioFeatures(i,1);
        updatedPosition[2] = ioFeatures(i,2);
        //## Insert position and normal into mesh
        ioMesh.set_point(vertexIt,updatedPosition);
    }

    //# Given the new positions, let ioMesh update its normals
    update_normals_safely(ioFeatures, ioMesh);

    //# Copy ioMesh's new normals into 'ioFeatures' last three columns.
    TriMesh::Point updatedNormal(0.0f,0.0f,0.0f);
    i = 0;
    vertexIt = ioMesh.vertices_begin();
    for ( ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
        //## Get normal
        updatedNormal = ioMesh.normal(vertexIt);
        //## Insert normal into ioFeatures
        ioFeatures(i,3) = updatedNormal[0];
        ioFeatures(i,4) = updatedNormal[1];
        ioFeatures(i,5) = updatedNormal[2];
    }

}

void update_normals_for_altered_positions(const Vec3Mat &inPositions,
                                        const FacesMat &inFaces,
                                        Vec3Mat &outNormals){
    /*
    GOAL
    This function takes the positions and faces, creates a mesh internally,
    and uses OpenMesh's functionality to recompute the vertex normals
    given those vertex positions.

    INPUT
    -inPositions
    -inFaces

    PARAMETERS

    OUTPUT
    -outNormals
    */

    //# Convert matrices to mesh
    TriMesh mesh;
    convert_matrices_to_mesh(inPositions, inFaces, mesh);

    //# let the mesh update its normals (note: can't use update_normals_safely() here
    //# because that requires having older values for the normals, whereas
    //# this function is called when there are no normals at all and they need to be
    //# determined.)
    mesh.request_face_normals();
    mesh.request_vertex_normals();
    mesh.update_normals();
    mesh.release_face_normals();

    //# Copy mesh's new normals into 'outNormals'.
    TriMesh::Point normal(0.0f,0.0f,0.0f);
    TriMesh::VertexIter vertexIt(mesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(mesh.vertices_end());
    for (unsigned int i = 0 ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
        //## Get normal
        normal = mesh.normal(vertexIt);
        //## Insert normal into ioFeatures
        outNormals(i,0) = normal[0];
        outNormals(i,1) = normal[1];
        outNormals(i,2) = normal[2];
    }

}


}//namespace registration
