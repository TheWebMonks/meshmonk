#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <nanoflann.hpp>

using namespace std;

const int NUM_FEATURES = 6;

typedef OpenMesh::TriMesh_ArrayKernelT<>  TriMesh;
typedef Eigen::Matrix< unsigned int, Eigen::Dynamic, Eigen::Dynamic> MatDynInt; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat; //matrix MxN of type float
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Vector3f Vec3Float;
typedef Eigen::Vector4f Vec4Float;
typedef Eigen::Matrix3f Mat3Float;
typedef Eigen::Matrix4f Mat4Float;
typedef Eigen::Matrix< float, Eigen::Dynamic, 3> Vec3Mat; //matrix Mx3 of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< float, 1, NUM_FEATURES> FeatureVec; //matrix Mx6 of type float
typedef Eigen::SparseMatrix<float, 0, int> SparseMat;
typedef Eigen::Triplet<float> Triplet;
typedef Eigen::SelfAdjointEigenSolver<Mat4Float> EigenVectorDecomposer;
typedef nanoflann::KDTreeEigenMatrixAdaptor<FeatureMat, NUM_FEATURES, nanoflann::metric_L2>  NanoKDTree6D;



void k_nearest_neighbours(const FeatureMat &inQueriedPoints,
                        const FeatureMat &inSourcePoints,
                        MatDynInt &outNeighbourIndices,
                        MatDynFloat &outNeighbourSquaredDistances,
                        const size_t paramK = 3,
                        const size_t paramLeafsize = 15){
    /*
    GOAL
    This function searches for the k nearest neighbours in 'inSourcePoints' for
    each element in the 'inQueriedPoints' set. It outputs the indices of each
    neighbour and the squared (!) distances between each element of
    'inQueriedPoints' and its neighbours.

    INPUT
    -inQueriedPoints:
    -inSourcePoints:

    PARAMETERS
    -paramK(= 3): number of nearest neighbours
    -paramLeafsize(= 15): should be between 5 and 50 or so

    OUTPUT
    -outNeighbourIndices
    -outNeighbourSquaredDistances.
    */

    //# Info and Initialization
    const size_t dimension = inSourcePoints.cols();
    const size_t numSourceElements = inSourcePoints.rows();
    const size_t numQueriedElements = inQueriedPoints.rows();
    outNeighbourIndices = MatDynInt::Zero(numQueriedElements,paramK);
    outNeighbourSquaredDistances = MatDynFloat::Zero(numQueriedElements,paramK);

    //# Construct kd-tree
    NanoKDTree6D kdTree(dimension, inSourcePoints, paramLeafsize);
    kdTree.index->buildIndex();

    //# Query the kd-tree
    //## Loop over the queried features
    //### Initialize variables we'll need during the loop
    unsigned int i = 0;
    unsigned int j = 0;
    std::vector<float> queriedFeature(dimension);
    std::vector<size_t> neighbourIndices(paramK);
    std::vector<float> neighbourSquaredDistances(paramK);
    nanoflann::KNNResultSet<float> knnResultSet(paramK);

    //### Execute loop
    for ( ; i < numQueriedElements ; ++i ) {
        //### Initiliaze the knnResultSet
        knnResultSet.init(&neighbourIndices[0], &neighbourSquaredDistances[0]);

        //### convert input features to 'queriedFeature' std::vector structure
        //### (required by nanoflann's kd-tree).
        for (j = 0 ; j < dimension ; ++j) {
            queriedFeature[j] = inQueriedPoints(i,j);
        }

        //### Query the kd-tree
//        size_t numNeighboursFound = kdTree.knnSearch(&queriedFeature[0], paramK, &neighbourIndices[0], &neighbourSquaredDistances[0]);
        kdTree.index->findNeighbors(knnResultSet, &queriedFeature[0],
                                    nanoflann::SearchParams(32, 0.0001 /*eps*/, true));

        //### Copy the result into the outputs by looping over the k nearest
        //### neighbours
        for (j = 0 ; j < paramK ; ++j) {
            outNeighbourIndices(i,j) = neighbourIndices[j];
            outNeighbourSquaredDistances(i,j) = neighbourSquaredDistances[j];
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

void load_obj_to_eigen_features(const string inObjFilename,
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
        cerr << "Read error \n";
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
    cerr << "Number of rows does not correspond with number of vertices when"
    << " calling eigen_features_to_openmesh()" << endl;
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
                                const string inObjFilename){
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
    cerr << "Number of rows does not correspond with number of vertices when"
    << " calling write_eigen_features_to_obj()" << endl;
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
    cerr << "Number of rows does not correspond with number of vertices when"
    << " calling eigen_features_to_openmesh()" << endl;
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

void wknn_affinity(const FeatureMat &inFeatures1,
                    const FeatureMat &inFeatures2,
                    SparseMat &outAffinity,
                    const size_t paramK = 3) {
    /*
    # GOAL
    For each element in inFeatures1, you're going to determine affinity weights
    which link it to elements in inFeatures2. This is based on the Euclidean
    distance of k nearest neighbours found for each element of inFeatures1.

    # INPUTS
    -inFeatures1
    -inFeatures2

    # PARAMETERS
    -paramK(=3):
    number of nearest neighbours

    # OUTPUT
    -outAffinity
    */

    //# Info and Initialization
    //## Number of elements in each feature matrix
    const size_t numVertices1 = inFeatures1.rows();
    const size_t numVertices2 = inFeatures2.rows();
    //## Initialize matrices to save k-nn results
    MatDynInt neighbourIndices = MatDynInt::Zero(numVertices1, paramK);
    MatDynFloat neighbourSquaredDistances = MatDynFloat::Zero(numVertices1, paramK);
    //## Initialize affinity matrix
    outAffinity = SparseMat(numVertices1, numVertices2);
    outAffinity.setZero();
    //## Initialize a vector of triplets which is used to insert elements into
    //## the affinity matrix later.
    const size_t numAffinityElements = numVertices1 * paramK;
    std::vector<Triplet> affinityElements(numAffinityElements, Triplet(0,0,0.0f));
    //## Parameters for the k-nn search
    const size_t maxLeafsize = 15;


    //# Determine the nearest neighbours
    k_nearest_neighbours(inFeatures1, inFeatures2, neighbourIndices,
                        neighbourSquaredDistances, paramK, maxLeafsize);

    //# Compute the affinity matrix
    //## Loop over the first feature set to determine their affinity with the
    //## second set.
    size_t i = 0;
    size_t j = 0;
    unsigned int counter = 0;
    for ( ; i < numVertices1 ; i++) {
        //### Loop over each found neighbour
        for ( j = 0 ; j < paramK ; j++) {
            //### Get index of neighbour and squared distance to it
            const int neighbourIndex = neighbourIndices(i,j);
            float distanceSquared = neighbourSquaredDistances(i,j);

            //### For numerical stability, check if the distance is very small
            const float eps1 = 0.000001f;
            if (distanceSquared < eps1) {distanceSquared = eps1;}

            //### Compute the affinity element as 1/distance*distance
            float affinityElement = 1.0f / distanceSquared;
            const float eps2 = 0.0001f;
            //### Check for numerical stability (the affinity elements will be
            //### normalized later, so dividing by a sum of tiny elements might
            //### go wrong.
            if (affinityElement < eps2) {affinityElement = eps2;}

            //### Write result to the affinity triplet list
            affinityElements[counter] = Triplet(i,neighbourIndex,affinityElement);
            counter++;
        }
    }

    //## Construct the sparse matrix with the computed element list
    outAffinity.setFromTriplets(affinityElements.begin(),
                                affinityElements.end());

    //# Normalize the rows of the affinity matrix
    normalize_sparse_matrix(outAffinity);
}//end wknn_affinity()

void fuse_affinities(SparseMat &ioAffinity1, const SparseMat &iAffinity2){
    /*
    # GOAL
    Fuse the two affinity matrices together. The result is inserted into
    ioAffinity1.

    # INPUTS
    -ioAffinity1
    -iAffinity2: dimensions should be transposed of ioAffinity1

    # PARAMETERS

    # OUTPUT
    -ioAffinity1

    # RETURNS
    """
    */
    //# Info and Initialization
    const size_t numRows1 = ioAffinity1.rows();
    const size_t numRows2 = iAffinity2.rows();
    const size_t numCols1 = ioAffinity1.rows();
    const size_t numCols2 = iAffinity2.rows();
    //## Safety check for input sizes
    if((numRows1 != numCols1) || (numCols1 != numRows2)) {
        cerr << "The sizes of the inputted matrices in fuse_affinities are wrong."
        << "Their sizes should be the transpose of each other!" << endl;
    }

    //# Fusing is done by simple averaging
    //## Sum matrices. Note: because Eigen's Sparse matrices require storage
    //## orders to match (row- or column-major) we need to construct a new
    //## temporary matrix of the tranpose of iAffinity2.
    ioAffinity1 += SparseMat(iAffinity2.transpose());
    //## Normalize result
    normalize_sparse_matrix(ioAffinity1);

}

void wkkn_correspondences(const FeatureMat &inFloatingFeatures,
                            const FeatureMat &inTargetFeatures,
                            FeatureMat &outCorrespondingFeatures,
                            const size_t paramK = 3,
                            const bool paramSymmetric = true) {
    /*
    # GOAL
    For each element in inFloatingFeatures, we're going to find corresponding
    features in the elements of inTargetFeatures. This is done in a symmetrical
    weighted nearest neighbour approach approach.
    The used weight is one over the distance squared.

    # INPUTS
    -inFloatingFeatures
    -inTargetFeatures

    # PARAMETERS
    -paramK(=3):
    number of nearest neighbours
    -paramSymmetric:
    If true, a push-pull approach is used to find correspondences. Not only
    will we look for correspondences for inFloatingFeatures in the
    inTargetFeatures set, we'll also do the opposite: find correspondences for
    the inTargetFeatures set in the inFloatingFeatures set. These findings are
    combined which creates an effect where the inTargetFeatures sort of attract
    the inFloatingFeatures towards itself.

    # OUTPUT
    -outCorrespondingFeatures
    */

    //# Info & Initialization
    const size_t numFloatingVertices = inFloatingFeatures.rows();
    const size_t numTargetVertices = inTargetFeatures.rows();

    SparseMat affinity = SparseMat(numFloatingVertices, numTargetVertices);
    SparseMat affinityPull;

    //# Compute the affinity for inFloatingFeatures towards inTargetFeatures
    wknn_affinity(inFloatingFeatures, inTargetFeatures, affinity, paramK);

    //# Ccompute symmetric correspondences (if required)
    if (paramSymmetric == true) {
        affinityPull = SparseMat(numTargetVertices, numFloatingVertices);
        //## Compute the affinity for inTargetFeatures towards inFloatingFeatures
        wknn_affinity(inTargetFeatures, inFloatingFeatures, affinityPull, paramK);

        //## Fuse the two affinities together.
        fuse_affinities(affinity, affinityPull);
    }

    //## Compute corresponding features as affinity matrix multiplied with the
    //## target features.
    outCorrespondingFeatures = affinity * inTargetFeatures;

}//end wkkn_correspondences()







void rigid_transformation(FeatureMat &ioFeatures,
                            const FeatureMat &inCorrespondingFeatures,
                            const VecDynFloat &inWeights,
                            const bool paramScaling = false) {
    /*
    # GOAL
    This function computes the rigid transformation between a set a features and
    a set of corresponding features. Each correspondence can be weighed between
    0.0 and 1.0.
    The features are automatically transformed, and the function returns the
    transformation matrix that was used.

    # INPUTS
    -ioFeatures
    -inCorrespondingFeatures
    -inWeights

    # PARAMETERS
    -paramScaling:
    Whether or not to allow scaling.

    # OUTPUTS
    -ioFeatures
    */

    //# Info & Initialization
    const size_t numVertices = ioFeatures.rows();
    const size_t numFeatures = ioFeatures.cols();
    MatDynFloat floatingPositions = MatDynFloat::Zero(3, numVertices);
    MatDynFloat correspondingPositions = MatDynFloat::Zero(3, numVertices);

    //## Tranpose the data if necessary
    if ((numVertices > numFeatures) && (numFeatures == NUM_FEATURES)) { //this should normally be the case
        floatingPositions = ioFeatures.leftCols(3).transpose();
        correspondingPositions = inCorrespondingFeatures.leftCols(3).transpose();
    }
    else {
        cerr<< "Warning: input of rigid transformation expects rows to correspond with elements, not features, and to have more elements than features per element." << endl;
    }

    //# Compute the tranformation in 10 steps.
    //## 1. Get the centroids of each set
    Vec3Float floatingCentroid = Vec3Float::Zero();
    Vec3Float correspondingCentroid = Vec3Float::Zero();
    float sumWeights = 0.0;
    //### Weigh and sum all features
    for (size_t i = 0 ; i < numVertices ; i++) {
        floatingCentroid += inWeights[i] * floatingPositions.col(i).segment(0,3);
        correspondingCentroid += inWeights[i] * correspondingPositions.col(i).segment(0,3);
        sumWeights += inWeights[i];
    }
    //### Divide by total weight
    floatingCentroid /= sumWeights;
    correspondingCentroid /= sumWeights;

    //## 2. Compute the Cross Variance matrix
    Mat3Float crossVarianceMatrix = Mat3Float::Zero();
    for(size_t i = 0 ; i < numVertices ; i++) {
        crossVarianceMatrix += inWeights[i] * floatingPositions.col(i) * correspondingPositions.col(i).transpose();
    }
    crossVarianceMatrix = crossVarianceMatrix / sumWeights - floatingCentroid*correspondingCentroid.transpose();

    //## 3. Compute the Anti-Symmetric matrix
    Mat3Float antiSymmetricMatrix = crossVarianceMatrix - crossVarianceMatrix.transpose();

    //## 4. Use the cyclic elements of the Anti-Symmetric matrix to construct delta
    Vec3Float delta = Vec3Float::Zero();
    delta[0] = antiSymmetricMatrix(1,2);
    delta[1] = antiSymmetricMatrix(2,0);
    delta[2] = antiSymmetricMatrix(0,1);

    //## 5. Compute Q
    Mat4Float Q = Mat4Float::Zero();
    Q(0,0) = crossVarianceMatrix.trace();
    Q.block<3,1>(1,0) = delta;
    Q.block<1,3>(0,1) = delta.transpose();
    Q.block<3,3>(1,1) = crossVarianceMatrix + crossVarianceMatrix.transpose()
                        - crossVarianceMatrix.trace() * Mat3Float::Identity();

    //## 6. Now compute the rotation quaternion by finding the eigenvector of Q
    //## of its largest eigenvalue.
    Vec4Float rotQuat = Vec4Float::Zero();
    EigenVectorDecomposer decomposer(Q);
    if (decomposer.info() != Eigen::Success) {
        cerr << "eigenvector decomposer on Q failed!" << std::endl;
        cerr << "Q : " << Q << std::endl;
    }
    size_t indexMaxVal = 0;
    float maxEigenValue = 0.0;
    for ( size_t i = 0 ; i < 4 ; i++ ){
        if ( decomposer.eigenvalues()[i] > maxEigenValue ){
            maxEigenValue = decomposer.eigenvalues()[i];
            indexMaxVal = i;
        }
    }
    rotQuat = decomposer.eigenvectors().col( indexMaxVal );

    //## 7. Construct the rotation matrix
    Mat3Float rotMatTemp = Mat3Float::Zero();
    //### diagonal elements
    rotMatTemp(0,0) = rotQuat[0] * rotQuat[0] + rotQuat[1] * rotQuat[1] - rotQuat[2] * rotQuat[2] - rotQuat[3] * rotQuat[3];
    rotMatTemp(1,1) = rotQuat[0] * rotQuat[0] + rotQuat[2] * rotQuat[2] - rotQuat[1] * rotQuat[1] - rotQuat[3] * rotQuat[3];
    rotMatTemp(2,2) = rotQuat[0] * rotQuat[0] + rotQuat[3] * rotQuat[3] - rotQuat[1] * rotQuat[1] - rotQuat[2] * rotQuat[2];
    //### remaining elements
    rotMatTemp(0,1) = 2.0 * (rotQuat[1] * rotQuat[2] - rotQuat[0] * rotQuat[3]);
    rotMatTemp(1,0) = 2.0 * (rotQuat[1] * rotQuat[2] + rotQuat[0] * rotQuat[3]);
    rotMatTemp(0,2) = 2.0 * (rotQuat[1] * rotQuat[3] + rotQuat[0] * rotQuat[2]);
    rotMatTemp(2,0) = 2.0 * (rotQuat[1] * rotQuat[3] - rotQuat[0] * rotQuat[2]);
    rotMatTemp(1,2) = 2.0 * (rotQuat[2] * rotQuat[3] - rotQuat[0] * rotQuat[1]);
    rotMatTemp(2,1) = 2.0 * (rotQuat[2] * rotQuat[3] + rotQuat[0] * rotQuat[1]);

    //## 8. Estimate scale (if required)
    float scaleFactor = 1.0; //>1 to grow ; <1 to shrink
    if (paramScaling == true){
        float numerator = 0.0;
        float denominator = 0.0;
        for (size_t i = 0 ; i < numVertices ; i++){
            //### Center and rotate the floating position
            Vec3Float newFloatingPos = rotMatTemp * (floatingPositions.block<3,1>(0,i) - floatingCentroid.segment(0, 3));
            //### Center the corresponding position
            Vec3Float newCorrespondingPos = correspondingPositions.block<3,1>(0,i) - correspondingCentroid.segment(0, 3);

            //### Increment numerator and denominator
            numerator += inWeights[i] * newCorrespondingPos.transpose() * newFloatingPos;
            denominator += inWeights[i] * newFloatingPos.transpose() * newFloatingPos;
        }
        scaleFactor = numerator / denominator;
    }


    //## 9. Compute the remaining translation necessary between the centroids
    Vec3Float translation = correspondingCentroid - scaleFactor * rotMatTemp * floatingCentroid;

    //## 10. Compute the entire transformation matrix.
    //### Initialize Matrices
    Mat4Float translationMatrix = Mat4Float::Identity();
    Mat4Float rotationMatrix = Mat4Float::Identity();
    Mat4Float transformationMatrix = Mat4Float::Identity();
    //### Convert to homogeneous transformation matrices
    translationMatrix.block<3, 1>(0,3) = translation;
    rotationMatrix.block<3, 3>(0, 0) = scaleFactor * rotMatTemp;
    //### Matrix transformations on data is performed from right to left.
    //### Translation should be performed before rotating, so translationMatrix
    //### stands right in the multiplication with rotationMatrix.
    transformationMatrix = rotationMatrix * translationMatrix;

    //# Apply the transformation
    //## initialize the position in a [x y z 1] representation
    Vec4Float position4d = Vec4Float::Ones();
    for (size_t i = 0 ; i < numVertices ; i++) {
        //## Extraxt position from feature matrix
        position4d.segment(0, 3) = floatingPositions.block<3,1>(0,i);
        //## Apply transformation and assign to 'position4d' variable again
        position4d = transformationMatrix * position4d;
        ioFeatures.block<1,3>(i,0) = position4d.segment(0, 3);
    }
}

int main()
{

    /*
    ############################################################################
    ##############################  INPUT  #####################################
    ############################################################################
    */
    //# IO variables
//    const string fuckedUpBunnyDir = "/home/jonatan/kuleuven-algorithms/examples/data/bunny_slightly_rotated.obj";
    const string fuckedUpBunnyDir = "/home/jonatan/kuleuven-algorithms/examples/data/fucked_up_bunny.obj";
    const string bunnyDir = "/home/jonatan/kuleuven-algorithms/examples/data/bunny90.obj";
    const string fuckedUpBunnyResultDir = "/home/jonatan/kuleuven-algorithms/examples/data/fucked_up_bunny_result.obj";

    //# Load meshes and convert to feature matrices
    TriMesh fuckedUpBunny;
    TriMesh bunny;
    FeatureMat floatingFeatures;
    FeatureMat targetFeatures;
    load_obj_to_eigen_features(fuckedUpBunnyDir, fuckedUpBunny, floatingFeatures);
    load_obj_to_eigen_features(bunnyDir, bunny, targetFeatures);
    //floatingFeatures = FeatureMat::Zero(3, NUM_FEATURES);
    //targetFeatures = FeatureMat::Zero(3, NUM_FEATURES);
    //floatingFeatures << 0.1f, 0.1f, 0.1f, 0.0f, 0.0f, 1.0f,
                        1.1f, 0.1f, 0.1f, 0.0f, 0.0f, 1.0f,
                        0.1f, 1.1f, 0.1f, 0.0f, 0.0f, 1.0f;
    //targetFeatures << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                      1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                      0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f;


    /*
    ############################################################################
    ##############################  RIGID ICP  #################################
    ############################################################################
    */
    //# Determine Correspondences
    //## Info & Initialization
    size_t numFloatingVertices = floatingFeatures.rows();
    size_t numTargetVertices = targetFeatures.rows();
    const size_t numNearestNeighbours = 3;
    FeatureMat correspondingFeatures = FeatureMat::Zero(numFloatingVertices, NUM_FEATURES);

    //## Compute symmetric wknn correspondences
    wkkn_correspondences(floatingFeatures, targetFeatures, correspondingFeatures, numNearestNeighbours, true);

    //# Inlier Detection
    VecDynFloat inlierWeights = VecDynFloat::Ones(numFloatingVertices);


    //# Compute the transformation
    rigid_transformation(floatingFeatures, correspondingFeatures, inlierWeights, false);

    /*
    ############################################################################
    ##########################  NON-RIGID ICP  #################################
    ############################################################################
    */
    //# Do non-rigid ICP


    /*
    ############################################################################
    ##############################  OUTPUT #####################################
    ############################################################################
    */
    //# Write result to file
    write_eigen_features_to_obj(floatingFeatures, fuckedUpBunny, fuckedUpBunnyResultDir);


    return 0;
}
