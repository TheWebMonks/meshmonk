#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <nanoflann.hpp>
#include <stdio.h>
#include "InlierDetector.hpp"
#include "global.hpp"

using namespace std;



typedef OpenMesh::TriMesh_ArrayKernelT<>  TriMesh;
typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> MatDynInt; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat; //matrix MxN of type float
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Vector3f Vec3Float;
typedef Eigen::Vector4f Vec4Float;
typedef Eigen::Matrix3f Mat3Float;
typedef Eigen::Matrix4f Mat4Float;
typedef Eigen::Matrix< float, Eigen::Dynamic, 3> Vec3Mat; //matrix Mx3 of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< float, 1, registration::NUM_FEATURES> FeatureVec; //matrix Mx6 of type float
typedef Eigen::SparseMatrix<float, 0, int> SparseMat;
typedef Eigen::Triplet<float> Triplet;
typedef Eigen::SelfAdjointEigenSolver<Mat4Float> EigenVectorDecomposer;


namespace registration {

template <typename VecMatType>
void k_nearest_neighbours(const VecMatType &inQueriedPoints,
                        const VecMatType &inSourcePoints,
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
    nanoflann::KDTreeEigenMatrixAdaptor<VecMatType> kdTree(dimension, inSourcePoints, paramLeafsize);
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

template <typename VecMatType>
void radius_nearest_neighbours(const VecMatType &inQueriedPoints,
                                const VecMatType &inSourcePoints,
                                MatDynInt &outNeighbourIndices,
                                MatDynFloat &outNeighbourSquaredDistances,
                                const float paramRadius = 3.0,
                                const size_t paramLeafsize = 15){
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
        cout << "Number of matches: " << nMatches << endl;

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


void fuse_affinities(SparseMat &ioAffinity1,
                    const SparseMat &inAffinity2){
    /*
    # GOAL
    Fuse the two affinity matrices together. The result is inserted into
    ioAffinity1.

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
    const size_t numRows2 = inAffinity2.rows();
    const size_t numCols1 = ioAffinity1.rows();
    const size_t numCols2 = inAffinity2.rows();
    //## Safety check for input sizes
    if((numRows1 != numCols1) || (numCols1 != numRows2)) {
        cerr << "The sizes of the inputted matrices in fuse_affinities are wrong."
        << "Their sizes should be the transpose of each other!" << endl;
    }

    //# Fusing is done by simple averaging
    //## Sum matrices. Note: because Eigen's Sparse matrices require storage
    //## orders to match (row- or column-major) we need to construct a new
    //## temporary matrix of the tranpose of inAffinity2.
    ioAffinity1 += SparseMat(inAffinity2.transpose());
    //## Normalize result
    normalize_sparse_matrix(ioAffinity1);

}


void affinity_to_correspondences(const FeatureMat &inTargetFeatures,
                                    const VecDynFloat &inTargetFlags,
                                    const MatDynFloat &inAffinity,
                                    FeatureMat &outCorrespondingFeatures,
                                    VecDynFloat &outCorrespondingFlags,
                                    const float paramFlagRoundingLimit = 0.9){
    /*
    # GOAL
    This function computes all the corresponding features and flags,
    when the affinity for a set of features and flags is given.

    # INPUTS
    -inTargetFeatures
    -inTargetFlags
    -inAffinity

    # PARAMETERS
    -paramFlagRoundingLimit:
    Flags are binary. Anything over this double is rounded up, whereas anything
    under it is rounded down. A suggested cut-off is 0.9, which means that if
    an element flagged zero contributes 10 percent or to the affinity, the
    corresponding element should be flagged a zero as well.

    # OUTPUT
    -outCorrespondingFeatures
    -outCorrespondingFlags
    */

    //# Info & Initialization
    const size_t numElements = inTargetFeatures.rows();

    //# Simple computation of corresponding features and flags
    outCorrespondingFeatures = inAffinity * inTargetFeatures;
    outCorrespondingFlags = inAffinity * inTargetFlags;

    //# Flag correction.
    //## Flags are binary. We will round them down if lower than the flag
    //## rounding limit (see explanation in parameter description).
    for (size_t i = 0 ; i < numElements ; i++) {
        if (outCorrespondingFlags[i] > paramFlagRoundingLimit){
            outCorrespondingFlags[i] = 1.0;
        }
        else {
            outCorrespondingFlags[i] = 0.0;
        }
    }

}


void wkkn_correspondences(const FeatureMat &inFloatingFeatures,
                            const FeatureMat &inTargetFeatures,
                            const VecDynFloat &inTargetFlags,
                            FeatureMat &outCorrespondingFeatures,
                            VecDynFloat &outCorrespondingFlags,
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
    -inTargetFlags

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
    -outCorrespondingFlags
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
    affinity_to_correspondences(inTargetFeatures, inTargetFlags, affinity,
                                outCorrespondingFeatures, outCorrespondingFlags,
                                0.9);
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





void vector_block_average(const Vec3Mat &inVectors,
                            const VecDynFloat &inWeights,
                            Vec3Float &outVector){
    /*
    GOAL
    This function computes the block average of input vectors (in other words,
    using uniform weights for each contribution).

    INPUT
    -inVectors:
    Values of the vectors of the vector field. The dimensions of this matrix
    are numVectors x 3
    -inWeights:
    The contribution of each vector field vector can be weighed.

    PARAMETERS

    OUTPUT
    -outVector:
    Computed through uniform averaging, taking additional weights (see
    fieldWeights) into account.
    */

    //# Info & Initialization
    size_t numVectors = inVectors.rows();
    size_t dimension = inVectors.cols(); //Assumed to be 3
    if (dimension != 3) {
        cerr << "Dimension is not equal to 3 in vector_block_average() call!" << endl;
    }
    outVector = Vec3Float::Zero();
    float sumWeights = 0.0;

    //# Sum all values
    for (size_t i = 0 ; i < numVectors ; i++) {
        outVector += inVectors.row(i);
        sumWeights += inWeights[i];
    }

    //# Normalize the vector
    outVector /= sumWeights;
}

template <typename VecType, typename VecMatType>
void gaussian_interpolate_scalar_field(const VecType &inQueriedPosition,
                            const VecDynFloat &inScalars,
                            const VecMatType &inVectorPositions,
                            const VecDynFloat &inVectorWeights,
                            float &outQueriedScalar,
                            const float paramSigma){
    /*
    GOAL
    This function computes the gaussian average of a scalar field for a queried
    position. (scalar field = collection of scalars at different vector positions).

    INPUT
    -inQueriedPosition:
    location you want to compute the gaussian average for
    -inScalars:
    Values of the scalars of the vector field.
    -inVectorPositions:
    Locations of the vectors of the vector field.
    -inVectorWeights:
    The contribution of each vector field vector can be weighed.

    PARAMETERS
    -paramSigma:
    The standard deviation of the Gaussian used.

    OUTPUT
    -outQueriedVector:
    Computed through gaussian averaging, taking additional weights (see
    inVectorWeights) into account.
    */

    //# Info & Initialization
    size_t numVectors = inVectorPositions.rows();
    size_t dimension = inVectorPositions.cols(); //Assumed to be 3
    if (dimension != 3) {
        cerr << "Dimension is not equal to 3 in scalar_gaussian_average() call!" << endl;
    }
    outQueriedScalar = 0.0;

    //# Compute Gaussian average;
    float sumWeights = 0.0;
    for (size_t i = 0 ; i < numVectors ; i++) {
        //## Compute distance
        VecType sourceVectorPosition = inVectorPositions.row(i);
        VecType distanceVector = sourceVectorPosition - inQueriedPosition;
        const float distanceSquared = distanceVector.squaredNorm();

        //## From the distance, we can compute the gaussian weight
        const float gaussianWeight = std::exp(-0.5 * distanceSquared / std::pow(paramSigma, 2.0));

        //## Let's take the user-defined weight into account
        const float combinedWeight = gaussianWeight * inVectorWeights[i];

        //## Weigh the current vector field vector and sum it
        outQueriedScalar += combinedWeight * inScalars[i];

        //## Add the weight to the total sum of weights for normalization after.
        sumWeights += combinedWeight;
    }

    //# Normalize the vector
    outQueriedScalar /= sumWeights;

}//end vector_gaussian_average()



template <typename VecType, typename VecMatType>
void gaussian_interpolate_vector_field(const VecType &inQueriedPosition,
                            const VecMatType &inVectors,
                            const VecMatType &inVectorPositions,
                            const VecDynFloat &inVectorWeights,
                            VecType &outQueriedVector,
                            const float paramSigma){
    /*
    GOAL
    This function computes the gaussian average of a vector field for a queried
    position. (vector field = collection of vectors at different vector positions).

    INPUT
    -inQueriedPosition:
    location you want to compute the gaussian average for
    -inVectors:
    Values of the vectors of the vector field. The dimensions of this matrix
    are numVectors x 3
    -inVectorPositions:
    Locations of the vectors of the vector field.
    -inVectorWeights:
    The contribution of each vector field vector can be weighed.

    PARAMETERS
    -paramSigma:
    The standard deviation of the Gaussian used.

    OUTPUT
    -outQueriedVector:
    Computed through gaussian averaging, taking additional weights (see
    inVectorWeights) into account.
    */

    //# Info & Initialization
    size_t numVectors = inVectors.rows();
    size_t dimension = inVectors.cols(); //Assumed to be 3
    if (dimension != 3) {
        cerr << "Dimension is not equal to 3 in vector_gaussian_average() call!" << endl;
    }
    outQueriedVector = VecType::Zero();

    //# Compute Gaussian average;
    float sumWeights = 0.0;
    for (size_t i = 0 ; i < numVectors ; i++) {
        //## Compute distance
        VecType sourceVectorPosition = inVectorPositions.row(i);
        VecType distanceVector = sourceVectorPosition - inQueriedPosition;
        const float distanceSquared = distanceVector.squaredNorm();

        //## From the distance, we can compute the gaussian weight
        const float gaussianWeight = std::exp(-0.5 * distanceSquared / std::pow(paramSigma, 2.0));

        //## Let's take the user-defined weight into account
        const float combinedWeight = gaussianWeight * inVectorWeights[i];

        //## Weigh the current vector field vector and sum it
        outQueriedVector += combinedWeight * inVectors.row(i);

        //## Add the weight to the total sum of weights for normalization after.
        sumWeights += combinedWeight;
    }

    //# Normalize the vector
    outQueriedVector /= sumWeights;

}//end vector_gaussian_average()


template <typename VecType, typename VecMatType>
void gaussian_smoothing_vector_field(const VecMatType &inVectors,
                                    const VecMatType &inVectorPositions,
                                    const VecDynFloat &inVectorWeights,
                                    VecMatType &outSmoothedVectors,
                                    const size_t paramNumNeighbours,
                                    const float paramSigma) {
    /*
    GOAL
    This function performs gaussian smoothing on an entire vector field.

    INPUT
    -inVectorPositions
    -inVectors:
    The vector field that should be smoothed.
    -inVectorWeights

    OUTPUT
    -regulatedFieldVectors:
    The resulting regulated displacement field.


    PARAMETERS
    -paramNumNeighbours:
    For the smoothing, the nearest neighbours for each floating positions have
    to be found. The number should be high enough so that every significant
    contribution (up to a distance of e.g. 3*gaussianSigma) is included. But low
    enough to keep computational speed high.
    -paramSigma:
    The value for sigma of the gaussian used for the smoothing.

    RETURNS
    */

    //# Info & Initialization
    const size_t numVectors = inVectorPositions.rows();
    const size_t numDimensions = inVectorPositions.cols();
    MatDynInt neighbourIndices = MatDynInt::Zero(numVectors, paramNumNeighbours);
    MatDynFloat neighbourSquaredDistances = MatDynFloat::Zero(numVectors, paramNumNeighbours);

    //# Determine for each field node the (closely) neighbouring nodes
    k_nearest_neighbours(inVectorPositions, inVectorPositions, neighbourIndices,
                        neighbourSquaredDistances, paramNumNeighbours, 15);

    //# Use the neighbouring field vectors to smooth each individual field vector
    VecType position = VecType::Zero(numDimensions);
    for (size_t i = 0 ; i < numVectors ; i++) {
        VecType position = inVectorPositions.row(i);

         //## For the current displacement, get all the neighbouring positions,
        //## vectors, and weights (needed for Gaussian smoothing!).
        VecMatType neighbourPositions = VecMatType::Zero(paramNumNeighbours, numDimensions);
        VecMatType neighbourVectors = VecMatType::Zero(paramNumNeighbours, numDimensions);
        VecDynFloat neighbourWeights = VecDynFloat::Zero(paramNumNeighbours);
        for (size_t j = 0 ; j < paramNumNeighbours ; j++) {
            const size_t neighbourIndex = neighbourIndices(i,j);
            neighbourPositions.row(j) = inVectorPositions.row(neighbourIndex);
            neighbourVectors.row(j) = inVectors.row(neighbourIndex);
            neighbourWeights[j] = inVectorWeights[neighbourIndex];
        }

        //## Gaussian averaging of neighbouring displacements
        VecType smoothedVector = VecType::Zero();
        gaussian_interpolate_vector_field(position, neighbourVectors, neighbourPositions,
                                            neighbourWeights, smoothedVector, paramSigma);

        outSmoothedVectors.row(i) = smoothedVector;
    }

}//end gaussian_smoothing_vector_field()

template <typename VecMatType>
void viscoelastic_transformation(VecMatType &ioFloatingPositions,
                                const VecMatType &inCorrespondingPositions,
                                const VecDynFloat &inFloatingWeights,
                                Vec3Mat &ioDisplacementField,
                                const size_t paramNumNeighbourDisplacements,
                                const float paramSigmaSmoothing = 1.0,
                                const size_t paramNumViscousSmoothingIterations = 1,
                                const size_t paramNumElasticSmoothingIterations = 1) {
    /*
    # GOAL
    This function computes the viscoelastic transformation between a set a
    features and a set of corresponding features. Each correspondence can be
    weighed between 0.0 and 1.0.
    The features are automatically transformed, and the function returns the
    displacement field that was used.

    # INPUTS
    -ioFloatingPositions
    -correspondingPositions
    -inFloatingWeights

    # OUTPUTS
    -toBeUpdatedDisplacementField:
    The current displacement field that should be updated.

    # PARAMETERS
    -paramNumNeighbourDisplacements:
    For the regularization, the nearest neighbours for each floating positions
    have to be found. The number should be high enough so that every significant
    contribution (up to a distance of e.g. 3*paramSigmaSmoothing) is included. But low
    enough to keep computational speed high.
    -paramSigmaSmoothing:
    The value for sigma of the gaussian used for the regularization.
    -paramNumViscousSmoothingIterations:
    Number of times the viscous deformation is smoothed.
    -paramNumElasticSmoothingIterations:
    Number of times the elastic deformation is smoothed

    # RETURNS
    */

    //# Info and Initialization
    const size_t numVertices = ioFloatingPositions.rows();

    //# Viscous Part
    //## The 'Force Field' is what drives the deformation: the difference between
    //## the floating vertices and their correspondences. By regulating it,
    //## viscous behaviour is obtained.
    //### Compute the Force Field
    Vec3Mat floatingPositions = ioFloatingPositions.leftCols(3);
    Vec3Mat forceField = inCorrespondingPositions.leftCols(3) - ioFloatingPositions.leftCols(3);


    //### Regulate the Force Field (Gaussian smoothing, iteratively)
    Vec3Mat regulatedForceField = Vec3Mat::Zero(numVertices,3);
    for (size_t i = 0 ; i < paramNumViscousSmoothingIterations ; i++) {
        //### Smooth the force field and save result in 'regulatedForceField'
        gaussian_smoothing_vector_field<Vec3Float, Vec3Mat> (forceField,
                                                            floatingPositions,
                                                            inFloatingWeights,
                                                            regulatedForceField,
                                                            paramNumNeighbourDisplacements,
                                                            paramSigmaSmoothing);

        //### Copy the result into forceField again (in case of more iterations)
        forceField = regulatedForceField;
    }


    //# Elastic Part
    //## Add the regulated Force Field to the current Displacement Field that has
    //## to be updated.
    Vec3Mat unregulatedDisplacementField = ioDisplacementField + regulatedForceField;

    //## Regulate the new Displacement Field (Gaussian smoothing, iteratively)
    for (size_t i = 0 ; i < paramNumElasticSmoothingIterations ; i++) {
        //### Smooth the force field and save result in 'ioDisplacementField'
        gaussian_smoothing_vector_field<Vec3Float, Vec3Mat> (unregulatedDisplacementField, Vec3Mat(ioFloatingPositions.leftCols(3)),
                                        inFloatingWeights, ioDisplacementField,
                                        paramNumNeighbourDisplacements, paramSigmaSmoothing);

        //### Copy the result into unregulatedDisplacementField again (in case
        //### more iterations need to be performed)
        unregulatedDisplacementField = ioDisplacementField;
    }

    //# Apply the transformation to the floating features
    for (size_t i = 0 ; i < numVertices ; i++) {
        ioFloatingPositions.row(i).head(3) += unregulatedDisplacementField.row(i);
    }

}//end viscoelastic_transformation()


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
//    floatingFeatures = FeatureMat::Zero(8, NUM_FEATURES);
//    targetFeatures = FeatureMat::Zero(8, NUM_FEATURES);
//    floatingFeatures << 0.1f, 0.1f, 0.1f, 0.0f, 0.0f, 1.0f,
//                        0.1f, 1.1f, 0.1f, 0.0f, 0.0f, 1.0f,
//                        1.1f, 0.1f, 0.1f, 0.0f, 0.0f, 1.0f,
//                        1.1f, 1.1f, 0.1f, 0.0f, 0.0f, 1.0f,
//                        0.1f, 0.1f, 1.1f, 0.0f, 0.0f, 1.0f,
//                        0.1f, 1.1f, 1.1f, 0.0f, 0.0f, 1.0f,
//                        1.1f, 0.1f, 1.1f, 0.0f, 0.0f, 1.0f,
//                        1.1f, 1.1f, 1.1f, 0.0f, 0.0f, 1.0f;
//    targetFeatures << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
//                      0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
//                      1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
//                      1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
//                      0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
//                      0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
//                      1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
//                      1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f;




     /*
    ############################################################################
    ##############################  TEST SHIZZLE  ##############################
    ############################################################################
    */
//    MatDynInt indices;
//    MatDynFloat distancesSquared;
//    radius_nearest_neighbours(floatingFeatures, targetFeatures, indices, distancesSquared, 1.0, 15);
//    cout << "Indices found: \n" << indices << endl;
//    cout << "distances squared found: \n" << distancesSquared << endl;

    /*
    ############################################################################
    ##############################  RIGID ICP  #################################
    ############################################################################
    */

    //# Info & Initialization
    //## Data and matrices
    size_t numFloatingVertices = floatingFeatures.rows();
    size_t numTargetVertices = targetFeatures.rows();
    VecDynFloat targetFlags = VecDynFloat::Ones(numTargetVertices);
    FeatureMat correspondingFeatures = FeatureMat::Zero(numFloatingVertices, NUM_FEATURES);
    VecDynFloat correspondingFlags = VecDynFloat::Ones(numFloatingVertices);
    //## Parameters
    const size_t numNearestNeighbours = 3;
    const size_t numRigidIterations = 20;

    //# Loop ICP
    for (size_t iteration = 0 ; iteration < numRigidIterations ; iteration++) {
        //# Determine Correspondences
        //## Compute symmetric wknn correspondences
        wkkn_correspondences(floatingFeatures, targetFeatures, targetFlags,
                            correspondingFeatures, correspondingFlags,
                            numNearestNeighbours, true);


        //# Inlier Detection
        VecDynFloat floatingWeights = VecDynFloat::Ones(numFloatingVertices);
        inlier_detection(floatingFeatures, correspondingFeatures,
                        correspondingFlags, floatingWeights, 3.0);

        //# Compute the transformation
        rigid_transformation(floatingFeatures, correspondingFeatures, floatingWeights, false);
    }
    /*
    ############################################################################
    ##########################  NON-RIGID ICP  #################################
    ############################################################################
    */
    const size_t numNonrigidIterations = 10;
    size_t smoothingIterations = numNonrigidIterations + 1;
    for (size_t i = 0 ; i < numNonrigidIterations ; i++) {
        //# Determine Correspondences
        //## Compute symmetric wknn correspondences
        wkkn_correspondences(floatingFeatures, targetFeatures, targetFlags,
                            correspondingFeatures, correspondingFlags,
                            numNearestNeighbours, true);


        //# Inlier Detection
        VecDynFloat floatingWeights = VecDynFloat::Ones(numFloatingVertices);
        inlier_detection(floatingFeatures, correspondingFeatures,
                        correspondingFlags, floatingWeights, 3.0);

        //# Visco-Elastic transformation
        Vec3Mat displacementField = Vec3Mat::Zero(numFloatingVertices, 3);
        viscoelastic_transformation(floatingFeatures,correspondingFeatures,
                                    floatingWeights, displacementField,
                                    10, 1.0, smoothingIterations, smoothingIterations);
        smoothingIterations--;
    }

    /*
    ############################################################################
    ##############################  OUTPUT #####################################
    ############################################################################
    */
    //# Write result to file
    write_eigen_features_to_obj(floatingFeatures, fuckedUpBunny, fuckedUpBunnyResultDir);


    return 0;
}

} //namespace registration
