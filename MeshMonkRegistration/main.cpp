#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <Eigen/Dense>
#include <nanoflann.hpp>

using namespace std;

const int NUM_FEATURES = 6;

typedef OpenMesh::TriMesh_ArrayKernelT<>  TriMesh;
typedef Eigen::Matrix< unsigned int, Eigen::Dynamic, Eigen::Dynamic> MatMxNi; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatMxNf; //matrix MxN of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, 3> MatMx3f; //matrix Mx3 of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, NUM_FEATURES> MatMx6f; //matrix Mx6 of type float
typedef nanoflann::KDTreeEigenMatrixAdaptor<MatMx6f, NUM_FEATURES, nanoflann::metric_L2>  NanoKDTree6D;



void k_nearest_neighbours(const MatMx6f &inQueriedPoints,
                        const MatMx6f &inSourcePoints,
                        MatMxNi &outNeighbourIndices,
                        MatMxNf &outNeighbourSquaredDistances,
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
    outNeighbourIndices = MatMxNi::Zero(numQueriedElements,paramK);
    outNeighbourSquaredDistances = MatMxNf::Zero(numQueriedElements,paramK);

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
                                MatMx6f &outputFeatures){
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
    outputFeatures = MatMx6f::Zero(numVertices,NUM_FEATURES);

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
                                MatMx6f &outFeatureMatrix){
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

void eigen_features_to_openmesh(const MatMx6f &inFeatures, TriMesh &outMesh){
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

void write_eigen_features_to_obj(const MatMx6f &inFeatures,
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

void update_normals_for_altered_positions(TriMesh &ioMesh, MatMx6f &ioFeatures){
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

int main()
{
    //# IO options
    OpenMesh::IO::Options readOptions;
    //# IO variables
    TriMesh fuckedUpBunny;
    TriMesh bunny;
    const string fuckedUpBunnyDir = "/home/jonatan/kuleuven-algorithms/examples/data/bunny_slightly_rotated.obj";
    const string bunnyDir = "/home/jonatan/kuleuven-algorithms/examples/data/bunny90.obj";
    const string fuckedUpBunnyResultDir = "/home/jonatan/kuleuven-algorithms/examples/data/fucked_up_bunny_result.obj";

    //# Load meshes and convert to feature matrices
    MatMx6f floatingFeatures;
    MatMx6f targetFeatures;
    load_obj_to_eigen_features(fuckedUpBunnyDir, fuckedUpBunny, floatingFeatures);
    load_obj_to_eigen_features(bunnyDir, bunny, targetFeatures);

    //# Do rigid ICP
    //## Find nearest neighbours
    const size_t k = 2;
    MatMxNi neighbourIndices;
    MatMxNf neighbourSquaredDistances;
    k_nearest_neighbours(floatingFeatures, targetFeatures, neighbourIndices, neighbourSquaredDistances, k, 15);






    //TEST: ASSIGN FIRST NEIGHBOURS AS NEW LOCATIONS
    for (unsigned int i = 0 ; i < floatingFeatures.rows() ; i++) {
        floatingFeatures.row(i) = targetFeatures.row(neighbourIndices(i,0));
    }

    //# Do non-rigid ICP

    //# Write result to file
    write_eigen_features_to_obj(floatingFeatures, fuckedUpBunny, fuckedUpBunnyResultDir);


    return 0;
}
