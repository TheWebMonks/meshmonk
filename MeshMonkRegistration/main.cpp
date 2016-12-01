#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <Eigen/Dense>
#include <nanoflann.hpp>

using namespace std;

typedef OpenMesh::TriMesh_ArrayKernelT<>  TriMesh;
typedef Eigen::Matrix< unsigned int, Eigen::Dynamic, Eigen::Dynamic> MatMxNi; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatMxNf; //matrix MxN of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, 3> MatMx3f; //matrix Mx3 of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, 6> MatMx6f; //matrix Mx6 of type float
typedef nanoflann::KDTreeEigenMatrixAdaptor<MatMx6f, 6, nanoflann::metric_L2>  NanoKDTree6D;

void k_nearest_neighbours(const MatMx6f &inQueriedFeatures,
                        const MatMx6f &inSourceFeatures,
                        MatMxNi &outNeighbourIndices,
                        MatMxNf &outNeighbourSquaredDistances,
                        const size_t paramK = 3,
                        const size_t paramLeafsize = 15){
    /*
    GOAL
    This function searches for the k nearest neighbours in 'inSourceFeatures' for
    each element in the 'inQueriedFeatures' set. It outputs the indices of each
    neighbour and the squared (!) distances between each element of
    'inQueriedFeatures' and its neighbours.

    INPUT
    -inQueriedFeatures:
    -inSourceFeatures:

    PARAMETERS
    -paramK(= 3): number of nearest neighbours
    -paramLeafsize(= 15): should be between 5 and 50 or so

    OUTPUT
    -outNeighbourIndices
    -outNeighbourSquaredDistances.
    */

    //# Info and Initialization
    const size_t dimension = inSourceFeatures.cols();
    const size_t numSourceElements = inSourceFeatures.rows();
    const size_t numQueriedElements = inQueriedFeatures.rows();
    cout << "dimension [" << dimension << "] - numSourceElements [" <<
    numSourceElements << "] - numQueriedElements [" << numQueriedElements << "]"
    << endl;
    cout << "inQueriedFeatures block: \n" << inQueriedFeatures.block<10,6>(0,0) << endl;
    cout << "inSourceFeatures block: \n" << inSourceFeatures.block<10,6>(0,0) << endl;

    //# Construct kd-tree
    NanoKDTree6D kdTree(dimension, inSourceFeatures, paramLeafsize);
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
    knnResultSet.init(&neighbourIndices[0], &neighbourSquaredDistances[0]);

    //### Execute loop
    for ( ; i < numQueriedElements ; ++i ) {
        //### Copy input features to 'queriedFeature'
        for (j = 0 ; j < dimension ; ++j) {
            queriedFeature[j] = inQueriedFeatures(i,j);
        }

        if (i < 3) {
            cout << "Step " << i << ": queried feature: " << queriedFeature[0] << endl;
            cout << "Step " << i << ": pre-knn neighbourIndices: " << neighbourIndices[0] << endl;
        }

        //### Query the kd-tree
        kdTree.index->findNeighbors(knnResultSet, &queriedFeature[0], nanoflann::SearchParams(32, 0.0001, true));
        if (i < 3) {
            cout << "Step " << i << ": post-knn neighbourIndices: " << neighbourIndices[0] << endl;
        }
        //### Copy the result into the outputs by looping over the k nearest
        //### neighbours
        for (j = 0 ; j < paramK ; ++j) {
            outNeighbourIndices(i,j) = neighbourIndices[j];
            outNeighbourSquaredDistances(i,j) = neighbourSquaredDistances[j];
        }
    }
}//end k_nearest_neighbours()

void openmesh_to_eigen_features(const TriMesh &inputMesh, MatMx6f &outputFeatures){
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
    const int numRows = outputFeatures.rows();
    if (numVertices != numRows) {
    cerr << "Number of vertices does not correspond with number of columns when"
    << " calling openmesh_to_numpy_features" << endl;
    cerr << "numVertices: " << numVertices << endl;
    cerr << "numRows: " << numRows << endl;
    }

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

void eigen_features_to_openmesh(const MatMx6f &inputFeatures, TriMesh &outputMesh){
    /*
    GOAL
    This function converts a feature representation of a mesh using Eigen
    matrices into a mesh representation using OpenMesh.

    INPUT
    -inputFeatures:
    a numVertices x 6 Eigen dense matrix where the first three columns are made
    up of the positions of the vertices, and the last three columns the normals
    of those vertices.

    PARAMETERS

    OUTPUT
    -outputMesh:
    this has to be a mesh of openmesh's TriMesh type
    */

    //# Info and Initialization
    const int numRows = inputFeatures.rows();
    const int numVertices = outputMesh.n_vertices();
     if (numRows != numVertices) {
    cerr << "Number of rows does not correspond with number of vertices when"
    << " calling eigen_features_to_openmesh()" << endl;
    }

    //# Put positions and normals back into mesh
    TriMesh::Point updatedPosition(0.0f,0.0f,0.0f);
    TriMesh::Normal updatedNormal(0.0f,0.0f,0.0f);
    unsigned int i = 0;
    TriMesh::VertexIter vertexIt(outputMesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(outputMesh.vertices_end());
    for ( ; vertexIt != vertexEnd ; ++i, ++vertexIt) {
        //## Get positions and normals from the feature matrix
        updatedPosition[0] = inputFeatures(i,0);
        updatedPosition[1] = inputFeatures(i,1);
        updatedPosition[2] = inputFeatures(i,2);
        updatedNormal[0] = inputFeatures(i,3);
        updatedNormal[1] = inputFeatures(i,4);
        updatedNormal[2] = inputFeatures(i,5);
        //## Insert position and normal into mesh
        outputMesh.set_point(vertexIt,updatedPosition);
        outputMesh.set_normal(vertexIt, updatedNormal);
    }
}

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
    // IO options
    OpenMesh::IO::Options readOptions;
    // IO variables
    TriMesh bunny;
    TriMesh fuckedUpBunny;
    const string bunnyDir = "/home/jonatan/kuleuven-algorithms/examples/data/bunny90.obj";
    const string fuckedUpBunnyDir = "/home/jonatan/kuleuven-algorithms/examples/data/bunny90.obj";
    const string fuckedUpBunnyResultDir = "/home/jonatan/kuleuven-algorithms/examples/data/fucked_up_bunny_result.obj";

    // Load meshes
    OpenMesh::IO::_OBJReader_(); //NOTE: this issue should only appear when STATICALLY linking OpenMesh
    if (!OpenMesh::IO::read_mesh(bunny,bunnyDir)){
        cerr << "Read error \n";
        exit(1);
    };
    if (!OpenMesh::IO::read_mesh(fuckedUpBunny,fuckedUpBunnyDir)){
        cerr << "Read error \n";
        exit(1);
    };

    // Update normals
    fuckedUpBunny.request_face_normals();
    fuckedUpBunny.request_vertex_normals();
    bunny.request_face_normals();
    bunny.request_vertex_normals();
    fuckedUpBunny.update_face_normals();
    bunny.update_face_normals();
    fuckedUpBunny.update_vertex_normals();
    bunny.update_vertex_normals();


    //# Convert to matrices (from Eigen library)
    const int numFloatingVertices = fuckedUpBunny.n_vertices();
    const int numTargetVertices = bunny.n_vertices();
    MatMx6f floatingFeatures = MatMx6f(numFloatingVertices,6);
    MatMx6f targetFeatures = MatMx6f(numTargetVertices,6);
    cout << "Number floating vertices: " << numFloatingVertices << endl;
    cout << "Number target vertices: " << numTargetVertices << endl;

    //## Call conversion function
    openmesh_to_eigen_features(fuckedUpBunny, floatingFeatures);
    openmesh_to_eigen_features(bunny, targetFeatures);

    //# Do rigid ICP
    //## Find nearest neighbours
    const size_t k = 3;
    MatMxNi neighbourIndices = MatMxNi::Zero(numFloatingVertices,k);
    MatMxNf neighbourSquaredDistances = MatMxNf::Zero(numFloatingVertices,k);
    k_nearest_neighbours(floatingFeatures, targetFeatures, neighbourIndices, neighbourSquaredDistances, k, 15);

    cout << "Post kNN: " << neighbourIndices.block<10,3>(0,0) << endl;





//	cout << "knnSearch(nn="<<k<<"): \n";
//	for (unsigned int i=0;i<k;i++) {
//		cout << "neighbourIndices["<<i<<"]=" << neighbourIndices[i] << " neighbourSquaredDistances=" << neighbourSquaredDistances[i] << endl;
//    }



    //# Do non-rigid ICP

    //# Convert matrices to OpenMesh structure
    eigen_features_to_openmesh(floatingFeatures, fuckedUpBunny);
    //# Write resulting meshes to file
    OpenMesh::IO::_OBJWriter_();
    if (!OpenMesh::IO::write_mesh(fuckedUpBunny, fuckedUpBunnyResultDir))
    {
        std::cerr << "write error\n";
        exit(1);
    }

    return 0;
}
