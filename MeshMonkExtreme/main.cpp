#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include <Eigen/Dense>
#include <stdio.h>
#include <RigidRegistration.hpp>
#include <NonrigidRegistration.hpp>
#include <InlierDetector.hpp>
#include <CorrespondenceFilter.hpp>
#include <SymmetricCorrespondenceFilter.hpp>
#include <RigidTransformer.hpp>
#include <ViscoElasticTransformer.hpp>
#include "global.hpp"
#include <helper_functions.hpp>

using namespace std;


typedef OpenMesh::DefaultTraits MyTraits;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  TriMesh;
typedef OpenMesh::Decimater::DecimaterT<TriMesh>    DecimaterType;
typedef OpenMesh::Decimater::ModQuadricT<TriMesh>::Handle HModQuadric;
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat; //matrix Mx3 of type unsigned int
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float

int main()
{

    /*
    ############################################################################
    ##############################  INPUT  #####################################
    ############################################################################
    */
    //# IO variables
//    const string fuckedUpBunnyDir = "/home/jonatan/meshmonk/examples/data/bunny_slightly_rotated.obj";
    const string fuckedUpBunnyDir = "/home/jonatan/projects/meshmonk/examples/data/fucked_up_bunny.obj";
    const string bunnyDir = "/home/jonatan/projects/meshmonk/examples/data/bunny90.obj";
    const string fuckedUpBunnyResultDir = "/home/jonatan/projects/meshmonk/examples/data/bunnyNonRigid.obj";
//    const string fuckedUpBunnyDir = "/home/jonatan/projects/meshmonk/examples/data/kul_gezichten/Outliers/Alspac/4707_template.obj";
//    const string bunnyDir = "/home/jonatan/projects/meshmonk/examples/data/kul_gezichten/Outliers/Alspac/4707_mevislab.obj";
//    const string fuckedUpBunnyResultDir = "/home/jonatan/projects/meshmonk/examples/data/kul_gezichten/Outliers/Alspac/4707_Monk.obj";

    //# Load meshes and convert to feature matrices
    FeatureMat floatingFeatures;
    FeatureMat targetFeatures;
    FacesMat floatingFaces;
    registration::import_data(fuckedUpBunnyDir, bunnyDir, floatingFeatures,
                              targetFeatures, floatingFaces);







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

//    std::cout << "test features: \n" << floatingFeatures.topRows(10) << std::endl;
//    std::cout << "test faces: \n" << floatingFaces.topRows(10) << std::endl;



    /*
    ############################################################################
    ##############################  DECIMATION  ################################
    ############################################################################
    */

    //# Load the bunny
    TriMesh fuckedUpBunny;
    OpenMesh::IO::_OBJReader_();
    if (!OpenMesh::IO::read_mesh(fuckedUpBunny,fuckedUpBunnyDir)){
        std::cerr << "Read error \n";
        exit(1);
    };

    //# Print Mesh property
    std::cout << "--------  Before processing " << std::endl;
    std::cout << "# Vertices " << fuckedUpBunny.n_vertices() << std::endl;
    std::cout << "# Edges " << fuckedUpBunny.n_edges() << std::endl;
    std::cout << "# Faces " << fuckedUpBunny.n_faces() << std::endl;

    //# block boundary vertices
    fuckedUpBunny.request_vertex_status();
    //## Get an iterator over all halfedges
    TriMesh::HalfedgeIter he_it, he_end=fuckedUpBunny.halfedges_end();
    //## If halfedge is boundary, lock the corresponding vertices
    for (he_it = fuckedUpBunny.halfedges_begin(); he_it != he_end ; ++he_it) {
      if (fuckedUpBunny.is_boundary(*he_it)) {
         fuckedUpBunny.status(fuckedUpBunny.to_vertex_handle(*he_it)).set_locked(true);
         fuckedUpBunny.status(fuckedUpBunny.from_vertex_handle(*he_it)).set_locked(true);
      }
    }
    //# Make sure mesh has necessary normals etc
    fuckedUpBunny.request_face_normals();
    fuckedUpBunny.update_face_normals();

    //# Set up the decimator
    DecimaterType decimater(fuckedUpBunny);  // a decimater object, connected to a mesh
    HModQuadric hModQuadric;      // use a quadric module
    bool addSucces = decimater.add( hModQuadric ); // register module at the decimater
    std::cout << "Adding quadric modifier to decimater : " << addSucces << std::endl;

    std::cout << decimater.module( hModQuadric ).name() << std::endl;
    decimater.module(hModQuadric).unset_max_err();

    //# Initialize the decimater
    bool rc = decimater.initialize();
    std::cout  << "Decimater Initialization: " << decimater.is_initialized() << std::endl;
    if (!rc){
        std::cerr << "  initializing failed!" << std::endl;
        std::cerr << "  maybe no priority module or more than one were defined!" << std::endl;
        return false;
    }


    std::cout << "Decimater Observer: " << decimater.observer() << std::endl;



    //# Run the decimater
//    rc = decimater.decimate_to(size_t(10));
    rc = decimater.decimate(size_t(10));
    std::cout << " Did we fucking decimate?? -> " << rc << std::endl;

    //# Collect garbage
    fuckedUpBunny.garbage_collection();

    //# Print Mesh property
    std::cout << "--------  After processing " << std::endl;
    std::cout << "# Vertices " << fuckedUpBunny.n_vertices() << std::endl;
    std::cout << "# Edges " << fuckedUpBunny.n_edges() << std::endl;
    std::cout << "# Faces " << fuckedUpBunny.n_faces() << std::endl;



    //# Write to file
    OpenMesh::IO::_OBJWriter_();
    if (!OpenMesh::IO::write_mesh(fuckedUpBunny, fuckedUpBunnyResultDir))
    {
        std::cerr << "write error\n";
        exit(1);
    }
//
//    //# convert to our features
//    registration::convert_mesh_to_matrices(fuckedUpBunny, floatingFeatures, floatingFaces);


    /*
    ############################################################################
    ##############################  RIGID ICP  #################################
    ############################################################################
    */

//    //# Info & Initialization
//    //## Data and matrices
//    size_t numFloatingVertices = floatingFeatures.rows();
//    size_t numTargetVertices = targetFeatures.rows();
//    VecDynFloat floatingFlags = VecDynFloat::Ones(numFloatingVertices);
//    VecDynFloat targetFlags = VecDynFloat::Ones(numTargetVertices);
//
//    //## Parameters
//    bool symmetric = true;
//    size_t numNearestNeighbours = 5;
//    size_t numRigidIterations = 10;
//    float kappa = 3.0f;
//
//    registration::RigidRegistration rigidRegistration;
//    rigidRegistration.set_input(&floatingFeatures, &targetFeatures, &floatingFlags, &targetFlags);
//    rigidRegistration.set_parameters(symmetric, numNearestNeighbours, kappa, numRigidIterations);
//    rigidRegistration.update();
//
//    /*
//    ############################################################################
//    ##########################  NON-RIGID ICP  #################################
//    ############################################################################
//    */
//    //## Initialization
//    size_t numNonrigidIterations = 20;
//    float sigmaSmoothing = 2.0f;
//    registration::NonrigidRegistration nonrigidRegistration;
//    nonrigidRegistration.set_input(&floatingFeatures, &targetFeatures, &floatingFaces, &floatingFlags, &targetFlags);
//    nonrigidRegistration.set_parameters(symmetric, numNearestNeighbours, kappa,
//                                        numNonrigidIterations, sigmaSmoothing,
//                                        100, 1,
//                                        100, 1);
//    nonrigidRegistration.update();

    /*
    ############################################################################
    ##############################  OUTPUT #####################################
    ############################################################################
    */
    //# Write result to file
//    registration::export_data(floatingFeatures,floatingFaces, fuckedUpBunnyResultDir);


    return 0;
}
