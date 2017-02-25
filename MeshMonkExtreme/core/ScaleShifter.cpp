#include "ScaleShifter.hpp"

namespace registration {



void ScaleShifter::set_input(const FeatureMat &inLowFeatures,
                       const VecDynInt &inLowOriginalIndices;
                       const VecDynInt &inHighOriginalIndices){
    _inLowFeatures = &inLowFeatures;
    _inLowOriginalIndices = &inLowOriginalIndices;
    _inHighOriginalIndices = &inHighOriginalIndices;

}//end set_input()


void ScaleShifter::set_output(FeatureMat &outHighFeatures){
    _outHighFeatures = &outHighFeatures;
}//end set_output()


void ScaleShifter::_find_corresponding_indices(){
    //## Put the indices into arrays

    //## Sort the arrays

    //## Find the intersect between the arrays
}

void ScaleShifter::update(){

    //# Build the list of corresponding indices
    //##    (The values in this list refer to the indides of the high
    //##    resolution mesh. The index of the list itself corresponds
    //##    with the index of the low resolution mesh. That gives us
    //##    a link between the vertices of the two meshes that correspond
    //##    to each other.)

    //## Loop over the low resolution vertices
    size_t numLowVertices = _inLowFeatures->rows();
    size_t numHighVertices = _outHighFeatures->rows();
    for (size_t i = 0 ; i < numLowVertices ; i++){

    }

















    //# Convert the input data to OpenMesh's mesh structure
//    TriMesh inMesh;
//    convert_matrices_to_mesh(*_inFeatures,
//                            *_inFaces,
//                            *_inFlags,
//                            mesh);
//
//
//
//
//    //# Downsample using OpenMesh library
//    //## Print Mesh property
//    std::cout << "--------  Before processing " << std::endl;
//    std::cout << "# Vertices " << mesh.n_vertices() << std::endl;
//    std::cout << "# Edges " << mesh.n_edges() << std::endl;
//    std::cout << "# Faces " << mesh.n_faces() << std::endl;
//    //## block boundary vertices
//    mesh.request_vertex_status();
//    //### Get an iterator over all halfedges
//    TriMesh::HalfedgeIter he_it, he_end=mesh.halfedges_end();
//    //### If halfedge is boundary, lock the corresponding vertices
//    for (he_it = mesh.halfedges_begin(); he_it != he_end ; ++he_it) {
//      if (mesh.is_boundary(*he_it)) {
//         mesh.status(mesh.to_vertex_handle(*he_it)).set_locked(true);
//         mesh.status(mesh.from_vertex_handle(*he_it)).set_locked(true);
//      }
//    }
//    //## Make sure mesh has necessary normals etc
//    mesh.request_face_normals();
//    mesh.update_face_normals();
//
//    //## Set up the decimator
//    DecimaterType decimater(mesh);  // a decimater object, connected to a mesh
//    HModQuadric hModQuadric;      // use a quadric module
//    bool addSucces = decimater.add( hModQuadric ); // register module at the decimater
//    std::cout << "Adding quadric modifier to decimater : " << addSucces << std::endl;
//
//    std::cout << decimater.module( hModQuadric ).name() << std::endl;
//    decimater.module(hModQuadric).unset_max_err();
//
//    //## Initialize the decimater
//    bool rc = decimater.initialize();
//    std::cout  << "Decimater Initialization - Success. " << std::endl;
//    if (!rc){
//        std::cerr << "  initializing failed!" << std::endl;
//        std::cerr << "  maybe no priority module or more than one were defined!" << std::endl;
//        return;
//    }
//    std::cout << "Decimater Observer: " << decimater.observer() << std::endl;
//
//
//
//    //## Run the decimater
//    size_t numVertices = mesh.n_vertices();
//    size_t numDecimations = size_t(round(_ScaleShifteratio * numVertices));
//    rc = decimater.decimate(numDecimations);
//    std::cout << " Decimation Succes? -> " << rc << std::endl;
//
//    //## Collect garbage
//    mesh.garbage_collection();
//
//    //## Print Mesh property
//    std::cout << "--------  After processing " << std::endl;
//    std::cout << "# Vertices " << mesh.n_vertices() << std::endl;
//    std::cout << "# Edges " << mesh.n_edges() << std::endl;
//    std::cout << "# Faces " << mesh.n_faces() << std::endl;
//
//    //# Convert the downsampled result to the output matrices
//    convert_mesh_to_matrices(mesh, *_outFeatures, *_outFaces, *_outFlags);

}


}//namespace registration
