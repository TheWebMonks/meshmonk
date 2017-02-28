#include "Downsampler.hpp"


namespace registration {



void Downsampler::set_input(const FeatureMat * const inFeatures,
                            const FacesMat * const inFaces,
                            const VecDynFloat * const inFlags){
    _inFeatures = inFeatures;
    _inFaces = inFaces;
    _inFlags = inFlags;

}//end set_input()


void Downsampler::set_output(FeatureMat &outFeatures,
                            FacesMat &outFaces,
                            VecDynFloat &outFlags,
                            VecDynInt &outOriginalIndices){

    _outFeatures = &outFeatures;
    _outFaces = &outFaces;
    _outFlags = &outFlags;
    _outOriginalIndices = &outOriginalIndices;
}//end set_output()


void Downsampler::set_output(FeatureMat &outFeatures,
                            FacesMat &outFaces,
                            VecDynFloat &outFlags){

    _outFeatures = &outFeatures;
    _outFaces = &outFaces;
    _outFlags = &outFlags;
    _outOriginalIndices = NULL;
}//end set_output()




void Downsampler::update(){
    //# Convert the input data to OpenMesh's mesh structure
    TriMesh mesh;
    convert_matrices_to_mesh(*_inFeatures,
                            *_inFaces,
                            *_inFlags,
                            mesh);

    //# Add the original indices as a custom property to each vertex
    //## This is so that we can extract the original index of each vertex
    //## after downsampling.
    OpenMesh::VPropHandleT<float> originalIndices;
    mesh.add_property(originalIndices, "originalIndices");
    mesh.property(originalIndices).set_persistent(true);
    TriMesh::VertexIter vertexIt(mesh.vertices_begin());
    TriMesh::VertexIter vertexEnd(mesh.vertices_end());
    for (size_t index = 0 ; vertexIt != vertexEnd ; index++, vertexIt++) {
        mesh.property(originalIndices, vertexIt) = index;
    }

    //# Downsample using OpenMesh library
    const size_t numOriginalVertices = mesh.n_vertices();
    const size_t numOriginalEdges = mesh.n_edges();
    const size_t numOriginalFaces = mesh.n_faces();
    //## block boundary vertices
    mesh.request_vertex_status();
    //### Get an iterator over all halfedges
    TriMesh::HalfedgeIter he_it, he_end=mesh.halfedges_end();
    //### If halfedge is boundary, lock the corresponding vertices
    for (he_it = mesh.halfedges_begin(); he_it != he_end ; ++he_it) {
      if (mesh.is_boundary(*he_it)) {
         mesh.status(mesh.to_vertex_handle(*he_it)).set_locked(true);
         mesh.status(mesh.from_vertex_handle(*he_it)).set_locked(true);
      }
    }
    //## Make sure mesh has necessary normals etc
    mesh.request_face_normals();
    mesh.update_face_normals();

    //## Set up the decimator
    DecimaterType decimater(mesh);  // a decimater object, connected to a mesh
    HModQuadric hModQuadric;      // use a quadric module
    bool addSucces = decimater.add( hModQuadric ); // register module at the decimater
    decimater.module(hModQuadric).unset_max_err();

    //## Initialize the decimater
    bool rc = decimater.initialize();
    if (!rc){
        std::cerr << "  initializing failed!" << std::endl;
        std::cerr << "  maybe no priority module or more than one were defined!" << std::endl;
        return;
    }



    //## Run the decimater
    size_t numVertices = mesh.n_vertices();
    size_t numDecimations = size_t(round(_downsampleRatio * numVertices));
    std::cout << "DEBUG - downsampleratio : " << _downsampleRatio << " || numdecimations : " << numDecimations << std::endl;
    if (numDecimations > 0 ){
        rc = decimater.decimate(numDecimations);
    }

    //## Collect garbage
    mesh.garbage_collection();

    //## Print Mesh property
    numVertices = mesh.n_vertices();
    const size_t numEdges = mesh.n_edges();
    const size_t numFaces = mesh.n_faces();
    if (rc){
        std::cout << "Downsampled mesh from " << numOriginalVertices << "v | " << numOriginalEdges << "e | " << numOriginalFaces
        << "f to " << numVertices << "v | " << numEdges << "e | " << numFaces << "f." << std::endl;
    }
    else{
        std::cerr << "DOWNSAMPLING FAILED !" << std::endl;
    }

    //# Convert the downsampled result to the output matrices
    convert_mesh_to_matrices(mesh, *_outFeatures, *_outFaces, *_outFlags);

    //# Extract the original indices.
    //## Get handle to the originalIndices property.
    if (_outOriginalIndices != NULL) {
        *_outOriginalIndices = VecDynInt::Zero(numVertices); //resizing
        bool propertyExist = mesh.get_property_handle(originalIndices, "originalIndices");
        if (!propertyExist)
        {
            std::cerr << "Tried to access the 'originalIndices' property of the mesh after downsampling - couldn't find handle\n";
            exit(1);
        }

        //## Loop over the vertices and save its original index
        vertexIt = mesh.vertices_begin();
        vertexEnd = mesh.vertices_end();
        //## Loop over every face
        for (size_t i = 0 ; vertexIt != vertexEnd ; i++, vertexIt++) {
                (*_outOriginalIndices)[i] = mesh.property(originalIndices, vertexIt);
        }
    }
}


}//namespace registration
