#include "DataImporter.hpp"

namespace registration{



void DataImporter::update(){
    std::cout << "Importing Data..." << std::endl;

    //# Safety check (do the input files exist?)
    if (!_safety_check()){
        std::cout<< "DataImporter can't update - input files not found!" << std::endl;
        return;
    }

    //# Load meshes from file
    //## Floating Mesh
    TriMesh floatingMesh;
    OpenMesh::IO::_OBJReader_();
    if (!OpenMesh::IO::read_mesh(floatingMesh,_inFloatingMeshPath)){
        std::cerr << "Read error \n";
        exit(1);
    };
    //## Target Mesh
    TriMesh targetMesh;
    if (!OpenMesh::IO::read_mesh(targetMesh,_inTargetMeshPath)){
        std::cerr << "Read error \n";
        exit(1);
    };

    //# Convert the meshes to matrices
    convert_openmesh_to_eigen(floatingMesh, _outFloatingFeatures, _outFloatingFaces);
    convert_openmesh_to_eigen(targetMesh, _outTargetFeatures);


    std::cout << "Imported Data" << std::endl;
}



bool DataImporter::_safety_check(){
    _inputfilesExist = true;
    //#Use ifstream to try to open the input files
    //## Floating Mesh
    std::ifstream infile(_inFloatingMeshPath);
    if (infile.good() != true){
        _inputfilesExist = false;
        std::cerr << "Floating mesh file does not exist" << std::endl;
    }
    infile.close();

    //## Target Mesh
    infile.open(_inTargetMeshPath);
    if (infile.good() != true){
        _inputfilesExist = false;
        std::cerr << "Target mesh file does not exist" << std::endl;
    }
    infile.close();

    return _inputfilesExist;
}




}//namespace registration




//def update(self):
//        print("Importing data....")
//        # Safety check
//        if not os.path.isfile(self.inFloatingMeshPath):
//            raise IOError('Floating mesh path {:s} does not exist.'.format(self.inFloatingMeshPath))
//
//        if not os.path.isfile(self.inTargetMeshPath):
//            raise IOError('Target mesh path {:s} does not exist.'.format(self.inTargetMeshPath))
//
//        # Load from file
//        openmesh.read_mesh(self._floatingMesh, self.inFloatingMeshPath)
//        openmesh.read_mesh(self._targetMesh, self.inTargetMeshPath)
//        # Obtain the floating and target mesh features (= positions and normals)
//        # and their faces
//        self.outFloatingFeatures, self.outFloatingFaces = openmesh_to_matrices(self._floatingMesh)
//        self.outTargetFeatures, self.outTargetFaces = openmesh_to_matrices(self._targetMesh)
