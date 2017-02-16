#ifndef DATAIMPORTER_HPP
#define DATAIMPORTER_HPP

#include <iostream>
#include <fstream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <Eigen/Dense>
#include "../global.hpp"
#include <helper_functions.hpp>

typedef OpenMesh::DefaultTraits MyTraits;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  TriMesh;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat;


namespace registration {

class DataImporter
{
    /*
    # GOAL
    This class takes filenames as input and prepares data that is required
    by the registration framework.

    # INPUTS
    -inFloatingMeshPath
    -inTargetMeshPath

    # PARAMETERS

    # OUTPUTS
    -outFloatingFeatures
    -outTargetFeatures
    -outFloatingFaces
    -outTargetFaces
    */
    public:

        void set_input(std::string floatingMeshPath, std::string targetMeshPath);
        void set_output(FeatureMat &floatingFeatures, FeatureMat &targetFeatures,
                        FacesMat &floatingFaces);
        void update();

    protected:

    private:
        //# Inputs
        std::string _inFloatingMeshPath;
        std::string _inTargetMeshPath;

        //# Outputs
        FeatureMat _outFloatingFeatures;
        FeatureMat _outTargetFeatures;
        FacesMat _outFloatingFaces;
        //FacesMat _outTargetFaces;

        //# User Parameters

        //# Internal Data structures

        //# Internal Parameters
        bool _inputfilesExist = false;
        //# Internal functions
        //## Safety check to see if input files exist
        bool _safety_check();
};

}//namespace registration

#endif // DATAIMPORTER_HPP
