#ifndef DOWNSAMPLER_HPP
#define DOWNSAMPLER_HPP


#include <Eigen/Dense>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
#include "../global.hpp"
#include <helper_functions.hpp>

typedef Eigen::Vector3f Vec3Float;
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::VectorXi VecDynInt;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat;
typedef OpenMesh::DefaultTraits MyTraits;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  TriMesh;
typedef OpenMesh::Decimater::DecimaterT<TriMesh>    DecimaterType;
typedef OpenMesh::Decimater::ModQuadricT<TriMesh>::Handle HModQuadric;

namespace registration{

class Downsampler
{
    public:

        void set_input(const FeatureMat * const inFeatures,
                       const FacesMat * const inFaces,
                       const VecDynFloat * const inFlags);
        void set_output(FeatureMat &outFeatures,
                        FacesMat &outFaces,
                        VecDynFloat &outFlags,
                        VecDynInt &outOriginalIndices);
        void set_output(FeatureMat &outFeatures,
                        FacesMat &outFaces,
                        VecDynFloat &outFlags);
        void set_parameters(float downsampleRatio = 0.8f){ _downsampleRatio = downsampleRatio;};
        void update();

    protected:

    private:
        //# Inputs
        const FeatureMat * _inFeatures = NULL;
        const FacesMat * _inFaces = NULL;
        const VecDynFloat * _inFlags = NULL;

        //# Outputs
        FeatureMat * _outFeatures = NULL;
        FacesMat * _outFaces = NULL;
        VecDynFloat * _outFlags = NULL;
        VecDynInt *_outOriginalIndices = NULL;


        //# User Parameters
        float _downsampleRatio = 0.8f; //must be between 0.0 and 1.0

        //# Internal Data structures

        //# Internal Parameters

        //# Internal functions
};

}//namespace registration

#endif // DOWNSAMPLER_HPP
