#ifndef SCALESHIFTER_HPP
#define SCALESHIFTER_HPP


#include <Eigen/Dense>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "../global.hpp"
#include <helper_functions.hpp>

typedef Eigen::Vector3f Vec3Float;
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat;
typedef OpenMesh::DefaultTraits MyTraits;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  TriMesh;

namespace registration{

class ScaleShifter
{
    public:

        void set_input(const FeatureMat &inLowFeatures,
                       const VecDynInt &inLowOriginalIndices;
                       const VecDynInt &inHighOriginalIndices);
        void set_output(FeatureMat &outHighFeatures);
        void update();

    protected:

    private:
        //# Inputs
        const FeatureMat * _inLowFeatures = NULL;
        const VecDynInt *_inLowOriginalIndices = NULL;
        const VecDynInt *_inHighOriginalIndices = NULL;

        //# Outputs
        FeatureMat * _outHighFeatures;


        //# User Parameters

        //# Internal Data structures

        //# Internal Parameters

        //# Internal functions
};

}//namespace registration

#endif // SCALESHIFTER_HPP
