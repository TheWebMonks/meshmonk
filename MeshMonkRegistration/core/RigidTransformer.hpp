#ifndef RIGIDTRANSFORMER_HPP
#define RIGIDTRANSFORMER_HPP

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <stdio.h>
#include <iostream>
#include "../global.hpp"

typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat; //matrix MxN of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Vector3f Vec3Float;
typedef Eigen::Vector4f Vec4Float;
typedef Eigen::Matrix3f Mat3Float;
typedef Eigen::Matrix4f Mat4Float;
typedef Eigen::SelfAdjointEigenSolver<Mat4Float> EigenVectorDecomposer;

namespace registration{

class RigidTransformer
{
    /*
    # GOAL
    This class computes the rigid transformation between a set a features and
    a set of corresponding features. Each correspondence can be weighed between
    0.0 and 1.0.

    # INPUTS
    -ioFeatures
    -inCorrespondingFeatures
    -inWeights

    # PARAMETERS
    -scaling:
    Whether or not to allow scaling.

    # OUTPUTS
    -ioFeatures
    */
    public:

        void set_input(const FeatureMat * const inCorrespondingFeatures, const VecDynFloat * const inWeights);
        void set_output(FeatureMat * const ioFeatures);
        void set_parameters(bool scaling);
        Mat4Float get_transformation() const {return _transformationMatrix;}
        void update();

    protected:

    private:
        //# Inputs
        FeatureMat * _ioFeatures = NULL;
        const FeatureMat * _inCorrespondingFeatures = NULL;
        const VecDynFloat * _inWeights = NULL;

        //# Outputs
        //_ioFeatures is used as both an input (to compute the transformation) and output

        //# User Parameters
        bool _scaling = false;

        //# Internal Data structures
        Mat4Float _transformationMatrix = Mat4Float::Identity();

        //# Internal Parameters
        size_t _numElements = 0;
        size_t _numFeatures = 0;

        //# Internal functions
        //## Function to update the transformation matrix and apply it to the floating positions
        void _update_transformation();
};

}//namespace registration

#endif // RIGIDTRANSFORMER_HPP
