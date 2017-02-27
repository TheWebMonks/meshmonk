#ifndef SCALESHIFTER_HPP
#define SCALESHIFTER_HPP


#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "../global.hpp"
#include <helper_functions.hpp>

typedef Eigen::Vector3f Vec3Float;
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::VectorXi VecDynInt;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat;
typedef OpenMesh::DefaultTraits MyTraits;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  TriMesh;

namespace registration{

class ScaleShifter
{
    public:

        void set_input(const FeatureMat &inLowFeatures,
                       const VecDynInt &inLowOriginalIndices,
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
        std::vector<std::pair<int,int>> _correspondingIndexPairs;
        std::vector<int> _newIndices;

        //# Internal Parameters
        size_t _numLowNodes = 0;
        size_t _numHighNodes = 0;
        size_t _numCorrespondingNodes = 0;
        size_t _numNewNodes = 0;

        //# Internal functions
        void _find_corresponding_and_new_indices();
};

}//namespace registration

#endif // SCALESHIFTER_HPP
