#ifndef SCALESHIFTER_HPP
#define SCALESHIFTER_HPP


#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "../global.hpp"
#include "helper_functions.hpp"
#include "NeighbourFinder.hpp"

typedef Eigen::Vector3f Vec3Float;
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::VectorXi VecDynInt;
typedef Eigen::Matrix< float, Eigen::Dynamic, 3> Vec3Mat;
typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> IntegerMat;
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> FloatMat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat;
typedef OpenMesh::DefaultTraits MyTraits;
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  TriMesh;

namespace registration{


/*
ScaleShifter class.

This class is tailored to 'shift scales' between pyramid layers of the nonrigid registration. This means it has to
transfer the features of the mesh of the previous scale to the mesh of the current scale.

This is done by finding matches between the elements of each mesh, and interpolating the features for the elements
for which no matches are found.

Because both meshes are a subsampled version of the same original mesh, matches can be efficiently found through the
original indices of each element. Therefor, the ScaleShifter expects the original indices for each mesh as input:
-inLowOriginalIndices for the original indices of the (lower sampled) mesh of the previous scale.
-inHighOriginalIndicies for the original indices of the (higher sampled) mesh of the current scale.
*/

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
        std::vector<std::pair<int,int> > _matchingIndexPairs;
        std::vector<int> _newIndices;

        //# Internal Parameters
        size_t _numLowNodes = 0;
        size_t _numHighNodes = 0;
        size_t _numMatchingNodes = 0;
        size_t _numNewNodes = 0;

        //# Internal functions
        void _find_matching_and_new_indices();
        void _interpolate_new_nodes();
        void _copy_matching_nodes();
};

}//namespace registration

#endif // SCALESHIFTER_HPP
