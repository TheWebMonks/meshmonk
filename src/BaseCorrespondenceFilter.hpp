#ifndef BASECORRESPONDENCEFILTER_HPP
#define BASECORRESPONDENCEFILTER_HPP


#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../global.hpp"
#include <iostream>

typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float
typedef Eigen::Matrix< float, 1, registration::NUM_FEATURES> FeatureVec; //matrix Mx6 of type float
typedef Eigen::Vector3f Vec3Float;
typedef Eigen::SparseMatrix<float, 0, int> SparseMat;
typedef Eigen::Triplet<float> Triplet;
typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> MatDynInt; //matrix MxN of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat;

namespace registration {

class BaseCorrespondenceFilter
{
    /*
    # GOAL
    This class serves as the base class for the correspondence filter classes.

    */

    public:
        BaseCorrespondenceFilter();
        virtual ~BaseCorrespondenceFilter();

        virtual void set_floating_input(const FeatureMat * const inFloatingFeatures,
                                const VecDynFloat * const inFloatingFlags){}
        virtual void set_target_input(const FeatureMat * const inTargetFeatures,
                            const VecDynFloat * const inTargetFlags){}
        void set_output(FeatureMat * const ioCorrespondingFeatures,
                        VecDynFloat * const ioCorrespondingFlags);
        SparseMat get_affinity() const {return _affinity;}
        virtual void set_parameters(const size_t numNeighbours,
                                    const float flagThreshold){}
        virtual void set_parameters(const size_t numNeighbours,
                                    const float flagThreshold,
                                    const bool equalizePushPull){}
        virtual void update(){}

    protected:

        //# Inputs
        const FeatureMat * _inFloatingFeatures = NULL;
        const VecDynFloat * _inFloatingFlags = NULL; //currently never used (only in the symmetric version)
        const FeatureMat * _inTargetFeatures = NULL;
        const VecDynFloat * _inTargetFlags = NULL;

        //# Outputs
        FeatureMat * _ioCorrespondingFeatures = NULL;
        VecDynFloat * _ioCorrespondingFlags = NULL;

        //# User Parameters
        size_t _numNeighbours = 3;
        float _flagThreshold = 0.99f;

        //# Internal Data structures
        SparseMat _affinity;

        //# Internal Parameters
        size_t _numFloatingElements = 0;
        size_t _numTargetElements = 0;

        //# Internal functions
        //## Function to update the sparse affinity matrix
        virtual void _update_affinity(){};
        //## Function to convert the sparse affinity weights into corresponding
        //## features and flags
        void _affinity_to_correspondences();

    private:


};

}//namespace registration
#endif // BASECORRESPONDENCEFILTER_HPP
