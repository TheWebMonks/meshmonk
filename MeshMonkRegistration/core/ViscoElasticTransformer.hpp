#ifndef VISCOELASTICTRANSFORMER_HPP
#define VISCOELASTICTRANSFORMER_HPP

#include <Eigen/Dense>
#include <NeighbourFinder.hpp>

typedef Eigen::Vector3f Vec3Float;
typedef Eigen::VectorXf VecDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, 3> Vec3Mat; //matrix Mx3 of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat;
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float

namespace registration{

class ViscoElasticTransformer
{
    public:

        void set_input(const FeatureMat * const inCorrespondingPositions, const VecDynFloat * const inWeights);
        void set_output(FeatureMat * const ioFloatingPositions);
        void set_parameters(size_t numNeighbours = 10, float sigma = 3.0,
                            size_t viscousIterations = 10, size_t elasticIterations = 10);
        Vec3Mat get_transformation() const {return _displacementField;}
        void update();

    protected:

    private:
        //# Inputs
        //##_ioFeatures is used as both an input (to compute the transformation) and output
        const FeatureMat * _inCorrespondingPositions = NULL;
        const VecDynFloat * _inWeights = NULL;

        //# Outputs
        FeatureMat * _ioFloatingPositions = NULL;


        //# User Parameters
        size_t _numNeighbours = 10;
        float _sigma = 3.0;
        size_t _viscousIterations = 0;
        size_t _elasticIterations = 0;

        //# Internal Data structures
        Vec3Mat _displacementField;
        Vec3Mat _oldDisplacementField;
        NeighbourFinder<FeatureMat> _neighbourFinder;
        MatDynFloat _smoothingWeights;

        //# Internal Parameters
        size_t _numElements = 0;
        bool _neighboursOutdated = true;
        bool _weightsOutdated = true;

        //# Internal functions
        //## Update the neighbour finder
        void _update_neighbours();
        //## Update the weights used for smoothing
        void _update_smoothing_weights();
        //## Update the displacement field in a viscous manner
        void _update_viscously();
        //## Update the displacement field in an elastic manner
        void _update_elastically();
        //## Function to update the transformation
        void _update_transformation();
        //## Function to apply the transformation
        void _apply_transformation();
};

}//namespace registration

#endif // VISCOELASTICTRANSFORMER_HPP
