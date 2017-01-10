#ifndef VISCOELASTICTRANSFORMER_HPP
#define VISCOELASTICTRANSFORMER_HPP

#include <Eigen/Dense>
#include <NeighbourFinder.hpp>

typedef Eigen::Matrix< float, Eigen::Dynamic, 3> Vec3Mat; //matrix Mx3 of type float

namespace registration{

class ViscoElasticTransformer
{
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
        size_t _numNeighbours = 10;
        size_t _sigma = 3.0;
        size_t _viscousIterations = 0;
        size_t _elasticIterations = 0;

        //# Internal Data structures
        Vec3Mat _displacementField;
        Vec3Mat _oldDisplacementField;
        NeighbourFinder _neighbourFinder;

        //# Internal Parameters
        size_t _numElements = 0;

        //# Internal functions
        //## Function to update the transformation (displacement field)
        void _update_transformation();
        //## Function to apply the transformation
        void _apply_transformation();
        //## Function to smooth a displacement field (iterative, gaussian smoothing)
        void _smooth_vector_field();
};

}//namespace registration

#endif // VISCOELASTICTRANSFORMER_HPP
