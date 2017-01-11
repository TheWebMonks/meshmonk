#include "ViscoElasticTransformer.hpp"

namespace registration {



void ViscoElasticTransformer::set_input(const FeatureMat * const inCorrespondingPositions,
                                        const VecDynFloat * const inWeights){
    _inCorrespondingPositions = inCorrespondingPositions;
    _inWeights = inWeights;
    _weightsOutdated = true; //if the user sets new _inWeights, we need to update our smoothing weights.
}//end set_input()


void ViscoElasticTransformer::set_output(FeatureMat * const ioFloatingPositions){
    _ioFloatingPositions = ioFloatingPositions;
    _neighboursOutdated = true; //if the user sets new floating positions, we need to update our neighbours, and hence our smoothing weights.
    _weightsOutdated = true;

}//end set_output()

void ViscoElasticTransformer::set_parameters(size_t numNeighbours, float sigma,
                                            size_t viscousIterations,
                                            size_t elasticIterations)
{
    if (_numNeighbours != numNeighbours) {
        _neighboursOutdated = true; //if number of requested neighbours changes, we need to update the neighbours and weights
        _weightsOutdated = true;
    }
    if (std::abs(_sigma - sigma) > 0.0001 * _sigma) {
        _weightsOutdated = true; //if sigma changes, we need to update the weights
    }
    _numNeighbours = numNeighbours;
    _sigma = sigma;
    _viscousIterations = viscousIterations;
    _elasticIterations = elasticIterations;
}


//## Update the neighbour finder
void ViscoElasticTransformer::_update_neighbours(){
    _neighbourFinder.set_source_points(_ioFloatingPositions);
    _neighbourFinder.set_queried_points(_ioFloatingPositions);
    _neighbourFinder.set_parameters(_numNeighbours);
    _neighbourFinder.update();
}//end _update_neighbours()


//## Update the weights used for smoothing
void ViscoElasticTransformer::_update_smoothing_weights(){
    /*
    The smoothing weights are the weights assigned to each vertex neighbour which
    will be used during the smoothing of the vector fields.

    The weight is a combination of the user inputted weights (_inWeights) and a
    gaussian weight based on the distance to each neighbour.

    Therefor, we initialize the smoothing weights as the squared distances to each neighbour.
    Next, we convert that distance to the gaussian weight.
    Then, we multiply that weight with the user inputted weights.
    And lastly, we normalize each row so that the sum of weights for each
    node's neighbours equals 1.0.
    */

    //# Initialize the smoothing weights as the squared distances to the neighbouring nodes.
    _smoothingWeights = _neighbourFinder.get_distances();

    //# Loop over each neighbour and compute its smoothing weight
    //## 1) compute gaussian weights based on the distance to each neighbour
    for (size_t i = 0 ; i < _numElements ; i++){
        float sumWeight = 0.0f;
        for (size_t j = 0 ; j < _numNeighbours ; j++){
            //## Get the distance to the neighbour
            float distanceSquared = _smoothingWeights(i,j); //smoothing weight still equals the squared distance here
            //## Compute the gaussian weight
            float gaussianWeight = std::exp(-0.5f * distanceSquared / std::pow(_sigma, 2.0f));
            //## Combine the gaussian weight with the user defined weight
            float combinedWeight = (*_inWeights)[i] * gaussianWeight;

            //## insert the combined weight into _smoothingWeights
            _smoothingWeights(i,j) = combinedWeight;
            sumWeight += combinedWeight;
        }
        //## normalize each row of weights
        _smoothingWeights.row(i) /= sumWeight;
    }
}//end _update_smoothing_weights()



void ViscoElasticTransformer::_update_viscously(){
    /*
    Viscosity is obtained by incrementing the displacement field with a regularized force field.

    The force field is the difference between the current floating and corresponding positions.
    A purely viscous transformation can be obtained by adding a regularized force fold to the
    total displacement field. If, however, that total displacement field is regularized as well,
    a more elastic behaviour is achieved.

    So here, we will:
    1) Determine the force field
    2) Regularize the force field
    3) Add it to the total displacement field
    */

    //# 1) Determine the force field (difference between current floating and corresponding
    //# positions).
    Vec3Mat forceField = _inCorrespondingPositions->leftCols(3) - _ioFloatingPositions->leftCols(3);

    //# 2) Regularize the force field through iterative weighted averaging.
    /*
    For each force field vector, we compute a weighted average of all its
    neighbouring vectors. So we'll need the indices of each neighbour (to
    retrieve the vector) and the smoothing weight (which is precomputed
    in _update_smoothing_weights().

    This smoothing is done iteratively for a number of iterations chosen by
    the user (_viscousIterations).
    */
    //## Initialize the regularized force field and get the neighbour indices
    Vec3Mat regularizedForceField = forceField;
    MatDynInt neighbourIndices = _neighbourFinder.get_indices();

    //## Start iterative loop
    for (size_t it = 0 ; it < _viscousIterations ; it++){
        for (size_t i = 0 ; i < _numElements ; i++) {
            //## For the current displacement, compute the weighted average of the neighbouring
            //## vectors.
            Vec3Float vectorAverage = Vec3Float::Zero();
            Vec3Float neighbourVector;
            for (size_t j = 0 ; j < _numNeighbours ; j++) {
                // get neighbour index
                size_t neighbourIndex = neighbourIndices(i,j);
                // get neighbour weight and vector
                float weight = _smoothingWeights(i,j);
                neighbourVector = forceField.row(neighbourIndex);

                // increment the weighted average with current weighted neighbour vector
                vectorAverage += weight * neighbourVector;
            }

            regularizedForceField.row(i) = vectorAverage;
        }
        forceField = regularizedForceField;
    }


    //# Elastic Part
    //#3) Add the regulated Force Field to the current Displacement Field
    _oldDisplacementField = _displacementField; //save the previous displcament field before overwriting it.
    _displacementField += regularizedForceField;

}



void ViscoElasticTransformer::_update_elastically(){

    //# Get the neighbour indices
    Vec3Mat unregulatedDisplacementField;
    MatDynInt neighbourIndices = _neighbourFinder.get_indices();

    //## Start iterative loop
    for (size_t it = 0 ; it < _elasticIterations ; it++){
        //## Copy the displacement field into a temporary variable.
        unregulatedDisplacementField = _displacementField;

        //## Loop over each unregularized displacement vector and smooth it.
        for (size_t i = 0 ; i < _numElements ; i++) {
            //## For the current displacement, compute the weighted average of the neighbouring
            //## vectors.
            Vec3Float vectorAverage = Vec3Float::Zero();
            Vec3Float neighbourVector;
            for (size_t j = 0 ; j < _numNeighbours ; j++) {
                // get neighbour index
                size_t neighbourIndex = neighbourIndices(i,j);
                // get neighbour weight and vector
                float weight = _smoothingWeights(i,j);
                neighbourVector = unregulatedDisplacementField.row(neighbourIndex);

                // increment the weighted average with current weighted neighbour vector
                vectorAverage += weight * neighbourVector;
            }

            _displacementField.row(i) = vectorAverage;
        }
    }
}
//## Function to update the transformation
void ViscoElasticTransformer::_update_transformation(){
    _update_viscously();
    _update_elastically();
}



void ViscoElasticTransformer::_apply_transformation(){
    //# Displace each current floating position by the difference between
    //# the old and new displacement fields.
    for (size_t i = 0 ; i < _numElements ; i++) {
        _ioFloatingPositions->row(i).head(3) += (_displacementField.row(i) - _oldDisplacementField.row(i));
    }
}



void ViscoElasticTransformer::update(){
    if (_neighboursOutdated == true) {
        _update_neighbours();
        _weightsOutdated = true;
    }
    if (_weightsOutdated == true) {
        _update_smoothing_weights();
    }

    //# update the transformation
    _update_transformation();
    //# apply the transformation
    _apply_transformation();
}


}//namespace registration
