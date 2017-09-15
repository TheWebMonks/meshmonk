#include "InlierDetector.hpp"

namespace registration {

void InlierDetector::set_input(const FeatureMat * const inFeatures,
                        const FeatureMat * const inCorrespondingFeatures,
                        const VecDynFloat * const inCorrespondingFlags)
{
    //# Set input
    _inFeatures = inFeatures;
    _inCorrespondingFeatures = inCorrespondingFeatures;
    _inCorrespondingFlags = inCorrespondingFlags;

    //# Update internal variables
    _numElements = _inFeatures->rows();
}

void InlierDetector::set_output(VecDynFloat * const ioProbability)
{
    _ioProbability = ioProbability;
}

void InlierDetector::set_parameters(const float kappa, const bool useOrientation)
{
    _kappa = kappa;
    _useOrientation = useOrientation;
}


void InlierDetector::_determine_neighbours(){
    Vec3Mat floatingPositions = _inFeatures->leftCols(3);
    _neighbourFinder.set_source_points(&floatingPositions);
    _neighbourFinder.set_queried_points(&floatingPositions);
    _neighbourFinder.set_parameters(_numNeighbours);
    _neighbourFinder.update();
}//end _determine_neighbours()


void InlierDetector::_update_smoothing_weights(){
    /*
    The smoothing weights are the weights assigned to each vertex neighbour which
    will be used during the smoothing of the inlier weights.

    The weight is based on the inverse of the squared distance to each neighbour.
    */

    //# Initialize the weights matrix (we'll overwrite these values later)
    //# and get the neighbour indices
    _smoothingWeights = _neighbourFinder.get_distances();
    MatDynInt neighbourIndices = _neighbourFinder.get_indices();

    //# Loop over each neighbour and compute its smoothing weight
    //## 1) compute gaussian weights based on the distance to each neighbour
    bool printedWarning = false;
    for (size_t i = 0 ; i < _numElements ; i++){
        float sumWeight = 0.0f;
        for (size_t j = 0 ; j < _numNeighbours ; j++){
            //## Get the distance to the neighbour
            const float distanceSquared = _smoothingWeights(i,j); //smoothing weight still equals the squared distance here
            //## Compute the weight
            float weight = 1.0f/distanceSquared;
            //## Crop the weight to [eps, 1.0]
            if (weight > 1.0f) {weight=1.0f;}
            if (weight < _minWeight) {weight=_minWeight;}

            //## insert the weight into _smoothingWeights
            _smoothingWeights(i,j) = weight;
            sumWeight += weight;
        }
        //## normalize each row of weights
        if (sumWeight > 0.000001f){
            _smoothingWeights.row(i) /= sumWeight;
        }
        else if (!printedWarning) {
            std::cout << "Sum of smoothing weights in ViscoElastic Transformer should never be smaller than epsilon." << std::endl;
            printedWarning = true;
        }
    }
}


void InlierDetector::_smooth_inlier_weights(){
    //# Get the neighbour indices
    VecDynFloat tempInlierWeights;
    MatDynInt neighbourIndices = _neighbourFinder.get_indices();

    //## Start iterative loop
    for (size_t it = 0 ; it < _numSmoothingPasses ; it++){
        //## Copy the inlier weights into a temporary variable.
        tempInlierWeights = (*_ioProbability);

        //## Loop over the nodes
        for (size_t i = 0 ; i < _numElements ; i++) {
            //## For the current node, compute the weighted average of the neighbouring inlier weights
            float inlierWeightAvg = 0.0f;
            float neighbourInlierWeight = 0.0f;
            float sumSmoothingWeights = 0.0f;
            for (size_t j = 0 ; j < _numNeighbours ; j++) {
                // get neighbour index
                size_t neighbourIndex = neighbourIndices(i,j);
                // get neighbour smoothing weight and inlier weight
                const float smoothingWeight = _smoothingWeights(i,j); //inlier weight is already incorporated in the smoothing weight.
                neighbourInlierWeight = tempInlierWeights(neighbourIndex);
                sumSmoothingWeights += smoothingWeight;

                // increment the weighted average with current weighted neighbour vector
                inlierWeightAvg += smoothingWeight * neighbourInlierWeight;
            }

            //## Determine the weighted average inlier weight by dividing by the sum of smoothing weights.
            inlierWeightAvg /= sumSmoothingWeights; //smoothing weights are already normalized, so this should be redundant.
            (*_ioProbability)[i] = inlierWeightAvg;
        }

        //# Multiply the resulting inlier weights with the deterministic corresponding flags again!
        (*_ioProbability) *= (*_inCorrespondingFlags);
    }
}//end _smooth_inlier_weights()

void InlierDetector::update() {

//    _parameterList["jos"] = 2.0f;
//    _parameterList.insert(std::make_pair("ljkewqr", 3));
//    _parameterList.insert(std::make_pair("beschrijving", "Wat een moeilijkheden"));
//    _parameterList.insert(std::make_pair("inlierKappa", int(3)));
//    std::cout << "Parameter list 'beschrijving':" << _parameterList["beschrijving"] << std::endl;
//
//    DictionaryType::iterator it = _parameterList.begin();
//    while(it != _parameterList.end())
//    {
//        std::cout<<it->first<<" :: "<<it->second<<std::endl;
//        it++;
//    }

    //# Flag based inlier/outlier classification
    //## If a floating node is attracted to a corresponding node with flag 0.0,
    //## its inlier weight should be made zero as well. So we can use the
    //## corresponding flags as a first way to determine inlier weights for the
    //## floating nodes.
    //## -> Initialize the probabilities as a copy of the flags
    *_ioProbability = *_inCorrespondingFlags;

    //# Distance based inlier/outlier classification
    const float numDistanceBasedIterations = 10;
    for (size_t it = 0 ; it < numDistanceBasedIterations ; it++) {
        //## Re-calculate the parameters sigma and lambda
        float sigmaa = 0.0;
        float lambdaa = 0.0;
        float sigmaNumerator = 0.0;
        float sigmaDenominator = 0.0;
        for (size_t i = 0 ; i < _numElements ; i++) {
            //### Compute distance (squared)
            FeatureVec difVector = _inCorrespondingFeatures->row(i) - _inFeatures->row(i);
            const float distanceSquared = difVector.squaredNorm();

            sigmaNumerator += (*_ioProbability)[i] * distanceSquared;
            sigmaDenominator += (*_ioProbability)[i];
        }
        //### sigma and lambda
        sigmaa = std::sqrt(sigmaNumerator/sigmaDenominator);
        if (sigmaa < _minimalSigma) {sigmaa = _minimalSigma;}
        else if (sigmaa > _maximalSigma) {sigmaa = _maximalSigma;}

        lambdaa = 1.0/(std::sqrt(2.0 * 3.14159) * sigmaa) * std::exp(-0.5 * _kappa * _kappa);

        //## Recalculate the distance-based probabilities
        for (size_t i = 0 ; i < _numElements ; i++) {
            //### Get squared distance
            FeatureVec difVector = _inCorrespondingFeatures->row(i) - _inFeatures->row(i);
            const float distanceSquared = difVector.squaredNorm();
            //### Compute probability
            float probability = 1.0/(std::sqrt(2.0 * 3.14159) * sigmaa) * std::exp(-0.5 * distanceSquared / std::pow(sigmaa, 2.0));
            probability /= (probability + lambdaa);
            (*_ioProbability)[i] *= probability;

        }
    }


    //#Gradient Based inlier/outlier classification
    if (_useOrientation){
        float averageOrientationInlierWeight = 0.0f; //simply to warn the user when this is too low, they probably have the normals flipped.
        for (size_t i = 0 ; i < _numElements ; i++) {
            const Vec3Float normal = _inFeatures->row(i).tail(3);
            const Vec3Float correspondingNormal = _inCorrespondingFeatures->row(i).tail(3);
            //## Dot product gives an idea of how well they point in the same
            //## direction. This gives a weight between -1.0 and +1.0
            const float dotProduct = normal.dot(correspondingNormal);
            //## Rescale this result so that it's continuous between 0.0 and +1.0
            float probability = dotProduct / 2.0 + 0.5;
            (*_ioProbability)[i] *= probability;
            averageOrientationInlierWeight += probability;
        }

        averageOrientationInlierWeight /= _numElements;
        if (averageOrientationInlierWeight < 0.5f) {
            std::cout << "Warning: very low inlier weights due to surface normals. Are you sure one of the surfaces doesn't have its vertex normals flipped?" <<std::endl;
        }
    }

}//end inlier_detection

}
