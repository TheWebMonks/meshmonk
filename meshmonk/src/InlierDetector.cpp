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
    std::cout << "inlier corresponding flags: " << _inCorrespondingFlags->topRows(5) << std::endl;
    *_ioProbability = *_inCorrespondingFlags;
    std::cout << "_ioProbability from flags : " << _ioProbability->topRows(5) << std::endl;


    //# Distance based inlier/outlier classification
    //## Re-calculate the parameters sigma and lambda
    float sigmaa = 0.0;
    float lambdaa = 0.0;
    float sigmaNumerator = 0.0;
    float sigmaDenominator = 0.0;
    //DEBUG
    std::cout << "===========================" << std::endl;
    std::cout << "Computing SIGMA." << std::endl;
    //END DEBUG
    for (size_t i = 0 ; i < _numElements ; i++) {
        //### Compute distance (squared)
        FeatureVec difVector = _inCorrespondingFeatures->row(i) - _inFeatures->row(i);
        const float distanceSquared = difVector.squaredNorm();

        sigmaNumerator += (*_ioProbability)[i] * distanceSquared;
        sigmaDenominator += (*_ioProbability)[i];

        //DEBUG
        if (i < 2){
            std::cout << "difVector            : " << difVector << std::endl;
            std::cout << "distanceSquared      : " << distanceSquared << std::endl;
            std::cout << "sigma numerator add  : " << (*_ioProbability)[i] * distanceSquared << std::endl;
            std::cout << "sigma numerator      : " << sigmaNumerator << std::endl;
            std::cout << "sigma denominator add: " << (*_ioProbability)[i] << std::endl;
            std::cout << "sigma denominator    : " << sigmaDenominator << std::endl;
        }
        //END DEBUG
    }
    //### sigma and lambda
    sigmaa = std::sqrt(sigmaNumerator/sigmaDenominator);
    std::cout << "sigma 1: " << sigmaa << std::endl;
    if (sigmaa < _minimalSigma) {sigmaa = _minimalSigma;}
    else if (sigmaa > _maximalSigma) {sigmaa = _maximalSigma;}
    std::cout << "sigma 2: " << sigmaa << std::endl;

    lambdaa = 1.0/(std::sqrt(2.0 * 3.14159) * sigmaa) * std::exp(-0.5 * _kappa * _kappa);
    std::cout << "lambdaa: " << lambdaa << std::endl;

    //DEBUG
    std::cout << "-----------------------" << std::endl;
    std::cout << "Computing Weights." << std::endl;
    //END DEBUG
    //## Recalculate the distance-based probabilities
    for (size_t i = 0 ; i < _numElements ; i++) {
        //### Get squared distance
        FeatureVec difVector = _inCorrespondingFeatures->row(i) - _inFeatures->row(i);
        const float distanceSquared = difVector.squaredNorm();
        //### Compute probability
        float probability = 1.0/(std::sqrt(2.0 * 3.14159) * sigmaa) * std::exp(-0.5 * distanceSquared / std::pow(sigmaa, 2.0));
        probability /= (probability + lambdaa);
        (*_ioProbability)[i] *= probability;

        //DEBUG
        if (i < 2){
            std::cout << "difVector      : " << difVector << std::endl;
            std::cout << "distanceSquared: " << distanceSquared << std::endl;
            std::cout << "probability 1  : " << 1.0/(std::sqrt(2.0 * 3.14159) * sigmaa) * std::exp(-0.5 * distanceSquared / std::pow(sigmaa, 2.0)) << std::endl;
            std::cout << "probability 2  : " << probability << std::endl;
            std::cout << "probability 3  : " << (*_ioProbability)[i] << std::endl;
        }
        //END DEBUG
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

    std::cout << "final probability node 0: " << (*_ioProbability)[0] << std::endl;
    std::cout << "final probability node 1: " << (*_ioProbability)[1] << std::endl;

//    //DEBUG
//    std::cout<< "Inlier Weights: \n"
//    << "[0]: " << (*_ioProbability)[0] << "  || "
//    << "[1]: " << (*_ioProbability)[1] << "  || "
//    << "[2]: " << (*_ioProbability)[2] << "  || "
//    << "[3]: " << (*_ioProbability)[3] << "  || "
//    << "[4]: " << (*_ioProbability)[4] << "  || "
//    << "[5]: " << (*_ioProbability)[5] << "  || "
//    << "[6]: " << (*_ioProbability)[6] << "  || "
//    << "[7]: " << (*_ioProbability)[7] << "  || "
//    << "[8]: " << (*_ioProbability)[8] << "  || "
//    << "[9]: " << (*_ioProbability)[9] << std::endl;
//    //END DEBUG

}//end inlier_detection

}
