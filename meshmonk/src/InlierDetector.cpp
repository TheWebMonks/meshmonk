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

void InlierDetector::set_parameters(const float kappa)
{
    _kappa = kappa;
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
    //## initialize the probabilities as a copy of the flags
    *_ioProbability = *_inCorrespondingFlags;


    //# Distance based inlier/outlier classification
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


    //#Gradient Based inlier/outlier classification
    for (size_t i = 0 ; i < _numElements ; i++) {
        const Vec3Float normal = _inFeatures->row(i).tail(3);
        const Vec3Float correspondingNormal = _inCorrespondingFeatures->row(i).tail(3);
        //## Dot product gives an idea of how well they point in the same
        //## direction. This gives a weight between -1.0 and +1.0
        const float dotProduct = normal.dot(correspondingNormal);
        //## Rescale this result so that it's continuous between 0.0 and +1.0
        float probability = dotProduct / 2.0 + 0.5;
        (*_ioProbability)[i] *= probability;
    }

}//end inlier_detection

}
