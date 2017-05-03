#include "ScaleShifter.hpp"

namespace registration {



void ScaleShifter::set_input(const FeatureMat &inLowFeatures,
                            const VecDynInt &inLowOriginalIndices,
                            const VecDynInt &inHighOriginalIndices){
    _inLowFeatures = &inLowFeatures;
    _inLowOriginalIndices = &inLowOriginalIndices;
    _inHighOriginalIndices = &inHighOriginalIndices;

    _numLowNodes = _inLowFeatures->rows();

    _numMatchingNodes = 0;
    _numNewNodes = 0;

    _matchingIndexPairs.clear();
    _newIndices.clear();
}//end set_input()


void ScaleShifter::set_output(FeatureMat &outHighFeatures){
    _outHighFeatures = &outHighFeatures;

    _numHighNodes = _outHighFeatures->rows();

    _numMatchingNodes = 0;
    _numNewNodes = 0;

    _matchingIndexPairs.clear();
    _newIndices.clear();
}//end set_output()


void ScaleShifter::_find_matching_and_new_indices(){
    //# Put index pairs of original and current indices into arrays
    //## Initialize the lists
    std::vector<std::pair<int,int>> lowIndexPairs;
    std::vector<std::pair<int,int>> highIndexPairs;
    //## Loop over the low sampled indices
    for (size_t i = 0 ; i < _numLowNodes ; i++){
        //## Get the original index
        int originalIndex = (*_inLowOriginalIndices)[i];
        //## Make a pair of the original index with its current index
        std::pair<int,int> indexPair(originalIndex,i);
        //## Save it into the list of index pairs.
        lowIndexPairs.push_back(indexPair);
    }

    //## Loop over the high sampled indices
    for (size_t i = 0 ; i < _numHighNodes ; i++){
        //## Get the original index
        int originalIndex = (*_inHighOriginalIndices)[i];
        //## Make a pair of the original index with its current index
        std::pair<int,int> indexPair(originalIndex,i);
        //## Save it into the list of index pairs.
        highIndexPairs.push_back(indexPair);
    }

    //# Sort the index pair lists by the original indices
    std::stable_sort(lowIndexPairs.begin(), lowIndexPairs.end(), [](auto &left, auto &right) {
        return left.first < right.first;
    }); //This construction makes sure we only sort using the first element of the std::pair
    std::stable_sort(highIndexPairs.begin(), highIndexPairs.end(), [](auto &left, auto &right) {
        return left.first < right.first;
    });


    //# Determine matching and new nodes of the high sampled mesh.
    //##    Determine the current index pairs between the matching nodes of the high and low
    //##    sampled mesh (using matches between the original indices). And determine which nodes
    //##    of the high sampled mesh have are new (have no matching original index in the low
    //##    sampled mesh).
    //## Loop over nodes of the high mesh
    size_t counterLow = 0;
    size_t counterHigh = 0;
    _matchingIndexPairs.clear();
    for ( ; counterHigh < _numHighNodes ; counterHigh++){
        //## Obtain original and current indices for high sampled mesh
        int highOriginalIndex = highIndexPairs[counterHigh].first;
        int highCurrentIndex = highIndexPairs[counterHigh].second;

        //## If the counter for the low sampled mesh becomes higher than its number
        //## of nodes, we ran out and all high samples mesh nodes are therefor new!
        if (counterLow >= _numLowNodes) {
            _newIndices.push_back(highCurrentIndex);
            continue;
        }
        else { //check if there is a match between the original indices of the high and low sampled node
            //## Obtain original and current indices for both meshes
            int lowOriginalIndex = lowIndexPairs[counterLow].first;
            int lowCurrentIndex = lowIndexPairs[counterLow].second;

            //## if the original indices match, we have found matching nodes!
            if (lowOriginalIndex == highOriginalIndex) {
                std::pair<int,int> matchingPair(highCurrentIndex, lowCurrentIndex);
                _matchingIndexPairs.push_back(matchingPair);
                counterLow++; //counterHigh is incremented in for-loop definition
            }
            //## if one index is smaller than the other, we have to increment that index
            else if (lowOriginalIndex > highOriginalIndex) {
                //## This node is new!
                _newIndices.push_back(highCurrentIndex);
                //## Move to the next pair of high indices
                continue; //skips the rest of the for-loop
            }
            else if (lowOriginalIndex < highOriginalIndex) {
                //## (!) This should never occur, it means that a node was found in the low sampled mesh that doesn't exist in the high sampled mesh.
                std::cerr << "the original indices in the low sampled mesh should be a subset of those of the high sampled mesh. Something went wrong?" << std::endl;
                counterLow++;
                continue;
            }
        }

    }

    //# Save the number of matching and new nodes
    _numMatchingNodes = _matchingIndexPairs.size();
    _numNewNodes = _newIndices.size();
    //## safety check
    if((_numMatchingNodes + _numNewNodes) != _numHighNodes){
        std::cerr << "Some nodes were missed as being new or matching nodes in ScaleShifter." << std::endl;
    }
}//end find_matching_and_new_indices()


void ScaleShifter::_interpolate_new_nodes(){
    /*
    The goal is here to interpolate the positions for the new nodes of
    the high sampled mesh. We'll do that by weighted k-nn. Therefor,
    we can reuse the CorrespondenceFilter class, which will gives us
    the 1/d_squared weighted average of k neighbouring positions.

    So we will set up a k-nn finder. The source points should be the nodes
    that match between the high and low sampled mesh, but with the feature
    values of the high sampled mesh (since these are the original features,
    unchanged by the registration process). The queried points are of course
    the new nodes of the high sampled mesh.

    After the nearest neighbours are found (using the old feature values of
    the high sampled mesh!), we then interpolate their feature values by
    weighted averaging of the features of the low sampled mesh.
    */

    //# Set up Source and Queried Points
    //## Construct matrix containing the features of the high sampled mesh, but
    //## only for those that have matching nodes in the low sampled mesh!
    FeatureMat matchingNodesOldFeatures = FeatureMat::Zero(_numMatchingNodes, NUM_FEATURES);
    for (size_t i = 0 ; i < _numMatchingNodes ; i++) {
        //## Get the index of the node of the high sampled mesh for which a
        //## matching node in the low sampled mesh was found. In other words,
        //## instead of using 'i' as the index, we are skipping all indices
        //## of the new nodes. We only want to copy the original features
        //## of those nodes of the high sampled mesh that are also present in
        //## the low sampled mesh.
        int matchingIndex = _matchingIndexPairs[i].first;
        matchingNodesOldFeatures.row(i) = _outHighFeatures->row(matchingIndex);
    }

    //# Set up Queried Points
    FeatureMat newNodesOldFeatures = FeatureMat::Zero(_numNewNodes, NUM_FEATURES);
    for (size_t i = 0 ; i < _numNewNodes ; i++) {
        //## Get the index of the new nodes of the high sampled mesh. We are skipping
        //## all the nodes for which a match was found in the low sampled mesh.
        int newIndex = _newIndices[i];
        newNodesOldFeatures.row(i) = _outHighFeatures->row(newIndex);
    }

    //# Set up a k-nn finder
    NeighbourFinder<FeatureMat> neighbourFinder;
    neighbourFinder.set_source_points(&matchingNodesOldFeatures);
    neighbourFinder.set_queried_points(&newNodesOldFeatures);
    size_t k = 3; //k = 3
    neighbourFinder.set_parameters(k);

    //# Get the indices and squared distances to the nearest neighbours
    neighbourFinder.update();
    const IntegerMat neighbourIndices = neighbourFinder.get_indices();
    const FloatMat neighbourSquaredDistances = neighbourFinder.get_distances();


    //# Set up the deformation field
    //## Initialization
    Vec3Mat deformationField = Vec3Mat::Zero(_numHighNodes, 3);
    Vec3Float deformation = Vec3Float::Zero();
    Vec3Float oldPosition = Vec3Float::Zero();
    Vec3Float newPosition = Vec3Float::Zero();
    //## Loop over the matching nodes
    for (size_t i = 0 ; i < _numMatchingNodes ; i++){
        //## Get the index pairs
        int highIndex = _matchingIndexPairs[i].first;
        int lowIndex = _matchingIndexPairs[i].second;
        //## Get the features and compute the differences
        oldPosition = (_outHighFeatures->row(highIndex)).head(3);
        newPosition = (_inLowFeatures->row(lowIndex)).head(3);
        deformation = newPosition - oldPosition;
        deformationField.row(highIndex) = deformation;
    }




    //# Compute the weighted average deformation for the new nodes
    //## Initialization
    Vec3Float newNodeOldNormal = Vec3Float::Zero();
    //## Loop over the indices of the new nodes
    for (size_t i = 0 ; i < _numNewNodes ; i++) {
        const int newNodeIndex = _newIndices[i];
        newNodeOldNormal = _outHighFeatures->row(newNodeIndex).tail(3);
        //## Initialize the deformation of the current new node to zero
        deformationField.row(newNodeIndex).setZero();
        float sumWeights = 0.0f;
        //## Loop over each found neighbour
        for ( size_t j = 0 ; j < k ; j++) {
            //### Get index of neighbour and squared distance to it
            const int neighbourIndex = neighbourIndices(i,j);
            float distanceSquared = neighbourSquaredDistances(i,j);

            //### The neighbourIndex now points to a row of matchingNodesOldFeatures.
            //### We need to find the matching index that points to the right row
            //### of the low sampled mesh.
            const int matchingHighNeighbourIndex = _matchingIndexPairs[neighbourIndex].first;

            //### For numerical stability, check if the distance is very small
            const float eps1 = 0.000001f;
            if (distanceSquared < eps1) {distanceSquared = eps1;}

            //### Compute the weight as 1/d_squared
            float weight = 1.0f / distanceSquared;

            //### Incorporate the orientation
            Vec3Float neighbourOldNormal = _outHighFeatures->row(matchingHighNeighbourIndex).tail(3);
            float dotProduct = newNodeOldNormal.dot(neighbourOldNormal);
            float orientationWeight = dotProduct / 2.0f + 0.5f;
            weight *= orientationWeight;

            //### Check for numerical stability (the weight will be
            //### normalized later, so dividing by a sum of tiny weights might
            //### go wrong.
            const float eps2 = 0.0001f;
            if (weight < eps2) {weight = eps2;}

            //### Weigh the matching features of the low sampled mesh
            deformationField.row(newNodeIndex) += weight * deformationField.row(matchingHighNeighbourIndex);
            sumWeights += weight;
        }
        //### Normalize
        deformationField.row(newNodeIndex) /= sumWeights;
    }


    //# Deform the new nodes
    for (size_t i = 0 ; i < _numNewNodes ; i++) {
        const int newNodeIndex = _newIndices[i];

        (_outHighFeatures->row(newNodeIndex)).head(3) += deformationField.row(newNodeIndex);
    }
}//end _interpolate_new_nodes()


void ScaleShifter::_copy_matching_nodes(){
    //# Copy the features of the lowly sampled mesh into the features
    //# of the matching nodes of the highly sampled mesh.
    for (size_t i = 0 ; i < _numMatchingNodes ; i++){
        //## Get the index pairs
        int highIndex = _matchingIndexPairs[i].first;
        int lowIndex = _matchingIndexPairs[i].second;
        //## copy the features of the low mesh into the features of the high mesh
        _outHighFeatures->row(highIndex) = _inLowFeatures->row(lowIndex);
    }
}//end _copy_matching_nodes()

void ScaleShifter::update(){

    //# Build the list of matching indices
    std::cout << "ScaleShifter: _find_matching_and_new_indices()" << std::endl;
    _find_matching_and_new_indices();



    //# Interpolate the features of new nodes
    std::cout << "ScaleShifter: _interpolate_new_nodes()" << std::endl;
    _interpolate_new_nodes();

    //# Copy the features of matching nodes.
    std::cout << "ScaleShifter: _copy_matching_nodes()" << std::endl;
    _copy_matching_nodes();

}//end update()


}//namespace registration
