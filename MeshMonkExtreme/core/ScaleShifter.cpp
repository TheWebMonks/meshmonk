#include "ScaleShifter.hpp"

namespace registration {



void ScaleShifter::set_input(const FeatureMat &inLowFeatures,
                       const VecDynInt &inLowOriginalIndices,
                       const VecDynInt &inHighOriginalIndices){
    _inLowFeatures = &inLowFeatures;
    _inLowOriginalIndices = &inLowOriginalIndices;
    _inHighOriginalIndices = &inHighOriginalIndices;

    _numLowNodes = _inLowFeatures->rows();
    _numHighNodes = _outHighFeatures->rows();

    _numCorrespondingNodes = 0;
    _numNewNodes = 0;

    _correspondingIndexPairs.clear();
    _newIndices.clear();
}//end set_input()


void ScaleShifter::set_output(FeatureMat &outHighFeatures){
    _outHighFeatures = &outHighFeatures;
}//end set_output()


void ScaleShifter::_find_corresponding_and_new_indices(){
    //# Put index pairs of original and current indices into arrays
    //## Initialize the lists
    std::vector<std::pair<int,int>> lowIndexPairs;
    std::vector<std::pair<int,int>> highIndexPairs;
    //## Loop over the low sampled indices
    for (int i = 0 ; i < _numLowNodes ; i++){
        //## Get the original index
        int originalIndex = (*_inLowOriginalIndices)[i];
        //## Make a pair of the original index with its current index
        std::pair<int,int> indexPair(originalIndex,i);
        //## Save it into the list of index pairs.
        lowIndexPairs.push_back(indexPair);
    }
    //## Loop over the high sampled indices
    for (int i = 0 ; i < _numHighNodes ; i++){
        //## Get the original index
        int originalIndex = (*_inHighOriginalIndices)[i];
        //## Make a pair of the original index with its current index
        std::pair<int,int> indexPair(originalIndex,i);
        //## Save it into the list of index pairs.
        highIndexPairs.push_back(indexPair);
    }


    //# Sort the index pair lists by the original indices
//    struct compare_first_only {
//        bool operator()(const std::pair<int, int>& p1, const std::pair<int, int>& p2) {
//            return p1.first < p2.first;
//        }
//    };
//    std::stable_sort(lowIndexPairs.begin(), lowIndexPairs.end(), compare_first_only());
    std::stable_sort(lowIndexPairs.begin(), lowIndexPairs.end(), [](auto &left, auto &right) {
        return left.first < right.first;
    }); //This construction makes sure we only sort using the first element of the std::pair


    //# Determine corresponding and new nodes of the high sampled mesh.
    //##    Determine the current index pairs between the corresponding nodes of the high and low
    //##    sampled mesh (using matches between the original indices). And determine which nodes
    //##    of the high sampled mesh have are new (have no corresponding original index in the low
    //##    sampled mesh).
    //## Loop over nodes of the high mesh
    int counterLow = 0;
    int counterHigh = 0;
    _correspondingIndexPairs.clear();
    for ( ; counterHigh < _numHighNodes ; counterHigh++){
        //## Obtain original and current indices for both meshes
        int lowOriginalIndex = lowIndexPairs[counterLow].first;
        int lowCurrentIndex = lowIndexPairs[counterLow].second;
        int highOriginalIndex = highIndexPairs[counterHigh].first;
        int highCurrentIndex = highIndexPairs[counterHigh].second;

        //## if the original indices match, we have found corresponding nodes!
        if (lowOriginalIndex == highOriginalIndex) {
            std::pair<int,int> correspondingPair(highCurrentIndex, lowCurrentIndex);
            _correspondingIndexPairs.push_back(correspondingPair);
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
            return;
        }
    }

    //# Save the number of corresponding and new nodes
    _numCorrespondingNodes = _correspondingIndexPairs.size();
    _numNewNodes = _newIndices.size();
    //## safety check
    if((_numCorrespondingNodes + _numNewNodes) != _numHighNodes){
        std::cerr << "Some nodes were missed as being new or corresponding nodes in ScaleShifter." << std::endl;
    }
}//end _find_corresponding_and_new_indices()


void ScaleShifter::update(){

    //# Build the list of corresponding indices
    _find_corresponding_and_new_indices();

    //# Copy the features of corresponding nodes
    //## Loop over the corresponding nodes
    for (int i = 0 ; i < _numCorrespondingNodes ; i++){
        //## Get the index pairs
        int highIndex = _correspondingIndexPairs[i].first;
        int lowIndex = _correspondingIndexPairs[i].second;
        //## copy the features of the low mesh into the features of the high mesh
        _outHighFeatures->row(highIndex) = _inLowFeatures->row(lowIndex);
    }

    //# Interpolate the features of new nodes




}//end update()


}//namespace registration
