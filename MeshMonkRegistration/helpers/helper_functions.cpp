
#include <helper_functions.hpp>

namespace registration {

void normalize_sparse_matrix(SparseMat &ioMat) {
    /*
    # GOAL
    Normalize each row of the inputted sparse matrix.

    # INPUTS
    -ioMat

    # PARAMETERS

    # OUTPUT
    -ioMat
    */

    //# Normalize the rows of the affinity matrix
    //## Initialize a vector that will contain the sum of each row
    const size_t numRows = ioMat.rows();
    std::vector<float> sumRows(numRows, 0.0f);

    //## Loop over the rows and compute the total sum of its elements
    for (size_t i = 0 ; i < ioMat.outerSize() ; i++) {
        for (SparseMat::InnerIterator innerIt(ioMat,i) ; innerIt ; ++innerIt) {
            //### Get index of current row we're in
            const unsigned int currentRowIndex = innerIt.row();
            //### Add current element to the sum of this row.
            const float currentElement = innerIt.value();
            sumRows[currentRowIndex] += currentElement;
        }
    }

    //## Loop over the rows and divide each row by the sum of its elements
    for (size_t i = 0 ; i < ioMat.outerSize() ; i++) {
        for (SparseMat::InnerIterator innerIt(ioMat,i) ; innerIt ; ++innerIt) {
            //### Get index and sum of the current row we're in
            const unsigned int currentRowIndex = innerIt.row();
            const float currentRowSum = sumRows[currentRowIndex];

            //### Get current element and divide it by the sum of its row.
            const float currentElement = innerIt.value();
            innerIt.valueRef() = currentElement / currentRowSum;
        }
    }
}//end normalize_sparse_matrix()

}//namespace registration
