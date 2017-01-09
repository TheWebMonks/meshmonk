#ifndef HELPER_FUNCTIONS_HPP_INCLUDED
#define HELPER_FUNCTIONS_HPP_INCLUDED

#include <Eigen/SparseCore>

typedef Eigen::SparseMatrix<float, 0, int> SparseMat;

namespace registration {

void normalize_sparse_matrix(SparseMat &ioMat);

}//namespace registration
#endif // HELPER_FUNCTIONS_HPP_INCLUDED
