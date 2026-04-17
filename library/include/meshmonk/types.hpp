#ifndef MESHMONK_TYPES_HPP
#define MESHMONK_TYPES_HPP

#include <Eigen/Dense>

namespace meshmonk {

using FeatureMat  = Eigen::Matrix<float, Eigen::Dynamic, 6>;
using FacesMat    = Eigen::Matrix<int,   Eigen::Dynamic, 3>;   // signed int32 for numpy/trimesh interop
using VecDynFloat = Eigen::VectorXf;
using Vec3Mat     = Eigen::Matrix<float, Eigen::Dynamic, 3>;   // put here, not in result.hpp, to avoid circular includes

inline constexpr int NUM_FEATURES = 6;

} // namespace meshmonk

#endif // MESHMONK_TYPES_HPP
