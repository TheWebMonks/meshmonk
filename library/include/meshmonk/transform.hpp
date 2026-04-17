#ifndef MESHMONK_TRANSFORM_HPP
#define MESHMONK_TRANSFORM_HPP

#include "types.hpp"
#include <Eigen/Dense>

namespace meshmonk {

// RigidTransform wraps a 4x4 homogeneous matrix representing a rigid (or similarity) transform.
// Semantic contract: apply(original_features) == aligned_features (up to float rounding)
struct RigidTransform {
    Eigen::Matrix4f matrix = Eigen::Matrix4f::Identity();

    // Returns (*this) * rhs — matrix multiplication order (apply rhs first, then *this)
    RigidTransform compose(const RigidTransform& rhs) const {
        return RigidTransform{matrix * rhs.matrix};
    }

    // Apply transform to feature matrix:
    //   positions (cols 0-2): R*p + t
    //   normals   (cols 3-5): R*n  (no translation for normals)
    FeatureMat apply(const FeatureMat& features) const {
        const int N = static_cast<int>(features.rows());
        Eigen::Matrix3f R = matrix.topLeftCorner<3, 3>();
        Eigen::Vector3f t = matrix.topRightCorner<3, 1>();

        FeatureMat result(N, 6);
        // Transform positions
        result.leftCols(3) = (features.leftCols(3) * R.transpose()).rowwise() + t.transpose();
        // Transform normals (rotation only, no translation)
        result.rightCols(3) = features.rightCols(3) * R.transpose();
        return result;
    }

    // Inverse of [sR|t] is [(1/s)R^T | -(1/s)R^T * t]
    // For pure rotation (s=1), this is simply [R^T | -R^T*t]
    RigidTransform inverse() const {
        Eigen::Matrix3f R = matrix.topLeftCorner<3, 3>();
        Eigen::Vector3f t = matrix.topRightCorner<3, 1>();
        // Compute scale as norm of first column of R (for similarity transforms)
        float s2 = R.col(0).squaredNorm();  // s^2
        Eigen::Matrix3f Rt = R.transpose();
        Eigen::Vector3f t_inv = -(Rt * t) / s2;

        Eigen::Matrix4f inv = Eigen::Matrix4f::Identity();
        inv.topLeftCorner<3, 3>() = Rt / s2;
        inv.topRightCorner<3, 1>() = t_inv;
        return RigidTransform{inv};
    }
};

} // namespace meshmonk

#endif // MESHMONK_TRANSFORM_HPP
