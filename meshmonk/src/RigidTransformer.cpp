#include "RigidTransformer.hpp"



namespace registration{


void RigidTransformer::set_input(const FeatureMat * const inCorrespondingFeatures, const VecDynFloat * const inWeights){
    _inCorrespondingFeatures = inCorrespondingFeatures;
    _inWeights = inWeights;
}
void RigidTransformer::set_output(FeatureMat * const ioFeatures){
    _ioFeatures = ioFeatures;
}
void RigidTransformer::set_parameters(bool scaling){
    _scaling = scaling;
}


void RigidTransformer::_update_transformation() {


    //# Info & Initialization
    _numElements = _ioFeatures->rows();
    _numFeatures = _ioFeatures->cols();
    MatDynFloat floatingPositions = MatDynFloat::Zero(3, _numElements);
    MatDynFloat correspondingPositions = MatDynFloat::Zero(3, _numElements);

    //## Tranpose the data if necessary
    if ((_numElements > _numFeatures) && (_numFeatures == NUM_FEATURES)) { //this should normally be the case
        floatingPositions = _ioFeatures->leftCols(3).transpose();
        correspondingPositions = _inCorrespondingFeatures->leftCols(3).transpose();
    }
    else {
        std::cerr << "Warning: input of rigid transformation expects rows to correspond with elements, not features, and to have more elements than features per element." << std::endl;
    }

    //# Compute the tranformation in 10 steps.
    //## 1. Get the centroids of each set
    Vec3Float floatingCentroid = Vec3Float::Zero();
    Vec3Float correspondingCentroid = Vec3Float::Zero();
    float sumWeights = 0.0;
    //### Weigh and sum all features
    for (size_t i = 0 ; i < _numElements ; i++) {
        floatingCentroid += (*_inWeights)[i] * floatingPositions.col(i).segment(0,3);
        correspondingCentroid += (*_inWeights)[i] * correspondingPositions.col(i).segment(0,3);
        sumWeights += (*_inWeights)[i];
    }
    //### Divide by total weight
    floatingCentroid /= sumWeights;
    correspondingCentroid /= sumWeights;

    //## 2. Compute the Cross Variance matrix
    Mat3Float crossVarianceMatrix = Mat3Float::Zero();
    for(size_t i = 0 ; i < _numElements ; i++) {
        crossVarianceMatrix += (*_inWeights)[i] * floatingPositions.col(i) * correspondingPositions.col(i).transpose();
    }
    crossVarianceMatrix = crossVarianceMatrix / sumWeights - floatingCentroid*correspondingCentroid.transpose();

    //## 3. Compute the Anti-Symmetric matrix
    Mat3Float antiSymmetricMatrix = crossVarianceMatrix - crossVarianceMatrix.transpose();

    //## 4. Use the cyclic elements of the Anti-Symmetric matrix to construct delta
    Vec3Float delta = Vec3Float::Zero();
    delta[0] = antiSymmetricMatrix(1,2);
    delta[1] = antiSymmetricMatrix(2,0);
    delta[2] = antiSymmetricMatrix(0,1);

    //## 5. Compute Q
    Mat4Float Q = Mat4Float::Zero();
    Q(0,0) = crossVarianceMatrix.trace();
    Q.block<3,1>(1,0) = delta;
    Q.block<1,3>(0,1) = delta.transpose();
    Q.block<3,3>(1,1) = crossVarianceMatrix + crossVarianceMatrix.transpose()
                        - crossVarianceMatrix.trace() * Mat3Float::Identity();

    //## 6. Now compute the rotation quaternion by finding the eigenvector of Q
    //## of its largest eigenvalue.
    Vec4Float rotQuat = Vec4Float::Zero();
    EigenVectorDecomposer decomposer(Q);
    if (decomposer.info() != Eigen::Success) {
        std::cerr << "eigenvector decomposer on Q failed!" << std::endl;
        std::cerr << "Q : " << Q << std::endl;
    }
    size_t indexMaxVal = 0;
    float maxEigenValue = 0.0;
    for ( size_t i = 0 ; i < 4 ; i++ ){
        if ( decomposer.eigenvalues()[i] > maxEigenValue ){
            maxEigenValue = decomposer.eigenvalues()[i];
            indexMaxVal = i;
        }
    }
    rotQuat = decomposer.eigenvectors().col( indexMaxVal );

    //## 7. Construct the rotation matrix
    Mat3Float rotMatTemp = Mat3Float::Zero();
    //### diagonal elements
    rotMatTemp(0,0) = rotQuat[0] * rotQuat[0] + rotQuat[1] * rotQuat[1] - rotQuat[2] * rotQuat[2] - rotQuat[3] * rotQuat[3];
    rotMatTemp(1,1) = rotQuat[0] * rotQuat[0] + rotQuat[2] * rotQuat[2] - rotQuat[1] * rotQuat[1] - rotQuat[3] * rotQuat[3];
    rotMatTemp(2,2) = rotQuat[0] * rotQuat[0] + rotQuat[3] * rotQuat[3] - rotQuat[1] * rotQuat[1] - rotQuat[2] * rotQuat[2];
    //### remaining elements
    rotMatTemp(0,1) = 2.0 * (rotQuat[1] * rotQuat[2] - rotQuat[0] * rotQuat[3]);
    rotMatTemp(1,0) = 2.0 * (rotQuat[1] * rotQuat[2] + rotQuat[0] * rotQuat[3]);
    rotMatTemp(0,2) = 2.0 * (rotQuat[1] * rotQuat[3] + rotQuat[0] * rotQuat[2]);
    rotMatTemp(2,0) = 2.0 * (rotQuat[1] * rotQuat[3] - rotQuat[0] * rotQuat[2]);
    rotMatTemp(1,2) = 2.0 * (rotQuat[2] * rotQuat[3] - rotQuat[0] * rotQuat[1]);
    rotMatTemp(2,1) = 2.0 * (rotQuat[2] * rotQuat[3] + rotQuat[0] * rotQuat[1]);

    //## 8. Estimate scale (if required)
    float scaleFactor = 1.0; //>1 to grow ; <1 to shrink
    if (_scaling == true){
        float numerator = 0.0;
        float denominator = 0.0;
        for (size_t i = 0 ; i < _numElements ; i++){
            //### Center and rotate the floating position
            Vec3Float newFloatingPos = rotMatTemp * (floatingPositions.block<3,1>(0,i) - floatingCentroid.segment(0, 3));
            //### Center the corresponding position
            Vec3Float newCorrespondingPos = correspondingPositions.block<3,1>(0,i) - correspondingCentroid.segment(0, 3);

            //### Increment numerator and denominator
            numerator += (*_inWeights)[i] * newCorrespondingPos.transpose() * newFloatingPos;
            denominator += (*_inWeights)[i] * newFloatingPos.transpose() * newFloatingPos;
        }
        scaleFactor = numerator / denominator;
    }


    //## 9. Compute the remaining translation necessary between the centroids
    Vec3Float translation = correspondingCentroid - scaleFactor * rotMatTemp * floatingCentroid;

    //## 10. Compute the entire transformation matrix.
    //### Initialize Matrices
    Mat4Float translationMatrix = Mat4Float::Identity();
    Mat4Float rotationMatrix = Mat4Float::Identity(); //to apply just the rotation to the surface normals later
    Mat4Float scaledRotationMatrix = Mat4Float::Identity();
    _transformationMatrix = Mat4Float::Identity();
    //### Convert to homogeneous transformation matrices
    translationMatrix.block<3, 1>(0,3) = translation;
    rotationMatrix.block<3, 3>(0, 0) = rotMatTemp;
    scaledRotationMatrix.block<3, 3>(0, 0) = scaleFactor * rotMatTemp;
    //### Matrix transformations on data is performed from right to left.
    //### Translation should be performed before rotating, so translationMatrix
    //### stands right in the multiplication with rotationMatrix.
    _transformationMatrix = scaledRotationMatrix * translationMatrix;

    //# Apply the transformation
    //## initialize a homogeneous vector in a [x y z 1] representation
    Vec4Float vector4d = Vec4Float::Ones();
    for (size_t i = 0 ; i < _numElements ; i++) {
        //## Transform the position
        //### Extract position from feature matrix
        vector4d.segment(0, 3) = floatingPositions.block<3,1>(0,i);
        //### Apply transformation to the position
        vector4d = _transformationMatrix * vector4d;
        _ioFeatures->block<1,3>(i,0) = vector4d.segment(0, 3);
        //## Rotate the vertex normals
        //### Extract normal from feature matrix
        vector4d.segment(0, 3) = _ioFeatures->block<1,3>(i,3);
        //### Apply rotation to the normal
        vector4d = rotationMatrix * vector4d;
        _ioFeatures->block<1,3>(i,3) = vector4d.segment(0, 3);
    }
}

void RigidTransformer::update() {
    _update_transformation();
}//end update

}//namespace registration
