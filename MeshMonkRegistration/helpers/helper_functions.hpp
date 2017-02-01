#ifndef HELPER_FUNCTIONS_HPP_INCLUDED
#define HELPER_FUNCTIONS_HPP_INCLUDED

#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <nanoflann.hpp>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../global.hpp"

typedef OpenMesh::TriMesh_ArrayKernelT<>  TriMesh;
typedef Eigen::SparseMatrix<float, 0, int> SparseMat;
typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic> MatDynInt; //matrix MxN of type unsigned int
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> MatDyn3Int; //matrix Mx3 of type unsigned int
typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic> MatDynFloat; //matrix MxN of type float
typedef Eigen::Matrix< float, Eigen::Dynamic, registration::NUM_FEATURES> FeatureMat; //matrix Mx6 of type float

namespace registration {

void fuse_affinities(SparseMat &ioAffinity1,
                    const SparseMat &inAffinity2);

void normalize_sparse_matrix(SparseMat &ioMat);

template <typename VecMatType>
void radius_nearest_neighbours(const VecMatType &inQueriedPoints,
                                const VecMatType &inSourcePoints,
                                MatDynInt &outNeighbourIndices,
                                MatDynFloat &outNeighbourSquaredDistances,
                                const float paramRadius = 3.0,
                                const size_t paramLeafsize = 15);






void convert_openmesh_to_eigen(const TriMesh &inMesh,
                                FeatureMat &outFeatures,
                                MatDyn3Int &outFaces);

void openmesh_to_eigen_features(const TriMesh &inputMesh,
                                FeatureMat &outputFeatures);


void load_obj_to_eigen_features(const std::string inObjFilename,
                                TriMesh &outMesh,
                                FeatureMat &outFeatureMatrix);


void eigen_features_to_openmesh(const FeatureMat &inFeatures,
                                TriMesh &outMesh);


void write_eigen_features_to_obj(const FeatureMat &inFeatures,
                                TriMesh &inMesh,
                                const std::string inObjFilename);


void update_normals_for_altered_positions(TriMesh &ioMesh,
                                        FeatureMat &ioFeatures);

}//namespace registration
#endif // HELPER_FUNCTIONS_HPP_INCLUDED
