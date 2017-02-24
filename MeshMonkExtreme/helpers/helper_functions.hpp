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
typedef Eigen::Matrix< int, Eigen::Dynamic, 3> FacesMat; //matrix Mx3 of type unsigned int
typedef Eigen::VectorXf VecDynFloat;
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






void convert_mesh_to_matrices(const TriMesh &inMesh,
                                FeatureMat &outFeatures,
                                FacesMat &outFaces);

void convert_mesh_to_matrices(const TriMesh &inMesh,
                                FeatureMat &outFeatures);

void convert_mesh_to_matrices(const TriMesh &inMesh,
                                FeatureMat &outFeatures,
                                FacesMat &outFaces,
                                VecDynFloat &outFlags);

void convert_matrices_to_mesh(const FeatureMat &inFeatures,
                                const FacesMat &inFaces,
                                TriMesh &outMesh);

void convert_matrices_to_mesh(const FeatureMat &inFeatures,
                                const FacesMat &inFaces,
                                const VecDynFloat &inFlags,
                                TriMesh &outMesh);


void load_obj_to_eigen(const std::string inObjFilename,
                                TriMesh &outMesh,
                                FeatureMat &outFeatureMatrix);


void write_eigen_to_obj(const FeatureMat &inFeatures,
                                TriMesh &inMesh,
                                const std::string inObjFilename);

bool import_data(const std::string inFloatingMeshPath,
                 const std::string inTargetMeshPath,
                 FeatureMat &outFloatingFeatures,
                 FeatureMat &outTargetFeatures,
                 FacesMat &outFloatingFaces);

bool export_data(FeatureMat &inResultFeatures,
                 FacesMat &inResultFaces,
                 const std::string inResultMeshPath);


void update_normals_for_altered_positions(TriMesh &ioMesh,
                                        FeatureMat &ioFeatures);

}//namespace registration
#endif // HELPER_FUNCTIONS_HPP_INCLUDED
