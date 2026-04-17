#include <iostream>
#include <string>
#include <vector>
#include <fstream> // Added for std::ofstream

// Argument parsing
#include "cxxopts.hpp"

// OpenMesh
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"

// Eigen
#include <Eigen/Dense>

// MeshMonk core library
// Assuming meshmonk.hpp and global.hpp are accessible via include paths
// set in the root CMakeLists.txt for meshmonk_shared
#include "meshmonk/meshmonk.hpp" // Relative to include path, e.g., project root
#include "meshmonk/global.hpp"   // Relative to include path

// Define a default mesh type for OpenMesh
struct CLIMeshTraits : public OpenMesh::DefaultTraits
{
  // Let OpenMesh store vertex normals
  VertexAttributes(OpenMesh::Attributes::Normal);
  // Let OpenMesh store face normals
  FaceAttributes(OpenMesh::Attributes::Normal);
};
typedef OpenMesh::TriMesh_ArrayKernelT<CLIMeshTraits> MyMesh;


// Placeholder for loading OBJ mesh
bool load_obj_mesh(const std::string& filename, 
                   Eigen::MatrixXd& V, 
                   Eigen::MatrixXi& F,
                   MyMesh& mesh) {
    // Read the mesh using OpenMesh
    if (!OpenMesh::IO::read_mesh(mesh, filename)) {
        std::cerr << "Error loading mesh from " << filename << std::endl;
        return false;
    }

    // Request and compute vertex normals if not already available
    mesh.request_vertex_normals();
    if (!mesh.has_vertex_normals()) {
        std::cerr << "Error: Could not request vertex normals for: " << filename << std::endl;
        // For MeshMonk, normals are often crucial. Let's consider this fatal for now.
        // If features are to be computed differently, this might change.
    }
    mesh.update_normals(); // Compute normals
    if (!mesh.has_vertex_normals()) {
        std::cerr << "Warning: Normals could not be computed for: " << filename << std::endl;
        // Proceeding without normals, FeatureMat will have zero normals.
    }

    // Populate V (vertices)
    V.resize(mesh.n_vertices(), 3);
    for (size_t i = 0; i < mesh.n_vertices(); ++i) {
        MyMesh::Point p = mesh.point(mesh.vertex_handle(i));
        V(i, 0) = p[0];
        V(i, 1) = p[1];
        V(i, 2) = p[2];
    }

    // Populate F (faces) - 0-indexed
    F.resize(mesh.n_faces(), 3);
    for (size_t i = 0; i < mesh.n_faces(); ++i) {
        MyMesh::FaceHandle fh = mesh.face_handle(i);
        int j = 0;
        for (MyMesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
            if (j < 3) { // Assuming triangular faces
                F(i, j) = fv_it->idx();
            }
            j++;
        }
        if (j != 3) {
            std::cerr << "Warning: Face " << i << " is not a triangle (or error in iteration)." << std::endl;
        }
    }
    return true;
}

// Placeholder for saving OBJ mesh
bool save_obj_mesh(const std::string& filename, 
                   const MyMesh& mesh_to_save) {
    // The V and F parameters are removed as mesh_to_save is the source of truth.
    std::cout << "Saving mesh to " << filename << std::endl;
    if (!OpenMesh::IO::write_mesh(mesh_to_save, filename)) {
        std::cerr << "Error saving mesh to " << filename << std::endl;
        return false;
    }
    std::cout << "Mesh saved successfully to " << filename << std::endl;
    return true;
}

// Placeholder for saving transformation matrix
bool save_transform_matrix(const std::string& filename, const Eigen::Matrix4f& transform) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file to save transform matrix: " << filename << std::endl;
        return false;
    }
    outfile << transform << std::endl;
    outfile.close();
    std::cout << "Transformation matrix saved to " << filename << std::endl;
    return true;
}


int main(int argc, char* argv[]) {
    cxxopts::Options options("meshmonk_cli", "Command-line interface for MeshMonk registration tasks.");

    options.add_options()
        ("h,help", "Print usage");

    options.add_options("Global")
        ("command", "The command to execute (pyramid_reg or rigid_reg)", cxxopts::value<std::string>())
        ("source", "Source mesh file (OBJ)", cxxopts::value<std::string>())
        ("target", "Target mesh file (OBJ)", cxxopts::value<std::string>())
        ("output", "Output mesh file (OBJ)", cxxopts::value<std::string>());
        
    options.parse_positional({"command", "source", "target", "output"});
    options.positional_help("<command> <source_mesh> <target_mesh> <output_mesh>");


    // Command-specific options
    // Pyramid Registration options
    options.add_options("pyramid_reg")
        ("num_iterations", "Number of iterations per level", cxxopts::value<size_t>()->default_value("90")) // MATLAB: 90
        ("smoothness", "Smoothness term weight (transformSigma)", cxxopts::value<float>()->default_value("3.0")) // MATLAB: 3.0
        ("num_pyramid_layers", "Number of pyramid layers", cxxopts::value<size_t>()->default_value("3")) // MATLAB: 3
        ("ds_float_start", "Downsample start % for float mesh", cxxopts::value<float>()->default_value("50.0")) // MATLAB: 50
        ("ds_target_start", "Downsample start % for target mesh", cxxopts::value<float>()->default_value("70.0")) // MATLAB: 70
        ("ds_float_end", "Downsample end % for float mesh", cxxopts::value<float>()->default_value("0.0")) // MATLAB: 0.0
        ("ds_target_end", "Downsample end % for target mesh", cxxopts::value<float>()->default_value("0.0")) // MATLAB: 0.0
        ("correspondences_symmetric", "Use symmetric correspondences", cxxopts::value<bool>()->default_value("true")) // MATLAB: true
        ("correspondences_num_neighbours", "Num neighbours for correspondences", cxxopts::value<size_t>()->default_value("3")) // MATLAB: 3
        ("correspondences_flag_threshold", "Flag threshold for correspondences", cxxopts::value<float>()->default_value("0.999")) // MATLAB: 0.999
        ("correspondences_equalize_push_pull", "Equalize push/pull for correspondences", cxxopts::value<bool>()->default_value("false")) // MATLAB: false
        ("inlier_kappa", "Kappa value for inlier detection", cxxopts::value<float>()->default_value("12.0")) // MATLAB: 12.0
        ("inlier_use_orientation", "Use orientation for inlier detection", cxxopts::value<bool>()->default_value("true")) // MATLAB: true
        ("transform_num_viscous_iterations_start", "Num viscous iterations (start)", cxxopts::value<size_t>()->default_value("90")) // MATLAB: 90 (numIterations)
        ("transform_num_viscous_iterations_end", "Num viscous iterations (end)", cxxopts::value<size_t>()->default_value("1")) // MATLAB: 1
        ("transform_num_elastic_iterations_start", "Num elastic iterations (start)", cxxopts::value<size_t>()->default_value("90")) // MATLAB: 90 (numIterations)
        ("transform_num_elastic_iterations_end", "Num elastic iterations (end)", cxxopts::value<size_t>()->default_value("1")); // MATLAB: 1

    // Rigid Registration options
    options.add_options("rigid_reg")
        ("transform_output", "Output file for the 4x4 transformation matrix (TXT)", cxxopts::value<std::string>()->default_value(""))
        ("rigid_num_iterations", "Number of iterations for rigid registration", cxxopts::value<size_t>()->default_value("80")) // MATLAB: 80
        ("rigid_correspondences_symmetric", "Use symmetric correspondences (rigid)", cxxopts::value<bool>()->default_value("true")) // MATLAB: true
        ("rigid_correspondences_num_neighbours", "Num neighbours for correspondences (rigid)", cxxopts::value<size_t>()->default_value("3")) // MATLAB: 3
        ("rigid_correspondences_flag_threshold", "Flag threshold for correspondences (rigid)", cxxopts::value<float>()->default_value("0.9")) // MATLAB: 0.9
        ("rigid_correspondences_equalize_push_pull", "Equalize push/pull for correspondences (rigid)", cxxopts::value<bool>()->default_value("false")) // MATLAB: false
        ("rigid_inlier_kappa", "Kappa value for inlier detection (rigid)", cxxopts::value<float>()->default_value("12.0")) // MATLAB: 12.0
        ("rigid_inlier_use_orientation", "Use orientation for inlier detection (rigid)", cxxopts::value<bool>()->default_value("true")) // MATLAB: true
        ("rigid_use_scaling", "Enable scaling in rigid registration", cxxopts::value<bool>()->default_value("false")); // MATLAB: false

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help({"", "Global", "pyramid_reg", "rigid_reg"}) << std::endl;
        return 0;
    }

    if (!result.count("command") || !result.count("source") || !result.count("target") || !result.count("output")) {
        std::cerr << "Error: Missing required arguments: command, source, target, output" << std::endl;
        std::cout << options.help({"", "Global"}) << std::endl;
        return 1;
    }

    std::string command = result["command"].as<std::string>();
    std::string source_file = result["source"].as<std::string>();
    std::string target_file = result["target"].as<std::string>();
    std::string output_file = result["output"].as<std::string>();

    std::cout << "Command: " << command << std::endl;
    std::cout << "Source: " << source_file << std::endl;
    std::cout << "Target: " << target_file << std::endl;
    std::cout << "Output: " << output_file << std::endl;

    // Eigen matrices for mesh data (double for OpenMesh compatibility, then cast to float for MeshMonk)
    Eigen::MatrixXd source_V_d, target_V_d;
    Eigen::MatrixXi source_F_i, target_F_i;
    
    // OpenMesh containers
    MyMesh source_mesh, target_mesh, result_mesh;

    // Load source and target meshes
    std::cout << "Loading source mesh: " << source_file << std::endl;
    if (!load_obj_mesh(source_file, source_V_d, source_F_i, source_mesh)) {
        std::cerr << "Failed to load source mesh." << std::endl;
        return 1;
    }
    std::cout << "Loading target mesh: " << target_file << std::endl;
    if (!load_obj_mesh(target_file, target_V_d, target_F_i, target_mesh)) {
        std::cerr << "Failed to load target mesh." << std::endl;
        return 1;
    }
    
    // Initialize result mesh with source data
    result_mesh = source_mesh;

    // Prepare data for MeshMonk
    // Faces (convert Eigen::MatrixXi to meshmonk::FacesMat)
    FacesMat source_faces_mat = source_F_i.cast<int>();
    FacesMat target_faces_mat = target_F_i.cast<int>();

    // Features (meshmonk::FeatureMat: N x 6 matrix of floats - pos_x, pos_y, pos_z, normal_x, normal_y, normal_z)
    FeatureMat source_features(source_mesh.n_vertices(), registration::NUM_FEATURES);
    for (size_t i = 0; i < source_mesh.n_vertices(); ++i) {
        MyMesh::Point p = source_mesh.point(source_mesh.vertex_handle(i));
        source_features(i, 0) = static_cast<float>(p[0]);
        source_features(i, 1) = static_cast<float>(p[1]);
        source_features(i, 2) = static_cast<float>(p[2]);
        if (source_mesh.has_vertex_normals()) {
            MyMesh::Normal n = source_mesh.normal(source_mesh.vertex_handle(i));
            source_features(i, 3) = static_cast<float>(n[0]);
            source_features(i, 4) = static_cast<float>(n[1]);
            source_features(i, 5) = static_cast<float>(n[2]);
        } else {
            source_features(i, 3) = 0.0f; source_features(i, 4) = 0.0f; source_features(i, 5) = 0.0f;
        }
    }

    FeatureMat target_features(target_mesh.n_vertices(), registration::NUM_FEATURES);
    for (size_t i = 0; i < target_mesh.n_vertices(); ++i) {
        MyMesh::Point p = target_mesh.point(target_mesh.vertex_handle(i));
        target_features(i, 0) = static_cast<float>(p[0]);
        target_features(i, 1) = static_cast<float>(p[1]);
        target_features(i, 2) = static_cast<float>(p[2]);
        if (target_mesh.has_vertex_normals()) {
            MyMesh::Normal n = target_mesh.normal(target_mesh.vertex_handle(i));
            target_features(i, 3) = static_cast<float>(n[0]);
            target_features(i, 4) = static_cast<float>(n[1]);
            target_features(i, 5) = static_cast<float>(n[2]);
        } else {
            target_features(i, 3) = 0.0f; target_features(i, 4) = 0.0f; target_features(i, 5) = 0.0f;
        }
    }

    // Flags (meshmonk::VecDynFloat: N x 1 vector of floats, typically initialized to 1.0)
    VecDynFloat source_flags = VecDynFloat::Ones(source_mesh.n_vertices());
    VecDynFloat target_flags = VecDynFloat::Ones(target_mesh.n_vertices());


    if (command == "pyramid_reg") {
        std::cout << "Executing pyramid registration..." << std::endl;
        
        meshmonk::pyramid_registration(
            source_features, target_features,
            source_faces_mat, target_faces_mat,
            source_flags, target_flags,
            result["num_iterations"].as<size_t>(),
            result["num_pyramid_layers"].as<size_t>(),
            result["ds_float_start"].as<float>(),
            result["ds_target_start"].as<float>(),
            result["ds_float_end"].as<float>(),
            result["ds_target_end"].as<float>(),
            result["correspondences_symmetric"].as<bool>(),
            result["correspondences_num_neighbours"].as<size_t>(),
            result["correspondences_flag_threshold"].as<float>(),
            result["correspondences_equalize_push_pull"].as<bool>(),
            result["inlier_kappa"].as<float>(),
            result["inlier_use_orientation"].as<bool>(),
            result["smoothness"].as<float>(), // transformSigma
            result["transform_num_viscous_iterations_start"].as<size_t>(),
            result["transform_num_viscous_iterations_end"].as<size_t>(),
            result["transform_num_elastic_iterations_start"].as<size_t>(),
            result["transform_num_elastic_iterations_end"].as<size_t>()
        );
        
        // Update result_mesh from the modified source_features (positions only)
        for (size_t i = 0; i < result_mesh.n_vertices(); ++i) {
            result_mesh.set_point(result_mesh.vertex_handle(i),
                                  MyMesh::Point(source_features(i, 0),
                                                source_features(i, 1),
                                                source_features(i, 2)));
        }
        result_mesh.request_vertex_normals(); // Request normals
        result_mesh.update_normals(); // Recompute normals after deformation

        std::cout << "Pyramid registration complete." << std::endl;

    } else if (command == "rigid_reg") {
        std::cout << "Executing rigid registration..." << std::endl;
        std::string transform_output_file = result["transform_output"].as<std::string>();
        Eigen::Matrix4f transform_matrix = Eigen::Matrix4f::Identity();

        meshmonk::rigid_registration(
            source_features, target_features,
            source_faces_mat, target_faces_mat,
            source_flags, target_flags,
            transform_matrix,
            result["rigid_num_iterations"].as<size_t>(),
            result["rigid_correspondences_symmetric"].as<bool>(),
            result["rigid_correspondences_num_neighbours"].as<size_t>(),
            result["rigid_correspondences_flag_threshold"].as<float>(),
            result["rigid_correspondences_equalize_push_pull"].as<bool>(),
            result["rigid_inlier_kappa"].as<float>(),
            result["rigid_inlier_use_orientation"].as<bool>(),
            result["rigid_use_scaling"].as<bool>()
        );
        
        // Update result_mesh from the modified source_features (positions only)
        for (size_t i = 0; i < result_mesh.n_vertices(); ++i) {
            result_mesh.set_point(result_mesh.vertex_handle(i),
                                  MyMesh::Point(source_features(i, 0),
                                                source_features(i, 1),
                                                source_features(i, 2)));
        }
        result_mesh.request_vertex_normals(); // Request normals
        result_mesh.update_normals(); // Recompute normals after transformation

        std::cout << "Rigid registration complete." << std::endl;

        if (!transform_output_file.empty()) {
            save_transform_matrix(transform_output_file, transform_matrix);
        }
    } else {
        std::cerr << "Unknown command: " << command << std::endl;
        std::cout << options.help({"", "Global", "pyramid_reg", "rigid_reg"}) << std::endl;
        return 1;
    }

    // Placeholder: Save the resulting mesh (transformed source)
    std::cout << "Saving output mesh..." << std::endl;
    if (!save_obj_mesh(output_file, result_mesh)) {
         std::cerr << "Failed to save output mesh to " << output_file << std::endl;
        return 1;
    }

    std::cout << "MeshMonk CLI execution finished." << std::endl;

    return 0;
}