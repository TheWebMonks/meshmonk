#include <Unittests/unittests_common.hh>

#include <gtest/gtest.h>

#include <vector>

namespace {

class OpenMeshTrimeshNavigation : public OpenMeshBase {
public:
    using VH = OpenMesh::VertexHandle;
    using FH = OpenMesh::FaceHandle;
    using EH = OpenMesh::EdgeHandle;
    using HEH = OpenMesh::HalfedgeHandle;

    // This function is called before each test is run.
    void SetUp() override {
        std::vector<VH> vh;

        //   3----2
        //   |   /|
        //   |  / |
        //   | /  |
        //   |/   |
        //   0----1

        vh.push_back(mesh_.add_vertex({0.0, 0.0, 0.0})); // vh[0]
        vh.push_back(mesh_.add_vertex({1.0, 0.0, 0.0})); // vh[1]
        vh.push_back(mesh_.add_vertex({1.0, 1.0, 0.0})); // vh[2]
        vh.push_back(mesh_.add_vertex({0.0, 1.0, 0.0})); // vh[3]
        mesh_.add_face(vh[0], vh[1], vh[2]);
        mesh_.add_face(vh[0], vh[2], vh[3]);
    }
};

TEST_F(OpenMeshTrimeshNavigation, EdgeHalfedgeDefault) {
    for (EH eh : mesh_.edges()) {
        EXPECT_EQ(mesh_.halfedge_handle(eh), mesh_.halfedge_handle(eh, 0));
    }
}

}
