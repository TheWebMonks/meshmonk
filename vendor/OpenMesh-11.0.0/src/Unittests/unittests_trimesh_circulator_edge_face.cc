#include "unittests_trimesh_circulator_edge.hh"

#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>

#include <vector>

namespace {

class OpenMeshTrimeshCirculatorEdgeFace : public OpenMeshTrimeshCirculatorEdge {
public:
    // Check that the _visited faces comply with the expected iteration order on _eh.
    void CheckIteration(EH _eh, const std::vector<FH>& _visited) {
        // Up to 2 elements, depending on whether incident faces exist.
        ASSERT_TRUE(mesh_.is_valid_handle(_eh));
        const auto eh = make_smart(_eh, mesh_);
        std::vector<FH> expected;
        for (int i = 0; i < 2; ++i) {
            const auto fh = eh.halfedge(i).face();
            if (fh.is_valid()) {
                expected.push_back(fh);
            }
        }

        // Elements must have been visited in the expected order.
        ASSERT_EQ(_visited.size(), expected.size());
        for (size_t i = 0; i < _visited.size(); ++i) {
            EXPECT_EQ(_visited[i], expected[i]);
        }
    }
};

// Mutable mesh.

TEST_F(OpenMeshTrimeshCirculatorEdgeFace, Iter) {
    for (auto eh : mesh_.edges()) {
        std::vector<FH> visited;
        for (auto it = mesh_.ef_iter(eh); it.is_valid(); ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeFace, BeginEnd) {
    for (auto eh : mesh_.edges()) {
        std::vector<FH> visited;
        for (auto it = mesh_.ef_begin(eh), end = mesh_.ef_end(eh); it != end; ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeFace, Range) {
    for (auto eh : mesh_.edges()) {
        std::vector<FH> visited;
        for (auto fh : mesh_.ef_range(eh)) {
            visited.push_back(fh);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeFace, SmartHandleRange) {
    for (auto eh : mesh_.edges()) {
        std::vector<FH> visited;
        auto smart_eh = make_smart(eh, mesh_);
        for (auto fh : smart_eh.faces()) {
            visited.push_back(fh);
        }
        CheckIteration(eh, visited);
    }
}

// const mesh.

TEST_F(OpenMeshTrimeshCirculatorEdgeFace, ConstIter) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<FH> visited;
        for (auto it = const_mesh.cef_iter(eh); it.is_valid(); ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeFace, ConstBeginEnd) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<FH> visited;
        for (auto it = const_mesh.cef_begin(eh), end = const_mesh.cef_end(eh); it != end; ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeFace, ConstRange) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<FH> visited;
        for (auto fh : const_mesh.ef_range(eh)) {
            visited.push_back(fh);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeFace, ConstSmartHandleRange) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<FH> visited;
        auto smart_eh = make_smart(eh, const_mesh);
        for (auto fh : smart_eh.faces()) {
            visited.push_back(fh);
        }
        CheckIteration(eh, visited);
    }
}

// Expected number of faces for interior, boundary, isolated edges.
TEST_F(OpenMeshTrimeshCirculatorEdgeFace, ExpectedNumber) {
    std::vector<std::pair<EH, int>> pairs = {
        { InteriorEdge, 2 },
        { BoundaryEdge, 1 },
        { IsolatedEdge, 0 },
    };
    for (const auto& pair : pairs) {
        const auto& eh = pair.first;
        const auto& expected_faces = pair.second;
        int visited_faces = 0;
        for (auto fh : mesh_.ef_range(eh)) {
            (void)fh; // Unused
            ++visited_faces;
        }
        EXPECT_EQ(visited_faces, expected_faces);
    }
}

}
