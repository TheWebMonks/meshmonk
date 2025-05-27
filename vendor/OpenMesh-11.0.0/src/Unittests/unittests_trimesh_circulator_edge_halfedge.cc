#include "unittests_trimesh_circulator_edge.hh"

#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>

#include <vector>

namespace {

class OpenMeshTrimeshCirculatorEdgeHalfedge : public OpenMeshTrimeshCirculatorEdge {
public:
    // Check that the _visited halfedges comply with the expected iteration order on _eh.
    void CheckIteration(EH _eh, const std::vector<HEH>& _visited) {
        ASSERT_TRUE(mesh_.is_valid_handle(_eh));
        const auto eh = make_smart(_eh, mesh_);
        std::vector<HEH> expected;
        for (int i = 0; i < 2; ++i) {
            expected.push_back(eh.halfedge(i));
        }

        // Elements must have been visited in the expected order.
        ASSERT_EQ(_visited.size(), expected.size());
        for (size_t i = 0; i < _visited.size(); ++i) {
            EXPECT_EQ(_visited[i], expected[i]);
        }
    }

    // Check that the _visited halfedges comply with the expected iteration order on _heh.edge(), starting at _heh.
    void CheckIteration(HEH _heh, const std::vector<HEH>& _visited) {
        // Always exactly 2 elements.
        ASSERT_TRUE(mesh_.is_valid_handle(_heh));
        std::vector<HEH> expected;
        expected.push_back(_heh);
        expected.push_back(mesh_.opposite_halfedge_handle(_heh));

        // Elements must have been visited in the expected order.
        ASSERT_EQ(_visited.size(), expected.size());
        for (size_t i = 0; i < _visited.size(); ++i) {
            EXPECT_EQ(_visited[i], expected[i]);
        }
    }
};

// Mutable mesh.

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, Iter) {
    for (auto eh : mesh_.edges()) {
        std::vector<HEH> visited;
        for (auto it = mesh_.eh_iter(eh); it.is_valid(); ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, BeginEnd) {
    for (auto eh : mesh_.edges()) {
        std::vector<HEH> visited;
        for (auto it = mesh_.eh_begin(eh), end = mesh_.eh_end(eh); it != end; ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, Range) {
    for (auto eh : mesh_.edges()) {
        std::vector<HEH> visited;
        for (auto vh : mesh_.eh_range(eh)) {
            visited.push_back(vh);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, SmartHandleRange) {
    for (auto eh : mesh_.edges()) {
        std::vector<HEH> visited;
        auto smart_eh = make_smart(eh, mesh_);
        for (auto vh : smart_eh.halfedges()) {
            visited.push_back(vh);
        }
        CheckIteration(eh, visited);
    }
}

// const mesh.

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, ConstIter) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<HEH> visited;
        for (auto it = const_mesh.ceh_iter(eh); it.is_valid(); ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, ConstBeginEnd) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<HEH> visited;
        for (auto it = const_mesh.ceh_begin(eh), end = const_mesh.ceh_end(eh); it != end; ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, ConstRange) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<HEH> visited;
        for (auto vh : const_mesh.eh_range(eh)) {
            visited.push_back(vh);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, ConstSmartHandleRange) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<HEH> visited;
        auto smart_eh = make_smart(eh, const_mesh);
        for (auto vh : smart_eh.halfedges()) {
            visited.push_back(vh);
        }
        CheckIteration(eh, visited);
    }
}

// Mutable mesh, with given start halfedge.

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, InitializedRange) {
    for (auto start_heh : mesh_.halfedges()) {
        std::vector<HEH> visited;
        for (auto vh : mesh_.eh_range(start_heh)) {
            visited.push_back(vh);
        }
        CheckIteration(start_heh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, InitializedSmartHandleRange) {
    for (auto start_heh : mesh_.halfedges()) {
        std::vector<HEH> visited;
        auto smart_start_heh = make_smart(start_heh, mesh_);
        auto smart_start_eh = smart_start_heh.edge();
        for (auto vh : smart_start_eh.halfedges(smart_start_heh)) {
            visited.push_back(vh);
        }
        CheckIteration(start_heh, visited);
    }
}

// const mesh, with given start halfedge.

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, ConstInitializedRange) {
    const auto& const_mesh = mesh_;
    for (auto start_heh : const_mesh.halfedges()) {
        std::vector<HEH> visited;
        for (auto vh : mesh_.eh_range(start_heh)) {
            visited.push_back(vh);
        }
        CheckIteration(start_heh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeHalfedge, ConstInitializedSmartHandleRange) {
    const auto& const_mesh = mesh_;
    for (auto start_heh : const_mesh.halfedges()) {
        std::vector<HEH> visited;
        auto smart_start_heh = make_smart(start_heh, const_mesh);
        auto smart_start_eh = smart_start_heh.edge();
        for (auto vh : smart_start_eh.halfedges(smart_start_heh)) {
            visited.push_back(vh);
        }
        CheckIteration(start_heh, visited);
    }
}

}
