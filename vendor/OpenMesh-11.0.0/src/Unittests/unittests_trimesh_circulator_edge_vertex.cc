#include "unittests_trimesh_circulator_edge.hh"

#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>

#include <vector>

namespace {

class OpenMeshTrimeshCirculatorEdgeVertex : public OpenMeshTrimeshCirculatorEdge {
public:
    // Check that the _visited vertices comply with the expected iteration order on _eh.
    void CheckIteration(EH _eh, const std::vector<VH>& _visited) {
        ASSERT_TRUE(mesh_.is_valid_handle(_eh));
        const auto eh = make_smart(_eh, mesh_);
        std::vector<VH> expected;
        for (int i = 0; i < 2; ++i) {
            expected.push_back(eh.vertex(i));
        }

        // Elements must have been visited in the expected order.
        ASSERT_EQ(_visited.size(), expected.size());
        for (size_t i = 0; i < _visited.size(); ++i) {
            EXPECT_EQ(_visited[i], expected[i]);
        }
    }
};

// Mutable mesh.

TEST_F(OpenMeshTrimeshCirculatorEdgeVertex, Iter) {
    for (auto eh : mesh_.edges()) {
        std::vector<VH> visited;
        for (auto it = mesh_.ev_iter(eh); it.is_valid(); ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeVertex, BeginEnd) {
    for (auto eh : mesh_.edges()) {
        std::vector<VH> visited;
        for (auto it = mesh_.ev_begin(eh), end = mesh_.ev_end(eh); it != end; ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeVertex, Range) {
    for (auto eh : mesh_.edges()) {
        std::vector<VH> visited;
        for (auto vh : mesh_.ev_range(eh)) {
            visited.push_back(vh);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeVertex, SmartHandleRange) {
    for (auto eh : mesh_.edges()) {
        std::vector<VH> visited;
        auto smart_eh = make_smart(eh, mesh_);
        for (auto vh : smart_eh.vertices()) {
            visited.push_back(vh);
        }
        CheckIteration(eh, visited);
    }
}

// const mesh.

TEST_F(OpenMeshTrimeshCirculatorEdgeVertex, ConstIter) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<VH> visited;
        for (auto it = const_mesh.cev_iter(eh); it.is_valid(); ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeVertex, ConstBeginEnd) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<VH> visited;
        for (auto it = const_mesh.cev_begin(eh), end = const_mesh.cev_end(eh); it != end; ++it) {
            visited.push_back(*it);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeVertex, ConstRange) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<VH> visited;
        for (auto vh : const_mesh.ev_range(eh)) {
            visited.push_back(vh);
        }
        CheckIteration(eh, visited);
    }
}

TEST_F(OpenMeshTrimeshCirculatorEdgeVertex, ConstSmartHandleRange) {
    const auto& const_mesh = mesh_;
    for (auto eh : const_mesh.edges()) {
        std::vector<VH> visited;
        auto smart_eh = make_smart(eh, const_mesh);
        for (auto vh : smart_eh.vertices()) {
            visited.push_back(vh);
        }
        CheckIteration(eh, visited);
    }
}

}
