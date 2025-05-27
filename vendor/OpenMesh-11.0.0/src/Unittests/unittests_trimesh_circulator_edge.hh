#include <Unittests/unittests_common.hh>

#include <vector>

namespace {

class OpenMeshTrimeshCirculatorEdge : public OpenMeshBase {
public:
    using VH = OpenMesh::VertexHandle;
    using FH = OpenMesh::FaceHandle;
    using EH = OpenMesh::EdgeHandle;
    using HEH = OpenMesh::HalfedgeHandle;

    // This function is called before each test is run.
    void SetUp() override {
        std::vector<VH> vh;

        //   3----2    5
        //   |   /|    |
        //   |  / |    |
        //   | /  |    |
        //   |/   |    |
        //   0----1    4

        // A quad consisting of two triangles.
        vh.push_back(mesh_.add_vertex({0.0, 0.0, 0.0})); // vh[0]
        vh.push_back(mesh_.add_vertex({1.0, 0.0, 0.0})); // vh[1]
        vh.push_back(mesh_.add_vertex({1.0, 1.0, 0.0})); // vh[2]
        vh.push_back(mesh_.add_vertex({0.0, 1.0, 0.0})); // vh[3]
        mesh_.add_face(vh[0], vh[1], vh[2]);
        mesh_.add_face(vh[0], vh[2], vh[3]);

        // An isolated edge.
        vh.push_back(mesh_.add_vertex({2.0, 0.0, 0.0})); // vh[4]
        vh.push_back(mesh_.add_vertex({2.0, 1.0, 0.0})); // vh[5]
        auto heh = mesh_.new_edge(vh[4], vh[5]);
        auto heh_opp = mesh_.opposite_halfedge_handle(heh);
        mesh_.set_halfedge_handle(vh[4], heh);
        mesh_.set_halfedge_handle(vh[5], heh_opp);

        InteriorEdge = Edge(0, 2);
        BoundaryEdge = Edge(0, 1);
        IsolatedEdge = Edge(4, 5);
    }

    // Helper function to quickly retrieve edges from their endpoint vertex IDs.
    EH Edge(int _vh0_idx, int _vh1_idx) {
        VH vh0(_vh0_idx);
        VH vh1(_vh1_idx);
        HEH heh = mesh_.find_halfedge(vh0, vh1);
        if (!heh.is_valid())
            return EH();
        return mesh_.edge_handle(heh);
    }

    // Handles of some interesting edges.
    EH InteriorEdge;
    EH BoundaryEdge;
    EH IsolatedEdge;
};

}
