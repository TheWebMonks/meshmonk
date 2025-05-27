#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>
#include <iostream>

namespace {
    
class OpenMeshPrevHalfedge : public OpenMeshBase {

    protected:

        // This function is called before each test is run
        virtual void SetUp() {
            
            // Do some initial stuff with the member data here...
        }

        // This function is called after all tests are through
        virtual void TearDown() {

            // Do some final stuff with the member data here...
        }

    // Member already defined in OpenMeshBase
    //Mesh mesh_;  
};

/*
 * ====================================================================
 * Define tests below
 * ====================================================================
 */

/* creates a Halfedge WITHOUT the PrevHalfedge Attribute
 */
TEST_F(OpenMeshPrevHalfedge, RemovePrevHalfedge) {

    mesh_.clear();

    struct MeshTraits : OpenMesh::DefaultTraits {
        HalfedgeAttributes(0);
    };
    using MeshT = OpenMesh::TriMesh_ArrayKernelT<MeshTraits>;


    // Check if Prev Halfedge is referenced
    EXPECT_EQ(12u, sizeof(MeshT::Halfedge) ) << "Wrong size of Halfedge";
    EXPECT_EQ(0u, MeshT::HasPrevHalfedge::my_bool )    << "Attribute HasPrevHalfedge is wrong";
}

/* creates a Halfedge WITH the PrevHalfedge Attribute
 */
TEST_F(OpenMeshPrevHalfedge, HavePrevHalfedge) {

    mesh_.clear();

    struct MeshTraits : OpenMesh::DefaultTraits {
        //HalfedgeAttributes(0);
    };
    using MeshT = OpenMesh::TriMesh_ArrayKernelT<MeshTraits>;


    // Check if Prev Halfedge is referenced
    EXPECT_EQ(16u, sizeof(MeshT::Halfedge) ) << "Wrong size of Halfedge";
    EXPECT_EQ(1u, MeshT::HasPrevHalfedge::my_bool )    << "Attribute HasPrevHalfedge is wrong";
}

}
