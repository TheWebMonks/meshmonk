
#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>
#include <OpenMesh/Tools/HoleFiller/HoleFillerT.hh>

namespace {

class OpenMeshHoleFiller_Triangle : public OpenMeshBase {

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

/*
 */
TEST_F(OpenMeshHoleFiller_Triangle,Triangle_Hole_Filling) {

  mesh_.clear();


  bool ok = OpenMesh::IO::read_mesh(mesh_, "cube_2holes.off");

  ASSERT_TRUE(ok);

  // Check setup
  EXPECT_EQ(1456u, mesh_.n_vertices() ) << "Wrong number of vertices";
  EXPECT_EQ(2864u, mesh_.n_faces() )    << "Wrong number of faces";


  // Initialize subdivider
  OpenMesh::HoleFiller::HoleFillerT<Mesh> filler(mesh_);


  // Execute the algorithm
  filler.fill_all_holes();

  if ( std::is_same<double,typename Mesh::Scalar>() ) {
      EXPECT_EQ(1504u, mesh_.n_vertices() ) << "Wrong number of vertices after smoothing?";
      EXPECT_EQ(3004u, mesh_.n_faces() )    << "Wrong number of faces after smoothing?";
  } else {
      EXPECT_EQ(1507u, mesh_.n_vertices() ) << "Wrong number of vertices after smoothing?";
      EXPECT_EQ(3010u, mesh_.n_faces() )    << "Wrong number of faces after smoothing?";
  }

}

}
