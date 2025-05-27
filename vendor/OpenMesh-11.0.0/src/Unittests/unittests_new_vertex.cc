#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>
#include <iostream>

namespace {
    
class OpenMeshNewVertexTriangleMesh : public OpenMeshBase {

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

class OpenMeshNewVertexPolyMesh : public OpenMeshBasePoly {

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

/* Takes a vertex position directly from the mesh and readds it as a new vertex position
 */
TEST_F(OpenMeshNewVertexTriangleMesh, CopyVertexinsideMeshTriangle) {

  mesh_.clear();

  // Add some vertices
  Mesh::VertexHandle vhandle[4];

  vhandle[0] = mesh_.add_vertex(Mesh::Point(0, 0, 0));
  vhandle[1] = mesh_.add_vertex(mesh_.point(Mesh::VertexHandle(0)));
  vhandle[2] = mesh_.add_vertex(mesh_.point(Mesh::VertexHandle(1)));
  vhandle[3] = mesh_.add_vertex(mesh_.point(Mesh::VertexHandle(2)));


  // Check setup
  EXPECT_EQ(4u, mesh_.n_vertices() ) << "Wrong number of vertices";
  EXPECT_EQ(0u, mesh_.n_faces() )    << "Wrong number of faces";

}

/* Takes a vertex position directly from the mesh and readds it as a new vertex position
 */
TEST_F(OpenMeshNewVertexPolyMesh, CopyVertexinsideMeshPoly) {

  mesh_.clear();

  // Add some vertices
  Mesh::VertexHandle vhandle[4];

  vhandle[0] = mesh_.add_vertex(Mesh::Point(0, 0, 0));
  vhandle[1] = mesh_.add_vertex(mesh_.point(Mesh::VertexHandle(0)));
  vhandle[2] = mesh_.add_vertex(mesh_.point(Mesh::VertexHandle(1)));
  vhandle[3] = mesh_.add_vertex(mesh_.point(Mesh::VertexHandle(2)));

  // Check setup
  EXPECT_EQ(4u, mesh_.n_vertices() ) << "Wrong number of vertices";
  EXPECT_EQ(0u, mesh_.n_faces() )    << "Wrong number of faces";

}

}
