# --- Example Executable ---
# Define the executable for example.cpp
add_executable(MeshMonkExample example.cpp)

# --- Link against meshmonk_shared ---
# This ensures it can find meshmonk_shared library and its headers
target_link_libraries(MeshMonkExample PRIVATE meshmonk_shared)

# Optional: If example.cpp needs direct access to Eigen or other specific includes
# not covered by meshmonk_shared's PUBLIC includes, add them here.
# target_include_directories(MeshMonkExample PRIVATE ${EIGEN3_INCLUDE_DIR})
# However, meshmonk_shared should export Eigen include directory correctly if it's a public dependency.
