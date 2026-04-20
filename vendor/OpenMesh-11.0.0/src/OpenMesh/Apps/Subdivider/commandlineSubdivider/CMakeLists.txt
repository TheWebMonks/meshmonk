include (VCICommon)

vci_add_executable (commandlineSubdivider ../subdivider.cc)

target_link_libraries (commandlineSubdivider
  OpenMeshCore
  OpenMeshTools
)

