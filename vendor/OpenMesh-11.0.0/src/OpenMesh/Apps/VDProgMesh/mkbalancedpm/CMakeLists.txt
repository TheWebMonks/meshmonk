include (VCICommon)

vci_add_executable (mkbalancedpm mkbalancedpm.cc)

target_link_libraries (mkbalancedpm
  OpenMeshCore
  OpenMeshTools
)

