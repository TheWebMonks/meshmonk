include (VCICommon)

vci_add_executable (Smoothing smooth.cc)

target_link_libraries (Smoothing
  OpenMeshCore
  OpenMeshTools
)

