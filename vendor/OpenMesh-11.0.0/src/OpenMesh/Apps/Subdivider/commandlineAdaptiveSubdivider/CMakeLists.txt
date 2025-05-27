include (VCICommon)

vci_add_executable (commandlineAdaptiveSubdivider ../adaptive_subdivider.cc)

target_link_libraries (commandlineAdaptiveSubdivider
  OpenMeshCore
  OpenMeshTools
)

