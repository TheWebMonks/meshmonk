include (VCICommon)

vci_add_executable (commandlineDecimater ../decimater.cc)

target_link_libraries (commandlineDecimater
  OpenMeshCore
  OpenMeshTools
)

