include (VCICommon)

vci_add_executable (Analyzer vdpmanalyzer.cc)

target_link_libraries (Analyzer
  OpenMeshCore
  OpenMeshTools
)

