# Specify the source files for the library
set(SOURCES
  PDSPTree.C
  AnalyzerCore.C
  dEdx_res.C
  Fit_Modified_Box.C
  MCS_Angle_Shapes.C
  MCS_Performance.C
  MCS_Tree_Maker.C
  PionKEScale.C
  PionQE.C
  Secondary_proton_ana.C
)

# Create the shared library
add_library(Analyzers SHARED ${SOURCES})

# Specify output directory for the shared library
set_target_properties(Analyzers PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib
)

# Link the library with dir1_lib and dir3_lib
target_link_libraries(Analyzers PRIVATE AnalyzerTools DataFormats ${ROOT_LIBRARIES})

# Specify include directories for this library
target_include_directories(Analyzers PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
