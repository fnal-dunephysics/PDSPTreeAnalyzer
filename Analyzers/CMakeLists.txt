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
  KaonTOF.C
  PionKI.C
  PionStudy.C
  Pion2Proton.C
  PionSkimTreeMaker.C
)

set(HDRFILES
  PDSPTree.h
  AnalyzerCore.h
  dEdx_res.h
  Fit_Modified_Box.h
  MCS_Angle_Shapes.h
  MCS_Performance.h
  MCS_Tree_Maker.h
  PionKEScale.h
  PionQE.h
  Secondary_proton_ana.h
  KaonTOF.h
  PionKI.h
  PionStudy.h
  Pion2Proton.h
  PionSkimTreeMaker.h
)

# Create the shared library
add_library(Analyzers SHARED ${SOURCES})

# Specify output directory for the shared library
set_target_properties(Analyzers PROPERTIES
  PUBLIC_HEADER "${HDRFILES}"
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib
)

#include_directories(${CMAKE_SOURCE_DIR}/DataFormats)
#include_directories(${CMAKE_SOURCE_DIR}/AnalyzerTools)

# Link the library with dir1_lib and dir3_lib
#target_link_libraries(Analyzers PRIVATE AnalyzerTools DataFormats ${ROOT_LIBRARIES})
target_link_libraries(Analyzers
PUBLIC
  AnalyzerTools
  DataFormats
  ${ROOT_LIBRARIES}
)

install(TARGETS Analyzers
    EXPORT PDSPTreeAnalyzer-targets
    PUBLIC_HEADER DESTINATION include/Analyzers COMPONENT Development)
