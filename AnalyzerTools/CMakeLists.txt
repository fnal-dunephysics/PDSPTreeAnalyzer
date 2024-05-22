# Specify the source files for the library
set(SOURCES
  BetheBloch.C
  Fittings.C
  GEANT4_XS.C
  MCCorrection.C
  TDirectoryHelper.C
)

set(HDRFILES
  BetheBloch.h
  Fittings.h
  GEANT4_XS.h
  MCCorrection.h
  TDirectoryHelper.h
)

# Create the shared library
add_library(AnalyzerTools SHARED ${SOURCES})

# Specify output directory for the shared library
set_target_properties(AnalyzerTools PROPERTIES
  PUBLIC_HEADER "${HDRFILES}"
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib
)

# Link the library
#target_link_libraries(AnalyzerTools PRIVATE DataFormats ${ROOT_LIBRARIES})
target_link_libraries(AnalyzerTools PUBLIC DataFormats ${ROOT_LIBRARIES})

install(TARGETS AnalyzerTools
    EXPORT PDSPTreeAnalyzer-targets
    PUBLIC_HEADER DESTINATION include/AnalyzerTools COMPONENT Development)
