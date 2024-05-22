# Specify the source files for the library
set(SOURCES
  BetheBloch.C
  Fittings.C
  GEANT4_XS.C
  MCCorrection.C
  TDirectoryHelper.C
)

# Create the shared library
add_library(AnalyzerTools SHARED ${SOURCES})

# Specify output directory for the shared library
set_target_properties(AnalyzerTools PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib
)

# Link the library with dir1_lib
target_link_libraries(AnalyzerTools PRIVATE DataFormats ${ROOT_LIBRARIES})

# Optionally, you can specify include directories for this library
target_include_directories(AnalyzerTools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
