# Specify the source files for the library
set(SOURCES
  Daughter.C
  Event.C
  MCSSegment.C
)

# Create the shared library
add_library(DataFormats SHARED ${SOURCES})

# Specify output directory for the shared library
set_target_properties(DataFormats PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib
)

# Link the library with ROOT
target_link_libraries(DataFormats PRIVATE ${ROOT_LIBRARIES})

# Optionally, you can specify include directories for this library
target_include_directories(DataFormats PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
