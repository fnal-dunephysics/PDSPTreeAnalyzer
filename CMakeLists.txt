cmake_minimum_required(VERSION 3.19.6 FATAL_ERROR)
project(PDSPTreeAnalyzer)

# require C++11
set (CMAKE_CXX_STANDARD 17)

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED COMPONENTS RIO MathMore)
if(NOT ROOT_FOUND)
  message(STATUS "ROOT package not found.")
  return()
endif()


###### Compiler set up
# Add ROOT include directories
include_directories(${ROOT_INCLUDE_DIRS})

# add sub directories
add_subdirectory(DataFormats)
add_subdirectory(AnalyzerTools)
add_subdirectory(Analyzers)
