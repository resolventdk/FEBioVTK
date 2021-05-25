# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( libvtk )

cmake_minimum_required(VERSION 3.9)

# Compiler settings
set(CMAKE_CXX, "g++")
set(CMAKE_CXX_FLAGS "-fPIC -g")

# FEBIO
SET(FEBIO_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../)
ADD_LIBRARY(febio STATIC IMPORTED)
SET_TARGET_PROPERTIES(febio PROPERTIES IMPORTED_LOCATION ${FEBIO_ROOT}/build/lib/libfecore.a)

# OPENMP
find_package(OpenMP)

# VTK
find_package(VTK 8.2 REQUIRED)
if ( NOT VTK_FOUND )
  message(STATUS "This project requires the VTK library, and will not be compiled.")
  return()  
endif()
include(${VTK_USE_FILE}) # include helper file
message(STATUS "Using VTK version: ${VTK_VERSION}")  # remember to put vtk-library dir in LD_LIBRARY_PATH for this cmake crap to work

# The plugin library
add_library(vtk_callback SHARED
  src/callback.cpp
  src/FEVTKExport.cpp
  )

target_link_libraries(vtk_callback ${VTK_LIBRARIES} febio OpenMP::OpenMP_CXX)
target_include_directories(vtk_callback PRIVATE ${FEBIO_ROOT})


