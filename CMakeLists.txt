cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project( libvtk )

# Set a default build type if none was specified
set(default_build_type "Release")

add_compile_options(-fPIC)

# FEBIO
if(NOT EXISTS ${FEBIO_LIB_DIR}/libfecore.a)
  message(SEND_ERROR "Could not find FEBio library (libfecore.a). Check FEBIO_LIB_DIR.")
  return()
endif()
ADD_LIBRARY(febio STATIC IMPORTED)
SET_TARGET_PROPERTIES(febio PROPERTIES IMPORTED_LOCATION ${FEBIO_LIB_DIR}/libfecore.a)

# OPENMP
find_package(OpenMP)
if ( NOT OPENMP_FOUND )
  message(STATUS "This project requires OpenMP, and will not be compiled.")
  return()  
endif()
add_compile_options(-fopenmp)

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

target_link_libraries(vtk_callback ${VTK_LIBRARIES} febio)
target_include_directories(vtk_callback PRIVATE ${FEBIO_ROOT})


