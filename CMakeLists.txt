# $Id: CMakeLists.txt 494 2011-07-26 10:55:57Z huber $

# Top-level CMakeLists.txt for the AGILE project

# The name of the project and the compilers
project (AGILE CXX)

# Set the minimum CMake version required
cmake_minimum_required (VERSION 2.6)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES -std=c++11")

ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 

###########################################################################
# Paths for out-of-source build.
###########################################################################
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${AGILE_BINARY_DIR}/lib)
message (STATUS "Library destination directory: ${AGILE_BINARY_DIR}/lib")

###########################################################################
# Compiler options for configurations.
###########################################################################
# If no configuration is set, it well default to 'Release'.
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  message (STATUS "Switching to default configuration 'Release'")
  set (CMAKE_BUILD_TYPE "Release")
endif()

# TODO: set reasonable CUDA_CXX_FLAGS

###########################################################################
# Default values for library switches.
###########################################################################
set (HAVE_MPI 0)

###########################################################################
# locate CUDA
###########################################################################
find_package (CUDA 2.2 REQUIRED)
set(CUDA_HOST_COMPILER ${CMAKE_C_COMPILER})
include_directories (${CUDA_INCLUDE_DIRS})
set(CUDA_LIBRARIES ${CUDA_LIBRARIES} cuda)
link_libraries (${CUDA_LIBRARIES})

set(ENABLE_GPU_DOUBLE ON CACHE BOOL "Set flag to enable double-precision data types on GPU")

if(ENABLE_GPU_DOUBLE)
  SET(CUDA_NVCC_FLAGS "-arch;sm_30 --expt-relaxed-constexpr")
endif(ENABLE_GPU_DOUBLE)

###########################################################################
# added by A. Huber
# for CUBLAS functionality
###########################################################################
link_libraries (${CUDA_CUBLAS_LIBRARIES})


###########################################################################
# added by H. Heigl
# for DICOM functionality
###########################################################################
ADD_DEFINITIONS(-DHAVE_CONFIG_H)

find_package (DCMTK REQUIRED)
if (DCMTK_FOUND)
  MESSAGE(STATUS "DCMTK_FOUND")
  include_directories(${DCMTK_INCLUDE_DIRS})
  link_libraries(${DCMTK_LIBRARIES})
endif()

###########################################################################
# try to locate the MPI library
###########################################################################
SET(ENABLE_MPI OFF CACHE BOOL "Enable/Disable MPI support.")
find_package (MPI)
if (MPI_FOUND AND ENABLE_MPI)
  set (HAVE_MPI 1)
  include_directories (${MPI_INCLUDE_PATH})
  link_libraries (${MPI_LIBRARIES})
endif()

###########################################################################
# the sub-directories to process
###########################################################################
add_subdirectory (doc)
add_subdirectory (include)
add_subdirectory (src)
add_subdirectory (test)
add_subdirectory (tutorial)
add_subdirectory (apps)
add_subdirectory (unittest)

# End of $Id: CMakeLists.txt 494 2011-07-26 10:55:57Z huber $.
