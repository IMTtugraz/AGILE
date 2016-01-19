# Install script for directory: /home2/GIT/AGILE/test

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home2/GIT/AGILE/build/test/dense_matrix_IO_Test/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/file_io_test/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gmres_test/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_crs_matrix/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_matrix_float/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_matrix_stdcomplexfloat/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_matrix_pitched/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_matrix_vs_gpu_matrix_pitched_float/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_matrix_vs_gpu_matrix_pitched_stdcomplexfloat/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_vector/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/lsqr_test/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/matlab_3d_matrix_I_O/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/matlab_fft_test/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/pcg_test/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/scalarproduct_test/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/vector_example/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/irgn/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_matrix_double/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_matrix_stdcomplexdouble/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_matrix_vs_gpu_matrix_pitched_double/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/test/gpu_matrix_vs_gpu_matrix_pitched_stdcomplexdouble/cmake_install.cmake")

endif()

