# Install script for directory: /home2/GIT/AGILE/tutorial

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
  include("/home2/GIT/AGILE/build/tutorial/tut01/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/tutorial/tut02/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/tutorial/tut03/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/tutorial/tut04/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/tutorial/tut05/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/tutorial/tut06/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/tutorial/tut07/cmake_install.cmake")
  include("/home2/GIT/AGILE/build/tutorial/tut08/cmake_install.cmake")

endif()

