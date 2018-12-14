# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

# generate the configuration file
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/gpu_config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/gpu_config.hpp)

# install this file
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/gpu_config.hpp
         DESTINATION include/agile/)

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.
