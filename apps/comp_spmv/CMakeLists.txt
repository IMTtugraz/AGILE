# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

# List of source files which have to be compiled for this app.
set (SOURCES comp_spmv.cpp ./spmv/SpMV_compute.cpp ./spmv/SpMV_gen.cpp)

# The resulting executable.
add_executable(comp_spmv.exe ${SOURCES})

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.
