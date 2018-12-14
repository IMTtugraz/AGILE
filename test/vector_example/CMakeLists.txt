# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

# List of source files which have to be compiled for this test.
set (SOURCES example.cu)

# The resulting executable.
cuda_add_executable(vector_example ${SOURCES})

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.

