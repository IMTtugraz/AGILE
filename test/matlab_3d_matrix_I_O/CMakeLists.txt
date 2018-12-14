# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

# The resulting executable.
add_executable (3D_matrix_I_O 3D_matrix_I_O.cpp)

# Copy Matlab stuff in build directory
configure_file ("${AGILE_SOURCE_DIR}/test/matlab_3d_matrix_I_O/img.bin"
                  "${AGILE_BINARY_DIR}/test/matlab_3d_matrix_I_O/img.bin" COPYONLY)

configure_file ("${AGILE_SOURCE_DIR}/test/matlab_3d_matrix_I_O/t2brainRawdata.mat"
                  "${AGILE_BINARY_DIR}/test/matlab_3d_matrix_I_O/t2brainRawdata.mat" COPYONLY)

configure_file ("${AGILE_SOURCE_DIR}/test/matlab_3d_matrix_I_O/matlab_test_3Dmatrix_I_O.m"
                  "${AGILE_BINARY_DIR}/test/matlab_3d_matrix_I_O/matlab_test_3Dmatrix_I_O.m" COPYONLY)

configure_file ("${AGILE_SOURCE_DIR}/test/matlab_3d_matrix_I_O/readBin2Matlab_3D.m"
                  "${AGILE_BINARY_DIR}/test/matlab_3d_matrix_I_O/readBin2Matlab_3D.m" COPYONLY)

configure_file ("${AGILE_SOURCE_DIR}/test/matlab_3d_matrix_I_O/writeMatlab2Bin_3D.m"
                  "${AGILE_BINARY_DIR}/test/matlab_3d_matrix_I_O/writeMatlab2Bin_3D.m" COPYONLY)

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.

