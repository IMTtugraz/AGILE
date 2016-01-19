# $Id: CMakeLists.txt 494 2011-07-26 10:55:57Z huber $

include_directories ("${AGILE_SOURCE_DIR}/include"
                     "${AGILE_BINARY_DIR}/include"
                    )

link_libraries(agile)
#link_libraries(agile_nvcc)

###########################################################################
# All test sub-directories.
###########################################################################
add_subdirectory (dense_matrix_IO_Test)
add_subdirectory (file_io_test)
add_subdirectory (gmres_test)
add_subdirectory (gpu_crs_matrix)
add_subdirectory (gpu_matrix_float)
add_subdirectory (gpu_matrix_stdcomplexfloat)
add_subdirectory (gpu_matrix_pitched)
add_subdirectory (gpu_matrix_vs_gpu_matrix_pitched_float)
add_subdirectory (gpu_matrix_vs_gpu_matrix_pitched_stdcomplexfloat)
add_subdirectory (gpu_vector)
add_subdirectory (lsqr_test)
add_subdirectory (matlab_3d_matrix_I_O)
add_subdirectory (matlab_fft_test)
add_subdirectory (pcg_test)
add_subdirectory (scalarproduct_test)
add_subdirectory (vector_example)
add_subdirectory (irgn)

if(ENABLE_GPU_DOUBLE)
add_subdirectory (gpu_matrix_double)
add_subdirectory (gpu_matrix_stdcomplexdouble)
add_subdirectory (gpu_matrix_vs_gpu_matrix_pitched_double)
add_subdirectory (gpu_matrix_vs_gpu_matrix_pitched_stdcomplexdouble)
endif(ENABLE_GPU_DOUBLE)
# End of $Id: CMakeLists.txt 494 2011-07-26 10:55:57Z huber $.

