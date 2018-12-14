# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

# The resulting executable.
add_executable (matlab_fft_test fft_test.cpp)
add_executable (matlab_fft_test_kspace fft_test_kspace.cpp)

# We need CuFFT for this test
cuda_add_cufft_to_target (matlab_fft_test)
cuda_add_cufft_to_target (matlab_fft_test_kspace)

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.

