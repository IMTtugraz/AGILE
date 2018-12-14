# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
	     	    ${CUDA_SDK_ROOT_DIR}/C/common/inc)

link_libraries("-L/${CUDA_SDK_ROOT_DIR}/C/lib")

# List of source files which have to be compiled for this test.
set (SOURCES1 nfft2d_tgv_recon_nogfx.cu nfft2d.cu ParseParam.cu)

# The resulting executable.
cuda_add_executable(nfft2d_tgv_recon_nogfx ${SOURCES1})
target_link_libraries(nfft2d_tgv_recon_nogfx ${CUDA_cufft_LIBRARY})

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.

