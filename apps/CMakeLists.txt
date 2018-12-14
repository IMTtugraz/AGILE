# $Id: CMakeLists.txt 506 2011-07-27 10:17:40Z freiberger $

include_directories ("${AGILE_SOURCE_DIR}/include"
                     "${AGILE_BINARY_DIR}/include"
                    )

link_libraries(agile)

###########################################################################
# All app sub-directories.
###########################################################################
add_subdirectory (comp_spmv)
add_subdirectory (interp2d)
add_subdirectory (spmvm)
add_subdirectory (gridding)
add_subdirectory (imt_irgn)

if(ENABLE_GPU_DOUBLE)
  add_subdirectory (tgv_radial_image_recon)
  endif(ENABLE_GPU_DOUBLE)

# End of $Id: CMakeLists.txt 506 2011-07-27 10:17:40Z freiberger $.
