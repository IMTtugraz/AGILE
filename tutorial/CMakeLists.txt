# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

# This is the common CMake control file for all tutorials.
include_directories ("${AGILE_SOURCE_DIR}/include"
                     "${AGILE_BINARY_DIR}/include"
                    )

link_libraries(agile)

###########################################################################
# The sub-directories to process. Please do not include tutxx in this list
# as this is only a template for a new tutorial.
###########################################################################
add_subdirectory (tut01)
add_subdirectory (tut02)
add_subdirectory (tut03)
add_subdirectory (tut04)
add_subdirectory (tut05)
add_subdirectory (tut06)
add_subdirectory (tut07)
add_subdirectory (tut08)

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.
