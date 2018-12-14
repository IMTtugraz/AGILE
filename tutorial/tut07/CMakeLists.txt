# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

# The name of the current tutorial. This is simply the last part of the
# current source directory (i.e. id should be one of tut01, tut02, tut03...).
get_filename_component (TUTORIAL ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# List of source files which have to be compiled for this tutorial.
set (SOURCES program.cpp)

# The resulting executable.
# TODO
#add_executable(${TUTORIAL} ${SOURCES})

# A custom target to run the tutorial.
# TODO
#add_custom_target(${TUTORIAL}-run ${TUTORIAL}
#                  DEPENDS ${TUTORIAL}
#                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.

