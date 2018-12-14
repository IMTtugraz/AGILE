# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

# List of source files which have to be compiled for this app.
#set (SOURCES source.cpp)

# The resulting executables.
add_executable(gridding_backward.exe ./src/gridding_backward.cpp)
add_executable(gridding_forward.exe ./src/gridding_forward.cpp)
add_executable(gridding_loop.exe ./src/gridding_loop.cpp)

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.
