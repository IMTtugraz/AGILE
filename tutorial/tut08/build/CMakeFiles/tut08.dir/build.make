# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dieheart/workspace/AGILE/tutorial/tut08

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dieheart/workspace/AGILE/tutorial/tut08/build

# Include any dependencies generated for this target.
include CMakeFiles/tut08.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/tut08.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tut08.dir/flags.make

CMakeFiles/tut08.dir/program.cpp.o: CMakeFiles/tut08.dir/flags.make
CMakeFiles/tut08.dir/program.cpp.o: ../program.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/dieheart/workspace/AGILE/tutorial/tut08/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/tut08.dir/program.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/tut08.dir/program.cpp.o -c /home/dieheart/workspace/AGILE/tutorial/tut08/program.cpp

CMakeFiles/tut08.dir/program.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tut08.dir/program.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/dieheart/workspace/AGILE/tutorial/tut08/program.cpp > CMakeFiles/tut08.dir/program.cpp.i

CMakeFiles/tut08.dir/program.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tut08.dir/program.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/dieheart/workspace/AGILE/tutorial/tut08/program.cpp -o CMakeFiles/tut08.dir/program.cpp.s

CMakeFiles/tut08.dir/program.cpp.o.requires:
.PHONY : CMakeFiles/tut08.dir/program.cpp.o.requires

CMakeFiles/tut08.dir/program.cpp.o.provides: CMakeFiles/tut08.dir/program.cpp.o.requires
	$(MAKE) -f CMakeFiles/tut08.dir/build.make CMakeFiles/tut08.dir/program.cpp.o.provides.build
.PHONY : CMakeFiles/tut08.dir/program.cpp.o.provides

CMakeFiles/tut08.dir/program.cpp.o.provides.build: CMakeFiles/tut08.dir/program.cpp.o

# Object files for target tut08
tut08_OBJECTS = \
"CMakeFiles/tut08.dir/program.cpp.o"

# External object files for target tut08
tut08_EXTERNAL_OBJECTS =

tut08: CMakeFiles/tut08.dir/program.cpp.o
tut08: CMakeFiles/tut08.dir/build.make
tut08: /usr/local/cuda/lib64/libcudart.so
tut08: /usr/local/cuda/lib64/libcublas.so
tut08: /usr/local/lib/libagile.a
tut08: CMakeFiles/tut08.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable tut08"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tut08.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tut08.dir/build: tut08
.PHONY : CMakeFiles/tut08.dir/build

CMakeFiles/tut08.dir/requires: CMakeFiles/tut08.dir/program.cpp.o.requires
.PHONY : CMakeFiles/tut08.dir/requires

CMakeFiles/tut08.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tut08.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tut08.dir/clean

CMakeFiles/tut08.dir/depend:
	cd /home/dieheart/workspace/AGILE/tutorial/tut08/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dieheart/workspace/AGILE/tutorial/tut08 /home/dieheart/workspace/AGILE/tutorial/tut08 /home/dieheart/workspace/AGILE/tutorial/tut08/build /home/dieheart/workspace/AGILE/tutorial/tut08/build /home/dieheart/workspace/AGILE/tutorial/tut08/build/CMakeFiles/tut08.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tut08.dir/depend

