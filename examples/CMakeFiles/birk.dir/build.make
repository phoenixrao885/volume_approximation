# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/vissarion/workspace/randgeom-code/examples

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vissarion/workspace/randgeom-code/examples

# Include any dependencies generated for this target.
include CMakeFiles/birk.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/birk.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/birk.dir/flags.make

CMakeFiles/birk.dir/polytope_generators/birk.o: CMakeFiles/birk.dir/flags.make
CMakeFiles/birk.dir/polytope_generators/birk.o: polytope_generators/birk.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vissarion/workspace/randgeom-code/examples/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/birk.dir/polytope_generators/birk.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/birk.dir/polytope_generators/birk.o -c /home/vissarion/workspace/randgeom-code/examples/polytope_generators/birk.cpp

CMakeFiles/birk.dir/polytope_generators/birk.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/birk.dir/polytope_generators/birk.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/vissarion/workspace/randgeom-code/examples/polytope_generators/birk.cpp > CMakeFiles/birk.dir/polytope_generators/birk.i

CMakeFiles/birk.dir/polytope_generators/birk.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/birk.dir/polytope_generators/birk.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/vissarion/workspace/randgeom-code/examples/polytope_generators/birk.cpp -o CMakeFiles/birk.dir/polytope_generators/birk.s

CMakeFiles/birk.dir/polytope_generators/birk.o.requires:
.PHONY : CMakeFiles/birk.dir/polytope_generators/birk.o.requires

CMakeFiles/birk.dir/polytope_generators/birk.o.provides: CMakeFiles/birk.dir/polytope_generators/birk.o.requires
	$(MAKE) -f CMakeFiles/birk.dir/build.make CMakeFiles/birk.dir/polytope_generators/birk.o.provides.build
.PHONY : CMakeFiles/birk.dir/polytope_generators/birk.o.provides

CMakeFiles/birk.dir/polytope_generators/birk.o.provides.build: CMakeFiles/birk.dir/polytope_generators/birk.o

# Object files for target birk
birk_OBJECTS = \
"CMakeFiles/birk.dir/polytope_generators/birk.o"

# External object files for target birk
birk_EXTERNAL_OBJECTS =

birk: CMakeFiles/birk.dir/polytope_generators/birk.o
birk: CMakeFiles/birk.dir/build.make
birk: /usr/lib/x86_64-linux-gnu/libmpfr.so
birk: /usr/lib/x86_64-linux-gnu/libgmp.so
birk: /home/vissarion/workspace/CGAL-4.4/lib/libCGAL_Core.so
birk: /home/vissarion/workspace/CGAL-4.4/lib/libCGAL.so
birk: /usr/lib/x86_64-linux-gnu/libboost_thread.so
birk: /usr/lib/x86_64-linux-gnu/libboost_system.so
birk: /usr/lib/x86_64-linux-gnu/libpthread.so
birk: /home/vissarion/workspace/CGAL-4.4/lib/libCGAL_Core.so
birk: /home/vissarion/workspace/CGAL-4.4/lib/libCGAL.so
birk: /usr/lib/x86_64-linux-gnu/libboost_thread.so
birk: /usr/lib/x86_64-linux-gnu/libboost_system.so
birk: /usr/lib/x86_64-linux-gnu/libpthread.so
birk: CMakeFiles/birk.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable birk"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/birk.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/birk.dir/build: birk
.PHONY : CMakeFiles/birk.dir/build

CMakeFiles/birk.dir/requires: CMakeFiles/birk.dir/polytope_generators/birk.o.requires
.PHONY : CMakeFiles/birk.dir/requires

CMakeFiles/birk.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/birk.dir/cmake_clean.cmake
.PHONY : CMakeFiles/birk.dir/clean

CMakeFiles/birk.dir/depend:
	cd /home/vissarion/workspace/randgeom-code/examples && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vissarion/workspace/randgeom-code/examples /home/vissarion/workspace/randgeom-code/examples /home/vissarion/workspace/randgeom-code/examples /home/vissarion/workspace/randgeom-code/examples /home/vissarion/workspace/randgeom-code/examples/CMakeFiles/birk.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/birk.dir/depend
