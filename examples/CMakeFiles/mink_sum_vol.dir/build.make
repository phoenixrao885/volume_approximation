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
include CMakeFiles/mink_sum_vol.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mink_sum_vol.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mink_sum_vol.dir/flags.make

CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o: CMakeFiles/mink_sum_vol.dir/flags.make
CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o: mink_sum_vol.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/vissarion/workspace/randgeom-code/examples/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o -c /home/vissarion/workspace/randgeom-code/examples/mink_sum_vol.cpp

CMakeFiles/mink_sum_vol.dir/mink_sum_vol.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mink_sum_vol.dir/mink_sum_vol.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/vissarion/workspace/randgeom-code/examples/mink_sum_vol.cpp > CMakeFiles/mink_sum_vol.dir/mink_sum_vol.i

CMakeFiles/mink_sum_vol.dir/mink_sum_vol.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mink_sum_vol.dir/mink_sum_vol.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/vissarion/workspace/randgeom-code/examples/mink_sum_vol.cpp -o CMakeFiles/mink_sum_vol.dir/mink_sum_vol.s

CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o.requires:
.PHONY : CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o.requires

CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o.provides: CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o.requires
	$(MAKE) -f CMakeFiles/mink_sum_vol.dir/build.make CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o.provides.build
.PHONY : CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o.provides

CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o.provides.build: CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o

# Object files for target mink_sum_vol
mink_sum_vol_OBJECTS = \
"CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o"

# External object files for target mink_sum_vol
mink_sum_vol_EXTERNAL_OBJECTS =

mink_sum_vol: CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o
mink_sum_vol: CMakeFiles/mink_sum_vol.dir/build.make
mink_sum_vol: /usr/lib/x86_64-linux-gnu/libmpfr.so
mink_sum_vol: /usr/lib/x86_64-linux-gnu/libgmp.so
mink_sum_vol: /home/vissarion/workspace/CGAL-4.2/lib/libCGAL_Core.so
mink_sum_vol: /home/vissarion/workspace/CGAL-4.2/lib/libCGAL.so
mink_sum_vol: /usr/lib/libboost_thread-mt.so
mink_sum_vol: /usr/lib/libboost_system-mt.so
mink_sum_vol: /usr/lib/x86_64-linux-gnu/libpthread.so
mink_sum_vol: /home/vissarion/workspace/CGAL-4.2/lib/libCGAL_Core.so
mink_sum_vol: /home/vissarion/workspace/CGAL-4.2/lib/libCGAL.so
mink_sum_vol: /usr/lib/libboost_thread-mt.so
mink_sum_vol: /usr/lib/libboost_system-mt.so
mink_sum_vol: /usr/lib/x86_64-linux-gnu/libpthread.so
mink_sum_vol: CMakeFiles/mink_sum_vol.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable mink_sum_vol"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mink_sum_vol.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mink_sum_vol.dir/build: mink_sum_vol
.PHONY : CMakeFiles/mink_sum_vol.dir/build

CMakeFiles/mink_sum_vol.dir/requires: CMakeFiles/mink_sum_vol.dir/mink_sum_vol.o.requires
.PHONY : CMakeFiles/mink_sum_vol.dir/requires

CMakeFiles/mink_sum_vol.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mink_sum_vol.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mink_sum_vol.dir/clean

CMakeFiles/mink_sum_vol.dir/depend:
	cd /home/vissarion/workspace/randgeom-code/examples && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vissarion/workspace/randgeom-code/examples /home/vissarion/workspace/randgeom-code/examples /home/vissarion/workspace/randgeom-code/examples /home/vissarion/workspace/randgeom-code/examples /home/vissarion/workspace/randgeom-code/examples/CMakeFiles/mink_sum_vol.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mink_sum_vol.dir/depend
