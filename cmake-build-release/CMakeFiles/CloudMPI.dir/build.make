# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = /home/prashant/CloudMPI

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/prashant/CloudMPI/cmake-build-release

# Include any dependencies generated for this target.
include CMakeFiles/CloudMPI.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CloudMPI.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CloudMPI.dir/flags.make

CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.o: CMakeFiles/CloudMPI.dir/flags.make
CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.o: ../src/cloudmpi.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/prashant/CloudMPI/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.o"
	/usr/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.o -c /home/prashant/CloudMPI/src/cloudmpi.cpp

CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.i"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prashant/CloudMPI/src/cloudmpi.cpp > CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.i

CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.s"
	/usr/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prashant/CloudMPI/src/cloudmpi.cpp -o CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.s

# Object files for target CloudMPI
CloudMPI_OBJECTS = \
"CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.o"

# External object files for target CloudMPI
CloudMPI_EXTERNAL_OBJECTS =

libCloudMPI.a: CMakeFiles/CloudMPI.dir/src/cloudmpi.cpp.o
libCloudMPI.a: CMakeFiles/CloudMPI.dir/build.make
libCloudMPI.a: CMakeFiles/CloudMPI.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/prashant/CloudMPI/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libCloudMPI.a"
	$(CMAKE_COMMAND) -P CMakeFiles/CloudMPI.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CloudMPI.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CloudMPI.dir/build: libCloudMPI.a

.PHONY : CMakeFiles/CloudMPI.dir/build

CMakeFiles/CloudMPI.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CloudMPI.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CloudMPI.dir/clean

CMakeFiles/CloudMPI.dir/depend:
	cd /home/prashant/CloudMPI/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/prashant/CloudMPI /home/prashant/CloudMPI /home/prashant/CloudMPI/cmake-build-release /home/prashant/CloudMPI/cmake-build-release /home/prashant/CloudMPI/cmake-build-release/CMakeFiles/CloudMPI.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CloudMPI.dir/depend

