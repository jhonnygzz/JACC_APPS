# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/xa2/Fall2024/miniBUDE

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/xa2/Fall2024/miniBUDE/buildKokkosOMP

# Utility rule file for NightlyCoverage.

# Include any custom commands dependencies for this target.
include kokkos/CMakeFiles/NightlyCoverage.dir/compiler_depend.make

# Include the progress variables for this target.
include kokkos/CMakeFiles/NightlyCoverage.dir/progress.make

kokkos/CMakeFiles/NightlyCoverage:
	cd /home/xa2/Fall2024/miniBUDE/buildKokkosOMP/kokkos && /usr/bin/ctest -D NightlyCoverage

NightlyCoverage: kokkos/CMakeFiles/NightlyCoverage
NightlyCoverage: kokkos/CMakeFiles/NightlyCoverage.dir/build.make
.PHONY : NightlyCoverage

# Rule to build all files generated by this target.
kokkos/CMakeFiles/NightlyCoverage.dir/build: NightlyCoverage
.PHONY : kokkos/CMakeFiles/NightlyCoverage.dir/build

kokkos/CMakeFiles/NightlyCoverage.dir/clean:
	cd /home/xa2/Fall2024/miniBUDE/buildKokkosOMP/kokkos && $(CMAKE_COMMAND) -P CMakeFiles/NightlyCoverage.dir/cmake_clean.cmake
.PHONY : kokkos/CMakeFiles/NightlyCoverage.dir/clean

kokkos/CMakeFiles/NightlyCoverage.dir/depend:
	cd /home/xa2/Fall2024/miniBUDE/buildKokkosOMP && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xa2/Fall2024/miniBUDE /home/xa2/Fall2024/miniBUDE/kokkosgithub /home/xa2/Fall2024/miniBUDE/buildKokkosOMP /home/xa2/Fall2024/miniBUDE/buildKokkosOMP/kokkos /home/xa2/Fall2024/miniBUDE/buildKokkosOMP/kokkos/CMakeFiles/NightlyCoverage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : kokkos/CMakeFiles/NightlyCoverage.dir/depend

