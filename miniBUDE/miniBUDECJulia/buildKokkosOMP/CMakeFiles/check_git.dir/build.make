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

# Utility rule file for check_git.

# Include any custom commands dependencies for this target.
include CMakeFiles/check_git.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/check_git.dir/progress.make

CMakeFiles/check_git: ../src/meta_vcs.h.in
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/xa2/Fall2024/miniBUDE/buildKokkosOMP/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Checking the git repository for changes..."
	/usr/bin/cmake -D_BUILD_TIME_CHECK_GIT=TRUE -DGIT_WORKING_DIR=/home/xa2/Fall2024/miniBUDE -DGIT_EXECUTABLE=/usr/bin/git -DGIT_STATE_FILE=/home/xa2/Fall2024/miniBUDE/buildKokkosOMP/git-state-hash -DPRE_CONFIGURE_FILE=/home/xa2/Fall2024/miniBUDE/src/meta_vcs.h.in -DPOST_CONFIGURE_FILE=/home/xa2/Fall2024/miniBUDE/buildKokkosOMP/generated/meta_vcs.h -DGIT_FAIL_IF_NONZERO_EXIT=FALSE -DGIT_IGNORE_UNTRACKED=FALSE -P /home/xa2/Fall2024/miniBUDE/cmake/git_watcher.cmake

check_git: CMakeFiles/check_git
check_git: CMakeFiles/check_git.dir/build.make
.PHONY : check_git

# Rule to build all files generated by this target.
CMakeFiles/check_git.dir/build: check_git
.PHONY : CMakeFiles/check_git.dir/build

CMakeFiles/check_git.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/check_git.dir/cmake_clean.cmake
.PHONY : CMakeFiles/check_git.dir/clean

CMakeFiles/check_git.dir/depend:
	cd /home/xa2/Fall2024/miniBUDE/buildKokkosOMP && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/xa2/Fall2024/miniBUDE /home/xa2/Fall2024/miniBUDE /home/xa2/Fall2024/miniBUDE/buildKokkosOMP /home/xa2/Fall2024/miniBUDE/buildKokkosOMP /home/xa2/Fall2024/miniBUDE/buildKokkosOMP/CMakeFiles/check_git.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/check_git.dir/depend
