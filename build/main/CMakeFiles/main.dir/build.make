# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /data/opt/brew/Cellar/cmake/3.10.1/bin/cmake

# The command to remove a file.
RM = /data/opt/brew/Cellar/cmake/3.10.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /data/hzchen/DNAPregel

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /data/hzchen/DNAPregel/build

# Include any dependencies generated for this target.
include main/CMakeFiles/main.dir/depend.make

# Include the progress variables for this target.
include main/CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include main/CMakeFiles/main.dir/flags.make

main/CMakeFiles/main.dir/main.cpp.o: main/CMakeFiles/main.dir/flags.make
main/CMakeFiles/main.dir/main.cpp.o: ../main/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/hzchen/DNAPregel/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object main/CMakeFiles/main.dir/main.cpp.o"
	cd /data/hzchen/DNAPregel/build/main && /data/opt/brew/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/main.cpp.o -c /data/hzchen/DNAPregel/main/main.cpp

main/CMakeFiles/main.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cpp.i"
	cd /data/hzchen/DNAPregel/build/main && /data/opt/brew/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/hzchen/DNAPregel/main/main.cpp > CMakeFiles/main.dir/main.cpp.i

main/CMakeFiles/main.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cpp.s"
	cd /data/hzchen/DNAPregel/build/main && /data/opt/brew/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/hzchen/DNAPregel/main/main.cpp -o CMakeFiles/main.dir/main.cpp.s

main/CMakeFiles/main.dir/main.cpp.o.requires:

.PHONY : main/CMakeFiles/main.dir/main.cpp.o.requires

main/CMakeFiles/main.dir/main.cpp.o.provides: main/CMakeFiles/main.dir/main.cpp.o.requires
	$(MAKE) -f main/CMakeFiles/main.dir/build.make main/CMakeFiles/main.dir/main.cpp.o.provides.build
.PHONY : main/CMakeFiles/main.dir/main.cpp.o.provides

main/CMakeFiles/main.dir/main.cpp.o.provides.build: main/CMakeFiles/main.dir/main.cpp.o


# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/main.cpp.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

../release/main: main/CMakeFiles/main.dir/main.cpp.o
../release/main: main/CMakeFiles/main.dir/build.make
../release/main: utils/iniparser/libDNAPregel-iniparser.a
../release/main: main/CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/data/hzchen/DNAPregel/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../release/main"
	cd /data/hzchen/DNAPregel/build/main && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
main/CMakeFiles/main.dir/build: ../release/main

.PHONY : main/CMakeFiles/main.dir/build

main/CMakeFiles/main.dir/requires: main/CMakeFiles/main.dir/main.cpp.o.requires

.PHONY : main/CMakeFiles/main.dir/requires

main/CMakeFiles/main.dir/clean:
	cd /data/hzchen/DNAPregel/build/main && $(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : main/CMakeFiles/main.dir/clean

main/CMakeFiles/main.dir/depend:
	cd /data/hzchen/DNAPregel/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/hzchen/DNAPregel /data/hzchen/DNAPregel/main /data/hzchen/DNAPregel/build /data/hzchen/DNAPregel/build/main /data/hzchen/DNAPregel/build/main/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : main/CMakeFiles/main.dir/depend

