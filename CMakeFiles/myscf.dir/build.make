# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.9.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.9.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/tarinihardikar/Desktop/Dartmouth/myscf

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/tarinihardikar/Desktop/Dartmouth/myscf

# Include any dependencies generated for this target.
include CMakeFiles/myscf.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/myscf.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/myscf.dir/flags.make

CMakeFiles/myscf.dir/plugin.cc.o: CMakeFiles/myscf.dir/flags.make
CMakeFiles/myscf.dir/plugin.cc.o: plugin.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/tarinihardikar/Desktop/Dartmouth/myscf/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/myscf.dir/plugin.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/myscf.dir/plugin.cc.o -c /Users/tarinihardikar/Desktop/Dartmouth/myscf/plugin.cc

CMakeFiles/myscf.dir/plugin.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/myscf.dir/plugin.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/tarinihardikar/Desktop/Dartmouth/myscf/plugin.cc > CMakeFiles/myscf.dir/plugin.cc.i

CMakeFiles/myscf.dir/plugin.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/myscf.dir/plugin.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/tarinihardikar/Desktop/Dartmouth/myscf/plugin.cc -o CMakeFiles/myscf.dir/plugin.cc.s

CMakeFiles/myscf.dir/plugin.cc.o.requires:

.PHONY : CMakeFiles/myscf.dir/plugin.cc.o.requires

CMakeFiles/myscf.dir/plugin.cc.o.provides: CMakeFiles/myscf.dir/plugin.cc.o.requires
	$(MAKE) -f CMakeFiles/myscf.dir/build.make CMakeFiles/myscf.dir/plugin.cc.o.provides.build
.PHONY : CMakeFiles/myscf.dir/plugin.cc.o.provides

CMakeFiles/myscf.dir/plugin.cc.o.provides.build: CMakeFiles/myscf.dir/plugin.cc.o


# Object files for target myscf
myscf_OBJECTS = \
"CMakeFiles/myscf.dir/plugin.cc.o"

# External object files for target myscf
myscf_EXTERNAL_OBJECTS =

myscf.so: CMakeFiles/myscf.dir/plugin.cc.o
myscf.so: CMakeFiles/myscf.dir/build.make
myscf.so: /anaconda/envs/psi4env/lib/python2.7/site-packages/psi4/core.so
myscf.so: CMakeFiles/myscf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/tarinihardikar/Desktop/Dartmouth/myscf/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module myscf.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/myscf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/myscf.dir/build: myscf.so

.PHONY : CMakeFiles/myscf.dir/build

CMakeFiles/myscf.dir/requires: CMakeFiles/myscf.dir/plugin.cc.o.requires

.PHONY : CMakeFiles/myscf.dir/requires

CMakeFiles/myscf.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/myscf.dir/cmake_clean.cmake
.PHONY : CMakeFiles/myscf.dir/clean

CMakeFiles/myscf.dir/depend:
	cd /Users/tarinihardikar/Desktop/Dartmouth/myscf && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/tarinihardikar/Desktop/Dartmouth/myscf /Users/tarinihardikar/Desktop/Dartmouth/myscf /Users/tarinihardikar/Desktop/Dartmouth/myscf /Users/tarinihardikar/Desktop/Dartmouth/myscf /Users/tarinihardikar/Desktop/Dartmouth/myscf/CMakeFiles/myscf.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/myscf.dir/depend
