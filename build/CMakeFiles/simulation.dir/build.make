# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 4.0

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
CMAKE_COMMAND = /opt/homebrew/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/yernurbaibolatov/Documents/Projects/NLDKit/build

# Include any dependencies generated for this target.
include CMakeFiles/simulation.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/simulation.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/simulation.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/simulation.dir/flags.make

CMakeFiles/simulation.dir/codegen:
.PHONY : CMakeFiles/simulation.dir/codegen

CMakeFiles/simulation.dir/src/main.cpp.o: CMakeFiles/simulation.dir/flags.make
CMakeFiles/simulation.dir/src/main.cpp.o: /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/main.cpp
CMakeFiles/simulation.dir/src/main.cpp.o: CMakeFiles/simulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/yernurbaibolatov/Documents/Projects/NLDKit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/simulation.dir/src/main.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/simulation.dir/src/main.cpp.o -MF CMakeFiles/simulation.dir/src/main.cpp.o.d -o CMakeFiles/simulation.dir/src/main.cpp.o -c /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/main.cpp

CMakeFiles/simulation.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/simulation.dir/src/main.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/main.cpp > CMakeFiles/simulation.dir/src/main.cpp.i

CMakeFiles/simulation.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/simulation.dir/src/main.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/main.cpp -o CMakeFiles/simulation.dir/src/main.cpp.s

CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.o: CMakeFiles/simulation.dir/flags.make
CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.o: /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/core/DynamicalSystem.cpp
CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.o: CMakeFiles/simulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/yernurbaibolatov/Documents/Projects/NLDKit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.o -MF CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.o.d -o CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.o -c /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/core/DynamicalSystem.cpp

CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/core/DynamicalSystem.cpp > CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.i

CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/core/DynamicalSystem.cpp -o CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.s

CMakeFiles/simulation.dir/src/core/Integrator.cpp.o: CMakeFiles/simulation.dir/flags.make
CMakeFiles/simulation.dir/src/core/Integrator.cpp.o: /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/core/Integrator.cpp
CMakeFiles/simulation.dir/src/core/Integrator.cpp.o: CMakeFiles/simulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/yernurbaibolatov/Documents/Projects/NLDKit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/simulation.dir/src/core/Integrator.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/simulation.dir/src/core/Integrator.cpp.o -MF CMakeFiles/simulation.dir/src/core/Integrator.cpp.o.d -o CMakeFiles/simulation.dir/src/core/Integrator.cpp.o -c /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/core/Integrator.cpp

CMakeFiles/simulation.dir/src/core/Integrator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/simulation.dir/src/core/Integrator.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/core/Integrator.cpp > CMakeFiles/simulation.dir/src/core/Integrator.cpp.i

CMakeFiles/simulation.dir/src/core/Integrator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/simulation.dir/src/core/Integrator.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp/src/core/Integrator.cpp -o CMakeFiles/simulation.dir/src/core/Integrator.cpp.s

# Object files for target simulation
simulation_OBJECTS = \
"CMakeFiles/simulation.dir/src/main.cpp.o" \
"CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.o" \
"CMakeFiles/simulation.dir/src/core/Integrator.cpp.o"

# External object files for target simulation
simulation_EXTERNAL_OBJECTS =

simulation: CMakeFiles/simulation.dir/src/main.cpp.o
simulation: CMakeFiles/simulation.dir/src/core/DynamicalSystem.cpp.o
simulation: CMakeFiles/simulation.dir/src/core/Integrator.cpp.o
simulation: CMakeFiles/simulation.dir/build.make
simulation: CMakeFiles/simulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/yernurbaibolatov/Documents/Projects/NLDKit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable simulation"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/simulation.dir/build: simulation
.PHONY : CMakeFiles/simulation.dir/build

CMakeFiles/simulation.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/simulation.dir/cmake_clean.cmake
.PHONY : CMakeFiles/simulation.dir/clean

CMakeFiles/simulation.dir/depend:
	cd /Users/yernurbaibolatov/Documents/Projects/NLDKit/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp /Users/yernurbaibolatov/Documents/Projects/NLDKit/cpp /Users/yernurbaibolatov/Documents/Projects/NLDKit/build /Users/yernurbaibolatov/Documents/Projects/NLDKit/build /Users/yernurbaibolatov/Documents/Projects/NLDKit/build/CMakeFiles/simulation.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/simulation.dir/depend

