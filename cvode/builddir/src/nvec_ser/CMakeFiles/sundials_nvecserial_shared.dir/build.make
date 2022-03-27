# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

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
CMAKE_COMMAND = /usr/bin/cmake.exe

# The command to remove a file.
RM = /usr/bin/cmake.exe -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/MM-PIHM/cvode

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/MM-PIHM/cvode/builddir

# Include any dependencies generated for this target.
include src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/depend.make

# Include the progress variables for this target.
include src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/progress.make

# Include the compile flags for this target's objects.
include src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/flags.make

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/flags.make
src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o: ../src/nvec_ser/nvector_serial.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/MM-PIHM/cvode/builddir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o"
	cd /home/MM-PIHM/cvode/builddir/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o   -c /home/MM-PIHM/cvode/src/nvec_ser/nvector_serial.c

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.i"
	cd /home/MM-PIHM/cvode/builddir/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/MM-PIHM/cvode/src/nvec_ser/nvector_serial.c > CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.i

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.s"
	cd /home/MM-PIHM/cvode/builddir/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/MM-PIHM/cvode/src/nvec_ser/nvector_serial.c -o CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.s

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o.requires:

.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o.requires

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o.provides: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o.requires
	$(MAKE) -f src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/build.make src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o.provides.build
.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o.provides

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o.provides.build: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o


src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/flags.make
src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o: ../src/sundials/sundials_math.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/MM-PIHM/cvode/builddir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o"
	cd /home/MM-PIHM/cvode/builddir/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o   -c /home/MM-PIHM/cvode/src/sundials/sundials_math.c

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.i"
	cd /home/MM-PIHM/cvode/builddir/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/MM-PIHM/cvode/src/sundials/sundials_math.c > CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.i

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.s"
	cd /home/MM-PIHM/cvode/builddir/src/nvec_ser && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/MM-PIHM/cvode/src/sundials/sundials_math.c -o CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.s

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o.requires:

.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o.requires

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o.provides: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o.requires
	$(MAKE) -f src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/build.make src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o.provides.build
.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o.provides

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o.provides.build: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o


# Object files for target sundials_nvecserial_shared
sundials_nvecserial_shared_OBJECTS = \
"CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o" \
"CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o"

# External object files for target sundials_nvecserial_shared
sundials_nvecserial_shared_EXTERNAL_OBJECTS =

src/nvec_ser/cygsundials_nvecserial-2.dll: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o
src/nvec_ser/cygsundials_nvecserial-2.dll: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o
src/nvec_ser/cygsundials_nvecserial-2.dll: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/build.make
src/nvec_ser/cygsundials_nvecserial-2.dll: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/MM-PIHM/cvode/builddir/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C shared library cygsundials_nvecserial-2.dll"
	cd /home/MM-PIHM/cvode/builddir/src/nvec_ser && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sundials_nvecserial_shared.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/build: src/nvec_ser/cygsundials_nvecserial-2.dll

.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/build

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/requires: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/nvector_serial.c.o.requires
src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/requires: src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/__/sundials/sundials_math.c.o.requires

.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/requires

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/clean:
	cd /home/MM-PIHM/cvode/builddir/src/nvec_ser && $(CMAKE_COMMAND) -P CMakeFiles/sundials_nvecserial_shared.dir/cmake_clean.cmake
.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/clean

src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/depend:
	cd /home/MM-PIHM/cvode/builddir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/MM-PIHM/cvode /home/MM-PIHM/cvode/src/nvec_ser /home/MM-PIHM/cvode/builddir /home/MM-PIHM/cvode/builddir/src/nvec_ser /home/MM-PIHM/cvode/builddir/src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/nvec_ser/CMakeFiles/sundials_nvecserial_shared.dir/depend

