# CXX Project Makefile

# Copyright (c) 2014, 2015  Lester Hedges <lester.hedges@gmail.com>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

################################ INFO #######################################

# This Makefile can be used to build a CXX project library along with its
# associated demos, unit tests, and documentation. For detailed information
# on using the Makefile run make without a target, i.e. simply run make at
# your command prompt.
#
# Makefile style adapted from http://clarkgrubb.com/make-file-style-guide
# Conventions:
#   - Environment and Makefile variables are in upper case, user
#     defined variables are in lower case.
#   - Variables are declared using the immediate assignment operator :=

############################### MACROS ########################################

define colorecho
	@tput setaf $1
	@echo $2
	@tput sgr0
endef

define boldcolorecho
	@tput bold
	@tput setaf $1
	@echo $2
	@tput sgr0
endef

define inlinecolorecho
	tput setaf $1; echo $2; tput sgr0
endef

############################## VARIABLES ######################################

# Set shell to bash.
SHELL := bash

# Suppress display of executed commands.
.SILENT:

# Default goal will print the help message.
.DEFAULT_GOAL := help

# Project name.
project := lsm

# Upper case project name (for use in library header file).
project_upper := `echo $(project) | tr a-z A-Z`

# C++ compiler.
CXX := g++

# Installation path.
PREFIX ?= /usr/local

# External libraries.
LIBS :=

# Path for source files.
src_dir := src

# Path for demo code.
demo_dir := demos

# Path for unit tests.
test_dir := tests

# Path for object files.
obj_dir := obj

# Path for the library.
lib_dir := lib

# Library header file.
library_header := $(src_dir)/$(project).h

# Generate library target name.
library := $(lib_dir)/lib$(project).a

# Install command.
install_cmd := install

# Install flags for executables.
iflags_exec := -m 0755

# Install flags for non-executable files.
iflags := -m 0644

# Git commit information.
commit := $(shell git describe --abbrev=4 --dirty --always --tags 2> /dev/null)

# C++ compiler flags for development build.
cxxflags_devel := -O0 -std=c++11 -g -Wall -Isrc -DCOMMIT=\"$(commit)\" $(OPTFLAGS)

# C++ compiler flags for release build.
cxxflags_release := -O3 -std=c++11 -funroll-loops -DNDEBUG -Isrc -DCOMMIT=\"$(commit)\" $(OPTFLAGS)

# Default to release build.
CXXFLAGS := $(cxxflags_release)

# The C++ header, source, object, and dependency files.
headers := $(wildcard $(src_dir)/*.h)
headers := $(filter-out $(library_header), $(headers))
sources := $(wildcard $(src_dir)/*.cpp)
temp := $(patsubst %.cpp,%.o,$(sources))
objects := $(subst $(src_dir),$(obj_dir),$(temp))
-include $(subst .o,.d,$(objects))

# Source files and executable names for demos.
demo_sources := $(wildcard $(demo_dir)/*.cpp)
demos := $(patsubst %.cpp,%,$(demo_sources))

# Source files and executable names for unit tests.
test_sources := $(wildcard $(test_dir)/*_tests.cpp)
tests := $(patsubst %.cpp,%,$(test_sources))

# Doxygen files.
dox_files := $(wildcard dox/*.dox)

# Check that unit tests exist.
is_tests := 1
ifeq ($(strip $(tests)),)
	is_tests := 0
endif

############################### TARGETS #######################################

# Print help message.
.PHONY: help
help:
	$(call boldcolorecho, 4, "About")
	@echo " This Makefile can be used to build the $(project) library along with its"
	@echo " associated demos, unit tests, and documentation."
	@echo
	@echo " If you would prefer to use CMake, e.g. if you aren't using make as your"
	@echo " generator, or are compiling on Windows, simply copy the contents of the"
	@echo " cmake directory into the root of the project folder then run cmake, i.e."
	@echo " (assuming that you are currently in the root directory)"
	@echo
	@echo " cp -r cmake/* ."
	@echo " cmake ."
	@echo
	@echo " Be warned that if you use make as your cmake generator then the original"
	@echo " Makefile will be overwritten. To recover it simply run: git co Makefile"
	@echo
	$(call boldcolorecho, 4, "Targets")
	@echo " help       -->  print this help message"
	@echo " build      -->  build library, demos, and tests (default=release)"
	@echo " devel      -->  build using development compiler flags (debug)"
	@echo " release    -->  build using release compiler flags (optmized)"
	@echo " doc        -->  generate source code documentation with doxygen"
	@echo " test       -->  run unit tests"
	@echo " valgrind   -->  run unit tests via valgrind"
	@echo " clean      -->  remove object and dependency files"
	@echo " clobber    -->  remove all files generated by make"
	@echo " install    -->  install library, demos, and documentation"
	@echo " uninstall  -->  uninstall library, demos, and documentation"
	@echo
	$(call boldcolorecho, 4, "Dependencies")
	@echo " Make sure that any external libraries are are accessible to your"
	@echo " shell environment. The LIBS variable should be used for any linker"
	@echo " flags, e.g. -lfftw3. Use the LDFLAGS environment variable for any"
	@echo " non-standard paths. This can be set in your shell profile, on the"
	@echo " command-line, or passed directly to make as an option. For example,"
	@echo " if using MacPorts on OS X"
	@echo
	@echo " export LDFLAGS='-L/opt/local/lib'"
	@echo
	$(call boldcolorecho, 4, "Tips")
	@echo " To set a different installation path run"
	@echo "     PREFIX=path make install"
	@echo
	@echo " Additional CXXFLAGS can be passed using OPTFLAGS, e.g."
	@echo "     OPTFLAGS=-Wall make devel"
	@echo
	@echo " Targets can be chained together, e.g."
	@echo "     make release doc test"

# Set development compilation flags and build.
devel: CXXFLAGS := $(cxxflags_devel)
devel: build

# Set release compilation flags and build.
release: CXXFLAGS := $(cxxflags_release)
release: build

# Print compiler flags.
devel release:
	$(call colorecho, 5, "--> CXXFLAGS: $(CXXFLAGS)")

# Save compiler flags to file if they differ from previous build.
# This target ensures that all object files are recompiled if the flags change.
.PHONY: force
.compiler_flags: force
	@echo '$(CXXFLAGS)' | cmp -s - $@ || echo '$(CXXFLAGS)' > $@

# Compile object files.
# Autodepenencies are handled using a recipe taken from
# http://scottmcpeak.com/autodepend/autodepend.html
$(obj_dir)/%.o: $(src_dir)/%.cpp .compiler_flags
	$(call colorecho, 2, "--> Building CXX object $*.o")
	$(CXX) $(CXXFLAGS) -c -o $(obj_dir)/$*.o $(src_dir)/$*.cpp
	$(CXX) -MM $(CXXFLAGS) $(src_dir)/$*.cpp > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$(obj_dir)/$*.o:|' < $*.d.tmp > $(obj_dir)/$*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
		sed -e 's/^ *//' -e 's/$$/:/' >> $(obj_dir)/$*.d
	@rm -f $*.d.tmp

# Build the library and demos.
.PHONY: build
build: $(obj_dir) $(library_header) $(library) $(demos) $(tests)

# Create output directory for object and dependency files.
$(obj_dir):
	mkdir $(obj_dir)

# Create library header file.
$(library_header): $(headers)
	$(call colorecho, 4, "--> Generating CXX library header $(library_header)")
	@echo -e "#ifndef _$(project_upper)_H\n#define _$(project_upper)_H\n" > $(library_header)
	@for i in $(headers);						\
		do h=`echo $$i | cut -d'/' -f 2`;		\
		echo "#include \"$$h\"";				\
	done | sort -g >> $(library_header)
	@echo -e "\n#endif" >> $(library_header)

# Build the static library.
$(library): $(objects)
	$(call colorecho, 1, "--> Linking CXX static library $(library)")
	mkdir -p $(lib_dir)
	ar rcs $@ $(objects)
	ranlib $@

# Compile demonstration code.
$(demos): %: %.cpp $(library_header) $(library)
	$(call colorecho, 1, "--> Linking CXX executable $@")
	$(CXX) $(CXXFLAGS) $@.cpp $(library) $(LIBS) $(LDFLAGS) -o $@

# Compile unit tests.
$(tests): %: %.cpp $(obj_dir) $(library_header) $(library)
	$(call colorecho, 1, "--> Linking CXX executable $@")
	$(CXX) $(CXXFLAGS) $@.cpp $(library) $(LIBS) $(LDFLAGS) -o $@

# Run unit tests
.PHONY: test
test: $(tests)
	@if [ $(is_tests) -eq 1 ]; then									\
		$(call inlinecolorecho, 6, "--> Running CXX unit tests");	\
		./tests/runtests;											\
	fi

# Compile unit tests and run via valgrind.
.PHONY: valgrind
valgrind:
	@VALGRIND="valgrind --log-file=/tmp/valgrind-%p.log" $(MAKE) test

# Build documentation using Doxygen.
doc: $(headers) $(sources) $(dox_files)
	$(call colorecho, 4, "--> Generating CXX source documentation with Doxygen")
	doxygen dox/Doxyfile

# Install the library and demos.
.PHONY: install
install: build doc
	$(call colorecho, 3, "--> Installing CXX static library $(library) to $(PREFIX)/lib")
	$(call colorecho, 3, "--> Installing CXX demos $(demos) to $(PREFIX)/share/$(project)-demos")
	$(call colorecho, 3, "--> Installing CXX Doxygen documentation to $(PREFIX)/share/doc/$(project)")
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/lib
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/include/$(project)
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/share/$(project)-demos
	$(install_cmd) -d $(iflags_exec) $(PREFIX)/share/doc/$(project)
	$(install_cmd) $(iflags) $(library) $(PREFIX)/lib
	$(install_cmd) $(iflags) $(headers) $(PREFIX)/include/$(project)
	$(install_cmd) $(iflags) $(library_header) $(PREFIX)/include/$(project)
	$(install_cmd) $(iflags) $(demo_sources) $(PREFIX)/share/$(project)-demos
	$(install_cmd) $(iflags_exec) $(demos) $(PREFIX)/share/$(project)-demos
	cp -r doc/html $(PREFIX)/share/doc/$(project)

# Uninstall the library and demos.
.PHONY: uninstall
uninstall:
	$(call colorecho, 3, "--> Uninstalling CXX static library $(library) from $(PREFIX)/lib")
	$(call colorecho, 3, "--> Uninstalling CXX demos $(demos) from $(PREFIX)/share/$(project)-demos")
	$(call colorecho, 3, "--> Uninstalling CXX Doxygen documentation from $(PREFIX)/share/doc/$(project)")
	rm -f $(PREFIX)/$(library)
	rm -rf $(PREFIX)/include/$(project)
	rm -rf $(PREFIX)/share/$(project)-demos
	rm -rf $(PREFIX)/share/doc/$(project)

# Clean up object and dependecy files.
.PHONY: clean
clean:
	$(call colorecho, 6, "--> Cleaning CXX object and dependency files")
	rm -rf $(obj_dir)

# Clean up everything produced by make.
.PHONY: clobber
clobber:
	$(call colorecho, 6, "--> Cleaning all output files")
	rm -rf $(obj_dir)
	rm -rf $(lib_dir)
	rm -rf doc
	rm -f $(demos)
	rm -f $(tests)
	rm -f $(library_header)
	rm -rf $(demo_dir)/*dSYM
	rm -rf $(test_dir)/*dSYM
	rm -f tests/tests.log
	rm -f .compiler_flags

.PHONY: sandwich
sandwich:
	if [ "$$(id -u)" != "0" ]; then                        \
		echo " What? Make it yourself."                   ;\
	else                                                   \
		echo "                      ____"                 ;\
		echo "          .----------'    '-."              ;\
		echo "         /  .      '     .   \\"            ;\
		echo "        /        '    .      /|"            ;\
		echo "       /      .             \ /"            ;\
		echo "      /  ' .       .     .  || |"           ;\
		echo "     /.___________    '    / //"            ;\
		echo "     |._          '------'| /|"             ;\
		echo "     '.............______.-' /"             ;\
		echo "     |-.                  | /"              ;\
		echo "     \`\"\"\"\"\"\"\"\"\"\"\"\"\"-.....-'"  ;\
	fi;                                                    \
