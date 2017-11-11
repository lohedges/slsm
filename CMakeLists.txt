CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
IF(COMMAND cmake_policy)
    CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

# Check for build type argument.
IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel." FORCE)
ELSE()
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel." FORCE)
ENDIF()

# Allow user to enable/disable parallelisation.
OPTION(PARALLEL "Enable parallelisation." ON)

# Get Git commit and branch information.
EXECUTE_PROCESS(
    COMMAND
    git describe --abbrev=4 --dirty --always --tags
    OUTPUT_VARIABLE GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
EXECUTE_PROCESS(
    COMMAND
    git rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE
    GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Make sure the submodules have been initialised.
IF(NOT EXISTS ${CMAKE_SOURCE_DIR}/external)
    MESSAGE("Unable to find submodules.")

    # Clones submodules.
    EXECUTE_PROCESS(
        COMMAND git submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
ENDIF()

# Fix NLopt cmake_minimum_required.
IF(EXISTS ${CMAKE_SOURCE_DIR}/external/nlopt/CMakeLists.txt)
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        EXECUTE_PROCESS(
            COMMAND sed -i "" "s#VERSION 3.0#VERSION 2.8#g" ${CMAKE_SOURCE_DIR}/external/nlopt/CMakeLists.txt
        )
    ELSE()
        EXECUTE_PROCESS(
            COMMAND sed -i "s#VERSION 3.0#VERSION 2.8#g" ${CMAKE_SOURCE_DIR}/external/nlopt/CMakeLists.txt
        )
    ENDIF()
ENDIF()

# Make sure demos folder exists in build directory.
if (NOT EXISTS ${CMAKE_BINARY_DIR}/demos)
    MESSAGE("Making directory for demonstration programs.")
    EXECUTE_PROCESS(
        COMMAND
        mkdir ${CMAKE_BINARY_DIR}/demos
    )
ENDIF()

# Make sure tests folder exists in build directory.
IF(NOT EXISTS ${CMAKE_BINARY_DIR}/tests)
    MESSAGE("Making directory for test programs.")
    EXECUTE_PROCESS(
        COMMAND
        mkdir ${CMAKE_BINARY_DIR}/tests
    )
ENDIF()

# Copy demo input files to build directory.
IF(NOT ${CMAKE_SOURCE_DIR} EQUAL ${CMAKE_BINARY_DIR})
    MESSAGE("Copying files to build directory.")
    EXECUTE_PROCESS(
        COMMAND
        cp -r ${CMAKE_SOURCE_DIR}/demos/shapes ${CMAKE_BINARY_DIR}/demos/shapes
    )
    EXECUTE_PROCESS(
        COMMAND
        cp ${CMAKE_SOURCE_DIR}/tests/runtests ${CMAKE_BINARY_DIR}/tests/
    )
ENDIF()

PROJECT(SLSM)

# Set library path.
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib CACHE INTERNAL "Single output directory for libraries.")

# Set the compiler flags.
IF("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    SET(CMAKE_CXX_FLAGS "-O0 -g -Wall -std=c++11")
ELSEIF("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    SET(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -std=c++11")
ENDIF()

# Add Git information.
ADD_DEFINITIONS(-DCOMMIT="${GIT_COMMIT}")
ADD_DEFINITIONS(-DBRANCH="${GIT_BRANCH}")

# Disable parallelistion.
IF(PARALLEL)
    ADD_DEFINITIONS(-DPARALLEL)
ENDIF()

# Output some useful info.
SITE_NAME(COMPUTER_NAME)
MESSAGE(STATUS "Computer name: " \"${COMPUTER_NAME}\")
MESSAGE(STATUS "Build type: " \"${CMAKE_BUILD_TYPE}\")
MESSAGE(STATUS "Compiler flags: " \"${CMAKE_CXX_FLAGS}\")

# Add NLOpt library.
SET(NLOPT_PYTHON OFF CACHE BOOL "Build NLopt Python bindings" FORCE)
SET(NLOPT_OCTAVE OFF CACHE BOOL "Build NLopt Octave bindings" FORCE)
SET(NLOPT_MATLAB OFF CACHE BOOL "Build NLopt Matlab bindings" FORCE)
SET(NLOPT_GUILE  OFF CACHE BOOL "Build NLopt Guile bindings" FORCE)
SET(NLOPT_SWIG   OFF CACHE BOOL "Use SWIG to build bindings" FORCE)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/external/nlopt)

# Add Pybind11.
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/external/pybind11)

# Search for Doxygen, add dox subdirectory if found.
# CMakeLists.txt in dox directory adds documentation dependencies and doc make target.
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
    OPTION(ENABLE_DOXYGEN "Enables building of documentation with doxygen" 1)
    IF(ENABLE_DOXYGEN)
        ADD_SUBDIRECTORY(dox)
    ENDIF(ENABLE_DOXYGEN)
ENDIF(DOXYGEN_FOUND)

# Make slsm.h library header file.
FILE(GLOB HEADER_FILES "${CMAKE_SOURCE_DIR}/src/*.h")
STRING(REPLACE "${CMAKE_SOURCE_DIR}/src/slsm.h;" "" HEADER_FILES "${HEADER_FILES}")
STRING(REPLACE "${CMAKE_SOURCE_DIR}/src/" "" HEADER_FILES "${HEADER_FILES}")
FILE(WRITE ${CMAKE_SOURCE_DIR}/src/slsm.h "#ifndef _SLSM_H\n#define _SLSM_H\n\n")
FOREACH(_FILENAME ${HEADER_FILES})
    FILE(APPEND ${CMAKE_SOURCE_DIR}/src/slsm.h "#include \"${_FILENAME}\"\n")
ENDFOREACH(_FILENAME ${HEADER_FILES})
FILE(APPEND ${CMAKE_SOURCE_DIR}/src/slsm.h "\n#endif\n")

# Update include directories.
INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/external/nlopt/api
)
AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/src SLSM_SRC)

# Create slsm library.
ADD_LIBRARY(slsm
    ${SLSM_SRC}
)

# Library should be lined against NLopt.
TARGET_LINK_LIBRARIES(slsm nlopt)

# Install.
FILE(GLOB _FILES "${CMAKE_SOURCE_DIR}/src/*.h")
INSTALL(FILES ${_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/slsm)
INSTALL(TARGETS slsm
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

# Uninstall target
CONFIGURE_FILE(
    "${CMAKE_SOURCE_DIR}/uninstall.cmake.in"
    "${CMAKE_SOURCE_DIR}/uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/uninstall.cmake)

# Generate a list of test source files.
FILE(GLOB TESTS RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/tests/*.cpp)

# Build tests.
FOREACH(TEST ${TESTS})
    STRING(REPLACE ".cpp" "" NAME ${TEST})
    STRING(REPLACE "tests/" "" NAME ${NAME})
    MESSAGE(STATUS "Found test: " ${NAME})
    ADD_EXECUTABLE(${NAME} ${TEST})
    TARGET_LINK_LIBRARIES(${NAME} slsm nlopt)
    SET_TARGET_PROPERTIES(${NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
    )
ENDFOREACH(TEST ${TESTS})

# Generate a list of demo source files.
FILE(GLOB DEMOS RELATIVE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/demos/*.cpp)

# Build demos.
FOREACH(DEMO ${DEMOS})
    STRING(REPLACE ".cpp" "" NAME ${DEMO})
    STRING(REPLACE "demos/" "" NAME ${NAME})
    MESSAGE(STATUS "Found demo: " ${NAME})
    ADD_EXECUTABLE(${NAME} ${DEMO})
    TARGET_LINK_LIBRARIES(${NAME} slsm nlopt)
    INSTALL(TARGETS ${NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/share/slsm-demo
    )
    SET_TARGET_PROPERTIES(${NAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/demos"
    )
ENDFOREACH(DEMO ${DEMOS})

# Build Python bindings.
PYBIND11_ADD_MODULE(
	pyslsm
    ${CMAKE_SOURCE_DIR}/python/bindings/pyslsm.cpp
    ${CMAKE_SOURCE_DIR}/python/bindings/bind_Boundary.cpp
    ${CMAKE_SOURCE_DIR}/python/bindings/bind_FastMarchingMethod.cpp
    ${CMAKE_SOURCE_DIR}/python/bindings/bind_Hole.cpp
    ${CMAKE_SOURCE_DIR}/python/bindings/bind_InputOutput.cpp
    ${CMAKE_SOURCE_DIR}/python/bindings/bind_LevelSet.cpp
    ${CMAKE_SOURCE_DIR}/python/bindings/bind_MersenneTwister.cpp
    ${CMAKE_SOURCE_DIR}/python/bindings/bind_Mesh.cpp
    ${CMAKE_SOURCE_DIR}/python/bindings/bind_Optimise.cpp
    ${CMAKE_SOURCE_DIR}/python/bindings/bind_Sensitivity.cpp
)

# Specify the path for python shared library.
SET_TARGET_PROPERTIES(
    pyslsm
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/python"
	COMPILE_DEFINITIONS PYBIND
)

# Link against NLopt.
TARGET_LINK_LIBRARIES(pyslsm PUBLIC nlopt)
