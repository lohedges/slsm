CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
IF(COMMAND cmake_policy)
    CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

# Check for build type argument.
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
ENDIF()

PROJECT(LSM)

# Output some useful info.
SITE_NAME(COMPUTER_NAME)
MESSAGE(STATUS "Computer name: " \"${COMPUTER_NAME}\")
MESSAGE(STATUS "Build type: " \"${CMAKE_BUILD_TYPE}\")

# Set paths.
SET(EXECUTABLE_OUTPUT_PATH ${LSM_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for executables.")
SET(LIBRARY_OUTPUT_PATH ${LSM_BINARY_DIR}/lib CACHE INTERNAL "Single output directory for libraries.")

# Set the compiler flags.
IF("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    SET(CMAKE_CXX_FLAGS "-O0 -std=c++11 -g -Wall")
ELSEIF("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    SET(CMAKE_CXX_FLAGS "-O3 -std=c++11 -funroll-loops -Wall")
ENDIF()

# Print compiler flags.
MESSAGE(STATUS "Compiler flags: " \"${CMAKE_CXX_FLAGS}\")

# Make sure Boost is installed.
FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${Boost_INCLUDE_DIRS})

# Search for Doxygen, add dox subdirectory if found.
# CMakeLists.txt in dox directory adds documentation dependencies and doc make target.
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
    OPTION(ENABLE_DOXYGEN "Enables building of documentation with doxygen" 1)
    IF(ENABLE_DOXYGEN)
        ADD_SUBDIRECTORY(dox)
    ENDIF(ENABLE_DOXYGEN)
ENDIF(DOXYGEN_FOUND)

# Make lsm.h library header file.
FILE(GLOB HEADER_FILES "${LSM_SOURCE_DIR}/src/*.h")
STRING(REPLACE "${LSM_SOURCE_DIR}/src/lsm.h;" "" HEADER_FILES "${HEADER_FILES}")
STRING(REPLACE "${LSM_SOURCE_DIR}/src/" "" HEADER_FILES "${HEADER_FILES}")
FILE(WRITE ${LSM_SOURCE_DIR}/src/lsm.h "#ifndef _LSM_H\n#define _LSM_H\n\n")
FOREACH(_FILENAME ${HEADER_FILES})
    FILE(APPEND ${LSM_SOURCE_DIR}/src/lsm.h "#include \"${_FILENAME}\"\n")
ENDFOREACH(_FILENAME ${HEADER_FILES})
FILE(APPEND ${LSM_SOURCE_DIR}/src/lsm.h "\n#endif\n")

# Update include directories.
INCLUDE_DIRECTORIES(src)
AUX_SOURCE_DIRECTORY(src LSM_SRC)

# Create lsm library.
ADD_LIBRARY(lsm
    ${LSM_SRC}
)

# Install.
FILE(GLOB _FILES "${LSM_SOURCE_DIR}/src/*.h")
INSTALL(FILES ${_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/lsm)
INSTALL(TARGETS lsm
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

# Uninstall target
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake)

# Set the test code path.
FILE(GLOB TESTS RELATIVE ${LSM_SOURCE_DIR} ${LSM_SOURCE_DIR}/tests/*.cpp)

# Build tests.
FOREACH(TEST ${TESTS})
    STRING(REPLACE ".cpp" "" NAME ${TEST})
    STRING(REPLACE "tests/" "" NAME ${NAME})
    MESSAGE(STATUS "Found test: " ${NAME})
    ADD_EXECUTABLE(${NAME} ${TEST})
    TARGET_LINK_LIBRARIES(${NAME} lsm)
ENDFOREACH(TEST ${TESTS})

# Set the demo code path.
FILE(GLOB DEMOS RELATIVE ${LSM_SOURCE_DIR} ${LSM_SOURCE_DIR}/demos/*.cpp)

# Build demos.
FOREACH(DEMO ${DEMOS})
    STRING(REPLACE ".cpp" "" NAME ${DEMO})
    STRING(REPLACE "demos/" "" NAME ${NAME})
    MESSAGE(STATUS "Found demo: " ${NAME})
    ADD_EXECUTABLE(${NAME} ${DEMO})
    TARGET_LINK_LIBRARIES(${NAME} lsm)
    INSTALL(TARGETS ${NAME} RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/share/lsm-demos)
ENDFOREACH(DEMO ${DEMOS})
