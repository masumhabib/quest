# Configures, makes and installs QUEST
# Copyright (c) 2014 K M Masum Habib
# Created on: June 02, 2014
#


cmake_minimum_required(VERSION 2.8)

project(quest CXX)

# Options
option(VERBOSE_MAKE "Verbose cmake for debugging" off)
option(LINK_STATIC "Link everything statically" off)
#set(CMAKE_BUILD_TYPE "Debug")

# To use locally compiled libraries (Armadillo and Boost), create a file
# that contains the local library paths and include it here. Some examples of
# such files can be found in cmake/local folder.
#include(cmake/local/Armadillo.cmake)

# Build and run testing
enable_testing()

# For custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")


include(cmake/Versioning.cmake)
get_version (QUEST_MAJOR QUEST_MINOR QUEST_PATCH QUEST_BUILD)

# Start configuration
set(QUEST_VERSION v${QUEST_MAJOR}.${QUEST_MINOR}.${QUEST_PATCH}${QUEST_BUILD})
set(QUEST_SHORT_VERSION ${QUEST_MAJOR}.${QUEST_MINOR}.${QUEST_PATCH})
message(STATUS "--------------------------------------------------------------")
message(STATUS "Configuring QUEST ${QUEST_VERSION}")
message(STATUS "System: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "--------------------------------------------------------------")

# Ensure that we have a C++11 compiler.
include(cmake/CXX11.cmake)
check_for_cxx11_compiler(HAS_CXX11)
if(NOT HAS_CXX11)
  message(FATAL_ERROR "No C++11 compiler available!")
endif(NOT HAS_CXX11)
enable_cxx11()

# C++ compiler flags
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Werror=return-type")

if (VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE      on)
endif()

#Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# Libraries
# -----------------------------------------------------------------------------
if(LINK_STATIC)
    set(LINK_SEARCH_START_STATIC    on)
endif()

# Armadillo
find_package(Armadillo REQUIRED)

# MPI mpich2
find_package(MPI REQUIRED)

# Python libs
find_package(PythonInterp 2.7 REQUIRED)
find_package(PythonLibs 2.7 REQUIRED)

# NumPy
include(cmake/FindNumPy.cmake)
find_package(NumPy REQUIRED)

# doxygen
find_package(Doxygen)

# Boost
# Shared linking of libboost_python is required for mpi::communicator
# to work both in python and C++.
set(Boost_USE_MULTITHREADED     off)
if(LINK_STATIC)
    set(Boost_USE_STATIC_LIBS       off)
endif()
find_package(Boost 1.55 REQUIRED python)
set(Boost_PYTHON_LIBRARIES_SHARED ${Boost_PYTHON_LIBRARIES})
unset(Boost_LIBRARIES)
unset(Boost_PYTHON_LIBRARIES)
if(LINK_STATIC)
    set(Boost_USE_STATIC_LIBS       on)
endif()
find_package(Boost 1.55 REQUIRED serialization random mpi mpi_python unit_test_framework)
include_directories(${Boost_INCLUDE_DIR})

# Prepare distribution directory
set(QUEST_BUILD_DIR ${PROJECT_BINARY_DIR}/..)
set(QUEST_DIST_DIR ${QUEST_BUILD_DIR}/${PROJECT_NAME}-${QUEST_SHORT_VERSION})
set(LIBRARY_OUTPUT_PATH ${QUEST_DIST_DIR}/lib)
set(QUEST_TEST_DIR ${QUEST_BUILD_DIR}/tests)

# installation
install(DIRECTORY ${QUEST_DIST_DIR} DESTINATION .)

add_subdirectory(lib)
add_subdirectory(pylib)
add_subdirectory(quester)
add_subdirectory(tests)

message(STATUS "--------------------------------------------------------------")

