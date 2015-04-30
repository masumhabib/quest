# Configures, makes and installs QMICAD
# Copyright (c) 2014 K M Masum Habib
# Created on: June 02, 2014
#


cmake_minimum_required(VERSION 2.8)

project(qmicad CXX)
# Version
set(QMICAD_MAJOR 0)
set(QMICAD_MINOR 10)
set(QMICAD_PATCH 0)

# Options
option(VERBOSE_MAKE "Verbose cmake for debugging" off)
option(LINK_STATIC "Link everything statically" off)

# For custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Add machine specific local usr path for locally compiled libraries
#include(cmake/machines/fir)

# optimization
set(CMAKE_CXX_FLAGS_RELEASE " -m64 -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
# Linker flags
#set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

# Start configuration
set(QMICAD_VERSION ${QMICAD_MAJOR}.${QMICAD_MINOR}.${QMICAD_PATCH})
message(STATUS "Configuring QMICAD ${QMICAD_VERSION}")
message(STATUS "System: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "--------------------------------------------------------------")

# Ensure that we have a C++11 compiler.
include(cmake/CXX11.cmake)
check_for_cxx11_compiler(HAS_CXX11)
if(NOT HAS_CXX11)
  message(FATAL_ERROR "No C++11 compiler available!")
endif(NOT HAS_CXX11)
enable_cxx11()


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

# Libraries
# -----------------------------------------------------------------------------
if(LINK_STATIC)
    set(LINK_SEARCH_START_STATIC    on)
endif()

# Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# MPI mpich2
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Python libs
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# NumPy
find_package(NumPy REQUIRED)
include_directories(${NUMPY_INCLUDE_DIRS})

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
find_package(Boost 1.55 REQUIRED serialization random mpi mpi_python)
include_directories(${Boost_INCLUDE_DIR})

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    set(QMICAD_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif(${CMAKE_BUILD_TYPE} MATCHES "Release")

if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(QMICAD_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

message(STATUS "Compiler flags: ${QMICAD_FLAGS}")
message(STATUS "Linker flags: ${CMAKE_SHARED_LINKER_FLAGS}")
message(STATUS "Installation location: ${CMAKE_INSTALL_PREFIX}/${QMICAD_INSTALL_DIR}")

# Prepare distribution directory
set(QMICAD_DIST_DIR ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-${QMICAD_VERSION})
set(LIBRARY_OUTPUT_PATH ${QMICAD_DIST_DIR}/lib)

# installation
install(DIRECTORY ${QMICAD_DIST_DIR} DESTINATION .)

# QMICAD C++ Library
add_subdirectory(lib)
add_subdirectory(pylib)
add_subdirectory(simulators/cpp/tmfsc)


