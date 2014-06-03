# Configures, makes and installs QMICAD
# Copyright (c) 2014 K M Masum Habib
# Created on: June 02, 2014
#


cmake_minimum_required(VERSION 2.8)

# Version
set(QMICAD_MAJOR 0)
set(QMICAD_MINOR 05)
set(QMICAD_PATCH 0)


# Add local usr path for locally compiled libraries
set(CMAKE_PREFIX_PATH $ENV{HOME}/usr)
set(Boost_INCLUDE_DIR $ENV{HOME}/usr/include)
set(Boost_LIBRARY_DIR $ENV{HOME}/usr/lib)

# Start configuration
set(QMICAD_VERSION ${QMICAD_MAJOR}.${QMICAD_MINOR}.${QMICAD_PATCH})
message(STATUS "Configuring QMICAD ${QMICAD_VERSION}")
project(qmicad)

# Options
# -----------------------------------------------------------------------------
set(CMAKE_VERBOSE_MAKEFILE      off)
set(Boost_USE_STATIC_LIBS       on)
set(LINK_SEARCH_START_STATIC    on)

# For custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE, Release)
    message(STATUS "Selecting Release build type by default. You might change it by using -DCMAKE_BUILD_TYPE=Release or -DCMAKE_BUILD_TYPE=Debug on your cmake command line.")
endif ()

# Compiler flags
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Libraries
# -----------------------------------------------------------------------------
# Armadillo
find_package(Armadillo REQUIRED)

# Boost
find_package(Boost 1.55 REQUIRED serialization python mpi mpi_python)

# MPI mpich2
find_package(MPI REQUIRED)

# Python libs
find_package(PythonLibs REQUIRED)

# MKL
find_package(MKL REQUIRED)

# Sources and headers
# -----------------------------------------------------------------------------
file(GLOB_RECURSE QMICAD_SOURCES "src/*.cpp")
file(GLOB_RECURSE QMICAD_HEADERS "include/*.h*")
set (QMICAD_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)

#set (QMICAD_INCLUDE_DIRS "")
#foreach (_headerFile ${QMICAD_HEADERS})
#    get_filename_component(_dir ${_headerFile} PATH)
#    list (APPEND HEADER_INCLUDE_DIRS ${_dir})
#endforeach()
#list(REMOVE_DUPLICATES QMICAD_INCLUDE_DIRS)

# includes
include_directories(${QMICAD_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})

# linking
add_library(qmicad SHARED ${QMICAD_SOURCES})
target_link_libraries (qmicad ${Boost_SERIALIZATION_LIBRARIES})
target_link_libraries (qmicad ${Boost_PYTHON_LIBRARIES})
target_link_libraries (qmicad ${Boost_MPI_LIBRARIES})
target_link_libraries (qmicad ${Boost_MPI_PYTHON_LIBRARIES})
target_link_libraries (qmicad ${MPI_LIBRARIES})
target_link_libraries (qmicad ${MKL_LIBRARIES})

#install(TARGETS qmicad DESTINATION lib/qmicad)
#install(FILES test.h DESTINATION include)
