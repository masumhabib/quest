# Configures, makes and installs QMICAD
# Copyright (c) 2014 K M Masum Habib
# Created on: June 02, 2014
#


cmake_minimum_required(VERSION 2.8)
project(qmicad)

# Add local usr path for locally compiled libraries
# -----------------------------------------------------------------------------
set(CMAKE_PREFIX_PATH ${HOME}/usr)
#INCLUDE_DIRECTORIES(/home/kh8hk/usr/includ)
#LINK_DIRECTORIES(/home/kh8hk/usr/lib)

# Version
# -----------------------------------------------------------------------------
set(QMICAD_VERSION_MAJOR "0")
set(QMICAD_VERSION_MINOR "05")
set(QMICAD_VERSION_PATCH "0")

set(QMICAD_VERSION ${QMICAD_VERSION_MAJOR}.${QMICAD_VERSION_MINOR}.${QMICAD_VERSION_PATCH})

# Options
# -----------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE, "Release")
    message(STATUS "Selecting Release build type by default. You might change it by using -DCMAKE_BUILD_TYPE=Release or -DCMAKE_BUILD_TYPE=Debug on your cmake command line.")
endif ()

# Libraries
# -----------------------------------------------------------------------------
# Armadillo
find_package(Armadillo REQUIRED)

# MPI mpich2
find_package(MPI REQUIRED)

# Python libs
find_package(PythonLibs REQUIRED)

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
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})

# linking
#set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})
#set(LIBS ${LIBS} ${PYTHON_LIBRARIES})
#set(LIBS ${LIBS} ${MPI_LIBRARIES})
add_library(qmicad SHARED ${QMICAD_SOURCES})

#install(TARGETS qmicad DESTINATION lib/qmicad)
#install(FILES test.h DESTINATION include)
