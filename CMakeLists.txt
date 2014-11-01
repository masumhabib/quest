# Configures, makes and installs QMICAD
# Copyright (c) 2014 K M Masum Habib
# Created on: June 02, 2014
#


cmake_minimum_required(VERSION 2.8)

# Version
set(QMICAD_MAJOR 0)
set(QMICAD_MINOR 08)
set(QMICAD_PATCH 1)

# Add local usr path for locally compiled libraries
set(CMAKE_PREFIX_PATH $ENV{HOME}/usr)
#set(Boost_INCLUDE_DIR $ENV{HOME}/usr/include)
#set(Boost_LIBRARY_DIR $ENV{HOME}/usr/lib)

# Start configuration
set(QMICAD_VERSION ${QMICAD_MAJOR}.${QMICAD_MINOR}.${QMICAD_PATCH})
message(STATUS "Configuring QMICAD ${QMICAD_VERSION}")
project(qmicad CXX)
message(STATUS "System: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")

# Options
# -----------------------------------------------------------------------------
set(CMAKE_VERBOSE_MAKEFILE      off)
#set(LINK_SEARCH_START_STATIC    on)

# For custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

#Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" )
endif()

# Compiler flags
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# optimization
#set(CMAKE_CXX_FLAGS_RELEASE " -m64 -O2 -MMD -MP -MF ")

# Libraries
# -----------------------------------------------------------------------------
# Armadillo
find_package(Armadillo REQUIRED)

# MPI mpich2
find_package(MPI REQUIRED)

# Python libs
find_package(PythonLibs REQUIRED)

# MKL
#find_package(MKL REQUIRED)

# doxygen
find_package(Doxygen)

# Boost
# Shared linking of libboost_python is required for mpi::communicator
# to work both in python and C++.
set(Boost_USE_MULTITHREADED     off)
#set(Boost_USE_STATIC_LIBS       off)
find_package(Boost 1.55 REQUIRED python)
set(Boost_PYTHON_LIBRARIES_SHARED ${Boost_PYTHON_LIBRARIES})
unset(Boost_LIBRARIES)
unset(Boost_PYTHON_LIBRARIES)
#set(Boost_USE_STATIC_LIBS       on)
find_package(Boost 1.55 REQUIRED serialization mpi mpi_python)

# Sources and headers
# -----------------------------------------------------------------------------

# Prepare distribution directory
set(QMICAD_DIST_DIR ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-${QMICAD_VERSION})
set(LIBRARY_OUTPUT_PATH ${QMICAD_DIST_DIR}/lib)

# version header
configure_file (
    "${CMAKE_SOURCE_DIR}/include/config.in.h"
    "${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include/config.h"
)

# includes
set(QMICAD_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)
list(REMOVE_DUPLICATES QMICAD_INCLUDE_DIRS)
list(APPEND QMICAD_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include)

include_directories(${QMICAD_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})

# -------- <QMICAD C++ Library> ------------------------------------------------
# sources
file(GLOB_RECURSE QMICAD_SOURCES "src/*.cpp")
file(GLOB_RECURSE QMICAD_HEADERS "include/*.h*")
add_library(qmicad SHARED ${QMICAD_SOURCES} ${QMICAD_HEADERS})


# linking
target_link_libraries (qmicad ${Boost_MPI_PYTHON_LIBRARIES})
target_link_libraries (qmicad ${Boost_MPI_LIBRARIES})
target_link_libraries (qmicad ${Boost_PYTHON_LIBRARIES})
target_link_libraries (qmicad -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})
target_link_libraries (qmicad ${Boost_SERIALIZATION_LIBRARIES})
#target_link_libraries (qmicad ${MPI_LIBRARIES})
target_link_libraries (qmicad ${ARMADILLO_LIBRARIES})
#target_link_libraries (qmicad ${MKL_LIBRARIES})
set_target_properties(qmicad PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

# Prepare qmicad package
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND cp -rp ${CMAKE_SOURCE_DIR}/include ${QMICAD_DIST_DIR}
    COMMAND cp -rp ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include/config.h ${QMICAD_DIST_DIR}/include
    COMMAND rm -f ${QMICAD_DIST_DIR}/include/config.in.h
)

# doxygen
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/qmicad.in.doxy ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/qmicad.doxy @ONLY)

    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/qmicad.doxy
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

    add_custom_command(
        TARGET doc  POST_BUILD
        COMMAND mkdir -p ${QMICAD_DIST_DIR}/doc
        COMMAND cp -rf ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/html ${QMICAD_DIST_DIR}/doc 
        COMMAND cd ${QMICAD_DIST_DIR}/doc && ln -sf html/index.html index.html
    )
endif(DOXYGEN_FOUND)

# -------- </QMICAD C++ Library> -----------------------------------------------

# -------- <QMICAD Python Library> ---------------------------------------------
# sources
set(PYQMICAD_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/python)
file(GLOB_RECURSE PYQMICAD_SOURCES "python/*.cpp")
file(GLOB_RECURSE PYQMICAD_HEADERS "python/*.h*")
include_directories(${PYQMICAD_INCLUDE_DIRS})
add_library(pyqmicad SHARED ${PYQMICAD_SOURCES} ${PYQMICAD_HEADERS})


# linking
target_link_libraries (pyqmicad qmicad)
set_target_properties(pyqmicad PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")
set_target_properties(pyqmicad PROPERTIES PREFIX "")
set_target_properties(pyqmicad PROPERTIES OUTPUT_NAME "qmicad")

# Prepare qmicad package
add_custom_command(
    TARGET pyqmicad POST_BUILD
    COMMAND mkdir -p  ${QMICAD_DIST_DIR}/python
    COMMAND cp    -rp ${CMAKE_SOURCE_DIR}/python/*.py ${QMICAD_DIST_DIR}/python
    COMMAND cp    -rp ${CMAKE_SOURCE_DIR}/python/simulators ${QMICAD_DIST_DIR}/python
    COMMAND cd        ${QMICAD_DIST_DIR} && ln -sf python qmicad
    COMMAND mv        ${QMICAD_DIST_DIR}/lib/qmicad.so  ${QMICAD_DIST_DIR}/python
    COMMAND cp    -rp ${QMICAD_DIST_DIR}/lib/libqmicad.so  ${QMICAD_DIST_DIR}/python
)

# -------- </QMICAD Python Library> ---------


# installation
install(DIRECTORY ${QMICAD_DIST_DIR} DESTINATION .)

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    set(QMICAD_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif(${CMAKE_BUILD_TYPE} MATCHES "Release")

if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(QMICAD_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

message(STATUS "Compiler flags: ${QMICAD_FLAGS}")
message(STATUS "Linker flags: ${CMAKE_SHARED_LINKER_FLAGS}")

message(STATUS "Installation location: ${CMAKE_INSTALL_PREFIX}/${QMICAD_INSTALL_DIR}")

