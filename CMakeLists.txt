# Configures, makes and installs QMICAD
# Copyright (c) 2014 K M Masum Habib
# Created on: June 02, 2014
#


cmake_minimum_required(VERSION 2.8)

# Version
set(QMICAD_MAJOR 0)
set(QMICAD_MINOR 06)
set(QMICAD_PATCH 0)


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
set(LINK_SEARCH_START_STATIC    on)

# For custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

#Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"x )
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
find_package(MKL REQUIRED)

# doxygen
find_package(Doxygen)

# Boost
# Shared linking of libboost_python is required for mpi::communicator
# to work both in python and C++.
set(Boost_USE_MULTITHREADED     off)
set(Boost_USE_STATIC_LIBS       off)
find_package(Boost 1.55 REQUIRED python)
set(Boost_PYTHON_LIBRARIES_SHARED ${Boost_PYTHON_LIBRARIES})
unset(Boost_LIBRARIES)
unset(Boost_PYTHON_LIBRARIES)
set(Boost_USE_STATIC_LIBS       on)
find_package(Boost 1.55 REQUIRED serialization mpi mpi_python)

# Sources and headers
# -----------------------------------------------------------------------------
# Prepare distribution directory
set(QMICAD_DIST_DIR ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-${QMICAD_VERSION})
message(STATUS "Preparing distribution directory: ${QMICAD_DIST_DIR}")
file(COPY ${CMAKE_SOURCE_DIR}/python/ DESTINATION ${QMICAD_DIST_DIR}/)
file(COPY ${CMAKE_SOURCE_DIR}/include DESTINATION ${QMICAD_DIST_DIR})
file(REMOVE "${QMICAD_DIST_DIR}/include/config.h.in")
file(MAKE_DIRECTORY "${QMICAD_DIST_DIR}/lib")

# version header
configure_file (
    "${CMAKE_SOURCE_DIR}/include/config.h.in"
    "${QMICAD_DIST_DIR}/include/config.h"
)

file(GLOB_RECURSE QMICAD_SOURCES "src/*.cpp")
set (QMICAD_INCLUDE_DIRS ${QMICAD_DIST_DIR}/include)
list(REMOVE_DUPLICATES QMICAD_INCLUDE_DIRS)

# includes
include_directories(${QMICAD_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})

# linking

add_library(qmicad SHARED ${QMICAD_SOURCES})
#target_link_libraries (qmicad ${Boost_LIBRARIES})
target_link_libraries (qmicad ${Boost_MPI_PYTHON_LIBRARIES})
target_link_libraries (qmicad ${Boost_MPI_LIBRARIES})
target_link_libraries (qmicad ${Boost_PYTHON_LIBRARIES})
target_link_libraries (qmicad -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})
target_link_libraries (qmicad ${Boost_SERIALIZATION_LIBRARIES})
target_link_libraries (qmicad ${MPI_LIBRARIES})
target_link_libraries (qmicad ${MKL_LIBRARIES})
set_target_properties(qmicad PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

# doxygen
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/qmicad.in.doxy ${CMAKE_CURRENT_BINARY_DIR}/doc/qmicad.doxy @ONLY)

    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/qmicad.doxy
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

# Prepare python package
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND cp ${CMAKE_BINARY_DIR}/lib${PROJECT_NAME}.so ${QMICAD_DIST_DIR}/lib 
    COMMAND cd ${QMICAD_DIST_DIR} && ln -sf lib/lib${PROJECT_NAME}.so qmicad.so
)

add_custom_command(
    TARGET doc  POST_BUILD
    COMMAND mkdir -p ${QMICAD_DIST_DIR}/doc
    COMMAND cp -rf ${CMAKE_BINARY_DIR}/doc/html ${QMICAD_DIST_DIR}/doc 
    COMMAND cd ${QMICAD_DIST_DIR}/doc && ln -sf html/index.html index.html
)

# installation
install(DIRECTORY ${QMICAD_DIST_DIR} DESTINATION lib)

# Symbolic links in python lib
#set(legacy_link   ${CMAKE_INSTALL_PREFIX}/lib/qmicad)
#set(legacy_target ${CMAKE_INSTALL_PREFIX}/${QMICAD_INSTALL_DIR})
#add_custom_command( OUTPUT ${legacy_link}
#    COMMAND ln -s ${legacy_target} ${legacy_link}
#    DEPENDS install ${legacy_target} 
#    COMMENT "Generating legacy symbolic link"
#)
#add_custom_target(install_legacy DEPENDS ${legacy_link} )

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    set(QMICAD_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif(${CMAKE_BUILD_TYPE} MATCHES "Release")

if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(QMICAD_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

message(STATUS "Compiler flags: ${QMICAD_FLAGS}")
message(STATUS "Linker flags: ${CMAKE_SHARED_LINKER_FLAGS}")

message(STATUS "Installation location: ${CMAKE_INSTALL_PREFIX}/${QMICAD_INSTALL_DIR}")

