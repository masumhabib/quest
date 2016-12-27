# Configures, makes and installs QUEST
# Copyright (c) 2014 K M Masum Habib
# Created on: June 02, 2014
#


cmake_minimum_required(VERSION 3.0)

project(quest CXX)

# Options
option(VERBOSE_MAKE "Verbose cmake for debugging" off)
option(LINK_STATIC "Link everything statically" on)

# Ensure that we have a C++11 compiler.
include(cmake/CXX11.cmake)
check_for_cxx11_compiler(HAS_CXX11)
if(NOT HAS_CXX11)
    message(FATAL_ERROR "No C++11 compiler available!")
endif(NOT HAS_CXX11)
enable_cxx11()
#set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++14")

# Version number
include(cmake/Versioning.cmake)
get_version (QUEST_MAJOR QUEST_MINOR QUEST_PATCH QUEST_BUILD)

set(QUEST_VERSION v${QUEST_MAJOR}.${QUEST_MINOR}.${QUEST_PATCH}${QUEST_BUILD})
set(QUEST_SHORT_VERSION ${QUEST_MAJOR}.${QUEST_MINOR}.${QUEST_PATCH})
message(STATUS "--------------------------------------------------------------")
message(STATUS "Configuring QUEST ${QUEST_VERSION}")
message(STATUS "System: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "--------------------------------------------------------------")

# Prefer static libs 
if (LINK_STATIC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLINK_STATIC")
    #set(LINK_SEARCH_START_STATIC    on)
    #set(LINK_SEARCH_END_STATIC    on)
    if (UNIX)
        set (CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    elseif (WIN32)
        set (CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif ()
endif ()

# For custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
include (cmake/FindNumPy.cmake)
include (cmake/boost_components.cmake)

# To use locally compiled libraries (Armadillo and Boost), create a file
# that contains the local library paths and include it here. Some examples of
# such files can be found in cmake/local folder.
#include(cmake/local/Armadillo.cmake)

# Path to local boost build directory 
set(boost_root ${CMAKE_SOURCE_DIR}/third_party/boost)
list(APPEND CMAKE_PREFIX_PATH ${boost_root})

# Path to local lapack build directory 
set(lapack_root ${CMAKE_SOURCE_DIR}/third_party/lapack)
list(APPEND CMAKE_PREFIX_PATH ${lapack_root})
#list(APPEND CMAKE_PREFIX_PATH ${lapack_root}/lib64/)

# Path to local arpack build directory 
set(arpack_root ${CMAKE_SOURCE_DIR}/third_party/arpack)
list(APPEND CMAKE_PREFIX_PATH ${arpack_root})

# Path to local superlu build directory 
set(superlu_root ${CMAKE_SOURCE_DIR}/third_party/superlu)
list(APPEND CMAKE_PREFIX_PATH ${superlu_root})

# Path to local armadillo build directory 
set(armadillo_root ${CMAKE_SOURCE_DIR}/third_party/armadillo)
list(APPEND CMAKE_PREFIX_PATH ${armadillo_root})

# Find linread algebra libraries: MKL, ACML, LAPACK/BLAS
include (cmake/FindLinearAlgebraLibs.cmake)

if (VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE      on)
endif()

# Build and run testing
enable_testing()


# C++ compiler flags
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Werror=return-type")
# Just to suppress deprecated warning from boost for using auto_ptr
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")


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

# MPI mpich2
find_package(MPI REQUIRED)

# Python libs
find_package(PythonInterp 3.0 REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)


# Armadillo
find_package(Armadillo REQUIRED)


# NumPy
find_package(NumPy REQUIRED)

# doxygen
find_package(Doxygen)

if(LINK_STATIC)
    set(Boost_USE_STATIC_LIBS       on)
endif()
find_package(Boost 1.55 REQUIRED python3 mpi_python ${boost_component_list})

include_directories(${Boost_INCLUDE_DIR})

# Prepare distribution directory
set(QUEST_BUILD_DIR ${PROJECT_BINARY_DIR}/..)
set(QUEST_DIST_DIR ${QUEST_BUILD_DIR}/${PROJECT_NAME}-${QUEST_SHORT_VERSION})
set(LIBRARY_OUTPUT_PATH ${QUEST_DIST_DIR}/lib)
set(QUEST_TEST_DIR ${QUEST_BUILD_DIR}/tests)

# installation
install(DIRECTORY ${QUEST_DIST_DIR} DESTINATION .)

add_subdirectory(lib)
add_subdirectory(qlog)
add_subdirectory(pylib)
add_subdirectory(pyengine)
add_subdirectory(quester)
add_subdirectory(tests)
add_subdirectory(doc)


message(STATUS "--------------------------------------------------------------")

