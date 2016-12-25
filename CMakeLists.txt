# Configures, makes and installs QUEST
# Copyright (c) 2014 K M Masum Habib
# Created on: June 02, 2014
#


cmake_minimum_required(VERSION 3.0)

project(quest CXX)

# Options
option(VERBOSE_MAKE "Verbose cmake for debugging" off)
option(LINK_STATIC "Link everything statically" off)
option(BUILD_BOOST "Build boost from source" off)
option(BUILD_BLAS "Build blas from source" off)
option(BUILD_LAPACK "Build lapack from source" off)
option(BUILD_ARMADILLO "Build boost from source" off)

# For custom cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# To use locally compiled libraries (Armadillo and Boost), create a file
# that contains the local library paths and include it here. Some examples of
# such files can be found in cmake/local folder.
#include(cmake/local/Armadillo.cmake)

# Path to local boost build directory 
set(boost_root ${CMAKE_SOURCE_DIR}/dependencies/boost)
list(APPEND CMAKE_PREFIX_PATH ${boost_root})
# ***NOTE:*** it is important to keep python and mpi at the front of this list
# otherwise, the build system wont build mpi.so required run Boost.MPI
# from python.
set (boost_component_list "python;mpi;unit_test_framework;filesystem;program_options;random;serialization;system")
#set (boost_component_list "mpi;unit_test_framework;filesystem;program_options;random;serialization;system")

# Path to local armadillo build directory 
set(armadillo_root ${CMAKE_SOURCE_DIR}/dependencies/armadillo)
list(APPEND CMAKE_PREFIX_PATH ${armadillo_root})

# Ensure that we have a C++11 compiler.
include(cmake/CXX11.cmake)
check_for_cxx11_compiler(HAS_CXX11)
if(NOT HAS_CXX11)
    message(FATAL_ERROR "No C++11 compiler available!")
endif(NOT HAS_CXX11)
enable_cxx11()

#set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++14")

if (VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE      on)
endif()

if(LINK_STATIC)
    set(LINK_SEARCH_START_STATIC    on)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DLINK_STATIC")
endif()

# MPI mpich2
find_package(MPI REQUIRED)

# Python libs
find_package(PythonInterp 3.0 REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)

# -----------------------------------------------------------------------------
# Build boost library from source
# -----------------------------------------------------------------------------
if (BUILD_BOOST)
    add_subdirectory (dependencies/boost)

# -----------------------------------------------------------------------------
# Build lapack library from source
# -----------------------------------------------------------------------------
elseif (BUILD_LAPACK)
    add_subdirectory (dependencies/lapack)

# -----------------------------------------------------------------------------
# Build armadillo library from source
# -----------------------------------------------------------------------------
elseif (BUILD_ARMADILLO)
    add_subdirectory (dependencies/armadillo)

# -----------------------------------------------------------------------------
# QUEST package
# -----------------------------------------------------------------------------
else ()

    #set(CMAKE_BUILD_TYPE "Debug")
    
    # Build and run testing
    enable_testing()
    
    
    
    include(cmake/Versioning.cmake)
    get_version (QUEST_MAJOR QUEST_MINOR QUEST_PATCH QUEST_BUILD)
    
    # Start configuration
    set(QUEST_VERSION v${QUEST_MAJOR}.${QUEST_MINOR}.${QUEST_PATCH}${QUEST_BUILD})
    set(QUEST_SHORT_VERSION ${QUEST_MAJOR}.${QUEST_MINOR}.${QUEST_PATCH})
    message(STATUS "--------------------------------------------------------------")
    message(STATUS "Configuring QUEST ${QUEST_VERSION}")
    message(STATUS "System: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
    message(STATUS "--------------------------------------------------------------")
    
    
    # C++ compiler flags
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
    
    # Armadillo
    find_package(Armadillo REQUIRED)
    
    
    # NumPy
    include(cmake/FindNumPy.cmake)
    find_package(NumPy REQUIRED)
    
    # doxygen
    find_package(Doxygen)
    
    if(LINK_STATIC)
        set(Boost_USE_STATIC_LIBS       on)
    endif()
    find_package(Boost 1.55 REQUIRED python3 mpi mpi_python filesystem program_options random serialization system unit_test_framework)
    
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
endif()

#add_subdirectory (dependencies)


message(STATUS "--------------------------------------------------------------")

