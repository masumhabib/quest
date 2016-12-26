# This script downloads and builds lapack from source.
#

cmake_minimum_required(VERSION 3.0)

# lapack settings
set (lapack_version "3.6.1")
set (lapack_sha1 "7220a9bc9c766cb513007f313b804d88f110c228")

# Path to our working location
set(lapack_root ${CMAKE_SOURCE_DIR})

set (build_dir "build")

# Assemble lapack name: lapack_X_YY_Z
set (lapack_name "lapack_${lapack_version}")
#string (REPLACE "." "_" lapack_name ${lapack_name})
message (STATUS "Configuring lapack-${lapack_version} ...")

#set (lapack_url "/tmp/${lapack_name}.tar.gz")
set (lapack_url "http://www.netlib.org/lapack/lapack-${lapack_version}.tgz")

set(build_path "${lapack_root}/${build_dir}")

include(ExternalProject)

ExternalProject_Add (lapack
    URL ${lapack_url}
    #URL_HASH SHA1=${lapack_sha1}
    PREFIX ${build_path}
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND cmake -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX= .
    #CONFIGURE_COMMAND ${armadillo_bootstrap_command} --with-icu --with-python=${PYTHON_EXECUTABLE} ${armadillo_libs} --prefix=${armadillo_root} 
    #BUILD_COMMAND ${armadillo_b2_command} variant=release debug-symbols=off threading=multi link=shared,static runtime-link=shared cxxflags=-std=c++14 cxxflags=-fPIC include="${PYTHON_INCLUDE_DIRS}" -j8 install 
    INSTALL_COMMAND make install DESTDIR=${lapack_root}
)

#ExternalProject_Add_Step (armadillo enable_mpi
#  COMMENT "Enabling armadillo.MPI ..."
#  COMMAND ${armadillo_enable_mpi} #&& ${armadillo_enable_python}
#  DEPENDEES configure 
#  DEPENDERS build
#  WORKING_DIRECTORY ${build_path}/src/armadillo
#)
#
#set(Boost_LIBRARY_DIR ${armadillo_root}/lib PARENT_SCOPE)
#set(Boost_INCLUDE_DIR ${armadillo_root}/include PARENT_SCOPE)
#list(APPEND CMAKE_PREFIX_PATH ${armadillo_root})
#set (CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} PARENT_SCOPE)


