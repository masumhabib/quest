# This script downloads and builds arpack from source.
#

cmake_minimum_required(VERSION 3.0)

# arpack settings
set (arpack_version "master")
set (arpack_sha1 "7220a9bc9c766cb513007f313b804d88f110c228")

# MPI for Boost.MPI
#find_package(MPI REQUIRED)

# Path to our working location
set (arpack_root ${CMAKE_SOURCE_DIR})
set (lapack_lib ${arpack_root}/../lapack/lib64)

set (build_dir "build")

# Assemble arpack name: arpack_X_YY_Z
set (arpack_name "arpack_${arpack_version}")
#string (REPLACE "." "_" arpack_name ${arpack_name})
message (STATUS "Configuring arpack-${arpack_version} ...")

#set (arpack_url "/tmp/${arpack_name}.tar.gz")
set (arpack_url "https://github.com/opencollab/arpack-ng/archive/${arpack_version}.zip")

set(build_path "${arpack_root}/${build_dir}")

include(ExternalProject)

ExternalProject_Add (arpack
    URL ${arpack_url}
    #URL_HASH SHA1=${arpack_sha1}
    PREFIX ${build_path}
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND cmake -DEXAMPLES=OFF -DMPI=ON -DCMAKE_C_FLAGS=-fPIC -DCMAKE_Fortran_FLAGS=-fPIC -DCMAKE_PREFIX_PATH=${lapack_lib} -DCMAKE_INSTALL_PREFIX= .
    INSTALL_COMMAND make install DESTDIR=${arpack_root}
)


