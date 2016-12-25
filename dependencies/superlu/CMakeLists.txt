# This script downloads and builds superlu from source.
#

cmake_minimum_required(VERSION 3.0)

# superlu settings
set (superlu_version "5.2.1")
set (superlu_sha1 "7220a9bc9c766cb513007f313b804d88f110c228")

# Path to our working location
#set(superlu_root ${CMAKE_SOURCE_DIR}/dependencies/superlu)

set (build_dir "build")

# Assemble superlu name: superlu_X_YY_Z
set (superlu_name "superlu_${superlu_version}")
#string (REPLACE "." "_" superlu_name ${superlu_name})
message (STATUS "Configuring superlu-${superlu_version} ...")

#set (superlu_url "/tmp/${superlu_name}.tar.gz")
set (superlu_url "http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_${superlu_version}.tar.gz")

set(build_path "${superlu_root}/${build_dir}")

include(ExternalProject)

ExternalProject_Add (superlu
    URL ${superlu_url}
    #URL_HASH SHA1=${superlu_sha1}
    PREFIX ${build_path}
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND cmake -Denable_blaslib=OFF -Denable_tests=OFF -DCMAKE_PREFIX_PATH=${lapack_root}/lib64/ -DCMAKE_INSTALL_PREFIX= .
    INSTALL_COMMAND make install DESTDIR=${superlu_root}
)


