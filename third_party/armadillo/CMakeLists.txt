# This script downloads and builds armadillo from source.
#

cmake_minimum_required(VERSION 3.0)

# armadillo settings
set (armadillo_version "7.500.2")
set (armadillo_sha1 "75b2266c57f8c51ba76979ee98e6338cb4c9f5fc")

# Path to our working location
set(armadillo_root ${CMAKE_SOURCE_DIR})
set (lapack_lib ${armadillo_root}/../lapack/lib64)
set (superlu_lib ${armadillo_root}/../superlu/lib64)
set (arpack_lib ${armadillo_root}/../arpack/lib64)

set (build_dir "build")

# Assemble armadillo name: armadillo_X_YY_Z
set (armadillo_name "armadillo_${armadillo_version}")
message (STATUS "Configuring armadillo-${armadillo_version} ...")

set (armadillo_url "https://sourceforge.net/projects/arma/files/armadillo-${armadillo_version}.tar.xz")

set(build_path "${armadillo_root}/${build_dir}")

include(ExternalProject)

message ("DBG: PREFIX ${CMAKE_PREFIX_PATH}")
ExternalProject_Add (armadillo
    URL ${armadillo_url}
    URL_HASH SHA1=${armadillo_sha1}
    PREFIX ${build_path}
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND cmake -DCMAKE_PREFIX_PATH=${lapack_lib}$<SEMICOLON>${arpack_lib}$<SEMICOLON>${superlu_lib} -DCMAKE_INSTALL_PREFIX= .
    INSTALL_COMMAND make install DESTDIR=${armadillo_root}
)



