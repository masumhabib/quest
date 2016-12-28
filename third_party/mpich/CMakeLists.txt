# This script downloads and builds mpich from source.
#

cmake_minimum_required(VERSION 3.0)

# mpich settings
set (mpich_version "3.2")
set (mpich_sha1 "5fae7df02f750b6966837b554f9985a9c574236f")

# Path to our working location
set(mpich_root ${CMAKE_SOURCE_DIR})

set (build_dir "build")

message (STATUS "Configuring mpich-${mpich_version} ...")

set (mpich_url "http://www.mpich.org/static/downloads/${mpich_version}/mpich-${mpich_version}.tar.gz")

set(build_path "${mpich_root}/${build_dir}")

include(ExternalProject)

set(ENV{CFLAGS} -fPIC)
set(ENV{CXXFLAGS} -fPIC)
set(ENV{FFLAGS} -fPIC)

ExternalProject_Add (mpich
    URL ${mpich_url}
    URL_HASH SHA1=${mpich_sha1}
    PREFIX ${build_path}
    #WORKING_DIRECTORY ${build_path}/src/mpich
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND export CFLAGS=-fPIC && export CXXFLAGS=-fPIC && export FFLAGS=-fPIC && ./configure --prefix=${mpich_root} 
    #BUILD_COMMAND make
    #INSTALL_COMMAND make install
)



