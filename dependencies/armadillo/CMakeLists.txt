# This script downloads and builds armadillo from source.
#
# Inspired from:
# https://github.com/arnaudgelas/ExternalProject/blob/master/External-Boost.cmake
#

cmake_minimum_required(VERSION 3.0)

# armadillo settings
set (armadillo_version "7.500.2")
set (armadillo_sha1 "75b2266c57f8c51ba76979ee98e6338cb4c9f5fc")

# Path to our working location
set(armadillo_root ${CMAKE_SOURCE_DIR}/dependencies/armadillo)

set (build_dir "build")

# Assemble armadillo name: armadillo_X_YY_Z
set (armadillo_name "armadillo_${armadillo_version}")
#string (REPLACE "." "_" armadillo_name ${armadillo_name})
message (STATUS "Configuring armadillo-${armadillo_version} ...")

#set (armadillo_url "/tmp/${armadillo_name}.tar.gz")
set (armadillo_url "https://sourceforge.net/projects/arma/files/armadillo-${armadillo_version}.tar.xz")

set(build_path "${armadillo_root}/${build_dir}")

#set (armadillo_bootstrap_command)
#set (armadillo_enable_mpi)
#if (WIN32)
#    set (armadillo_bootstrap_command bootstrap.bat)
#    set (armadillo_b2_command b2.exe)
#	set (armadillo_enable_mpi echo "using mpi $<SEMICOLON>" >> project-config.jam)
#else ()
#    set (armadillo_bootstrap_command ./bootstrap.sh)
#    set (armadillo_b2_command ./b2)
#	set (armadillo_enable_mpi echo "using mpi $<SEMICOLON>" >> project-config.jam)
#endif()
#
include(ExternalProject)
#
#set (armadillo_libs)
#foreach (armadillo_component ${armadillo_component_list})
#    if (${armadillo_component} STREQUAL "unit_test_framework")
#        set (armadillo_component "test")
#    endif()
#    set (armadillo_libs ${armadillo_libs} --with-libraries=${armadillo_component})
#endforeach()


ExternalProject_Add (armadillo
    URL ${armadillo_url}
    URL_HASH SHA1=${armadillo_sha1}
    PREFIX ${build_path}
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND cmake -DCMAKE_INSTALL_PREFIX= .
    #CONFIGURE_COMMAND ${armadillo_bootstrap_command} --with-icu --with-python=${PYTHON_EXECUTABLE} ${armadillo_libs} --prefix=${armadillo_root} 
    #BUILD_COMMAND ${armadillo_b2_command} variant=release debug-symbols=off threading=multi link=shared,static runtime-link=shared cxxflags=-std=c++14 cxxflags=-fPIC include="${PYTHON_INCLUDE_DIRS}" -j8 install 
    INSTALL_COMMAND make install DESTDIR=${armadillo_root}
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


