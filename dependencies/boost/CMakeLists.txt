# This script downloads and builds boost from source.
#
# Inspired from:
# https://github.com/arnaudgelas/ExternalProject/blob/master/External-Boost.cmake
#

cmake_minimum_required(VERSION 3.0)

# Boost settings
#set (boost_version "1.60.0")
#set (boost_sha1 "ac74db1324e6507a309c5b139aea41af624c8110")
set (boost_version "1.62.0")
set (boost_sha1 "34a745901533cef1a5c066199b50fe83368b961b")

# Path to our working location
set(boost_root ${CMAKE_SOURCE_DIR}/dependencies/boost)

set (build_dir "build")

# Assemble boost name: boost_X_YY_Z
set (boost_name "boost_${boost_version}")
string (REPLACE "." "_" boost_name ${boost_name})
message (STATUS "Configuring boost-${boost_version} ...")

set (boost_url "/tmp/${boost_name}.tar.gz")
#set (boost_url "https://sourceforge.net/projects/boost/files/boost/${boost_version}/${boost_name}.tar.gz")

set(build_path "${boost_root}/${build_dir}")

set (boost_bootstrap_command)
set (boost_enable_mpi)
if (WIN32)
    set (boost_bootstrap_command bootstrap.bat)
    set (boost_b2_command b2.exe)
	set (boost_enable_mpi echo "using mpi $<SEMICOLON>" >> project-config.jam)
else ()
    set (boost_bootstrap_command ./bootstrap.sh)
    set (boost_b2_command ./b2)
	set (boost_enable_mpi echo "using mpi $<SEMICOLON>" >> project-config.jam)
endif()

include(ExternalProject)

set (boost_libs)
foreach (boost_component ${boost_component_list})
    if (${boost_component} STREQUAL "unit_test_framework")
        set (boost_component "test")
    endif()
    set (boost_libs ${boost_libs} --with-libraries=${boost_component})
endforeach()

ExternalProject_Add (boost
    URL ${boost_url}
    URL_HASH SHA1=${boost_sha1}
    PREFIX ${build_path}
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ${boost_bootstrap_command} ${boost_libs} --with-python=/usr/bin/python3 --prefix=${boost_root} 
    BUILD_COMMAND  ${boost_b2_command} threading=multi link=static runtime-link=static variant=release cxxflags=-std=c++14 -j8 --disable-icu install 
    INSTALL_COMMAND ""
)

ExternalProject_Add_Step(boost enable_mpi
  COMMENT "Enabling Boost.MPI ..."
  COMMAND ${boost_enable_mpi}
  DEPENDEES configure 
  DEPENDERS build
  WORKING_DIRECTORY ${build_path}/src/boost
)

set(Boost_LIBRARY_DIR ${boost_root}/lib)
set(Boost_INCLUDE_DIR ${boost_root}/include)


