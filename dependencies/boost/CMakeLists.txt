# This script downloads and builds boost from source.
#
# Inspired from:
# https://github.com/arnaudgelas/ExternalProject/blob/master/External-Boost.cmake
#

# Boost settings
set (boost_version "1.62.0")
set (boost_sha1 "34a745901533cef1a5c066199b50fe83368b961b")

# Path to our working location
set(boost_root ${CMAKE_SOURCE_DIR}/dependencies/boost)

set (build_dir "build")

# Assemble boost name: boost_X_YY_Z
set (boost_name "boost_${boost_version}")
string (REPLACE "." "_" boost_name ${boost_name})
message (STATUS "Building boost-${boost_version} from source ...")

set (boost_url "https://sourceforge.net/projects/boost/files/boost/${boost_version}/${boost_name}.tar.gz")

set(build_path "${boost_root}/${build_dir}")

set (boost_bootstrap_command)
if (WIN32)
    set (boost_bootstrap_command bootstrap.bat)
    set (boost_b2_command b2.exe)
else ()
    set (boost_bootstrap_command ./bootstrap.sh)
    set (boost_b2_command ./b2)
endif()

include(ExternalProject)

ExternalProject_Add(boost
    URL ${boost_url}
    URL_HASH SHA1=${boost_sha1}
    PREFIX ${build_path}
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ${boost_bootstrap_command} 
        --with-libraries=filesystem
    BUILD_COMMAND  ${boost_b2_command}
        --disable-icu
        --prefix=${boost_root}
        --threading=single
        --link=shared,static
        --variant=release
        -j8
    INSTALL_COMMAND ${boost_b2_command} --prefix=${boost_root} install
)


set(Boost_LIBRARY_DIR ${boost_root}/lib)
set(Boost_INCLUDE_DIR ${boost_root}/include)


