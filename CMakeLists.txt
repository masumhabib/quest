# Configures, makes and installs QMICAD
# Copyright (c) 2014 K M Masum Habib
# Created on: June 02, 2014
#


cmake_minimum_required(VERSION 2.8)
project(qmicad)

# Version
# -----------------------------------------------------------------------------
set(QMICAD_VERSION_MAJOR "0")
set(QMICAD_VERSION_MINOR "05")
set(QMICAD_VERSION_PATCH "0")

set(QMICAD_VERSION ${QMICAD_VERSION_MAJOR}.${QMICAD_VERSION_MINOR}.${QMICAD_VERSION_PATCH})


# Options
# -----------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE, "Release")
    message(STATUS "Selecting Release build type by default. You might change it by using -DCMAKE_BUILD_TYPE=Release or -DCMAKE_BUILD_TYPE=Debug on your cmake command line.")
endif ()

# Sources and headers
# -----------------------------------------------------------------------------
add_subdirectory(src)
add_subdirectory(doc)


#add_library(qmicad SHARED test.c)

install(TARGETS qmicad DESTINATION lib/qmicad)
#install(FILES test.h DESTINATION include)
