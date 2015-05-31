# Copyright (c) 2014 K M Masum Habib

# TMFSC Simulator C++ Library
set(TMFSC_ROOT_DIR ${CMAKE_SOURCE_DIR}/simulators/cpp/tmfsc)

# includes
set(TMFSC_INCLUDE_DIRS ${TMFSC_ROOT_DIR}/lib/include)
list(REMOVE_DUPLICATES TMFSC_INCLUDE_DIRS)

include_directories(${QMICAD_INCLUDE_DIRS})
include_directories(${TMFSC_INCLUDE_DIRS})

set(PYTMFSC_INCLUDE_DIRS ${TMFSC_ROOT_DIR}/pylib)
include_directories(${PYQMICAD_INCLUDE_DIRS})
include_directories(${PYTMFSC_INCLUDE_DIRS})

# sources
file(GLOB_RECURSE TMFSC_SOURCES "lib/src/*.cpp")
file(GLOB_RECURSE PYTMFSC_SOURCES "pylib/*.cpp")
add_library(tmfsc SHARED ${TMFSC_SOURCES} ${PYTMFSC_SOURCES})

# linking
target_link_libraries (tmfsc qmicad)
target_link_libraries (tmfsc pyqmicad)
target_link_libraries (tmfsc ${PYTHON_LIBRARIES})
target_link_libraries (tmfsc ${Boost_MPI_PYTHON_LIBRARIES})
target_link_libraries (tmfsc ${Boost_PYTHON_LIBRARIES})
target_link_libraries (tmfsc -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})
set_target_properties (tmfsc PROPERTIES PREFIX "")
set_target_properties (tmfsc PROPERTIES OUTPUT_NAME "tmfsc")

# Prepare qmicad package
add_custom_command(
    TARGET tmfsc POST_BUILD
    COMMAND cp        ${QMICAD_DIST_DIR}/lib/tmfsc${CMAKE_SHARED_LIBRARY_SUFFIX}  ${QMICAD_DIST_DIR}/qmicad/tmfsc.so
)

