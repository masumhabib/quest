# Copyright (c) 2015 K M Masum Habib<masum.habib@gmail.com>

# includes
set(PYTMFSC_INCLUDE_DIRS ${TMFSC_ROOT_DIR}/pylib)
list(REMOVE_DUPLICATES PYTMFSC_INCLUDE_DIRS)

include_directories(${PYTMFSC_INCLUDE_DIRS})
include_directories(${TMFSC_INCLUDE_DIRS})

include_directories(${PYQMICAD_INCLUDE_DIRS})
include_directories(${QMICAD_INCLUDE_DIRS})

# sources
file(GLOB_RECURSE PYTMFSC_SOURCES "${TMFSC_ROOT_DIR}/pylib/*.cpp")
add_library(pytmfsc SHARED ${PYTMFSC_SOURCES})

# linking
target_link_libraries (pytmfsc tmfsc)
target_link_libraries (pytmfsc pyqmicad)
target_link_libraries (pytmfsc ${PYTHON_LIBRARIES})
target_link_libraries (pytmfsc ${Boost_MPI_PYTHON_LIBRARIES})
target_link_libraries (pytmfsc ${Boost_PYTHON_LIBRARIES})
target_link_libraries (pytmfsc -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})
set_target_properties (pytmfsc PROPERTIES PREFIX "")
set_target_properties (pytmfsc PROPERTIES OUTPUT_NAME "tmfsc")

# Prepare qmicad package
add_custom_command(
    TARGET pytmfsc POST_BUILD
    COMMAND cp        ${QMICAD_DIST_DIR}/lib/tmfsc${CMAKE_SHARED_LIBRARY_SUFFIX}  ${QMICAD_DIST_DIR}/qmicad/tmfsc.so
    COMMAND mkdir -p  ${QMICAD_DIST_DIR}/simulators/tmfsc
    COMMAND cp    -rp ${TMFSC_ROOT_DIR}/simulators/python/* ${QMICAD_DIST_DIR}/simulators/tmfsc
    COMMAND cp    -rp ${TMFSC_ROOT_DIR}/tests/* ${QMICAD_DIST_DIR}/tests
)


