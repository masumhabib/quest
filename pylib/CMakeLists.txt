# Copyright (c) 2014 K M Masum Habib

# -------- <QMICAD Python Library> ---------------------------------------------
# sources
set(PYQMICAD_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/pylib)
set(PYQMICAD_INCLUDE_DIRS ${PYQMICAD_INCLUDE_DIRS} PARENT_SCOPE)
file(GLOB_RECURSE PYQMICAD_SOURCES "./*.cpp")
file(GLOB_RECURSE PYQMICAD_HEADERS "./*.h*")

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${QMICAD_INCLUDE_DIRS})
include_directories(${PYQMICAD_INCLUDE_DIRS})
add_library(pyqmicad SHARED ${PYQMICAD_SOURCES})


# linking
target_link_libraries (pyqmicad qmicad)
target_link_libraries (pyqmicad ${PYTHON_LIBRARIES})
target_link_libraries (pyqmicad ${Boost_MPI_PYTHON_LIBRARIES})
target_link_libraries (pyqmicad ${Boost_PYTHON_LIBRARIES})
target_link_libraries (pyqmicad -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})
#set_target_properties (pyqmicad PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")
set_target_properties (pyqmicad PROPERTIES PREFIX "")
set_target_properties (pyqmicad PROPERTIES OUTPUT_NAME "qmicad")

# Prepare qmicad python package
add_custom_target(pyqmicad_package ALL
    COMMAND mkdir -p  ${QMICAD_DIST_DIR}/qmicad
    COMMAND mkdir -p  ${QMICAD_DIST_DIR}/qmicad/simulators
    COMMAND cp    -rp ${CMAKE_SOURCE_DIR}/pylib/python/*.py ${QMICAD_DIST_DIR}/qmicad
    COMMAND cp    -rp ${CMAKE_SOURCE_DIR}/simulators/python/* ${QMICAD_DIST_DIR}/qmicad/simulators
    COMMAND cp        ${QMICAD_DIST_DIR}/lib/qmicad${CMAKE_SHARED_LIBRARY_SUFFIX}  ${QMICAD_DIST_DIR}/qmicad/qmicad.so
    COMMAND cp    -rp ${QMICAD_DIST_DIR}/lib/libqmicad${CMAKE_SHARED_LIBRARY_SUFFIX} ${QMICAD_DIST_DIR}/qmicad

    DEPENDS qmicad pyqmicad 
)

add_custom_target(qmicad_utils ALL
    COMMAND mkdir -p  ${QMICAD_DIST_DIR}/utils
    COMMAND cp    -rp ${CMAKE_SOURCE_DIR}/utils/matlab ${QMICAD_DIST_DIR}/utils
    DEPENDS qmicad pyqmicad pyqmicad_package
)

# -------- </QMICAD Python Library> ---------

