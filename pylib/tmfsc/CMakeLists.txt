# Copyright (c) 2015 K M Masum Habib<masum.habib@gmail.com>

# includes
set(PYTMFSC_INCLUDE_DIRS ${TMFSC_ROOT_DIR}/pylib)
list(REMOVE_DUPLICATES PYTMFSC_INCLUDE_DIRS)

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PYTMFSC_INCLUDE_DIRS})
include_directories(${TMFSC_INCLUDE_DIRS})

include_directories(${PYQUEST_INCLUDE_DIRS})
include_directories(${QUEST_INCLUDE_DIRS})

# sources
file(GLOB_RECURSE PYTMFSC_SOURCES "${TMFSC_ROOT_DIR}/pylib/*.cpp")
add_library(pytmfsc SHARED ${PYTMFSC_SOURCES})

# linking
target_link_libraries (pytmfsc tmfsc)
target_link_libraries (pytmfsc pyquest)
target_link_libraries (pytmfsc ${PYTHON_LIBRARIES})
target_link_libraries (pytmfsc ${Boost_MPI_PYTHON3_LIBRARIES})
target_link_libraries (pytmfsc ${Boost_PYTHON3_LIBRARIES})
target_link_libraries (pytmfsc -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})
set_target_properties (pytmfsc PROPERTIES PREFIX "")
set_target_properties (pytmfsc PROPERTIES OUTPUT_NAME "tmfsc")

# Prepare quest package
add_custom_command(
    TARGET pytmfsc POST_BUILD
    COMMAND cp        ${QUEST_DIST_DIR}/lib/tmfsc${CMAKE_SHARED_LIBRARY_SUFFIX}  ${QUEST_DIST_DIR}/quest/tmfsc.so
    COMMAND mkdir -p  ${QUEST_DIST_DIR}/simulators/tmfsc
    COMMAND cp    -rp ${TMFSC_ROOT_DIR}/simulators/python/* ${QUEST_DIST_DIR}/simulators/tmfsc
    COMMAND cp    -rp ${TMFSC_ROOT_DIR}/tests/* ${QUEST_DIST_DIR}/tests
)


