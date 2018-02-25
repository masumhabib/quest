# Copyright (c) 2014 K M Masum Habib

# -------- <QUEST Python Library> ---------------------------------------------
# sources
set(PYQUEST_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/pylib)
set(PYQUEST_INCLUDE_DIRS ${PYQUEST_INCLUDE_DIRS} PARENT_SCOPE)
file(GLOB_RECURSE PYQUEST_SOURCES "./*.cpp")
file(GLOB_RECURSE PYQUEST_HEADERS "./*.h*")

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${QUEST_INCLUDE_DIRS})
include_directories(${PYQUEST_INCLUDE_DIRS})
add_library(pyquest SHARED ${PYQUEST_SOURCES})


# linking
target_link_libraries (pyquest quest)
target_link_libraries (pyquest ${PYTHON_LIBRARIES})
target_link_libraries (pyquest ${Boost_MPI_PYTHON_LIBRARIES})
target_link_libraries (pyquest ${Boost_PYTHON_LIBRARIES})
target_link_libraries (pyquest -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})
#set_target_properties (pyquest PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")
set_target_properties (pyquest PROPERTIES PREFIX "")
set_target_properties (pyquest PROPERTIES OUTPUT_NAME "quest")

# Prepare quest python package
add_custom_target(pyquest_package ALL
    COMMAND mkdir -p  ${QUEST_DIST_DIR}/quest
    COMMAND mkdir -p  ${QUEST_DIST_DIR}/quest/simulators
    COMMAND cp    -rp ${CMAKE_SOURCE_DIR}/pylib/python/*.py ${QUEST_DIST_DIR}/quest
    COMMAND cp    -rp ${CMAKE_SOURCE_DIR}/simulators/python/* ${QUEST_DIST_DIR}/quest/simulators
    COMMAND cp        ${QUEST_DIST_DIR}/lib/quest${CMAKE_SHARED_LIBRARY_SUFFIX}  ${QUEST_DIST_DIR}/quest/quest.so
    COMMAND cp    -rp ${QUEST_DIST_DIR}/lib/libquest${CMAKE_SHARED_LIBRARY_SUFFIX} ${QUEST_DIST_DIR}/quest

    DEPENDS quest pyquest 
)

add_custom_target(quest_utils ALL
    COMMAND mkdir -p  ${QUEST_DIST_DIR}/utils
    COMMAND cp    -rp ${CMAKE_SOURCE_DIR}/utils/matlab ${QUEST_DIST_DIR}/utils
    DEPENDS quest pyquest pyquest_package
)

# -------- </QUEST Python Library> ---------

