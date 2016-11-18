##
# Build system for the python engine.
#
#


include_directories(${QUEST_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${PYQUEST_INCLUDE_DIRS})

# sources
set(PYENGINE_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/pyengine)
list(REMOVE_DUPLICATES PYENGINE_INCLUDE_DIRS)
set(PYENGINE_INCLUDE_DIRS ${PYENGINE_INCLUDE_DIRS} PARENT_SCOPE)
file(GLOB_RECURSE PYENGINE_SOURCES ${CMAKE_SOURCE_DIR}/pyengine/*.cpp)

add_library (pyengine STATIC ${PYENGINE_SOURCES})

# linking
#target_link_libraries (pyengine qlog)
#target_link_libraries (pyengine ${ARMADILLO_LIBRARIES} -fopenmp)
target_link_libraries (pyengine ${MPI_LIBRARIES})
target_link_libraries (pyengine ${Boost_MPI_LIBRARIES})
target_link_libraries (pyengine ${Boost_FILESYSTEM_LIBRARIES})
target_link_libraries (pyengine ${Boost_SYSTEM_LIBRARIES})
#target_link_libraries (pyengine ${Boost_SERIALIZATION_LIBRARIES})
#target_link_libraries (pyengine ${Boost_RANDOM_LIBRARIES}) 
target_link_libraries (pyengine ${PYTHON_LIBRARIES})
target_link_libraries (pyengine ${Boost_PYTHON_LIBRARIES})
target_link_libraries (pyengine -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})


