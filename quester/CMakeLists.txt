## The QUESTER executable
#
#

include_directories(${QUEST_INCLUDE_DIRS})
include_directories(${QLOG_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${PYQUEST_INCLUDE_DIRS})
include_directories(${PYENGINE_INCLUDE_DIRS})

# sources
file(GLOB_RECURSE QUEST_SOURCES ${CMAKE_SOURCE_DIR}/quester/*.cpp)
file(GLOB_RECURSE QUEST_HEADERS ${CMAKE_SOURCE_DIR}/quester/*.h*)
add_executable(quester ${QUEST_SOURCES})

# linking
target_link_libraries (quester pyengine)
#target_link_libraries (quester ${ARMADILLO_LIBRARIES} -fopenmp)
#target_link_libraries (quester ${MPI_LIBRARIES})
#target_link_libraries (quester ${Boost_MPI_LIBRARIES})
target_link_libraries (quester ${Boost_PROGRAM_OPTIONS_LIBRARIES})
target_link_libraries (quester ${Boost_FILESYSTEM_LIBRARIES})
target_link_libraries (quester ${Boost_SYSTEM_LIBRARIES})
#target_link_libraries (quester ${Boost_SERIALIZATION_LIBRARIES})
#target_link_libraries (quester ${Boost_RANDOM_LIBRARIES}) 
#target_link_libraries (quester ${PYTHON_LIBRARIES})
#target_link_libraries (quester ${Boost_PYTHON_LIBRARIES})
#target_link_libraries (quester -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})

set_target_properties(quester PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY ${QUEST_DIST_DIR})

