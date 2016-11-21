# Copyright (c) 2014 K M Masum Habib

# -------- <QUEST C++ Library> ------------------------------------------------
# version header
configure_file (
    "${CMAKE_SOURCE_DIR}/lib/include/config.in.h"
    "${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include/config.h"
)

# includes
set(QUEST_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/include)
list(REMOVE_DUPLICATES QUEST_INCLUDE_DIRS)
list(APPEND QUEST_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include)

include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${QUEST_INCLUDE_DIRS})
set(QUEST_INCLUDE_DIRS ${QUEST_INCLUDE_DIRS} PARENT_SCOPE)

# sources
file(GLOB_RECURSE QUEST_SOURCES ${CMAKE_SOURCE_DIR}/lib/src/*.cpp)
file(GLOB_RECURSE QUEST_HEADERS ${CMAKE_SOURCE_DIR}/lib/include/*.h*)
add_library(quest SHARED ${QUEST_SOURCES})

# linking
target_link_libraries (quest ${MPI_LIBRARIES})
target_link_libraries (quest ${Boost_MPI_LIBRARIES})
target_link_libraries (quest ${Boost_SERIALIZATION_LIBRARIES})
target_link_libraries (quest ${Boost_RANDOM_LIBRARIES}) 
target_link_libraries (quest ${ARMADILLO_LIBRARIES})
#target_link_libraries (quest ${ARMADILLO_LIBRARIES} -fopenmp)
#set_target_properties(quest PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

# Prepare quest package
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND cp -rp ${CMAKE_SOURCE_DIR}/lib/include ${QUEST_DIST_DIR}
    COMMAND cp -rp ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include/config.h ${QUEST_DIST_DIR}/include
    COMMAND rm -f ${QUEST_DIST_DIR}/include/config.in.h
)


# -------- </QUEST C++ Library> -----------------------------------------------

