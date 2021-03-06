# Copyright (c) 2014 K M Masum Habib

# -------- <QUEST C++ Library> ------------------------------------------------
# version header
#configure_file (
#    "${CMAKE_SOURCE_DIR}/lib/include/config.in.h"
#    "${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include/config.h"
#)

# includes
set(QLOG_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/qlog/include)
list(REMOVE_DUPLICATES QLOG_INCLUDE_DIRS)
#list(APPEND QUEST_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include)

include_directories(${QLOG_INCLUDE_DIRS})
#include_directories(${ARMADILLO_INCLUDE_DIRS})
#include_directories(${MPI_INCLUDE_PATH})
set(QLOG_INCLUDE_DIRS ${QLOG_INCLUDE_DIRS} PARENT_SCOPE)

# sources
file(GLOB_RECURSE QLOG_SOURCES ${CMAKE_SOURCE_DIR}/qlog/src/*.cpp)
file(GLOB_RECURSE QLOG_HEADERS ${CMAKE_SOURCE_DIR}/qlog/src/*.h*)
add_library(qlog STATIC ${QLOG_SOURCES})

# linking
#target_link_libraries (quest ${MPI_LIBRARIES})
#target_link_libraries (quest ${Boost_MPI_LIBRARIES})
#target_link_libraries (quest ${Boost_SERIALIZATION_LIBRARIES})
#target_link_libraries (quest ${Boost_RANDOM_LIBRARIES}) 
#target_link_libraries (quest ${ARMADILLO_LIBRARIES})
#target_link_libraries (quest ${ARMADILLO_LIBRARIES} -fopenmp)
#set_target_properties(quest PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

# Prepare quest package
#add_custom_command(
#    TARGET ${PROJECT_NAME} POST_BUILD
#    COMMAND cp -rp ${CMAKE_SOURCE_DIR}/lib/include ${QUEST_DIST_DIR}
#    COMMAND cp -rp ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include/config.h ${QUEST_DIST_DIR}/include
#    COMMAND rm -f ${QUEST_DIST_DIR}/include/config.in.h
#)

# doxygen
#if(DOXYGEN_FOUND)
#    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/quest.in.doxy ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/quest.doxy @ONLY)
#
#    add_custom_target(doc
#        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/quest.doxy
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir
#        COMMENT "Generating API documentation with Doxygen" VERBATIM
#    )
#
#    add_custom_command(
#        TARGET doc  POST_BUILD
#        COMMAND mkdir -p ${QUEST_DIST_DIR}/doc
#        COMMAND cp -rf ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/html ${QUEST_DIST_DIR}/doc 
#        COMMAND cd ${QUEST_DIST_DIR}/doc && ln -sf html/index.html index.html
#    )
#endif(DOXYGEN_FOUND)

# -------- </QUEST C++ Library> -----------------------------------------------

