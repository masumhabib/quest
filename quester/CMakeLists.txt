# Copyright (c) 2014 K M Masum Habib

# -------- <QUEST Executable> ------------------------------------------------

include_directories(${QUEST_INCLUDE_DIRS})
include_directories(${QLOG_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${PYQUEST_INCLUDE_DIRS})

# sources
file(GLOB_RECURSE QUEST_SOURCES ${CMAKE_SOURCE_DIR}/quester/*.cpp)
file(GLOB_RECURSE QUEST_HEADERS ${CMAKE_SOURCE_DIR}/quester/*.h*)
add_executable(quester ${QUEST_SOURCES})

# linking
target_link_libraries (quester qlog)
#target_link_libraries (quester ${ARMADILLO_LIBRARIES} -fopenmp)
target_link_libraries (quester ${MPI_LIBRARIES})
target_link_libraries (quester ${Boost_MPI_LIBRARIES})
target_link_libraries (quester ${Boost_PROGRAM_OPTIONS_LIBRARIES})
target_link_libraries (quester ${Boost_FILESYSTEM_LIBRARIES})
target_link_libraries (quester ${Boost_SYSTEM_LIBRARIES})
#target_link_libraries (quester ${Boost_SERIALIZATION_LIBRARIES})
#target_link_libraries (quester ${Boost_RANDOM_LIBRARIES}) 
target_link_libraries (quester ${PYTHON_LIBRARIES})
target_link_libraries (quester ${Boost_PYTHON3_LIBRARIES})
target_link_libraries (quester -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})

set_target_properties(quester PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY ${QUEST_DIST_DIR})
#set_target_properties (quester PROPERTIES OUTPUT_NAME "quest")

#set_target_properties(quest PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

## Prepare quest package
#add_custom_command(
#    TARGET ${PROJECT_NAME} POST_BUILD
#    COMMAND cp -rp ${CMAKE_SOURCE_DIR}/lib/include ${QUEST_DIST_DIR}
#    COMMAND cp -rp ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include/config.h ${QUEST_DIST_DIR}/include
#    COMMAND rm -f ${QUEST_DIST_DIR}/include/config.in.h
#)
#
## doxygen
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

# -------- </QUEST Executable> -----------------------------------------------

