# Copyright (c) 2014 K M Masum Habib

# -------- <QMICAD C++ Library> ------------------------------------------------
# version header
configure_file (
    "${CMAKE_SOURCE_DIR}/lib/include/config.in.h"
    "${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include/config.h"
)

# includes
set(QMICAD_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/lib/include)
list(REMOVE_DUPLICATES QMICAD_INCLUDE_DIRS)
list(APPEND QMICAD_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include)

include_directories(${QMICAD_INCLUDE_DIRS})
set(QMICAD_INCLUDE_DIRS ${QMICAD_INCLUDE_DIRS} PARENT_SCOPE)

# sources
file(GLOB_RECURSE QMICAD_SOURCES "src/*.cpp")
file(GLOB_RECURSE QMICAD_HEADERS "include/*.h*")
add_library(qmicad SHARED ${QMICAD_SOURCES})

# linking
target_link_libraries (qmicad ${MPI_LIBRARIES})
target_link_libraries (qmicad ${Boost_MPI_LIBRARIES})
target_link_libraries (qmicad ${Boost_SERIALIZATION_LIBRARIES})
target_link_libraries (qmicad ${Boost_RANDOM_LIBRARIES}) 
target_link_libraries (qmicad ${ARMADILLO_LIBRARIES})
#set_target_properties(qmicad PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

# Prepare qmicad package
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND cp -rp ${CMAKE_SOURCE_DIR}/lib/include ${QMICAD_DIST_DIR}
    COMMAND cp -rp ${CMAKE_BINARY_DIR}/CMakeFiles/${PROJECT_NAME}.dir/include/config.h ${QMICAD_DIST_DIR}/include
    COMMAND rm -f ${QMICAD_DIST_DIR}/include/config.in.h
)

# doxygen
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/qmicad.in.doxy ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/qmicad.doxy @ONLY)

    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/qmicad.doxy
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

    add_custom_command(
        TARGET doc  POST_BUILD
        COMMAND mkdir -p ${QMICAD_DIST_DIR}/doc
        COMMAND cp -rf ${CMAKE_BINARY_DIR}/CMakeFiles/doc.dir/html ${QMICAD_DIST_DIR}/doc 
        COMMAND cd ${QMICAD_DIST_DIR}/doc && ln -sf html/index.html index.html
    )
endif(DOXYGEN_FOUND)

# -------- </QMICAD C++ Library> -----------------------------------------------
