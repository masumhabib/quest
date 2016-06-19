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

include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})
include_directories(${QMICAD_INCLUDE_DIRS})
set(QMICAD_INCLUDE_DIRS ${QMICAD_INCLUDE_DIRS} PARENT_SCOPE)

# sources
file(GLOB_RECURSE QMICAD_SOURCES ${CMAKE_SOURCE_DIR}/lib/src/*.cpp)
file(GLOB_RECURSE QMICAD_HEADERS ${CMAKE_SOURCE_DIR}/lib/include/*.h*)
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

# Tests
file(GLOB QMICAD_TEST_SRCS ${CMAKE_SOURCE_DIR}/lib/tests/*.cpp)

#Run through each source
foreach(testSrc ${QMICAD_TEST_SRCS})
        #Extract the filename without an extension (NAME_WE)
        get_filename_component(testName ${testSrc} NAME_WE)

        #Add compile target
        add_executable(${testName} ${testSrc})

        #link to Boost libraries AND your targets and dependencies
        target_link_libraries(${testName} ${Boost_LIBRARIES} qmicad)

        #I like to move testing binaries into a testBin directory
        set_target_properties(${testName} PROPERTIES 
            RUNTIME_OUTPUT_DIRECTORY ${QMICAD_DIST_DIR}/tests)

        #Finally add it to test execution - 
        #Notice the WORKING_DIRECTORY and COMMAND
        add_test(NAME ${testName} 
                 WORKING_DIRECTORY ${QMICAD_DIST_DIR}/tests 
                 COMMAND ${QMICAD_DIST_DIR}/tests/${testName} )
endforeach(testSrc)


