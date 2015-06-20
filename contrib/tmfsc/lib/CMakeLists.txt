# Copyright (c) 2015 K M Masum Habib<masum.habib@gmail.com>

# includes
set(TMFSC_INCLUDE_DIRS ${TMFSC_ROOT_DIR}/lib/include)
list(REMOVE_DUPLICATES TMFSC_INCLUDE_DIRS)
set(TMFSC_INCLUDE_DIRS ${TMFSC_INCLUDE_DIRS} PARENT_SCOPE)

include_directories(${QMICAD_INCLUDE_DIRS})
include_directories(${TMFSC_INCLUDE_DIRS})

# sources
file(GLOB_RECURSE TMFSC_SOURCES "${TMFSC_ROOT_DIR}/lib/src/*.cpp")
add_library(tmfsc SHARED ${TMFSC_SOURCES})

# linking
target_link_libraries (tmfsc qmicad)

#I like to keep test files in a separate source directory called test
file(GLOB TMFSC_TEST_SRCS ${TMFSC_ROOT_DIR}/lib/tests/*.cpp)

#Run through each source
foreach(testSrc ${TMFSC_TEST_SRCS})
        #Extract the filename without an extension (NAME_WE)
        get_filename_component(testName ${testSrc} NAME_WE)

        #Add compile target
        add_executable(${testName} ${testSrc})

        #link to Boost libraries AND your targets and dependencies
        target_link_libraries(${testName} ${Boost_LIBRARIES} tmfsc)

        #I like to move testing binaries into a testBin directory
        set_target_properties(${testName} PROPERTIES 
            RUNTIME_OUTPUT_DIRECTORY ${QMICAD_DIST_DIR}/tests)

        #Finally add it to test execution - 
        #Notice the WORKING_DIRECTORY and COMMAND
        add_test(NAME ${testName} 
                 WORKING_DIRECTORY ${QMICAD_DIST_DIR}/tests 
                 COMMAND ${QMICAD_DIST_DIR}/tests/${testName} )
endforeach(testSrc)


