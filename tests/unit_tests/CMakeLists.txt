# Copyright (c) 2014 K M Masum Habib

# -------- <QMICAD C++ Library> ------------------------------------------------
# Tests
file(GLOB QMICAD_TEST_SRCS ${CMAKE_SOURCE_DIR}/lib/tests/*.cpp)

#Run through each source
foreach(testSrc ${QMICAD_TEST_SRCS})
        #Extract the filename without an extension (NAME_WE)
        get_filename_component(testName ${testSrc} NAME_WE)

        #Add compile target
        add_executable(${testName} ${testSrc})

        #link to Boost libraries AND your targets and dependencies
        target_link_libraries(${testName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARIES} qmicad)

        #I like to move testing binaries into a testBin directory
        set_target_properties(${testName} PROPERTIES 
            RUNTIME_OUTPUT_DIRECTORY ${QMICAD_DIST_DIR}/tests)

        #Finally add it to test execution - 
        #Notice the WORKING_DIRECTORY and COMMAND
        add_test(NAME ${testName} 
                 WORKING_DIRECTORY ${QMICAD_DIST_DIR}/tests 
                 COMMAND ${QMICAD_DIST_DIR}/tests/${testName} )
endforeach(testSrc)


