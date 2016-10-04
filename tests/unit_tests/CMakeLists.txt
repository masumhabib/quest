# Copyright (c) 2014 K M Masum Habib

# Tests
file(GLOB QUEST_TEST_SRCS test_*.cpp)
include_directories(${QUEST_INCLUDE_DIRS})

#Run through each source
foreach(testSrc ${QUEST_TEST_SRCS})
        get_filename_component(testName ${testSrc} NAME_WE)

        message (STATUS "Adding unit test: ${testName}")
        add_executable(${testName} ${testSrc})

        target_link_libraries(${testName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARIES} quest)

        set_target_properties(${testName} PROPERTIES 
            RUNTIME_OUTPUT_DIRECTORY ${QUEST_BUILD_DIR}/tests/unit_tests)

        add_test(NAME ${testName} 
            WORKING_DIRECTORY ${QUEST_BUILD_DIR}/tests/unit_tests
            COMMAND ${QUEST_BUILD_DIR}/tests/unit_tests/${testName} 
            --color_output=yes )
endforeach(testSrc)


