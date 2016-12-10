# Copyright (c) 2014 K M Masum Habib

# Tests
file(GLOB PYENGINE_TEST_SRCS test_*.cpp)
include_directories(${PYENGINE_INCLUDE_DIRS})

#Run through each source
foreach(testSrc ${PYENGINE_TEST_SRCS})
        get_filename_component(testName ${testSrc} NAME_WE)

        message (STATUS "Adding unit test: ${testName}")
        add_executable(${testName} ${testSrc})

        target_link_libraries(${testName} pyengine)
        target_link_libraries(${testName} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARIES})
        target_link_libraries (${testName} ${PYTHON_LIBRARIES})
        target_link_libraries (${testName} ${Boost_PYTHON3_LIBRARIES})
        target_link_libraries (${testName} -L${Boost_LIBRARY_DIRS} ${Boost_PYTHON_LIBRARIES_SHARED})

        set_target_properties(${testName} PROPERTIES 
            RUNTIME_OUTPUT_DIRECTORY ${QUEST_BUILD_DIR}/tests/pyengine/unit_tests)

        add_test(NAME ${testName} 
            WORKING_DIRECTORY ${QUEST_BUILD_DIR}/tests/pyengine/unit_tests
            COMMAND ${QUEST_BUILD_DIR}/tests/pyengine/unit_tests/${testName} 
            --color_output=yes )
endforeach(testSrc)


