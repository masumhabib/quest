# Copyright (c) 2014 K M Masum Habib

# Tests

set(QUEST_MINI_REGRESSION_DIR ${QUEST_TEST_DIR}/mini_regression)
message(STATUS "Adding regression test: mini_regression")


add_custom_target(mini_regression ALL
    COMMAND mkdir -p ${QUEST_MINI_REGRESSION_DIR} 
    COMMAND cp    -rp *.py ${QUEST_MINI_REGRESSION_DIR}
    COMMAND cp    -rp *.m ${QUEST_MINI_REGRESSION_DIR}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/mini_regression
    DEPENDS pyquest
)

add_custom_command(TARGET mini_regression POST_BUILD
    WORKING_DIRECTORY ${QUEST_MINI_REGRESSION_DIR}
    COMMAND ln -sf ${QUEST_DIST_DIR}/quest .
)

