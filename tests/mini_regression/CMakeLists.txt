# Copyright (c) 2014 K M Masum Habib

# Tests

set(QMICAD_MINI_REGRESSION_DIR ${QMICAD_TEST_DIR}/mini_regression)
message(STATUS "Adding regression test: mini_regression")


add_custom_target(mini_regression ALL
    COMMAND mkdir -p ${QMICAD_MINI_REGRESSION_DIR} 
    COMMAND cp    -rp *.py ${QMICAD_MINI_REGRESSION_DIR}
    COMMAND cp    -rp *.m ${QMICAD_MINI_REGRESSION_DIR}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/mini_regression
    DEPENDS pyqmicad
)

