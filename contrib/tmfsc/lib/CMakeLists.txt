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
