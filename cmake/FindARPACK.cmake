# The original form of this file was a part of Armadillo 
# library: http://arma.sourceforge.net. 
#
# Original copyright notice:
#
# Copyright (C) 2008-2016 National ICT Australia (NICTA)
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
# -------------------------------------------------------------------
# 
# Written by Conrad Sanderson - http://conradsanderson.id.au
# Written by Ryan Curtin
# Written by Clement Creusot
#
# Adopted and revised by the QUEST project.


# - Try to find ARPACK
# Once done this will define
#
#  ARPACK_FOUND        - system has ARPACK
#  ARPACK_LIBRARY      - Link this to use ARPACK


find_library(ARPACK_LIBRARY
  NAMES arpack
  PATHS ${CMAKE_SYSTEM_LIBRARY_PATH} /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib /opt/local/lib64 /opt/local/lib
  )


IF (ARPACK_LIBRARY)
  SET(ARPACK_FOUND YES)
ELSE ()
  SET(ARPACK_FOUND NO)
ENDIF ()


IF (ARPACK_FOUND)
  IF (NOT ARPACK_FIND_QUIETLY)
     MESSAGE(STATUS "Found ARPACK: ${ARPACK_LIBRARY}")
  ENDIF (NOT ARPACK_FIND_QUIETLY)
ELSE (ARPACK_FOUND)
  IF (ARPACK_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find ARPACK")
  ENDIF (ARPACK_FIND_REQUIRED)
ENDIF (ARPACK_FOUND)
