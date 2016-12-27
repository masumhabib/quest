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


# - Find AMD's ACML library (no includes) which provides optimised BLAS and LAPACK functions
# This module defines
#  ACML_LIBRARIES, the libraries needed to use ACML.
#  ACML_FOUND, If false, do not try to use ACML.
# also defined, but not for general use are
#  ACML_LIBRARY, where to find the ACML library.

SET(ACML_NAMES ${ACML_NAMES} acml)
FIND_LIBRARY(ACML_LIBRARY
  NAMES ${ACML_NAMES}
  PATHS /usr/lib64 /usr/lib /usr/*/lib64 /usr/*/lib /usr/*/gfortran64/lib/ /usr/*/gfortran32/lib/ /usr/local/lib64 /usr/local/lib /opt/lib64 /opt/lib /opt/*/lib64 /opt/*/lib /opt/*/gfortran64/lib/ /opt/*/gfortran32/lib/
  )

IF (ACML_LIBRARY)
  SET(ACML_LIBRARIES ${ACML_LIBRARY})
  SET(ACML_FOUND "YES")
ELSE (ACML_LIBRARY)
  SET(ACML_FOUND "NO")
ENDIF (ACML_LIBRARY)


IF (ACML_FOUND)
   IF (NOT ACML_FIND_QUIETLY)
      MESSAGE(STATUS "Found ACML: ${ACML_LIBRARIES}")
   ENDIF (NOT ACML_FIND_QUIETLY)
ELSE (ACML_FOUND)
   IF (ACML_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find ACML")
   ENDIF (ACML_FIND_REQUIRED)
ENDIF (ACML_FOUND)

# Deprecated declarations.
GET_FILENAME_COMPONENT (NATIVE_ACML_LIB_PATH ${ACML_LIBRARY} PATH)

MARK_AS_ADVANCED(
  ACML_LIBRARY
  )
