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


# - Find a BLAS library (no includes)
# This module defines
#  BLAS_LIBRARIES, the libraries needed to use BLAS.
#  BLAS_FOUND, If false, do not try to use BLAS.
# also defined, but not for general use are
#  BLAS_LIBRARY, where to find the BLAS library.

SET(BLAS_NAMES ${BLAS_NAMES} blas)
FIND_LIBRARY(BLAS_LIBRARY
  NAMES ${BLAS_NAMES}
  PATHS /usr/lib64/atlas /usr/lib/atlas /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib
  )

IF (BLAS_LIBRARY)
  SET(BLAS_LIBRARIES ${BLAS_LIBRARY})
  SET(BLAS_FOUND "YES")
ELSE (BLAS_LIBRARY)
  SET(BLAS_FOUND "NO")
ENDIF (BLAS_LIBRARY)


IF (BLAS_FOUND)
   IF (NOT BLAS_FIND_QUIETLY)
      MESSAGE(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
   ENDIF (NOT BLAS_FIND_QUIETLY)
ELSE (BLAS_FOUND)
   IF (BLAS_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find BLAS")
   ENDIF (BLAS_FIND_REQUIRED)
ENDIF (BLAS_FOUND)

# Deprecated declarations.
GET_FILENAME_COMPONENT (NATIVE_BLAS_LIB_PATH ${BLAS_LIBRARY} PATH)

MARK_AS_ADVANCED(
  BLAS_LIBRARY
  )
