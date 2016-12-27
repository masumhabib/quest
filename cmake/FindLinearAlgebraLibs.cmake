# Find linear algebra libraries: BLAS, LAPACK, ACML, MKL
#
# Sets:
#       QUEST_LAPACK_FOUND
#       QUEST_LAPACL_LIBRARIES
#       QUEST_ARPACK_FOUND
#       QUEST_ARPACK_LIBRARIES
#       QUEST_SUPERLU_FOUND
#       QUEST_SUPERLU_LIBRARIES
#

set (QUEST_LAPACK_FOUND false)

# find lapack and blas
if (APPLE)
    set (APPLE_LAPACK_FOUND true)
elseif (UNIX)
    include (FindMKL)
    include (FindACMLMP)
    include (FindACML)
    include (FindLAPACK)
    include (FindOpenBLAS)
    include (FindATLAS)
    include (FindBLAS)

    if (MKL_FOUND)
        set (QUEST_LAPACK_LIBRARIES ${MKL_LIBRARIES})
        set (QUEST_LAPACK_FOUND true)
    elseif (ACMLMP_FOUND)
        set (QUEST_LAPACK_LIBRARIES ${ACMLMP_LIBRARIES})
        set (QUEST_LAPACK_FOUND true)
    elseif (ACML_FOUND)
        set (QUEST_LAPACK_LIBRARIES ${ACML_LIBRARIES})
        set (QUEST_LAPACK_FOUND true)
    elseif (APPLE_LAPACK_FOUND)
        set (QUEST_LAPACK_LIBRARIES "-framework Accelerate")
        set (QUEST_LAPACK_FOUND true)
    elseif (LAPACK_FOUND)
        set (QUEST_LAPACK_LIBRARIES ${LAPACK_LIBRARIES})

        if (OpenBLAS_FOUND)
            list (APPEND QUEST_LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
        elseif (ATLAS_FOUND)
            list (APPEND QUEST_LAPACK_LIBRARIES ${ATLAS_LIBRARIES})
        elseif (BLAS_FOUND)
            list (APPEND QUEST_LAPACK_LIBRARIES ${BLAS_LIBRARIES})
        endif ()
    else ()
        message (FATAL_ERROR "No LAPACK/BLAS libraries were found, aborting.")
    endif ()
else ()
endif ()

include (FindARPACK)
if (ARPACK_FOUND)
    set (QUEST_ARPACK_FOUND true)
    set (QUEST_ARPACK_LIBRARIES ${ARPACK_LIBRARY})
endif ()

include (FindSuperLU5)
if (SUPERLU_FOUND)
    set (QUEST_SUPERLU_FOUND true)
    set (QUEST_SUPERLU_LIBRARIES ${SUPERLU_LIBRARY})
endif ()



