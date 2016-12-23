#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "armadillo" for configuration ""
set_property(TARGET armadillo APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(armadillo PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_NOCONFIG "/opt/intel/mkl/lib/intel64/libmkl_rt.so;/usr/lib64/libarpack.so;/usr/lib64/libsuperlu.so"
  IMPORTED_LOCATION_NOCONFIG "/usr/lib64/libarmadillo.so.7.500.2"
  IMPORTED_SONAME_NOCONFIG "libarmadillo.so.7"
  )

list(APPEND _IMPORT_CHECK_TARGETS armadillo )
list(APPEND _IMPORT_CHECK_FILES_FOR_armadillo "/usr/lib64/libarmadillo.so.7.500.2" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
