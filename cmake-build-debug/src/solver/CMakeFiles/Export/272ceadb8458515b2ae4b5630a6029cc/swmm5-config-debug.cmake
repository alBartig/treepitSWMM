#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "swmm5" for configuration "Debug"
set_property(TARGET swmm5 APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(swmm5 PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/lib/libswmm5.dll.a"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/libswmm5.dll"
  )

list(APPEND _cmake_import_check_targets swmm5 )
list(APPEND _cmake_import_check_files_for_swmm5 "${_IMPORT_PREFIX}/lib/libswmm5.dll.a" "${_IMPORT_PREFIX}/bin/libswmm5.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
