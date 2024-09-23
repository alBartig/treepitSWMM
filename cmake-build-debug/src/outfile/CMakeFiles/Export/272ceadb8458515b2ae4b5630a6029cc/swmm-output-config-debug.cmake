#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "swmm-output" for configuration "Debug"
set_property(TARGET swmm-output APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(swmm-output PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/lib/libswmm-output.dll.a"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/libswmm-output.dll"
  )

list(APPEND _cmake_import_check_targets swmm-output )
list(APPEND _cmake_import_check_files_for_swmm-output "${_IMPORT_PREFIX}/lib/libswmm-output.dll.a" "${_IMPORT_PREFIX}/bin/libswmm-output.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
