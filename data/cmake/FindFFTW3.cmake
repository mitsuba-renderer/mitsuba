# - Find FFTW3
# This module searches for the FFTW3 library.
# 
# The module defines the following variables:
#  FFTW3_INCLUDE_DIRS - where to find fftw3.h, etc.
#  FFTW3_LIBRARIES    - the libraries needed to use FFTW3.
#  FFTW3_FOUND        - if false, do not try to use FFTW3.
#
# Variables used by this module, they can change the default behavior and need
# to be set before calling find_package:
#
#  FFTW3_PRECISION - String which determines the variant of the library to find
#                    based on the floating point type. Its value is one of
#                    "DOUBLE" (default), "SINGLE", "LONG DOUBLE"
#  FFTW3_ROOT_DIR  - The preferred installation prefix for searching for FFTW3.
#                    Set if this module has problems finding the proper
#                    FFTW3 instalation.
#

# ============================================================================
# Originally created by Robert Osfield. [OpenSceneGraph]
# ============================================================================

# Search FFTW3_ROOT_DIR is set
set (_FFTW3_SEARCHES "")
if (FFTW3_ROOT_DIR)
  set(_FFTW3_SEARCH_ROOT PATHS ${FFTW3_ROOT_DIR} NO_DEFAULT_PATH)
  list(APPEND _FFTW3_SEARCHES _FFTW3_SEARCH_ROOT)
endif()

# Normal search
set(_FFTW3_SEARCH_NORMAL
  PATHS "/usr/local"
        "/opt/local/"
        "/opt/local/Library/Frameworks"
        "/opt")
list(APPEND _FFTW3_SEARCHES _FFTW3_SEARCH_NORMAL)

# Default precision is double
if (NOT FFTW3_PRECISION)
  set (FFTW3_PRECISION "DOUBLE")
else()
  string (TOUPPER "${FFTW3_PRECISION}" FFTW3_PRECISION)
endif()

# Suffix based on the precision
if (FFTW3_PRECISION STREQUAL "DOUBLE")
  set(_FFTW3_SUFFIX "")
elseif (FFTW3_PRECISION STREQUAL "SINGLE")
  set(_FFTW3_SUFFIX "f")
elseif (FFTW3_PRECISION STREQUAL "LONG DOUBLE")
  set(_FFTW3_SUFFIX "l")
else()
  message(FATAL_ERROR "Unknown FFTW3 precision: ${FFTW3_PRECISION}")
endif()

# Names of the library, including two forms of version naming
set(_FFTW3_NAMES "fftw3${_FFTW3_SUFFIX}" "fftw${_FFTW3_SUFFIX}-3")
if (WIN32)
  list(APPEND _FFTW3_NAMES "libfftw3${_FFTW3_SUFFIX}" "libfftw${_FFTW3_SUFFIX}-3")
endif()
foreach(version_minor 4 3 2 1)
  list(APPEND _FFTW3_NAMES "fftw${_FFTW3_SUFFIX}-3.${version_minor}"
                           "fftw3${_FFTW3_SUFFIX}-${version_minor}")
  if (WIN32)
    list(APPEND _FFTW3_NAMES "libfftw${_FFTW3_SUFFIX}-3.${version_minor}"
                             "libfftw3${_FFTW3_SUFFIX}-${version_minor}")
  endif()
endforeach()

# The separate threads and OpenMP library have a suffix
if (UNIX AND NOT APPLE)
  set(_FFTW3_THREADS_NAMES "")
  set(_FFTW3_OMP_NAMES "")
  foreach(name ${_FFTW3_NAMES})
    list(APPEND _FFTW3_THREADS_NAMES "${name}_threads")
    list(APPEND _FFTW3_OMP_NAMES     "${name}_omp")
  endforeach()
endif()

# Try each search configuration
foreach (search ${_FFTW3_SEARCHES})
  find_path (FFTW3_INCLUDE_DIR NAMES "fftw3.h"       ${${search}} PATH_SUFFIXES "include")
  find_library (FFTW3_LIBRARY  NAMES ${_FFTW3_NAMES} ${${search}} PATH_SUFFIXES "include")
  if (UNIX AND NOT APPLE)
    find_library (FFTW3_THREADS_LIBRARY
      NAMES ${_FFTW3_THREADS_NAMES} ${${search}} PATH_SUFFIXES "include")
    find_library (FFTW3_OMP_LIBRARY
      NAMES ${_FFTW3_OMP_NAMES} ${${search}} PATH_SUFFIXES "include")
  endif()
endforeach()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY)
if (UNIX AND NOT APPLE)
  mark_as_advanced(FFTW3_THREADS_LIBRARY FFTW3_OMP_LIBRARY)
endif()

# TODO: Define an option to look for an appropriate threads variant (pthreads or OMP)
set(_FFTW3_LIB_VARS FFTW3_LIBRARY)
if (UNIX AND NOT APPLE)
  list(APPEND _FFTW3_LIB_VARS FFTW3_THREADS_LIBRARY)
endif()

# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW3 DEFAULT_MSG FFTW3_INCLUDE_DIR ${_FFTW3_LIB_VARS})

if (FFTW3_FOUND)
  # TODO: create a compile test to find if FFTW_DLL is needed on Windows
  set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
  set(FFTW3_LIBRARIES "")
  foreach(var ${_FFTW3_LIB_VARS})
    list(APPEND FFTW3_LIBRARIES ${${var}})
  endforeach()
endif()
