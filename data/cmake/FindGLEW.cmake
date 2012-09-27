# - Find GLEW
# Find the native GLEW includes and libraries.
# This module defines the following read-only variables:
#  GLEW_INCLUDE_DIRS - where to find GL/glew.h, etc.
#  GLEW_LIBRARIES    - libraries to link against to use GLEW.
#  GLEW_DEFINITIONS  - compiler definitions necessary for using GLEW.
#  GLEW_FOUND        - if false, do not try to use GLEW.
#
# These variables alter the behavior of the module when defined before calling
# find_package(GLEW)
#  GLEW_MX       - Set to a value which evaluates to true to look for the
#                  multi-context version of GLEW
#  GLEW_ROOT_DIR - Base location of the GLEW (e.g. where the files were
#                  unzipped.)
# 
#=============================================================================
# Originally from:
# http://code.google.com/p/nvidia-texture-tools/source/browse/trunk/cmake/FindGLEW.cmake
#=============================================================================

# Additional modules
include(FindReleaseAndDebug)
include(FindPackageHandleStandardArgs)
include(CheckCSourceCompiles)

# First of all, GLEW depends on OpenGL
find_package(OpenGL REQUIRED)

if (NOT WIN32)
  set (GLEW_generic_include_path "/usr/include" "/usr/local/include"
    "/sw/include" "/opt/local/include")
  set (GLEW_generic_lib_path
    "/usr/lib" "/usr/local/lib" "/sw/lib" "/opt/local/lib")
else()
  set (GLEW_generic_include_path "")
  set (GLEW_generic_lib_path "")
endif()


# Build the names for the library
if (NOT DEFINED CMAKE_C_COMPILER_ID)
  message (AUTHOR_WARNING
    "The variable to check for the compiler ID has changed!")
endif()
if (MSVC OR (WIN32 AND CMAKE_C_COMPILER_ID MATCHES "Intel"))
  set (GLEW_NAMES "glew32" "glew")
else()
  set (GLEW_NAMES "GLEW" "glew")
endif()

# Assumes that GLEW_MX means that the library has the "mx" suffix
if (GLEW_MX)
  set (GLEW_NAMES_MX "")
  foreach(name ${GLEW_NAMES})
    set(GLEW_NAMES_MX ${GLEW_NAMES_MX} "${name}mx")
  endforeach()
  set (GLEW_NAMES ${GLEW_NAMES_MX})
  unset (GLEW_NAMES_MX)
endif()


# Finds the include files directory
if(NOT APPLE)
  set(GLEW_BASE_DIR "GL")
else()
  set(GLEW_BASE_DIR "OpenGL")
endif()
find_path(GLEW_INCLUDE_DIR ${GLEW_BASE_DIR}/glew.h
  HINTS ${GLEW_ROOT_DIR}/include
  PATHS ${GLEW_generic_include_path}
  DOC "The directory where ${GLEW_BASE_DIR}/glew.h resides"
)
unset(GLEW_BASE_DIR)
if (GLEW_INCLUDE_DIR)
  mark_as_advanced (GLEW_INCLUDE_DIR)
endif()

# Tries to find the required libraries
FIND_RELEASE_AND_DEBUG(GLEW NAMES ${GLEW_NAMES} DEFAULT_SUFFIXES
  PATHS ${GLEW_ROOT_DIR}/lib ${GLEW_generic_lib_path})
  
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLEW
  DEFAULT_MSG GLEW_INCLUDE_DIR GLEW_LIBRARY)
  
if (GLEW_FOUND)
  set (GLEW_LIBRARIES ${OPENGL_LIBRARIES} ${GLEW_LIBRARY})
  set (GLEW_INCLUDE_DIRS ${GLEW_INCLUDE_DIR})
  if (GLEW_MX)
    set (GLEW_DEFINITIONS -DGLEW_MX)
  else()
    set (GLEW_DEFINITIONS "")
  endif()
endif()
  
# On Windows, try to check if the library is static
if (GLEW_FOUND AND WIN32)
  set (CMAKE_REQUIRED_DEFINITIONS_OLD ${CMAKE_REQUIRED_DEFINITIONS})
  set (CMAKE_REQUIRED_INCLUDES_OLD    ${CMAKE_REQUIRED_INCLUDES})
  set (CMAKE_REQUIRED_LIBRARIES_OLD   ${CMAKE_REQUIRED_LIBRARIES})
  set (CMAKE_REQUIRED_DEFINITIONS ${GLEW_DEFINITIONS})
  set (CMAKE_REQUIRED_INCLUDES    ${GLEW_INCLUDE_DIRS})
  set (CMAKE_REQUIRED_LIBRARIES   ${GLEW_LIBRARIES})
  CHECK_C_SOURCE_COMPILES("
#include <GL/glew.h>
#include <stdio.h>
int main(int argc, char **argv) {
    puts(glewGetString(GLEW_VERSION));
    return 0;
}
" GLEW_WIN_CFLAGS)
  # If the test failed, check if it is static
  if (NOT GLEW_WIN_CFLAGS)
  CHECK_C_SOURCE_COMPILES("
#define GLEW_STATIC
#include <GL/glew.h>
#include <stdio.h>
int main(int argc, char **argv) {
    puts(glewGetString(GLEW_VERSION));
    return 0;
}
" GLEW_WIN_CFLAGS_STATIC)
    if (GLEW_WIN_CFLAGS_STATIC)
      set (GLEW_DEFINITIONS ${GLEW_DEFINITIONS} -DGLEW_STATIC)
    else()
      message(WARNING "Could not determine the compile flags for GLEW.")
    endif()
  endif()

  set (CMAKE_REQUIRED_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS_OLD})
  set (CMAKE_REQUIRED_INCLUDES    ${CMAKE_REQUIRED_INCLUDES_OLD})
  set (CMAKE_REQUIRED_LIBRARIES   ${CMAKE_REQUIRED_LIBRARIES_OLD})
  unset (CMAKE_REQUIRED_DEFINITIONS_OLD)
  unset (CMAKE_REQUIRED_INCLUDES_OLD)
  unset (CMAKE_REQUIRED_LIBRARIES_OLD)
endif()
