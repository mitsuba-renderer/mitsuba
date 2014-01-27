# - Find COLLADA
# This module searches for the COLLADA library, by default using version 1.4
# of the COLLADA DOM Schema.
# 
# The module defines the following variables:
#  COLLADA_INCLUDE_DIRS - where to find dae.h, etc.
#  COLLADA_LIBRARIES    - the libraries needed to use COLLADA.
#  COLLADA_DEFINITIONS  - preprocessor definitions to use with COLLADA.
#  COLLADA_NAMESPACE    - boolean to indicate whether this version contains the
#                         namespaces and new functions introduced with
#                         COLLADA-DOM 2.4.
# COLLADA_NAMESPACE_141 - boolean to indicate whether the namespace
#                         ColladaDOM141 is supported.
# COLLADA_NAMESPACE_150 - boolean to indicate whether the namespace
#                         ColladaDOM150 is supported.
#  COLLADA_FOUND        - if false, do not try to use COLLADA.
#
# Variables used by this module, they can change the default behavior and need
# to be set before calling find_package:
#
#  COLLADADOM_15    - if set to a value which evaluates to true, the module
#                     will look for the COLLADA DOM Schema version 1.5
#                     components instead of the default (1.4.)
#  COLLADA_ROOT_DIR - The preferred installation prefix for searching for
#                     COLLADA. This corresponds to
#                     ./configure --prefix=$COLLADA_ROOT_DIR
#                     Set if this module has problems finding the proper
#                     COLLADA installation.
#

# ============================================================================
# Originally created by Robert Osfield. [OpenSceneGraph]
# ============================================================================


# Standarnd issue macros
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)
include (FindReleaseAndDebug)
include (CheckCXXSourceCompiles)

if (COLLADA_15)
  set (COLLADADOM_VERSION    "15")
  set (COLLADADOM_VERSION_PT "1.5")
else()
  set (COLLADADOM_VERSION    "14")
  set (COLLADADOM_VERSION_PT "1.4")
endif()
set(COLLADA_VERSIONS 24 "2.4" 23 22 21 20 2)


# Default search paths
set (COLLADA_generic_include_path
  ~/Library/Frameworks
  /Library/Frameworks
  /opt/local/Library/Frameworks #macports
  /usr/local/include
  /usr/local/include/colladadom
  /usr/include/
  /usr/include/colladadom
  /sw/include # Fink
  /opt/local/include # DarwinPorts
  /opt/csw/include # Blastwave
  /opt/include
  /usr/freeware/include
)
set (COLLADA_generic_lib_path
  /opt/local/Library/Frameworks #macports
  /usr/local/lib
  /usr/lib
  /sw/lib
  /opt/local/lib
  /opt/csw/lib
  /opt/lib
  /usr/freeware/lib
)


# Macro to assemble a helper state variable
macro (SET_STATE_VAR varname)
  set (tmp_lst ColladaDOM | ${COLLADA_INCLUDE_DIR} | ${COLLADA_LIBRARY})
  set (${varname} "${tmp_lst}")
  unset (tmp_lst)
endmacro ()

# Macro to search for an include directory
macro (PREFIX_FIND_INCLUDE_DIR prefix includefile libpath_var)
  string (TOUPPER ${prefix}_INCLUDE_DIR tmp_varname)
  find_path(${tmp_varname} ${includefile}
    HINTS ${${libpath_var}}
    PATHS ${COLLADA_generic_include_path}
    PATH_SUFFIXES
      "collada-dom2.4" "collada_dom2.4"
	    "include/collada-dom2.4" "include/collada_dom2.4"
      "collada-dom" "collada_dom" "include"
	    "include/collada-dom" "include/collada_dom"
  )
  if (${tmp_varname})
    mark_as_advanced (${tmp_varname})
  endif ()
  unset (tmp_varname)
endmacro ()

# Macro to test for namespace support. Factorized just to keep the code cleaner
macro (COLLADA_CHECK_NAMESPACE)
  
  # COLLADA-DOM 2.4 needs special treatment
  set(_COLLADA_OLD_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS})
  set(_COLLADA_OLD_INCLUDES    ${CMAKE_REQUIRED_INCLUDES})
  set(_COLLADA_OLD_LIBRARIES   ${CMAKE_REQUIRED_LIBRARIES})

  # Flags for the test source:
  #  (1 << 0) - COLLADA-DOM 1.4.1
  #  (1 << 1) - COLLADA-DOM 1.5.0
  set(_COLLADA_TEST_SRC "
#define NO_BOOST
#ifdef COLLADA_DOM_SUPPORT150
# undef COLLADA_DOM_SUPPORT150
#endif
#ifdef COLLADA_DOM_SUPPORT141
# undef COLLADA_DOM_SUPPORT141
#endif
#ifdef COLLADA_DOM_NAMESPACE
# undef COLLADA_DOM_NAMESPACE
#endif
#ifndef COLLADA_FLAGS
# define COLLADA_FLAGS 0
#endif

#define COLLADA_DOM_NAMESPACE
#if (COLLADA_FLAGS & (1 << 0)) != 0
# define COLLADA_DOM_SUPPORT141
#endif
#if (COLLADA_FLAGS & (1 << 1)) != 0
# define COLLADA_DOM_SUPPORT150
#endif
#include <dae.h>
int main(int argc, char** argv) {
    DAE* dae = 0;
    int result = 0;
#ifdef COLLADA_DOM_SUPPORT141
    ColladaDOM141::domCOLLADA* root1 = dae->getRoot141(argv[0]);
    result += ((int*)root1)[argc];
#endif
#ifdef COLLADA_DOM_SUPPORT150
    ColladaDOM150::domCOLLADA* root2 = dae->getRoot150(argv[0]);
    result += ((int*)root2)[argc];
#endif
    return result;
}
  ")
  
  # Run the tests
  set(CMAKE_REQUIRED_INCLUDES  ${COLLADA_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${COLLADA_LIBRARIES})
  
  set(CMAKE_REQUIRED_DEFINITIONS "-DCOLLADA_FLAGS=3")
  CHECK_CXX_SOURCE_COMPILES("${_COLLADA_TEST_SRC}" HAVE_COLLADA_DOM_141_150)
  # Perhaps it was compiled with only one of the namespaces
  if (NOT HAVE_COLLADA_DOM_141_150)
    set(CMAKE_REQUIRED_DEFINITIONS "-DCOLLADA_FLAGS=1")
    CHECK_CXX_SOURCE_COMPILES("${_COLLADA_TEST_SRC}" HAVE_COLLADA_DOM_141)
    set(CMAKE_REQUIRED_DEFINITIONS "-DCOLLADA_FLAGS=2")
    CHECK_CXX_SOURCE_COMPILES("${_COLLADA_TEST_SRC}" HAVE_COLLADA_DOM_150)
  endif()
  
  set(CMAKE_REQUIRED_DEFINITIONS ${_COLLADA_OLD_DEFINITIONS})
  set(CMAKE_REQUIRED_INCLUDES    ${_COLLADA_OLD_INCLUDES})
  set(CMAKE_REQUIRED_LIBRARIES   ${_COLLADA_OLD_LIBRARIES})
  unset(_COLLADA_TEST_SRC)
  unset(_COLLADA_OLD_INCLUDES)
  unset(_COLLADA_OLD_LIBRARIES)
endmacro()


# Encode the current state of the external variables into a string
SET_STATE_VAR (COLLADA_CURRENT_STATE)

# If the state has changed, clear the cached variables
if (COLLADA_CACHED_STATE AND
    NOT COLLADA_CACHED_STATE STREQUAL COLLADA_CURRENT_STATE)
  foreach (libvar ${COLLADA_CACHED_VARS})
    unset (${libvar} CACHE)
  endforeach ()
endif ()


# Where to look for the libraries
if(APPLE)
  set(COLLADA_BUILDNAME "mac")
elseif(MINGW)
  set(COLLADA_BUILDNAME "mingw")
elseif(MSVC)
  math(EXPR COLLADA_BUILDNAME "(${MSVC_VERSION} - 600) / 100")
  set(COLLADA_BUILDNAME "vc${COLLADA_BUILDNAME}")
endif()

set(Collada_library_paths
  ${COLLADA_ROOT_DIR}
  ${COLLADA_ROOT_DIR}/lib
  ${COLLADA_ROOT_DIR}/build/${COLLADA_BUILDNAME}-${COLLADADOM_VERSION_PT}
)

# Construct the possible names for the DOM library
set(COLLADADOM_NAMES
  collada${COLLADADOM_VERSION}dom
  Collada${COLLADADOM_VERSION}Dom
  collada-dom
)
if (MSVC)
  list(APPEND COLLADADOM_NAMES
    libCollada${COLLADADOM_VERSION}Dom
    libcollada${COLLADADOM_VERSION}dom
    libcollada-dom
  )
endif()

# Version suffixes for MSVC
if (MSVC)
  math(EXPR VC_SUFFIX_A "(${MSVC_VERSION} - 600) / 10")
  math(EXPR VC_SUFFIX_B "${VC_SUFFIX_A} / 10")
  set(VC_SUFFIXES "vc${VC_SUFFIX_A}-mt" "vc${VC_SUFFIX_B}-mt"
                  "vc${VC_SUFFIX_A}"    "vc${VC_SUFFIX_B}")
  unset(VC_SUFFIX_A)
  unset(VC_SUFFIX_B)
endif()
set(COLLADADOM_NAMES_BASE ${COLLADADOM_NAMES})
foreach(name ${COLLADADOM_NAMES_BASE})
  foreach(version ${COLLADA_VERSIONS})
    list(APPEND COLLADADOM_NAMES "${name}${version}" "${name}${version}-dp")
    if(MSVC)
      foreach(vc_suffix ${VC_SUFFIXES})
        list(APPEND COLLADADOM_NAMES "${name}${version}-${vc_suffix}"
          "${name}${version}-dp-${vc_suffix}")
      endforeach()
    endif()
  endforeach()
endforeach()
unset(COLLADADOM_NAMES_BASE)
unset(VC_SUFFIXES)

# Locate the header files
PREFIX_FIND_INCLUDE_DIR (COLLADA dae.h COLLADA_ROOT_DIR)

# Locate the actual library
FIND_RELEASE_AND_DEBUG(COLLADA NAMES ${COLLADADOM_NAMES} DEFAULT_SUFFIXES
  PATHS ${Collada_library_paths} ${COLLADA_generic_lib_path})
  
# Create the list of variables that might need to be cleared.
# The libraries and include path are not cleared since they might be manually
# set by the user. Besides they will be processed again if the user deletes
# them from the CMake cache. Only the internal cache variables need to be
# specifically cleared so that the tests run again.
set (COLLADA_CACHED_VARS
  HAVE_COLLADA_DOM_141_150 HAVE_COLLADA_DOM_141 HAVE_COLLADA_DOM_150
  CACHE INTERNAL "Variables set by FindCOLLADA.cmake" FORCE)

# Store the current state so that variables might be cleared if required
set (COLLADA_CACHED_STATE ${COLLADA_CURRENT_STATE}
  CACHE INTERNAL "State last seen by FindCOLLADA.cmake" FORCE)
  
# Use the standard function to handle COLLADA_FOUND
FIND_PACKAGE_HANDLE_STANDARD_ARGS (COLLADA DEFAULT_MSG
  COLLADA_INCLUDE_DIR COLLADA_LIBRARY)
  
# Set the uncached variables for the appropriate version
if (COLLADA_FOUND)
  set(COLLADA_INCLUDE_DIRS
    ${COLLADA_INCLUDE_DIR}
    ${COLLADA_INCLUDE_DIR}/${COLLADADOM_VERSION_PT})
  set(COLLADA_LIBRARIES ${COLLADA_LIBRARY})
  
  # Check if this version of COLLADA-DOM contains namespaces introduced in v2.4
  COLLADA_CHECK_NAMESPACE()
  
  # Set the results and definitions
  if (HAVE_COLLADA_DOM_141_150 OR HAVE_COLLADA_DOM_141 OR HAVE_COLLADA_DOM_150)
    set(COLLADA_NAMESPACE ON)
    # TODO Check if double support was actually enabled
    set(COLLADA_DEFINITIONS "-DCOLLADA_DOM_NAMESPACE -DCOLLADA_DOM_DAEFLOAT_IS64")
    if (HAVE_COLLADA_DOM_141_150 OR HAVE_COLLADA_DOM_141)
      set(COLLADA_NAMESPACE_141 ON)
      set(COLLADA_DEFINITIONS "${COLLADA_DEFINITIONS} -DCOLLADA_DOM_SUPPORT141")
    endif()
    if (HAVE_COLLADA_DOM_141_150 OR HAVE_COLLADA_DOM_150)
      set(COLLADA_NAMESPACE_150 ON)
      set(COLLADA_DEFINITIONS "${COLLADA_DEFINITIONS} -DCOLLADA_DOM_SUPPORT150")
    endif()
  else()
    set(COLLADA_NAMESPACE OFF)
    set(COLLADA_NAMESPACE_141 OFF)
    set(COLLADA_NAMESPACE_150 OFF)
    set(COLLADA_DEFINITIONS "")
  endif()
endif()
