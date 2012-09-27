# Experimental CMake file for Mitusba
# Tested only on Windows and Linux (Ubuntu 10.10)
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

# Internal variable to know whether this is the first time CMake runs
if (NOT DEFINED MTS_CMAKE_INIT)
  set(MTS_CMAKE_INIT ON CACHE INTERNAL "Is this the initial CMake run?")
else()
  set(MTS_CMAKE_INIT OFF CACHE INTERNAL "Is this the initial CMake run?")
endif()

# Allow to override the default project name "mitsuba"
if (NOT DEFINED MTS_PROJECT_NAME)
  set(MTS_PROJECT_NAME "mitsuba")
endif()
project(${MTS_PROJECT_NAME})

# Tell cmake where to find the additional modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/cmake")
# Make sure the cmake-provided modules use the versions they expect
if(NOT CMAKE_VERSION VERSION_LESS "2.8.4")
  cmake_policy(SET CMP0017 NEW)
endif()

# Enable folders for projects in Visual Studio
if (CMAKE_GENERATOR MATCHES "Visual Studio")
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

# Set CMAKE_BUILD_TYPE to Release by default
if (MTS_CMAKE_INIT AND DEFINED CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug, Release, RelWithDebInfo, MinSizeRel." FORCE)
endif()

# Load the required modules
include (MitsubaUtil)
include (CheckCXXSourceCompiles)
include (CMakeDependentOption)

# Read version (MTS_VERSION) from include/mitsuba/core/version.h
file(STRINGS "include/mitsuba/core/version.h" MITSUBA_H REGEX "^#define MTS_VERSION \"[^\"]*\"$")
string(REGEX REPLACE "^.*MTS_VERSION \"([0-9]+).*$" "\\1" MTS_VERSION_MAJOR "${MITSUBA_H}")
string(REGEX REPLACE "^.*MTS_VERSION \"[0-9]+\\.([0-9]+).*$" "\\1" MTS_VERSION_MINOR  "${MITSUBA_H}")
string(REGEX REPLACE "^.*MTS_VERSION \"[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" MTS_VERSION_PATCH "${MITSUBA_H}")
set(MTS_VERSION "${MTS_VERSION_MAJOR}.${MTS_VERSION_MINOR}.${MTS_VERSION_PATCH}")
set(MITSUBA_H)
if("${MTS_VERSION_MAJOR}" MATCHES "[0-9]+" AND
   "${MTS_VERSION_MINOR}" MATCHES "[0-9]+" AND
   "${MTS_VERSION_PATCH}" MATCHES "[0-9]+")
  message(STATUS "mitsuba ${MTS_VERSION}")
else()
  message(FATAL_ERROR "The mitsuba version could not be determined!")
endif()


########################### External libraries ################################
include (MitsubaExternal)


###############################################################################
#         CONFIGURATION AND DEFAULT DEFINITIONS & INCLUDES                    #
###############################################################################

# Image format definitions
if (PNG_FOUND)
  add_definitions(-DMTS_HAS_LIBPNG=1)
endif()
if (JPEG_FOUND)
  add_definitions(-DMTS_HAS_LIBJPEG=1)
endif()
if (OPENEXR_FOUND)
  add_definitions(-DMTS_HAS_OPENEXR=1)
endif()

# Top level configuration definitions
option(MTS_DEBUG "Enable assertions etc. Usually a good idea." ON)
if (MTS_DEBUG)
  add_definitions(-DMTS_DEBUG)
endif()
option(MTS_KD_DEBUG "Enable additional checks in the kd-Tree.
This is quite slow and mainly useful to track down bugs when they are suspected."
OFF)
if (MTS_KD_DEBUG)
  add_definitions(-DMTS_KD_DEBUG)
endif()
option(MTS_KD_CONSERVE_MEMORY
  "Use less memory for storing geometry (at the cost of speed)." OFF)
if (MTS_KD_CONSERVE_MEMORY)
  add_definitions(-DMTS_KD_CONSERVE_MEMORY)
endif()
option(MTS_SINGLE_PRECISION
  "Do all computation in single precision. This is usually sufficient." ON)
if (MTS_SINGLE_PRECISION)
  add_definitions(-DSINGLE_PRECISION)
endif()

set(MTS_SPECTRUM_SAMPLES 3 CACHE STRING
  "Number of spectral samples used to render. The default is 3 (RGB-mode).
For high-quality spectral rendering, this should be set to 30 or higher.")
if(NOT "${MTS_SPECTRUM_SAMPLES}" MATCHES "^[1-9][0-9]*$" OR
    MTS_SPECTRUM_SAMPLES LESS 3 OR MTS_SPECTRUM_SAMPLES GREATER 2048)
  message(FATAL_ERROR
    "Invalid number of spectrum samples: ${MTS_SPECTRUM_SAMPLES}. Valid values: [3,2048]")
else()
  add_definitions(-DSPECTRUM_SAMPLES=${MTS_SPECTRUM_SAMPLES})
endif()

CMAKE_DEPENDENT_OPTION (MTS_DOUBLE_PRECISION
  "Do all computation in double precision." ON
  "NOT MTS_SINGLE_PRECISION" OFF)
if (MTS_DOUBLE_PRECISION)
  add_definitions(-DDOUBLE_PRECISION)
endif()

CMAKE_DEPENDENT_OPTION (MTS_SSE "Activate optimized SSE routines." ON
  "NOT MTS_DOUBLE_PRECISION" OFF)
if (MTS_SSE)
  add_definitions(-DMTS_SSE)
endif ()

CMAKE_DEPENDENT_OPTION (MTS_HAS_COHERENT_RT
  "Include coherent ray tracing support." ON
  "MTS_SSE" OFF)
if (MTS_HAS_COHERENT_RT)
  add_definitions(-DMTS_HAS_COHERENT_RT)
endif()

CMAKE_DEPENDENT_OPTION (MTS_DEBUG_FP
  "Generated NaNs will cause floating point exceptions, which can be caught in a debugger (very slow!)" OFF
  "NOT MTS_DOUBLE_PRECISION" OFF)
if (MTS_DEBUG_FP)
  add_definitions(-DMTS_DEBUG_FP)
endif()


# Options to disable MSVC STL debug + security features (slow..!)
if (MSVC OR (WIN32 AND CMAKE_C_COMPILER_ID MATCHES "Intel"))
  # _SECURE_SCL already defaults to 0 in release mode in MSVC 2010
  if(MSVC_VERSION LESS 1600)
    option(MTS_NO_CHECKED_ITERATORS  "Disable checked iterators in MSVC" OFF)
    option(MTS_NO_ITERATOR_DEBUGGING "Disable iterator debugging in MSVC" OFF)
  else()
    set(MTS_NO_CHECKED_ITERATORS  OFF)
    set(MTS_NO_ITERATOR_DEBUGGING OFF)
  endif()
  option(MTS_NO_BUFFER_CHECKS "Disable the buffer security checks in MSVC" ON)
  
  if (MTS_NO_ITERATOR_DEBUGGING)
    add_definitions (-D_HAS_ITERATOR_DEBUGGING=0)
  endif()
  if (MTS_NO_CHECKED_ITERATORS OR MTS_NO_ITERATOR_DEBUGGING)
    add_definitions (-D_SECURE_SCL=0 -D_SCL_SECURE_NO_WARNINGS)
    message (WARNING "The secure iterators were manually disabled. There might be incompatibility problems.")
  endif ()
  if (MTS_NO_BUFFER_CHECKS)
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} /GS-")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /GS-")
  endif()
endif()

# Platform-specific definitions
if (WIN32 AND CMAKE_SIZEOF_VOID_P EQUAL 8)
  add_definitions(-DWIN64)
endif()


# Main mitsuba include directory
include_directories("include")

# Includes for the common libraries
include_directories(${Boost_INCLUDE_DIRS} ${Eigen_INCLUDE_DIR})

# If we are using the system OpenEXR, add its headers which half.h requires
if (OPENEXR_FOUND)
  include_directories(${ILMBASE_INCLUDE_DIRS})
endif()



###############################################################################
#                   CORE MODULES, APPLICATIONS AND PLUGINS                    #
###############################################################################

# ===== Prerequisite resources =====

# Process the XML schemas
add_subdirectory(data/schema)
# Add the IOR database
add_subdirectory(data/ior)
# Microfacet precomputed data
add_subdirectory(data/microfacet)


# ===== Build the support libraries ====

# Core support library
add_subdirectory(src/libcore)
# Rendering-related APIs
add_subdirectory(src/librender)
# Hardware acceleration
add_subdirectory(src/libhw)
# Bidirectional support library
add_subdirectory(src/libbidir)
# Python binding library
if (BUILD_PYTHON)
  add_subdirectory(src/libpython)
elseif(NOT PYTHON_FOUND)
  message(STATUS "Python was not found. The bindings will not be built.")
endif()


# Additional files to add to main executables
if(APPLE)
  set(MTS_DARWIN_STUB "${CMAKE_CURRENT_SOURCE_DIR}/src/mitsuba/darwin_stub.mm")
endif()


# ===== Build the applications =====

# Build the command-line binaries
add_subdirectory(src/mitsuba)

# Build the COLLADA converter
if (COLLADA_FOUND)
  add_subdirectory(src/converter)
else()
  message(STATUS "Collada DOM was not found. The importer will not be built.")
endif()

# Build the Qt-based GUI binaries
if (BUILD_GUI)
  add_subdirectory(src/mtsgui)
elseif(NOT QT4_FOUND)
  message(STATUS "Qt4 was not found. The mitsuba gui will not be built.")
endif()


# ===== Build the plugins =====

# Utilities
add_subdirectory(src/utils)
# Surface scattering models
add_subdirectory(src/bsdfs)
# Phase functions
add_subdirectory(src/phase)
# Intersection shapes
add_subdirectory(src/shapes)
# Sample generators
add_subdirectory(src/samplers)
# Reconstruction filters
add_subdirectory(src/rfilters)
# Film implementations
add_subdirectory(src/films)
# Sensors
add_subdirectory(src/sensors)
# Emitters
add_subdirectory(src/emitters)
# Participating media
add_subdirectory(src/medium)
# Volumetric data sources
add_subdirectory(src/volume)
# Sub-surface integrators
add_subdirectory(src/subsurface)
# Texture types
add_subdirectory(src/textures)
# Integrators
add_subdirectory(src/integrators)
# Testcases
add_subdirectory(src/tests)


# ===== Packaging =====

# Use a subdirectory to enforce that packaging runs after all other targets
add_subdirectory(data/cmake/packaging)
