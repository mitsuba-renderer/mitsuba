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


# Find the external libraries and setup the paths
include (MitsubaExternal)

# Setup the build options, include paths and compile definitions
include (MitsubaBuildOptions)


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
else()
  set(MTS_DARWIN_STUB "")
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
