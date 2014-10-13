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

# Remove Debug from CMAKE_CONFIGURATION_TYPES as the dependencies do not contain the necessary debug libraries
if(MSVC AND CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_CONFIGURATION_TYPES Release MinSizeRel RelWithDebInfo)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
    "Remove Debug from available configuration types"
     FORCE)
endif()


# make config 'RelWithDebInfo' debuggable by disabling optimizations
if(MSVC)
  # /MD link against release libs, /Od disable optimizations, /Ob0
  if(CMAKE_C_FLAGS_RELWITHDEBINFO MATCHES "/O2")
    string(REGEX REPLACE "/O2" "/Od" CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  endif()
  if(CMAKE_CXX_FLAGS_RELWITHDEBINFO MATCHES "/O2")
    string(REGEX REPLACE "/O2" "/Od" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  endif()
  if(CMAKE_C_FLAGS_RELWITHDEBINFO MATCHES "/Ob1")
    string(REGEX REPLACE "/Ob1" "/Ob0" CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  endif()
  if(CMAKE_CXX_FLAGS_RELWITHDEBINFO MATCHES "/Ob1")
    string(REGEX REPLACE "/Ob1" "/Ob0" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  endif()
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} /RTC1")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /RTC1")
  list(REMOVE_DUPLICATES CMAKE_C_FLAGS_RELWITHDEBINFO)
  list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_RELWITHDEBINFO)
endif()

# Set CMAKE_BUILD_TYPE to Release by default
if (MTS_CMAKE_INIT AND DEFINED CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug, Release, RelWithDebInfo, MinSizeRel." FORCE)
endif()

# Load the required modules
include (MitsubaUtil)
include (MtsGetVersionInfo)
include (CheckCXXSourceCompiles)
include (CMakeDependentOption)

# Read the version information
MTS_GET_VERSION_INFO()
if (MTS_HAS_VALID_REV)
  message(STATUS "mitsuba ${MTS_VERSION}-hg${MTS_REV_ID} (${MTS_DATE})")
else()
  message(STATUS "mitsuba ${MTS_VERSION} (${MTS_DATE})")
endif()

# Setup the build options
include (MitsubaBuildOptions)

# Find the external libraries and setup the paths
include (MitsubaExternal)

# Main mitsuba include directory
include_directories("include")

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
  set(MTS_WINDOWS_STUB "")
elseif(WIN32)
  set(MTS_DARWIN_STUB "")
  set(MTS_WINDOWS_STUB "${CMAKE_CURRENT_SOURCE_DIR}/data/windows/wmain_stub.cpp")
else()
  set(MTS_DARWIN_STUB "")
  set(MTS_WINDOWS_STUB "")
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
