###############################################################################
#                EXTERNAL LIBRARIES DETECTION                                 #
###############################################################################

if (NOT DEFINED MTS_VERSION)
  message(FATAL_ERROR "This file has to be included from the main build file.")
endif()

# Set up CMake to use the Mitsuba bundled libraries. Set the variable
# "MTS_NO_DEPENDENCIES" to a value which evaluates to TRUE to avoid
# using the Mitsuba dependencies even if they are present.
set(MTS_DEPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dependencies")
if((MSVC OR APPLE) AND NOT MTS_NO_DEPENDENCIES AND
   IS_DIRECTORY "${MTS_DEPS_DIR}")
  set(MTS_DEPENDENCIES ON)
  set(CMAKE_PROGRAM_PATH "${MTS_DEPS_DIR}/bin")
  set(CMAKE_INCLUDE_PATH "${MTS_DEPS_DIR}/include")
  set(Boost_NO_SYSTEM_PATHS TRUE)
  
  if (MSVC)
    if(CMAKE_CL_64)
      set(MTS_ARCH "x64")
    else()
      set(MTS_ARCH "i386")
    endif()
    list(APPEND CMAKE_INCLUDE_PATH "${MTS_DEPS_DIR}/qt/include")
    set(CMAKE_LIBRARY_PATH "${MTS_DEPS_DIR}/lib/${MTS_ARCH}/"
      "${MTS_DEPS_DIR}/qt/${MTS_ARCH}/lib/")
    set(QT_BINARY_DIR "${MTS_DEPS_DIR}/qt/${MTS_ARCH}/bin")
  elseif(APPLE)
    set(CMAKE_LIBRARY_PATH   "${MTS_DEPS_DIR}/lib")
    set(CMAKE_FRAMEWORK_PATH "${MTS_DEPS_DIR}/frameworks")
    set(QT_BINARY_DIR        "${MTS_DEPS_DIR}/bin")
    # Create a shell script to set the paths for dyld
    file(WRITE "${PROJECT_BINARY_DIR}/binaries/mitsuba_dyld.sh"
"#!/bin/sh
# DYLD paths for the mitsuba dependencies. Created automatically by CMake.
export DYLD_FALLBACK_FRAMEWORK_PATH=\"${MTS_DEPS_DIR}/frameworks\":$DYLD_FALLBACK_FRAMEWORK_PATH
export DYLD_FALLBACK_LIBRARY_PATH=\"${MTS_DEPS_DIR}/lib\":$DYLD_FALLBACK_LIBRARY_PATH
")
  endif()
else()
  set(MTS_DEPENDENCIES OFF)
  unset(MTS_DEPS_DIR)
endif()


# Qt4 (optional)
find_package(Qt4 4.7 COMPONENTS
  QtCore QtGui QtXml QtXmlPatterns QtNetwork QtOpenGL)
CMAKE_DEPENDENT_OPTION(BUILD_GUI "Built the Qt4-based mitsuba GUI." ON
  "QT4_FOUND" OFF)
  
# System threading library, used for custom options
set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)


###########################################################################
# Boost 
find_package(Boost 1.44 REQUIRED COMPONENTS "filesystem" "system" "thread")
# As of CMake 2.8.2, FindBoost doesn't honor the "REQUIRED" flag
if (NOT Boost_FOUND)
  set(BOOST_ROOT "" CACHE PATH
    "Preferred installation prefix for searching for Boost.")
  message(FATAL_ERROR
    "Boost is missing. The required modules are math, filesystem and system.")
endif()
mark_as_advanced(Boost_LIB_DIAGNOSTIC_DEFINITIONS)

# Check if spirit works: the version of Clang in Ubuntu 11.04 does not support
# the system version of Boost Spirit
set(CMAKE_REQUIRED_INCLUDES ${Boost_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${Boost_LIBRARIES})
CHECK_CXX_SOURCE_COMPILES("
#include <boost/spirit/include/qi.hpp>
int main (int argc, char **argv) {
    return 0;
}
" BOOST_SPIRIT_WORKS)


# Try to figure out if this boost distro has Boost::python. If we include
# python in the main boost components list above, CMake will abort if it
# is not found. So we resort to checking for the boost_python library's
# existance to get a soft failure
if (APPLE AND MTS_DEPENDENCIES)
  set(mts_boost_python_names boost_python boost_python26 boost_python27
    boost_python32 boost_python33)
else()
  set(mts_boost_python_names boost_python)
endif()
find_library (mts_boost_python_lib NAMES ${mts_boost_python_names}
              HINTS ${Boost_LIBRARY_DIRS} NO_DEFAULT_PATH)
mark_as_advanced (mts_boost_python_lib)
if (NOT mts_boost_python_lib AND Boost_SYSTEM_LIBRARY_RELEASE)
    get_filename_component (mts_boost_PYTHON_rel
                            ${Boost_SYSTEM_LIBRARY_RELEASE} NAME
                           )
    string (REGEX REPLACE "^(lib)?(.+)_system(.+)$" "\\2_python\\3"
            mts_boost_PYTHON_rel ${mts_boost_PYTHON_rel}
           )
    find_library (mts_boost_PYTHON_LIBRARY_RELEASE
                  NAMES ${mts_boost_PYTHON_rel} lib${mts_boost_PYTHON_rel}
                  HINTS ${Boost_LIBRARY_DIRS}
                  NO_DEFAULT_PATH
                 )
    mark_as_advanced (mts_boost_PYTHON_LIBRARY_RELEASE)
endif ()
if (NOT mts_boost_python_lib AND Boost_SYSTEM_LIBRARY_DEBUG)
    get_filename_component (mts_boost_PYTHON_dbg
                            ${Boost_SYSTEM_LIBRARY_DEBUG} NAME
                           )
    string (REGEX REPLACE "^(lib)?(.+)_system(.+)$" "\\2_python\\3"
            mts_boost_PYTHON_dbg ${mts_boost_PYTHON_dbg}
           )
    find_library (mts_boost_PYTHON_LIBRARY_DEBUG
                  NAMES ${mts_boost_PYTHON_dbg} lib${mts_boost_PYTHON_dbg}
                  HINTS ${Boost_LIBRARY_DIRS}
                  NO_DEFAULT_PATH
                 )
    mark_as_advanced (mts_boost_PYTHON_LIBRARY_DEBUG)
endif ()
if (mts_boost_python_lib OR
    mts_boost_PYTHON_LIBRARY_RELEASE OR mts_boost_PYTHON_LIBRARY_DEBUG)
    set (mts_boost_PYTHON_FOUND ON)
else ()
    set (mts_boost_PYTHON_FOUND OFF)
endif ()

###########################################################################

find_package(Eigen 3.0 REQUIRED)

find_package(JPEG 6 REQUIRED)
find_package(ZLIB 1.2 REQUIRED)
find_package(PNG 1.2 REQUIRED)
add_definitions(${PNG_DEFINITIONS})

find_package(IlmBase)
find_package(OpenEXR)
if (OPENEXR_FOUND AND WIN32)
  set(CMAKE_REQUIRED_INCLUDES ${ILMBASE_INCLUDE_DIRS} ${OPENEXR_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${ILMBASE_LIBRARIES} ${OPENEXR_LIBRARIES})

  CHECK_CXX_SOURCE_COMPILES("
#define OPENEXR_DLL
#include <OpenEXR/half.h>
#include <OpenEXR/ImfRgbaFile.h>
int main(int argc, char **argv) {
    half x = 1.5f;
    Imf::RgbaInputFile file(static_cast<const char*>(0));
    file.readPixels(0,0);
    return x > 0 ? 0 : 1;
}
" OPENEXR_IS_DLL)

  unset (CMAKE_REQUIRED_INCLUDES)
  unset (CMAKE_REQUIRED_LIBRARIES)
  
  if (OPENEXR_IS_DLL)
    add_definitions(-DOPENEXR_DLL)
  endif()
endif()

# XERCES_ROOT_DIR
find_package(Xerces 3.0 REQUIRED)


# ColladaDOM (optional)
find_package(COLLADA)
if (COLLADA_FOUND)
  add_definitions(-DMTS_HAS_COLLADA=1)
endif()

find_package(OpenGL REQUIRED)
set (GLEW_MX ON)
find_package(GLEW REQUIRED)
if (GLEW_FOUND)
  set (GLEW_STATE_VARS ${GLEW_INCLUDE_DIRS} ${GLEW_LIBRARIES})
  if (NOT GLEW_TEST_STATE)
    set (GLEW_TEST_STATE "${GLEW_STATE_VARS}" CACHE INTERNAL "GLEW State")
  endif ()
  if (NOT GLEW_TEST_STATE STREQUAL "${GLEW_STATE_VARS}")
    set (GLEW_TEST_STATE "${GLEW_STATE_VARS}" CACHE INTERNAL "GLEW State" FORCE)
    unset (GLEW_VERSION_IS_OK CACHE)
  endif ()
  set (CMAKE_REQUIRED_INCLUDES  ${GLEW_INCLUDE_DIRS})
  set (CMAKE_REQUIRED_LIBRARIES ${GLEW_LIBRARIES})
  CHECK_CXX_SOURCE_COMPILES("
#if defined(__APPLE__)
#include <OpenGL/glew.h>
#else
#include <GL/glew.h>
#endif
int main (int argc, char **argv) {
    int i = GL_VERTEX_ATTRIB_ARRAY_UNIFIED_NV;
    return 0;
}
" GLEW_VERSION_IS_OK)
  if (NOT GLEW_VERSION_IS_OK)
	  #    message (SEND_ERROR "The version of GLEW seems to be outdated!") XXX
  endif ()
endif ()

# Try to get OpenMP support
find_package(OpenMP)
CMAKE_DEPENDENT_OPTION(MTS_OPENMP "Enable OpenMP support" ON
  "OPENMP_FOUND" OFF)
if (MTS_OPENMP)
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else ()
  add_definitions (-DMTS_NO_OPENMP)
endif()

# Linux requires X11
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  find_package(X11 REQUIRED)
  if (NOT X11_xf86vmode_FOUND)
    message(FATAL_ERROR "X11 vmode was not found.")
  endif()
endif()

# Mac OS X Frameworks
if (APPLE)
  find_library(COCOA_LIBRARY Cocoa)
  find_library(BWTOOLKIT_LIBRARY BWToolkitFramework)
  find_path(BWTOOLKIT_INCLUDE_DIR BWToolkitFramework/BWToolkitFramework.h)
  find_library (SECURITY_LIBRARY Security 
    NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH)
  find_path(SECURITY_INCLUDE_DIR "Authorization.h"
    NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH)
  mark_as_advanced (COCOA_LIBRARY)
  mark_as_advanced (BWTOOLKIT_LIBRARY BWTOOLKIT_INCLUDE_DIR)
  mark_as_advanced (SECURITY_LIBRARY  SECURITY_INCLUDE_DIR)
endif()


# The Python libraries.
find_package (PythonLibs "2.6")
CMAKE_DEPENDENT_OPTION(BUILD_PYTHON "Build the Python bindings." ON
  "PYTHONLIBS_FOUND;mts_boost_PYTHON_FOUND" OFF)
if (PYTHONLIBS_FOUND AND mts_boost_PYTHON_FOUND)
  set (PYTHON_FOUND TRUE)
endif ()


# Includes for the common libraries
include_directories(${Boost_INCLUDE_DIRS} ${Eigen_INCLUDE_DIR})

# If we are using the system OpenEXR, add its headers which half.h requires
if (OPENEXR_FOUND)
  include_directories(${ILMBASE_INCLUDE_DIRS})
endif()

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
