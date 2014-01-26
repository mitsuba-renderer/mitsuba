###############################################################################
#         CONFIGURATION AND DEFAULT DEFINITIONS & INCLUDES                    #
###############################################################################

if (NOT DEFINED MTS_VERSION)
  message(FATAL_ERROR "This file has to be included from the main build file.")
endif()

# Default initial compiler flags which may be modified by advanced users
if (MTS_CMAKE_INIT)
  set(MTS_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(MTS_CXX_FLAGS "-fvisibility=hidden -pipe -march=nocona -ffast-math -Wall -Winvalid-pch")
  endif()
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(MTS_CXX_FLAGS "${MTS_CXX_FLAGS} -mfpmath=sse")
  endif()
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(MTS_CXX_FLAGS "${MTS_CXX_FLAGS} -ftemplate-depth=512")
  endif()
  if (MTS_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${MTS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}" CACHE
      STRING "Flags used by the compiler during all build types." FORCE)
    set(MTS_CXX_FLAGS)
  endif()
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
if (MSVC AND MTS_SSE AND NOT CMAKE_CL_64)
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} /arch:SSE2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
endif()
