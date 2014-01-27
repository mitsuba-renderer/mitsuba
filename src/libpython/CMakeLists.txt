# mitsuba Python bindings
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

# The Mitsuba provided dependencies on OS X and Windows add a version suffix to
# boost-python to provide support for multiple version in the same
# distribution. CMake does not know how to handle this, requiring a little hack
set(mts_python "python")
if((APPLE OR WIN32) AND MTS_DEPENDENCIES)
  # Try to guess the suffix from the library version
  if (PYTHON_LIBRARY MATCHES ".+python([23])\\.?([0-9])\\.[^.]+$")
    set(mts_python "python${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
  endif()
endif()

# Call the FindBoost module again to set up the components. At this point it is
# assumed that the call will succeed
find_package(Boost 1.44 REQUIRED
  COMPONENTS "filesystem" "system" ${mts_python})

include_directories(
  ${PYTHON_INCLUDE_PATH} ${Boost_INCLUDE_DIRS} ${XERCES_INCLUDE_DIRS})

# Headers
set(HDRS
  base.h
)

# Common sources
set(SRCS
  core.cpp
  render.cpp
)

add_definitions(-DMTS_BUILD_MODULE=MTS_MODULE_PYTHON)

# core.cpp and render.cpp are too large
if (MSVC)
  set_property (SOURCE core.cpp render.cpp APPEND PROPERTY COMPILE_FLAGS "/bigobj")
endif ()


add_library (PyMitsuba MODULE ${HDRS} ${SRCS})
target_link_libraries (PyMitsuba mitsuba-core mitsuba-render
  ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${XERCES_LIBRARIES})
  
# Set the output directory using the macro included from MitsubaUtil.cmake
SET_OUTPATH_CFG (PyMitsuba LIBRARY_OUTPUT_DIRECTORY
  "${PROJECT_BINARY_DIR}/binaries/@CFGNAME@/${MTS_PYTHON_DEST}"
)
 
set_target_properties (PyMitsuba PROPERTIES OUTPUT_NAME mitsuba PREFIX "")
if (APPLE)
  # Python seems to only look for plugins that end in .so, not .dylib
  set_target_properties (PyMitsuba PROPERTIES SUFFIX ".so")
elseif (WIN32)
  # On Windows, Python looks for .pyd instead of .dll
  set_target_properties (PyMitsuba PROPERTIES SUFFIX ".pyd")
endif()


install(TARGETS PyMitsuba
  RUNTIME DESTINATION ${MTS_PYTHON_DEST} COMPONENT Runtime
  LIBRARY DESTINATION ${MTS_PYTHON_DEST} COMPONENT Runtime)
