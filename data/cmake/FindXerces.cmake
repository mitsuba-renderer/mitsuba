################################################################################
#
# CMake script for finding XERCES.
# If the optional XERCES_ROOT_DIR environment variable exists, header files and
# libraries will be searched in the XERCES_ROOT_DIR/include and XERCES_ROOT_DIR/libs
# directories, respectively. Otherwise the default CMake search process will be
# used.
#
# This script creates the following variables:
#  XERCES_FOUND: Boolean that indicates if the package was found
#  XERCES_INCLUDE_DIRS: Paths to the necessary header files
#  XERCES_LIBRARIES: Package libraries
#
# http://svn.mech.kuleuven.be/websvn/orocos/trunk/rtt/config/FindXerces.cmake
################################################################################

include(FindPackageHandleStandardArgs)

# Get hint from environment variable (if any)
if(NOT XERCES_ROOT_DIR AND DEFINED ENV{XERCES_ROOT_DIR})
  set(XERCES_ROOT_DIR "$ENV{XERCES_ROOT_DIR}" CACHE PATH "XERCES base directory location (optional, used for nonstandard installation paths)")
  mark_as_advanced(XERCES_ROOT_DIR)
endif()

# Search path for nonstandard locations
if(XERCES_ROOT_DIR)
  set(XERCES_INCLUDE_PATH PATHS "${XERCES_ROOT_DIR}/include" NO_DEFAULT_PATH)
  set(XERCES_LIBRARY_PATH PATHS "${XERCES_ROOT_DIR}/lib" NO_DEFAULT_PATH)
endif()

# Find headers and libraries
find_path(XERCES_INCLUDE_DIR       NAMES xercesc/dom/DOM.hpp ${XERCES_INCLUDE_PATH})
find_library(XERCES_C_LIBRARY      NAMES xerces-c xerces-c_3 ${XERCES_LIBRARY_PATH})
#find_library(XERCES_DEPDOM_LIBRARY NAMES xerces-depdom       ${XERCES_LIBRARY_PATH})

# Set Xerces_FOUND honoring the QUIET and REQUIRED arguments
find_package_handle_standard_args(Xerces DEFAULT_MSG 
  XERCES_C_LIBRARY 
 #XERCES_DEPDOM_LIBRARY 
  XERCES_INCLUDE_DIR)

# Output variables
if(XERCES_FOUND)
  # Include dirs
  set(XERCES_INCLUDE_DIRS ${XERCES_INCLUDE_DIR})

  # Libraries
  set(XERCES_LIBRARIES ${XERCES_C_LIBRARY}
                       ${XERCES_DEPDOM_LIBRARY})

endif()

# Advanced options for not cluttering the cmake UIs
mark_as_advanced(XERCES_INCLUDE_DIR XERCES_C_LIBRARY XERCES_DEPDOM_LIBRARY)
