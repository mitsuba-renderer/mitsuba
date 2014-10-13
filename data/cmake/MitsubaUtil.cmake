# Macros to build the mitsuba targets. They are to be used by the CMake scripts
# only, otherwise they don't make any sense at all.

include (CMakeParseArguments)
include (CMakeDependentOption)
include (PCHTargets)

# Function to check that the assumed configurations exist
function (mts_check_configurations)
  if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(WARNING "The build type is not set. Set the value of CMAKE_BUILD_TYPE.")
    return ()
  endif ()
  if (CMAKE_BUILD_TYPE)
    set (configs ${ARGN})
    list (FIND configs ${CMAKE_BUILD_TYPE} idx)
    if (idx LESS 0)
      message (AUTHOR_WARNING "Unexpected configuration '${CMAKE_BUILD_TYPE}' Check the value of CMAKE_BUILD_TYPE.")
    endif ()
    return ()
  endif ()
  set (configs "")
  foreach (cfg ${CMAKE_CONFIGURATION_TYPES})
    set (configs ${configs} ${cfg})
  endforeach ()
  set (myconfigs "")
  foreach (cfg ${ARGN})
    set (myconfigs ${myconfigs} ${cfg})
    list (FIND configs ${cfg} idx)
    if (idx LESS 0)
      message (AUTHOR_WARNING "The assumed configuration '${cfg}' is not available.")
    endif ()
  endforeach ()
  foreach (cfg ${configs})
    list (FIND myconfigs ${cfg} idx)
    if (idx LESS 0)
      message (AUTHOR_WARNING "Unexpected configuration '${cfg}' found.")
    endif ()
  endforeach ()
endfunction()

# Check the standard configurations
# NO Debug! Does not work under windows as dependencies do not include debug libraries
mts_check_configurations (Release MinSizeRel RelWithDebInfo)


# Option to enable interprocedural optimizations
option(MTS_LTCG "Enable interprocedural optimizations on Release targets" ON)
mark_as_advanced (MTS_LTCG)

# Macro to enable interprocedural optimizations on a target
macro (mts_target_ltcg target)
  if (MTS_LTCG)
    set_target_properties (${target} PROPERTIES
      INTERPROCEDURAL_OPTIMIZATION_RELEASE 1
    )
  endif()
endmacro()


# Macro to enable parallel compilation on MSVC
macro (mts_msvc_mp target)
  if (MSVC_IDE AND MSVC_VERSION GREATER 1400)
    set_property (TARGET ${target} APPEND
      PROPERTY COMPILE_FLAGS " /MP")
  endif ()
endmacro()



# Option for precompiled headers
if (PCH_MSVC OR PCH_GCC)
  set (MTS_PCH_DEFAULT ON)
else ()
  set (MTS_PCH_DEFAULT OFF)
endif ()
CMAKE_DEPENDENT_OPTION (MTS_USE_PCH "Use precompiled headers (PCH)."
  ${MTS_PCH_DEFAULT} "PCH_SUPPORTED" OFF)
CMAKE_DEPENDENT_OPTION (MTS_USE_PCH_ALL_PLUGINS
  "Use PCH for all plugins irrespective of their number of source files."
  OFF "MTS_USE_PCH" OFF)
  
# Default project-wide header to be precompiled
set (MTS_DEFAULT_PCH "${PROJECT_SOURCE_DIR}/data/pch/mitsuba_precompiled.hpp")



# Function to configure the output path according to the configurations
# The output path configure expression, contained in "path_cfgstr" should
# have the placeholder @CFGNAME@ for proper substitution
function (SET_OUTPATH_CFG target_name property_suffix path_cfgstr)
  if (CMAKE_CONFIGURATION_TYPES)
    foreach (CFGNAME ${CMAKE_CONFIGURATION_TYPES})
      string (TOUPPER ${CFGNAME} cfg_upper)
      string (CONFIGURE "${path_cfgstr}" outpath @ONLY)
      set_target_properties (${target_name} PROPERTIES
        ${property_suffix}_${cfg_upper} "${outpath}")
    endforeach ()
  else ()
    set (CFGNAME ".")
    string (CONFIGURE "${path_cfgstr}" outpath @ONLY)
    set_target_properties (${target_name} PROPERTIES
        ${property_suffix} "${outpath}")
  endif ()
endfunction ()



# Function to set up the default windows resource file:
#  target_filename - where to write the configured file
#  name        - base name WITHOUT extension, eg "mitsuba"
#  extension   - file extension including the . ".exe"
#  description - File description to be presented to users
#  ICON <iconfile> - icon to be shown for executables
function(mts_win_resource target_filename name ext description)
  CMAKE_PARSE_ARGUMENTS(_res "" "ICON" "" ${ARGN})
  if (NOT WIN32)
    message(AUTHOR_WARNING "This is not a Windows build!")
  elseif(NOT MTS_VERSION)
    message(AUTHOR_WARNING "The mitsuba version variable is not set")
  endif()
  set(RC_FILE "${PROJECT_SOURCE_DIR}/data/windows/mitsuba_res.rc.in")
  
  #TODO Set up the VS_FF_PRERELEASE and VS_FF_PRIVATEBUILD flags adequately
  string(TOLOWER "${ext}" ext_lower)
  if(ext_lower STREQUAL ".dll")
    set(RC_FILETYPE "VFT_DLL")
  elseif(ext_lower STREQUAL ".exe")
    set(RC_FILETYPE "VFT_APP")
  elseif(ext_lower STREQUAL ".lib")
    set(RC_FILETYPE "VFT_STATIC_LIB")
  else()
    message(AUTHOR_WARNING "Unknown windows file type: ${ext_lower}")
    set(RC_FILETYPE "VFT_UNKNOWN")
  endif()
  
  if(_res_ICON)
    get_filename_component(_res_dir "${target_filename}" PATH)
    file(RELATIVE_PATH RC_ICON "${_res_dir}" "${_res_ICON}")
  else()
    set(RC_ICON "")
  endif()
  
  set(RC_DESCRIPTION "${description}")
  if (MTS_HAS_VALID_REV)
    set(RC_VERSION "${MTS_VERSION}-${MTS_VERSION_BUILD}hg${MTS_REV_ID}")
  else()
    set(RC_VERSION "${MTS_VERSION}")
  endif()
  set(RC_VERSION_COMMA "${MTS_VERSION_MAJOR},${MTS_VERSION_MINOR},${MTS_VERSION_PATCH},0")
  set(RC_FILENAME "${name}${ext}")
  set(RC_NAME "${name}")
  # MTS_DATE has the format YYYY.MM.DD
  string(SUBSTRING "${MTS_DATE}" 0 4 RC_YEAR)
  
  configure_file("${RC_FILE}" "${target_filename}" ESCAPE_QUOTES @ONLY)
endfunction()



# Constant with the bundle name for Mitsuba
set(MTS_BUNDLE_NAME "Mitsuba.app")
set(MTS_BUNDLE_RESOURCES "${MTS_BUNDLE_NAME}/Contents/Resources")

# Flag to use either simple or traditional Unix paths (e.g. share/mitsuba/...)
CMAKE_DEPENDENT_OPTION(MTS_SIMPLE_PATHS
  "Use a simple, Windows-like dir structure instead of the typical Unix one."
  ON "NOT WIN32; NOT APPLE" OFF)


# Constant with the destination for the Python bindings
if (WIN32 OR MTS_SIMPLE_PATHS)
  set(MTS_PYTHON_DEST "python")
elseif (APPLE)
  set(MTS_PYTHON_DEST "${MTS_BUNDLE_NAME}/python")
else()
  set(MTS_PYTHON_DEST "share/mitsuba/python")
endif()


# Constant with the destination for the core libraries
if (WIN32 OR MTS_SIMPLE_PATHS)
  set(MTS_LIB_DEST ".")
elseif (APPLE)
  set(MTS_LIB_DEST "${MTS_BUNDLE_NAME}/Contents/Frameworks")
else ()
  set(MTS_LIB_DEST "lib")
endif ()

# Macro to add a build target for a mitsuba core library.
#
# Usage:
#
# add_mts_corelib (name source1 [source2 ...]
#                   [LINK_LIBRARIES external_lib1 ...] )
#
# The plugin name is taken from the first argument. Additional libraries
# (for example, libpng) may be specified after the optionl LINK_LIBRARIES
# keyword.
#
# Each time this macro adds a target, it adds a new element to the
# variable "MTS_CORELIBS", which will contain all the generated core libs

macro (add_mts_corelib _corelib_name)
  CMAKE_PARSE_ARGUMENTS(_corelib "" "" "LINK_LIBRARIES" ${ARGN})
  set (_corelib_srcs ${_corelib_UNPARSED_ARGUMENTS})
  set(MTS_CORELIBS ${MTS_CORELIBS} ${_corelib_name} PARENT_SCOPE)
  
  if (WIN32)
    set(_corelib_res "${CMAKE_CURRENT_BINARY_DIR}/${_corelib_name}_res.rc")
    set(_corelib_description "Mitsuba core library: ${_corelib_name}")
    mts_win_resource("${_corelib_res}"
      "lib${_corelib_name}" ".dll" "${_corelib_description}")
    list(APPEND _corelib_srcs "${_corelib_res}")
  endif()
  
  if (MTS_USE_PCH)
    pch_add_library (${_corelib_name} SHARED
      PCH_HEADER "${MTS_DEFAULT_PCH}" ${_corelib_srcs})
  else ()
    add_library (${_corelib_name} SHARED ${_corelib_srcs})
  endif ()
  target_link_libraries (${_corelib_name} ${_corelib_LINK_LIBRARIES})
  if (WIN32)
    set_target_properties (${_corelib_name} PROPERTIES 
      PREFIX "lib"
      VERSION "${MTS_VERSION}")
  endif()
  if (WIN32)
    set (_corelib_property_suffix "RUNTIME_OUTPUT_DIRECTORY")
  else ()
    set (_corelib_property_suffix "LIBRARY_OUTPUT_DIRECTORY")
  endif ()
  SET_OUTPATH_CFG (${_corelib_name} ${_corelib_property_suffix}
    "${PROJECT_BINARY_DIR}/binaries/@CFGNAME@/${MTS_LIB_DEST}"
  )
  mts_target_ltcg (${_corelib_name})
  mts_msvc_mp (${_corelib_name})
  install(TARGETS ${_corelib_name}
    RUNTIME DESTINATION "${MTS_LIB_DEST}" COMPONENT Runtime
    LIBRARY DESTINATION "${MTS_LIB_DEST}" COMPONENT Runtime
    ARCHIVE DESTINATION "sdk/lib" COMPONENT Developer
  )
endmacro()



# Constant with the destination for the plugins
if (WIN32 OR MTS_SIMPLE_PATHS)
  set(MTS_PLUGIN_DEST "plugins")
elseif (APPLE)
  set(MTS_PLUGIN_DEST "${MTS_BUNDLE_NAME}/plugins")
else()
  set(MTS_PLUGIN_DEST "share/mitsuba/plugins")
endif()
  
# Macro to add a build target for a mitsuba plugin (based on the OIIO one).
#
# Usage:
#
# add_mts_plugin (name source1 [source2 ...]
#                 [LINK_LIBRARIES external_lib1 ...]
#                 [MTS_HW] [MTS_BIDIR]
#                 [NO_MTS_PCH]
#                 [TYPE plugin_type] )
#
# The plugin name is taken from the first argument and the
# source is automatically linked against mitsuba libraries.  Additional
# libraries (for example, libpng) may be specified after the optional
# LINK_LIBRARIES keyword.  NO_MTS_PCH makes the target not to use the 
# default mitsuba precompiled header (on supported platforms).
#
# By default the plugins are linked against mitsuba-core and mitsuba-render.
# When MTS_HW is set, the plugin will be linked against with mitsuba-hw. When
# MTS_BIDIR is specified, the plugin will also be linked against mitsuba-bidir.
# 
# The plugin type (i.e. camera, bsdf, luminaire) may be specified
# after the TYPE keyword. Currently doing this modifies the IDE project name
# in order to have a nicer organization.

macro (add_mts_plugin _plugin_name)
  CMAKE_PARSE_ARGUMENTS(_plugin "MTS_HW;MTS_BIDIR;NO_MTS_PCH"
    "TYPE" "LINK_LIBRARIES" ${ARGN})
  set (_plugin_srcs ${_plugin_UNPARSED_ARGUMENTS})
  
  if (WIN32)
    set(_plugin_res "${CMAKE_CURRENT_BINARY_DIR}/${_plugin_name}_res.rc")
    if (_plugin_TYPE)
      set(_plugin_description "Mitsuba ${_plugin_TYPE} plugin: ${_plugin_name}")
    else()
      set(_plugin_description "Mitsuba plugin: ${_plugin_name}")
    endif()
    mts_win_resource("${_plugin_res}"
      "${_plugin_name}" ".dll" "${_plugin_description}")
    list(APPEND _plugin_srcs "${_plugin_res}")
  endif()
  
  # Use the PCH only with the plugins with 2 or more source files
  set(_plugin_cxx_count "0")
  foreach(_src ${_plugin_srcs})
    if(_src MATCHES ".+\\.[cC][^.]*$")
      math(EXPR _plugin_cxx_count "${_plugin_cxx_count} + 1")
    endif()
  endforeach()
  
  if (NOT _plugin_NO_MTS_PCH AND MTS_USE_PCH AND
      (MTS_USE_PCH_ALL_PLUGINS OR _plugin_cxx_count GREATER 1))
    pch_add_library (${_plugin_name} MODULE
      PCH_HEADER "${MTS_DEFAULT_PCH}" ${_plugin_srcs})
  else ()
    add_library (${_plugin_name} MODULE ${_plugin_srcs})
  endif ()
  
  set(_plugin_core_libraries "mitsuba-core" "mitsuba-render")
  if (_plugin_MTS_HW)
    list(APPEND _plugin_core_libraries "mitsuba-hw")
  endif()
  if (_plugin_MTS_BIDIR)
    list(APPEND _plugin_core_libraries "mitsuba-bidir")
  endif()
  target_link_libraries (${_plugin_name} 
    ${_plugin_core_libraries} ${_plugin_LINK_LIBRARIES})
  
  set_target_properties (${_plugin_name} PROPERTIES PREFIX "")
  if (APPLE)
    set_target_properties (${_plugin_name} PROPERTIES SUFFIX ".dylib")
  endif ()
  if (WIN32)
    set_target_properties (${_plugin_name} PROPERTIES VERSION "${MTS_VERSION}")
  endif()
  set (_plugin_FOLDER "plugins")
  if (_plugin_TYPE)
    if (CMAKE_GENERATOR MATCHES "Visual Studio")
      set (_plugin_FOLDER "plugins/${_plugin_TYPE}")
    else()
      set_target_properties (${_plugin_name} PROPERTIES
        PROJECT_LABEL "${_plugin_TYPE}-${_plugin_name}")
    endif()
  endif()
  set_target_properties (${_plugin_name} PROPERTIES
    FOLDER ${_plugin_FOLDER})
  unset (_plugin_FOLDER)
  SET_OUTPATH_CFG (${_plugin_name} LIBRARY_OUTPUT_DIRECTORY
    "${PROJECT_BINARY_DIR}/binaries/@CFGNAME@/${MTS_PLUGIN_DEST}"
  )
  mts_target_ltcg (${_plugin_name})
  mts_msvc_mp (${_plugin_name})
  install(TARGETS ${_plugin_name}
    RUNTIME DESTINATION ${MTS_PLUGIN_DEST} COMPONENT Runtime
    LIBRARY DESTINATION ${MTS_PLUGIN_DEST} COMPONENT Runtime)
endmacro ()



# Constant with the executables destination
if (WIN32 OR MTS_SIMPLE_PATHS)
  set (MTS_EXE_DEST ".")
elseif (APPLE)
  set (MTS_EXE_DEST "${MTS_BUNDLE_NAME}/Contents/MacOS")
else()
  set (MTS_EXE_DEST "bin")
endif()

# Macro to add a build target for a mitsuba application.
#
# Usage:
#
# add_mts_exe (name [WIN32] source1 [source2 ...]
#                   [LINK_LIBRARIES external_lib1 ...]
#                   [RES_ICON filename]
#                   [RES_DESCRIPTION "Description string"]
#                   [NO_INSTALL]
#                   [MTS_HW] [MTS_BIDIR]
#                   [NO_MTS_PCH | PCH pch_header] )
#
# The executable name is taken from the first argument. The target gets
# automatically linked against mitusba's core libraries, as defined
# in the "MTS_CORELIBS" variable Additional libraries
# (for example, libpng) may be specified after the optionl LINK_LIBRARIES
# keyword.
#
# By default the executables are linked against mitsuba-core and mitsuba-render.
# When MTS_HW is set, the executable will be linked against with mitsuba-hw.
# When MTS_BIDIR is specified, the executable will also be linked against 
# mitsuba-bidir.
#
# The optional keyword WIN32, if presents, gets passed to add_executable(...)
# to produce a Windows executable using winmain, thus it won't have a
# console. The NO_INSTALL keyword causes the target not to be installed.
# NO_MTS_PCH makes the target not to use the default mitsuba precompiled header
# (on supported platforms). PCH specifies a custom precompiler header to use
# which is more suitable for the application.
#
# The optional RES_ICON parameter specified an icon to be bundled into the
# executable. This only works on Windows builds. The optional RES_DESCRIPTION
# parameters sets a specific executable description to be used in the Windows
# builds; other platforms simply ignore this value as with RES_ICON.

macro (add_mts_exe _exe_name)
  CMAKE_PARSE_ARGUMENTS(_exe "WIN32;NO_INSTALL;MTS_HW;MTS_BIDIR;NO_MTS_PCH"
    "PCH;RES_ICON;RES_DESCRIPTION" "LINK_LIBRARIES" ${ARGN})
  set (_exe_srcs ${_exe_UNPARSED_ARGUMENTS})
  if (_exe_WIN32)
    set(_exe_TYPE WIN32)
  endif()
  
  if (WIN32)
    set(_exe_res "${CMAKE_CURRENT_BINARY_DIR}/${_exe_name}_res.rc")
    set(_exe_res_args "${_exe_res}" "${_exe_name}" ".exe")
    if (_exe_RES_DESCRIPTION)
      set(_exe_description "${_exe_RES_DESCRIPTION}")
    else()
      set(_exe_description "Mitsuba application: ${_exe_name}")
    endif()
    list(APPEND _exe_res_args "${_exe_description}")
    if (_exe_RES_ICON)
      list(APPEND _exe_res_args "ICON" "${_exe_RES_ICON}")
    endif()
    mts_win_resource(${_exe_res_args})
    list(APPEND _exe_srcs "${_exe_res}")
  endif()
  
  if (MTS_USE_PCH AND (NOT _exe_NO_MTS_PCH OR _exe_PCH))
    set (_exe_pch_header "${MTS_DEFAULT_PCH}")
    if (_exe_PCH)
      set (_exe_pch_header "${_exe_PCH}")
      if (_exe_NO_MTS_PCH)
        message (AUTHOR_WARNING
	  "'NO_MTS_PCH' ignored due to 'PCH ${_exe_PCH}'.")
      endif ()
    endif ()
    pch_add_executable (${_exe_name} ${_exe_TYPE}
      PCH_HEADER "${_exe_pch_header}" ${_exe_srcs})
  else ()
    add_executable (${_exe_name} ${_exe_TYPE} ${_exe_srcs})
  endif ()

  set(_exe_core_libraries "mitsuba-core" "mitsuba-render")
  if (_exe_MTS_HW)
    list(APPEND _exe_core_libraries "mitsuba-hw")
  endif()
  if (_exe_MTS_BIDIR)
    list(APPEND _exe_core_libraries "mitsuba-bidir")
  endif()
  target_link_libraries (${_exe_name} ${_exe_core_libraries} ${_exe_LINK_LIBRARIES})
  if (WIN32)
    set_target_properties (${_exe_name} PROPERTIES VERSION "${MTS_VERSION}")
  endif()
  set_target_properties (${_exe_name} PROPERTIES FOLDER "apps")
  SET_OUTPATH_CFG (${_exe_name} RUNTIME_OUTPUT_DIRECTORY
    "${PROJECT_BINARY_DIR}/binaries/@CFGNAME@/${MTS_EXE_DEST}"
  )
  mts_target_ltcg (${_exe_name})
  mts_msvc_mp (${_exe_name})
  if (NOT _exe_NO_INSTALL)
    install(TARGETS ${_exe_name}
      RUNTIME DESTINATION ${MTS_EXE_DEST} COMPONENT Runtime)
  endif()
endmacro()



# Constant with the headers destination
if (WIN32 OR MTS_SIMPLE_PATHS)
  set (MTS_HEADER_DEST "sdk/include")
elseif (APPLE)
  set (MTS_HEADER_DEST "${MTS_BUNDLE_NAME}/Headers")
else()
  set (MTS_HEADER_DEST "include")
endif()

# Macro to install header files. The FOLDER option specifies a subdirectory
# on which the given headers will be installed. Usage:
#  mts_install_headers(header1 header2 ... [FOLDER subdir])
macro (mts_install_headers)
  CMAKE_PARSE_ARGUMENTS(_hdrs "" "FOLDER" "" ${ARGN})
  set (_hdrs_files ${_hdrs_UNPARSED_ARGUMENTS})
  if (NOT _hdrs_FOLDER)
    set (_hdrs_FOLDER ".")
  endif ()
  install (FILES ${_hdrs_files}
    PERMISSIONS "OWNER_READ" "GROUP_READ" "WORLD_READ"
    DESTINATION "${MTS_HEADER_DEST}/${_hdrs_FOLDER}"
    COMPONENT Developer)
endmacro ()



# Function to get a list of paths contained in variables which match the regex
# "LIBRAR(Y|IES)$". This is intended to be used with the FIXUP_BUNDLE macro.
# Example usage:
#  mts_library_paths (paths)
#  message (${paths})
function (mts_library_paths outvar_name)
  set (library_paths "")
  get_cmake_property (variables VARIABLES)
  foreach (var ${variables})
    if (var MATCHES "LIBRAR(Y|IES)$")
      foreach (path ${${var}})
        get_filename_component (libpath "${path}" PATH)
        if (libpath AND EXISTS "${libpath}")
          get_filename_component (libpath "${libpath}" REALPATH)
          set (library_paths ${library_paths} ${libpath})
        endif ()
      endforeach ()
    endif ()
  endforeach ()
  list (REMOVE_DUPLICATES library_paths)

  # Try to add the "bin" directories that might exist as well
  set (bin_paths "")
  if (WIN32)
    foreach (libpath ${library_paths})
      if (EXISTS "${libpath}/../bin")
        get_filename_component (binpath "${libpath}/../bin" REALPATH)
        set (bin_paths ${bin_paths} ${binpath})
      endif ()
      if (EXISTS "${libpath}/bin")
        get_filename_component (binpath "${libpath}/bin" REALPATH)
        set (bin_paths ${bin_paths} ${binpath})
      endif ()
    endforeach ()
  endif ()
  list (REMOVE_DUPLICATES bin_paths)
  set (${outvar_name} ${bin_paths} ${library_paths} PARENT_SCOPE)
endfunction ()
