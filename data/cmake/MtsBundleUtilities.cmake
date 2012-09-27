# - Functions to help assemble a standalone bundle application (with @rpath).
# A collection of CMake utility functions useful for dealing with .app
# bundles on the Mac and bundle-like directories on any OS.
# This module is based on the "BundleUtilities" module included with
# CMake 2.8.5, the difference is that for the Mac this module assumes
# that all dependencies (both dynamic libraries and frameworks) are installed
# under <bundle>/Contents/Frameworks. Thus the fixup stage will use @rpath
# instead of @executable_path.
#
# The following functions are provided by this module:
#   mts_fixup_bundle
#   mts_clear_bundle_keys (privately used by mts_fixup_bundle)
#   mts_fixup_bundle_item (privately used by mts_fixup_bundle)
# Requires CMake 2.6 or greater because it uses function, break and
# PARENT_SCOPE. Also depends on GetPrerequisites.cmake and the standard
# BundleUtilities.cmake.
#
#  MTS_FIXUP_BUNDLE(<app> <libs> <dirs>)
# Fix up a bundle in-place and make it standalone, such that it can be
# drag-n-drop copied to another machine and run on that machine as long as all
# of the system libraries are compatible.
#
# If you pass plugins to fixup_bundle as the libs parameter, you should install
# them or copy them into the bundle before calling fixup_bundle. The "libs"
# parameter is a list of libraries that must be fixed up, but that cannot be
# determined by otool output analysis. (i.e., plugins)
#
# Gather all the keys for all the executables and libraries in a bundle, and
# then, for each key, copy each prerequisite into the bundle. Then fix each one
# up according to its own list of prerequisites.
#
# Then clear all the keys and call verify_app on the final bundle to ensure
# that it is truly standalone.
#
#  MTS_FIXUP_BUNDLE_ITEM(<resolved_embedded_item> <exepath> <dirs>)
# Get the direct/non-system prerequisites of the resolved embedded item. For
# each prerequisite, change the way it is referenced to the value of the
# _EMBEDDED_ITEM keyed variable for that prerequisite. (Most likely changing to
# an "@executable_path" style reference.)
#
# This function requires that the resolved_embedded_item be "inside" the bundle
# already. In other words, if you pass plugins to fixup_bundle as the libs
# parameter, you should install them or copy them into the bundle before
# calling fixup_bundle. The "libs" parameter is a list of libraries that must
# be fixed up, but that cannot be determined by otool output analysis. (i.e.,
# plugins)
#
# Also, change the id of the item being fixed up to its own _EMBEDDED_ITEM
# value.
#
# Accumulate changes in a local variable and make *one* call to
# install_name_tool at the end of the function with all the changes at once.
#
# If the BU_CHMOD_BUNDLE_ITEMS variable is set then bundle items will be
# marked writable before install_name_tool tries to change them.

#=============================================================================
# Copyright 2008-2009 Kitware, Inc.
#
# CMake - Cross Platform Makefile Generator
# Copyright 2000-2009 Kitware, Inc., Insight Software Consortium
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# 
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# ----------------------------------------------------------------------------
# 
# The above copyright and license notice applies to distributions of
# CMake in source and binary form.  Some source files contain additional
# notices of original copyright by their contributors; see each source
# for details.  Third-party software packages supplied with CMake under
# compatible licenses provide their own copyright notices documented in
# corresponding subdirectories.
# 
# ----------------------------------------------------------------------------
# 
# CMake was initially developed by Kitware with the following sponsorship:
# 
#  * National Library of Medicine at the National Institutes of Health
#    as part of the Insight Segmentation and Registration Toolkit (ITK).
# 
#  * US National Labs (Los Alamos, Livermore, Sandia) ASC Parallel
#    Visualization Initiative.
# 
#  * National Alliance for Medical Image Computing (NAMIC) is funded by the
#    National Institutes of Health through the NIH Roadmap for Medical
#    Research, Grant U54 EB005149.
# 
#=============================================================================


# The module still depends on the standard BundleUtilities and GetPrerequites modules
include(BundleUtilities)
include(GetPrerequisites)


function(mts_clear_bundle_keys keys_var)
  foreach(key ${${keys_var}})
    set(${key}_ITEM PARENT_SCOPE)
    set(${key}_RESOLVED_ITEM PARENT_SCOPE)
    set(${key}_DEFAULT_EMBEDDED_PATH PARENT_SCOPE)
    set(${key}_EMBEDDED_ITEM PARENT_SCOPE)
    set(${key}_RESOLVED_EMBEDDED_ITEM PARENT_SCOPE)
    set(${key}_COPYFLAG PARENT_SCOPE)
    set(${key}_EMBEDDED_RPATH PARENT_SCOPE)
    set(${key}_RPATH_DIR PARENT_SCOPE)
  endforeach(key)
  set(${keys_var} PARENT_SCOPE)
endfunction(mts_clear_bundle_keys)


function(mts_fixup_bundle_item resolved_embedded_item exepath dirs)
  # This item's key is "ikey":
  #
  get_item_key("${resolved_embedded_item}" ikey)

  # Ensure the item is "inside the .app bundle" -- it should not be fixed up if
  # it is not in the .app bundle... Otherwise, we'll modify files in the build
  # tree, or in other varied locations around the file system, with our call to
  # install_name_tool. Make sure that doesn't happen here:
  #
  get_dotapp_dir("${exepath}" exe_dotapp_dir)
  string(LENGTH "${exe_dotapp_dir}/" exe_dotapp_dir_length)
  string(LENGTH "${resolved_embedded_item}" resolved_embedded_item_length)
  set(path_too_short 0)
  set(is_embedded 0)
  if(${resolved_embedded_item_length} LESS ${exe_dotapp_dir_length})
    set(path_too_short 1)
  endif()
  if(NOT path_too_short)
    string(SUBSTRING "${resolved_embedded_item}" 0 ${exe_dotapp_dir_length} item_substring)
    if("${exe_dotapp_dir}/" STREQUAL "${item_substring}")
      set(is_embedded 1)
    endif()
  endif()
  if(NOT is_embedded)
    message("  exe_dotapp_dir/='${exe_dotapp_dir}/'")
    message("  item_substring='${item_substring}'")
    message("  resolved_embedded_item='${resolved_embedded_item}'")
    message("")
    message("Install or copy the item into the bundle before calling mts_fixup_bundle.")
    message("Or maybe there's a typo or incorrect path in one of the args to mts_fixup_bundle?")
    message("")
    message(FATAL_ERROR "cannot fixup an item that is not in the bundle...")
  endif()

  set(prereqs "")
  get_prerequisites("${resolved_embedded_item}" prereqs 1 0 "${exepath}" "${dirs}")

  set(changes "")
  set(rpath_dirlst "")

  foreach(pr ${prereqs})
    # Each referenced item's key is "rkey" in the loop:
    #
    get_item_key("${pr}" rkey)

    if(NOT "${${rkey}_EMBEDDED_ITEM}" STREQUAL "")
      set(changes ${changes} "-change" "${pr}" "${${rkey}_EMBEDDED_ITEM}")
    else(NOT "${${rkey}_EMBEDDED_ITEM}" STREQUAL "")
      message("warning: unexpected reference to '${pr}'")
    endif(NOT "${${rkey}_EMBEDDED_ITEM}" STREQUAL "")
    
    if(NOT "${${rkey}_EMBEDDED_RPATH}" STREQUAL "" AND
       NOT "${${rkey}_RPATH_DIR}" STREQUAL "")
      # Assumes the rpath dir is already an absolute path
      list(APPEND rpath_dirlst "${${rkey}_RPATH_DIR}")
    endif()
  endforeach(pr)

  if(BU_CHMOD_BUNDLE_ITEMS)
    execute_process(COMMAND chmod u+w "${resolved_embedded_item}")
  endif()
  
  # Assumes that the only required location
  if(rpath_dirlst)
    list(REMOVE_DUPLICATES rpath_dirlst)
    list(LENGTH rpath_dirlst rpath_len)
    if (rpath_len GREATER 1)
      message(WARNING "Only one location for rpath is supported, ${rpath_len} provided.")
    endif ()
    list(GET rpath_dirlst 0 rpath_dir)
    
    # Determine how to get from the current component directory to rpath_dir
    get_filename_component(resolved_embedded_path "${resolved_embedded_item}" PATH)
    file(RELATIVE_PATH rpath_relative "${resolved_embedded_path}" "${rpath_dir}")
    if("${rpath_relative}" STREQUAL "")
      set(embedded_rpath "@loader_path/.")
    else()
      set(embedded_rpath "@loader_path/${rpath_relative}")
    endif()
    list(APPEND changes "-add_rpath" "${embedded_rpath}")
  endif()

  # Change this item's id and all of its references in one call
  # to install_name_tool:
  #
  execute_process(COMMAND install_name_tool
    ${changes} -id "${${ikey}_EMBEDDED_ITEM}" "${resolved_embedded_item}"
  )
endfunction(mts_fixup_bundle_item)


function(mts_fixup_bundle app libs dirs)
  message(STATUS "mts_fixup_bundle")
  message(STATUS "  app='${app}'")
  message(STATUS "  libs='${libs}'")
  message(STATUS "  dirs='${dirs}'")

  get_bundle_and_executable("${app}" bundle executable valid)
  if(valid)
    get_filename_component(exepath "${executable}" PATH)

    message(STATUS "mts_fixup_bundle: preparing...")
    get_bundle_keys("${app}" "${libs}" "${dirs}" keys)

    message(STATUS "mts_fixup_bundle: copying...")
    list(LENGTH keys n)
    math(EXPR n ${n}*2)

    set(i 0)
    foreach(key ${keys})
      math(EXPR i ${i}+1)
      if(${${key}_COPYFLAG})
        message(STATUS "${i}/${n}: copying '${${key}_RESOLVED_ITEM}'")
      else(${${key}_COPYFLAG})
        message(STATUS "${i}/${n}: *NOT* copying '${${key}_RESOLVED_ITEM}'")
      endif(${${key}_COPYFLAG})
      
      # Modify the embedded flag to use rpath. Assumes all frameworks and
      # libraries were copied to <bundle>/Contents/Frameworks
      if (APPLE AND NOT "${${key}_EMBEDDED_ITEM}" STREQUAL "" AND
          "${${key}_ITEM}" MATCHES "[^/]+(\\.framework/|\\.dylib$)")
        file(RELATIVE_PATH irelpath "${bundle}/Contents/Frameworks" "${${key}_RESOLVED_EMBEDDED_ITEM}")
        # Check if the item is in fact in the "Frameworks" directory
        string(SUBSTRING "${irelpath}" 0 3 irelpath_start)
        if (NOT "${irelpath_start}" STREQUAL "../")
          # Replace the embedded key for one using rpath
          set (${key}_EMBEDDED_ITEM "@rpath/${irelpath}")
          
          # Flag this key as using RPATH
          set (${key}_EMBEDDED_RPATH 1)
          set (${key}_RPATH_DIR "${bundle}/Contents/Frameworks")
        endif ()
      endif ()

      set(show_status 0)
      if(show_status)
        message(STATUS "key='${key}'")
        message(STATUS "item='${${key}_ITEM}'")
        message(STATUS "resolved_item='${${key}_RESOLVED_ITEM}'")
        message(STATUS "default_embedded_path='${${key}_DEFAULT_EMBEDDED_PATH}'")
        message(STATUS "embedded_item='${${key}_EMBEDDED_ITEM}'")
        message(STATUS "resolved_embedded_item='${${key}_RESOLVED_EMBEDDED_ITEM}'")
        message(STATUS "copyflag='${${key}_COPYFLAG}'")
        message(STATUS "embedded_rpath='${${key}_EMBEDDED_RPATH}'")
        message(STATUS "rpath_dir='${${key}_RPATH_DIR}'")
        message(STATUS "")
      endif(show_status)

      if(${${key}_COPYFLAG})
        set(item "${${key}_ITEM}")
        if(item MATCHES "[^/]+\\.framework/")
          copy_resolved_framework_into_bundle("${${key}_RESOLVED_ITEM}"
            "${${key}_RESOLVED_EMBEDDED_ITEM}")
        else()
          copy_resolved_item_into_bundle("${${key}_RESOLVED_ITEM}"
            "${${key}_RESOLVED_EMBEDDED_ITEM}")
        endif()
      endif(${${key}_COPYFLAG})
    endforeach(key)

    message(STATUS "mts_fixup_bundle: fixing...")
    foreach(key ${keys})
      math(EXPR i ${i}+1)
      if(APPLE)
        message(STATUS "${i}/${n}: fixing up '${${key}_RESOLVED_EMBEDDED_ITEM}'")
        mts_fixup_bundle_item("${${key}_RESOLVED_EMBEDDED_ITEM}" "${exepath}" "${dirs}")
      else(APPLE)
        message(STATUS "${i}/${n}: fix-up not required on this platform '${${key}_RESOLVED_EMBEDDED_ITEM}'")
      endif(APPLE)
    endforeach(key)

    message(STATUS "mts_fixup_bundle: cleaning up...")
    mts_clear_bundle_keys(keys)

    message(STATUS "mts_fixup_bundle: verifying...")
    verify_app("${app}")
  else(valid)
    message(SEND_ERROR "error: mts_fixup_bundle: not a valid bundle")
  endif(valid)

  message(STATUS "mts_fixup_bundle: done")
endfunction(mts_fixup_bundle)
