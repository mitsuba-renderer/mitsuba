# Dependencies installation, bundle fixup and CPack
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

include (CMakeDependentOption)


# Offer an option to automatically bundle qt4image through BundleUtilities
include(CMakeDependentOption)
if(APPLE OR WIN32)
  set(MTS_FIXUP_BUNDLE_DEFAULT ON)
else()
  set(MTS_FIXUP_BUNDLE_DEFAULT OFF)
endif()
option(MTS_FIXUP_BUNDLE "Automatically bundle all mitsuba dependencies."
  ${MTS_FIXUP_BUNDLE_DEFAULT})


if (APPLE)
  set(MTS_QT_PLUGIN_DEST "${MTS_BUNDLE_NAME}/Contents/qt4plugins")
  set(MTS_QTCONF_DEST    "${MTS_BUNDLE_RESOURCES}")
elseif (WIN32)
  set(MTS_QT_PLUGIN_DEST "qt4plugins")
  set(MTS_QTCONF_DEST    ".")
else ()
  set(MTS_QT_PLUGIN_DEST "bin/qt4plugins")
  set(MTS_QTCONF_DEST    "bin")
endif()



# Install the needed Qt plugins
if (MTS_FIXUP_BUNDLE AND QT4_FOUND AND BUILD_GUI AND
    NOT(QT_CONFIG MATCHES "static"))
  set (MTS_BUNDLE_QTPLUGINS TRUE)
else ()
  set (MTS_BUNDLE_QTPLUGINS FALSE)
endif ()

# Imageformats, excluding svg,mng and gif (most likely not actually used by mitsuba)
if(MTS_BUNDLE_QTPLUGINS AND EXISTS "${QT_PLUGINS_DIR}/imageformats")
  if(NOT WIN32 OR NOT (QT_QTGUI_LIBRARY_RELEASE AND QT_QTGUI_LIBRARY_DEBUG))
    install(DIRECTORY "${QT_PLUGINS_DIR}/imageformats"
      DESTINATION "${MTS_QT_PLUGIN_DEST}"
      COMPONENT "Runtime"
      REGEX ".+(d4\\.dll|_debug\\.(dylib|so))$" EXCLUDE
      REGEX "svg|mng|gif" EXCLUDE)
  else()
    # There are different plugins for release and debug
    install(DIRECTORY "${QT_PLUGINS_DIR}/imageformats"
      DESTINATION "${MTS_QT_PLUGIN_DEST}"
      COMPONENT "Runtime" CONFIGURATIONS Debug
      FILES_MATCHING
      REGEX ".+(d4\\.dll|_debug\\.(dylib|so))$"
      REGEX "svg|mng|gif|\\.pdb$|\\.lib$" EXCLUDE)
    
    install(DIRECTORY "${QT_PLUGINS_DIR}/imageformats"
      DESTINATION "${MTS_QT_PLUGIN_DEST}"
      COMPONENT "Runtime" CONFIGURATIONS Release MinSizeRel RelWithDebInfo
      REGEX ".+(d4\\.dll|_debug\\.(dylib|so))$" EXCLUDE
      REGEX "svg|mng|gif|\\.pdb$|\\.lib$" EXCLUDE)
  endif()
endif()

# Bearer management (introduced with Qt 4.7)
if(MTS_BUNDLE_QTPLUGINS AND EXISTS "${QT_PLUGINS_DIR}/bearer")
  if(NOT WIN32 OR NOT (QT_QTNETWORK_LIBRARY_RELEASE AND QT_QTNETWORK_LIBRARY_DEBUG))
    install(DIRECTORY "${QT_PLUGINS_DIR}/bearer"
      DESTINATION "${MTS_QT_PLUGIN_DEST}"
      COMPONENT "Runtime"
      REGEX ".+(d4\\.dll|_debug\\.(dylib|so))$" EXCLUDE)
  else()
    # There are different plugins for release and debug
    install(DIRECTORY "${QT_PLUGINS_DIR}/bearer"
      DESTINATION "${MTS_QT_PLUGIN_DEST}"
      COMPONENT "Runtime" CONFIGURATIONS Debug
      FILES_MATCHING
      REGEX ".+(d4\\.dll|_debug\\.(dylib|so))$"
      REGEX "\\.pdb$|\\.lib$" EXCLUDE)
    
    install(DIRECTORY "${QT_PLUGINS_DIR}/bearer"
      DESTINATION "${MTS_QT_PLUGIN_DEST}"
      COMPONENT "Runtime" CONFIGURATIONS Release MinSizeRel RelWithDebInfo
      REGEX ".+(d4\\.dll|_debug\\.(dylib|so))$" EXCLUDE
      REGEX "\\.pdb$|\\.lib$" EXCLUDE)
  endif()
endif()

# The Cocoa version needs extra resources
set(QTCOCOA_MENU_NIB "${QT_QTGUI_LIBRARY}/Resources/qt_menu.nib")
if(MTS_FIXUP_BUNDLE AND QT4_FOUND AND QT_MAC_USE_COCOA
   AND EXISTS "${QTCOCOA_MENU_NIB}")
  install(CODE "
    message(STATUS \"Copying qt_menu.nib into the Bundle resources\")
    file(MAKE_DIRECTORY \"\${CMAKE_INSTALL_PREFIX}/${MTS_BUNDLE_RESOURCES}/qt_menu.nib\")
    execute_process(COMMAND \${CMAKE_COMMAND} -E copy_directory
      \"${QTCOCOA_MENU_NIB}\" 
      \"\${CMAKE_INSTALL_PREFIX}/${MTS_BUNDLE_RESOURCES}/qt_menu.nib\")
    " COMPONENT "Runtime")
endif()

if (MTS_BUNDLE_QTPLUGINS) 
  install(CODE "
message(STATUS \"Writing qt.conf\")
file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${MTS_QTCONF_DEST}/qt.conf\"
  \"[Paths]
Plugins = qt4plugins
\")
" COMPONENT "Runtime")
endif ()



if (MTS_FIXUP_BUNDLE)
  # Construct the list of directories with the currently set up libraries as
  # a way to automate "fixup_bundle" as much as possible
  mts_library_paths (MTS_LIBPATHS)
  
  configure_file ("MtsFixupBundle.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/MtsFixupBundle.cmake" @ONLY)
  install (SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/MtsFixupBundle.cmake"
    COMPONENT Runtime)

endif () # MTS_FIXUP_BUNDLE
