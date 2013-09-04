# qtgui (mtsgui)
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

option(MTS_GUI_SOFTWARE_FALLBACK
  "Forces the GUI to use a software fallback mode, which is considerably slower and removes the realtime preview."
  OFF)
mark_as_advanced(MTS_GUI_SOFTWARE_FALLBACK)
if (MTS_GUI_SOFTWARE_FALLBACK)
  add_definitions(-DMTS_GUI_SOFTWARE_FALLBACK)
endif()


# Set up the Qt4 macros and variables
if (WIN32)
  set(QT_USE_QTMAIN 0)
  set(MTS_EXE_TYPE "WIN32")
endif()
include(${QT_USE_FILE})

# Add all the sources
set (UI_FILES
  aboutdlg.ui
  acknowledgmentdlg.ui
  addserverdlg.ui
  importdlg.ui
  loaddlg.ui
  locateresourcedlg.ui
  mainwindow.ui
  previewsettingsdlg.ui
  programsettingsdlg.ui
  rendersettingsdlg.ui
  sceneinfodlg.ui
  updatedlg.ui
)
set (RC_FILES
  resources.qrc
)
set (SRCS
  aboutdlg.cpp
  acknowledgmentdlg.cpp
  addserverdlg.cpp
  glwidget.cpp
  importdlg.cpp
  loaddlg.cpp
  locateresourcedlg.cpp
  logwidget.cpp
  main.cpp
  mainwindow.cpp
  preview.cpp
  previewsettingsdlg.cpp
  programsettingsdlg.cpp
  rendersettingsdlg.cpp
  save.cpp
  sceneimporter.cpp
  sceneinfodlg.cpp
  sceneloader.cpp
  server.cpp
  symlinks_auth.cpp
  tabbar.cpp
  updatedlg.cpp
  upgrade.cpp
  xmltreemodel.cpp
)

# Headers implementing the UI files
set (HDRS_UI
  aboutdlg.h
  acknowledgmentdlg.h
  addserverdlg.h
  importdlg.h
  loaddlg.h
  locateresourcedlg.h
  mainwindow.h
  previewsettingsdlg.h
  programsettingsdlg.h
  rendersettingsdlg.h
  sceneinfodlg.h
  updatedlg.h
)

# Headers implementing QObjects
set (HDRS_QOBJECT
  glwidget.h
  logwidget.h
  preview.h
  server.h
  sceneimporter.h
  tabbar.h
  upgrade.h
  xmltreemodel.h
)

# Normal headers
set (HDRS
  common.h
  save.h
  sceneloader.h
)

# Fast SSE Tonemapper
if (MTS_SSE)
  set(SSE_SRCS simdtonemap.h simdtonemap.cpp)
  list(APPEND SRCS ${SSE_SRCS})
  # Because we reuse source files between mtsgui and this test case, on OSX
  # there are problems when using PCHs: simdtonemap.cpp gets a file-level
  # dependency on the mtsgui later on. This causes a conflict with the
  # target-level PCH on the test case
  if (APPLE)
    set(test_simdtonemap_PCH "NO_MTS_PCH")
  endif()
  add_mts_plugin(test_simdtonemap ${test_simdtonemap_PCH} 
    test_simdtonemap.cpp ${SSE_SRCS} TYPE testcase)
endif()

# Add the cocoa-specific files
if (APPLE)
  list (APPEND SRCS
    previewsettingsdlg_cocoa.cpp previewsettingsdlg_cocoa_impl.mm)
  list (APPEND HDRS_UI
    previewsettingsdlg_cocoa.h)
endif()


# Process the UI files through uic and store the generated headers
QT4_WRAP_UI(UI_HDRS ${UI_FILES})
source_group("UI Files" FILES ${UI_FILES})

# Add the directory where the generated headers are
include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

# Process the resource files through rcc and store the generated cxx files
QT4_ADD_RESOURCES(RC_SRCS ${RC_FILES})
source_group("Qt Resources" FILES ${RC_FILES})

# Add our own qtmain and use the static file dialogs
if (WIN32)
  add_definitions(-D MTS_CUSTOM_QTMAIN -DMTSGUI_STATIC_QFILEDIALOG=1)
  list(APPEND SRCS qtmain_win.cpp)
endif()

# Generate the moc files
QT4_WRAP_CPP(MOC_HDRS ${HDRS_UI} ${HDRS_QOBJECT})

# Create a source group for the generated files
source_group("Qt Generated" FILES ${MOC_HDRS} ${UI_HDRS} ${RC_SRCS})

include_directories(${XERCES_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
add_definitions(${GLEW_DEFINITIONS})

set(LIBS ${QT_LIBRARIES} ${GLEW_LIBRARIES})
# This static library is built by the mtsimport target (src/converter)
if (COLLADA_FOUND)
  if (NOT TARGET mtsconverter_lib)
    message(FATAL_ERROR "Cannot find the mtsconverter_lib target")
  endif()
  list(APPEND LIBS mtsconverter_lib)
endif()

if (APPLE)
  include_directories(${BWTOOLKIT_INCLUDE_DIR} ${SECURITY_INCLUDE_DIR})
  list(APPEND LIBS ${BWTOOLKIT_LIBRARY} ${SECURITY_LIBRARY})
endif()

add_mts_exe(mtsgui ${MTS_EXE_TYPE} ${SRCS} ${HDRS}
  ${HDRS_UI} ${HDRS_QOBJECT}
  ${UI_HDRS} ${RC_SRCS} ${MOC_HDRS}
  PCH "${CMAKE_CURRENT_SOURCE_DIR}/pch/mtsgui_precompiled.hpp"
  RES_ICON "${PROJECT_SOURCE_DIR}/data/windows/mitsuba.ico"
  RES_DESCRIPTION "Mitsuba interactive Qt-based frontend"
  MTS_HW LINK_LIBRARIES ${LIBS})
  
 
# Custom targets to copy the bundle files and the PreviewSettings objects
# into the build directory. Also adds the special symlinks_install tool
if (APPLE)
  add_executable (symlinks_install "symlinks_install.c")
  SET_OUTPATH_CFG (symlinks_install RUNTIME_OUTPUT_DIRECTORY
    "${PROJECT_BINARY_DIR}/binaries/@CFGNAME@/${MTS_EXE_DEST}")
  install(TARGETS symlinks_install
      RUNTIME DESTINATION ${MTS_EXE_DEST} COMPONENT Runtime)

  set (CONTENTS_DIR  "${PROJECT_BINARY_DIR}/binaries/${CMAKE_CFG_INTDIR}/${MTS_BUNDLE_NAME}/Contents")
  set (RESOURCES_DIR "${PROJECT_BINARY_DIR}/binaries/${CMAKE_CFG_INTDIR}/${MTS_BUNDLE_RESOURCES}")
  set (DATA_DIR      "${PROJECT_SOURCE_DIR}/data/darwin")
 
  # CMake script to poplutate the bundle resources
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/BundleComponents.cmake"
"cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

# Configure the Info.plist to add the version
if (NOT DEFINED MTS_VERSION)
  set (MTS_VERSION \"${MTS_VERSION}\")  
endif ()
configure_file(\"${DATA_DIR}/Info.plist.in\" \"${CONTENTS_DIR}/Info.plist\" @ONLY)
# Just copy the PkgInfo file
configure_file(\"${DATA_DIR}/PkgInfo\" \"${CONTENTS_DIR}/PkgInfo\" COPYONLY)
")

  add_custom_command (
    OUTPUT "${CONTENTS_DIR}/Info.plist" "${CONTENTS_DIR}/PkgInfo"
    COMMAND "${CMAKE_COMMAND}" -P
      "${CMAKE_CURRENT_BINARY_DIR}/BundleComponents.cmake"
    DEPENDS "${DATA_DIR}/Info.plist.in" "${DATA_DIR}/PkgInfo"
  )
  add_custom_command (
    OUTPUT "${RESOURCES_DIR}/mitsuba.icns"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "${DATA_DIR}/Resources/mitsuba.icns"
	    "${RESOURCES_DIR}/mitsuba.icns"
    DEPENDS "${DATA_DIR}/Resources/mitsuba.icns"
  )
  add_custom_command (
    OUTPUT "${RESOURCES_DIR}/PreviewSettings.nib"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "${DATA_DIR}/PreviewSettings.nib"
	    "${RESOURCES_DIR}/PreviewSettings.nib"
    DEPENDS "${DATA_DIR}/PreviewSettings.nib"
  )
  
  add_custom_target (mtsgui_BundleResources DEPENDS
    "${CONTENTS_DIR}/Info.plist" "${CONTENTS_DIR}/PkgInfo"
    "${RESOURCES_DIR}/mitsuba.icns" "${RESOURCES_DIR}/PreviewSettings.nib"
  )
  
  # Actual installation
  install(FILES
    "${CONTENTS_DIR}/Info.plist"
    "${CONTENTS_DIR}/PkgInfo"
    DESTINATION "${MTS_BUNDLE_NAME}/Contents"
    COMPONENT "Runtime")
  install(FILES
    "${RESOURCES_DIR}/mitsuba.icns"
    DESTINATION "${MTS_BUNDLE_RESOURCES}"
    COMPONENT "Runtime")
  install(FILES
    "${RESOURCES_DIR}/PreviewSettings.nib"
    DESTINATION "${MTS_BUNDLE_RESOURCES}"
    COMPONENT "Runtime")

  # Add the dependencies to mtsgui to build the custom targets
  add_dependencies(mtsgui mtsgui_BundleResources)
endif ()
