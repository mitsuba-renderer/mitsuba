# IOR Database
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

if (WIN32 OR MTS_SIMPLE_PATHS)
  set(IOR_DESTINATION "data/ior")
elseif (APPLE)
  set(IOR_DESTINATION "${MTS_BUNDLE_NAME}/data/ior")
else()
  set(IOR_DESTINATION "share/mitsuba/data/ior")
endif()

set (IOR_DIR "${PROJECT_BINARY_DIR}/binaries/${CMAKE_CFG_INTDIR}/${IOR_DESTINATION}")
  
# Actual installation
install(DIRECTORY "."
  DESTINATION "${IOR_DESTINATION}"
  FILE_PERMISSIONS "OWNER_READ" "GROUP_READ" "WORLD_READ"
  COMPONENT Runtime
  FILES_MATCHING PATTERN "*.spd"
)

# Custom target wich creates a phony file, just to flag that the database
# has already been copied
add_custom_command(OUTPUT "${IOR_DIR}/.iordatabase.flag"
  COMMAND "${CMAKE_COMMAND}" -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}" "${IOR_DIR}"
  COMMAND "${CMAKE_COMMAND}" -E remove -f "${IOR_DIR}/CMakeLists.txt"
  COMMAND "${CMAKE_COMMAND}" -E touch "${IOR_DIR}/.iordatabase.flag"
  COMMENT "Copying IOR database directory"
)
add_custom_target (ior_database DEPENDS "${IOR_DIR}/.iordatabase.flag")
set_target_properties (ior_database PROPERTIES FOLDER "data")
