# Scene XML Schema and XSL transformer
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

if (WIN32 OR MTS_SIMPLE_PATHS)
  set(SCHEMA_DESTINATION "data/schema")
elseif (APPLE)
  set(SCHEMA_DESTINATION "${MTS_BUNDLE_NAME}/data/schema")
else()
  set(SCHEMA_DESTINATION "share/mitsuba/data/schema")
endif()
  
# Copy the XML components into the output directory, mainly for easier debugging
set (SCHEMA_DIR "${PROJECT_BINARY_DIR}/binaries/${CMAKE_CFG_INTDIR}/${SCHEMA_DESTINATION}")
add_custom_command (
  OUTPUT "${SCHEMA_DIR}/scene.xsd"
         "${SCHEMA_DIR}/upgrade_0.3.0.xsl"
         "${SCHEMA_DIR}/upgrade_0.4.0.xsl"
         "${SCHEMA_DIR}/upgrade_0.5.0.xsl"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/scene.xsd" "${SCHEMA_DIR}/scene.xsd"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/upgrade_0.3.0.xsl" "${SCHEMA_DIR}/upgrade_0.3.0.xsl"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/upgrade_0.4.0.xsl" "${SCHEMA_DIR}/upgrade_0.4.0.xsl"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/upgrade_0.5.0.xsl" "${SCHEMA_DIR}/upgrade_0.5.0.xsl"
  MAIN_DEPENDENCY "scene.xsd"
  DEPENDS "scene.xsd" "upgrade_0.3.0.xsl" "upgrade_0.4.0.xsl" "upgrade_0.5.0.xsl"
  COMMENT "Copying XML data: scene schema and upgrade XSL files"
)
add_custom_target (scene_schema
  DEPENDS "${SCHEMA_DIR}/scene.xsd"
          "${SCHEMA_DIR}/upgrade_0.3.0.xsl" "${SCHEMA_DIR}/upgrade_0.4.0.xsl" "${SCHEMA_DIR}/upgrade_0.5.0.xsl"
)
set_target_properties (scene_schema PROPERTIES FOLDER "data")

# Actual installation
install(FILES "scene.xsd" "upgrade_0.3.0.xsl" "upgrade_0.4.0.xsl" "upgrade_0.5.0.xsl"
  DESTINATION ${SCHEMA_DESTINATION}
  PERMISSIONS "OWNER_READ" "GROUP_READ" "WORLD_READ"
  COMPONENT Runtime
)
