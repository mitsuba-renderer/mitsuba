# Microfacet precomputed data
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

if (WIN32 OR MTS_SIMPLE_PATHS)
  set(MICROFACET_DESTINATION "data/microfacet")
elseif (APPLE)
  set(MICROFACET_DESTINATION "${MTS_BUNDLE_NAME}/data/microfacet")
else()
  set(MICROFACET_DESTINATION "share/mitsuba/data/microfacet")
endif()
  
# Copy the microfacet data into the output directory, mainly for easier debugging
set (MICROFACET_DIR "${PROJECT_BINARY_DIR}/binaries/${CMAKE_CFG_INTDIR}/${MICROFACET_DESTINATION}")
add_custom_command (
  OUTPUT "${MICROFACET_DIR}/beckmann.dat"
         "${MICROFACET_DIR}/ggx.dat"
         "${MICROFACET_DIR}/phong.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/beckmann.dat" "${MICROFACET_DIR}/beckmann.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/ggx.dat" "${MICROFACET_DIR}/ggx.dat"
  COMMAND "${CMAKE_COMMAND}" -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/phong.dat" "${MICROFACET_DIR}/phong.dat"
  DEPENDS "beckmann.dat" "ggx.dat" "phong.dat"
  COMMENT "Copying microfacet precomputed data"
)
add_custom_target (microfacet_data
  DEPENDS "${MICROFACET_DIR}/beckmann.dat"
          "${MICROFACET_DIR}/ggx.dat"
          "${MICROFACET_DIR}/phong.dat"
)
set_target_properties (microfacet_data PROPERTIES FOLDER "data")

# Actual installation
install(FILES "beckmann.dat" "ggx.dat" "phong.dat"
  DESTINATION ${MICROFACET_DESTINATION}
  PERMISSIONS "OWNER_READ" "GROUP_READ" "WORLD_READ"
  COMPONENT Runtime
)
