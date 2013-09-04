# main programs
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

include_directories(${XERCES_INCLUDE_DIRS})

add_mts_exe(mtssrv  mtssrv.cpp ${MTS_WINDOWS_STUB}
  RES_DESCRIPTION "Mitsuba compute node application")
add_mts_exe(mitsuba mitsuba.cpp ${MTS_WINDOWS_STUB}
  RES_DESCRIPTION "Mitsuba command line interface frontend")
add_mts_exe(mtsutil mtsutil.cpp ${MTS_DARWIN_STUB} ${MTS_WINDOWS_STUB}
  RES_DESCRIPTION "Mitsuba utility launcher")
