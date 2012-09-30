# mitsuba-hw
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/mitsuba/hw)
set(HDRS
  ${INCLUDE_DIR}/basicshader.h
  ${INCLUDE_DIR}/device.h
  ${INCLUDE_DIR}/font.h
  ${INCLUDE_DIR}/glgeometry.h
  ${INCLUDE_DIR}/glprogram.h
  ${INCLUDE_DIR}/glrenderer.h
  ${INCLUDE_DIR}/glsync.h
  ${INCLUDE_DIR}/gltexture.h
  ${INCLUDE_DIR}/gpugeometry.h
  ${INCLUDE_DIR}/gpuprogram.h
  ${INCLUDE_DIR}/gpusync.h
  ${INCLUDE_DIR}/gputexture.h
  ${INCLUDE_DIR}/renderer.h
  ${INCLUDE_DIR}/shadow.h
  ${INCLUDE_DIR}/session.h
  ${INCLUDE_DIR}/viewer.h
  ${INCLUDE_DIR}/vpl.h
)

set(SRCS
  basicshader.cpp
  device.cpp
  font.cpp
  glgeometry.cpp
  glprogram.cpp
  glrenderer.cpp
  glsync.cpp
  gltexture.cpp
  gpugeometry.cpp
  gpuprogram.cpp
  gpusync.cpp
  gputexture.cpp
  renderer.cpp
  shadow.cpp
  session.cpp
  viewer.cpp
  vpl.cpp
)

if (WIN32)
  list(APPEND HDRS
    ${INCLUDE_DIR}/wgldevice.h
    ${INCLUDE_DIR}/wglrenderer.h
    ${INCLUDE_DIR}/wglsession.h
  )
  list(APPEND SRCS wgldevice.cpp wglsession.cpp wglrenderer.cpp)
elseif(APPLE)
  list(APPEND HDRS
    ${INCLUDE_DIR}/nsgldevice.h
    ${INCLUDE_DIR}/nsglkeys.h
    ${INCLUDE_DIR}/nsglrenderer.h
    ${INCLUDE_DIR}/nsglsession.h
  )
  list(APPEND SRCS nsgldevice.mm nsglsession.mm nsglrenderer.mm)
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  list(APPEND HDRS
    ${INCLUDE_DIR}/glxdevice.h
    ${INCLUDE_DIR}/glxrenderer.h
	${INCLUDE_DIR}/x11device.h
    ${INCLUDE_DIR}/x11session.h
  )
  list(APPEND SRCS glxdevice.cpp glxrenderer.cpp x11device.cpp x11session.cpp)
endif()

add_definitions(-DMTS_BUILD_MODULE=MTS_MODULE_HW)
add_definitions(${GLEW_DEFINITIONS})

include_directories(${OPENGL_INCLUDE_DIR} ${GLEW_INCLUDE_DIRS})

set(LIBS mitsuba-core mitsuba-render ${OPENGL_LIBRARIES} ${GLEW_LIBRARIES})
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  include_directories(${X11_INCLUDE_DIR})
  list(APPEND LIBS ${X11_LIBRARIES} ${X11_Xxf86vm_LIB})
endif()

add_mts_corelib(mitsuba-hw ${HDRS} ${SRCS} LINK_LIBRARIES ${LIBS})

mts_install_headers (${HDRS} FOLDER "mitsuba/hw")
