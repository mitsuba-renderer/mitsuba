# mitsuba-render
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/mitsuba/render)
set(HDRS
  ${INCLUDE_DIR}/bsdf.h
  ${INCLUDE_DIR}/common.h
  ${INCLUDE_DIR}/emitter.h
  ${INCLUDE_DIR}/film.h
  ${INCLUDE_DIR}/fwd.h
  ${INCLUDE_DIR}/gatherproc.h
  ${INCLUDE_DIR}/gkdtree.h
  ${INCLUDE_DIR}/imageblock.h
  ${INCLUDE_DIR}/imageproc.h
  ${INCLUDE_DIR}/integrator.h
  ${INCLUDE_DIR}/irrcache.h
  ${INCLUDE_DIR}/medium.h
  ${INCLUDE_DIR}/mipmap.h
  ${INCLUDE_DIR}/noise.h
  ${INCLUDE_DIR}/particleproc.h
  ${INCLUDE_DIR}/phase.h
  ${INCLUDE_DIR}/photon.h
  ${INCLUDE_DIR}/photonmap.h
  ${INCLUDE_DIR}/range.h
  ${INCLUDE_DIR}/records.inl
  ${INCLUDE_DIR}/rectwu.h
  ${INCLUDE_DIR}/renderjob.h
  ${INCLUDE_DIR}/renderproc.h
  ${INCLUDE_DIR}/renderqueue.h
  ${INCLUDE_DIR}/sahkdtree2.h
  ${INCLUDE_DIR}/sahkdtree3.h
  ${INCLUDE_DIR}/sampler.h
  ${INCLUDE_DIR}/scene.h
  ${INCLUDE_DIR}/scenehandler.h
  ${INCLUDE_DIR}/sensor.h
  ${INCLUDE_DIR}/shader.h
  ${INCLUDE_DIR}/shape.h
  ${INCLUDE_DIR}/skdtree.h
  ${INCLUDE_DIR}/spiral.h
  ${INCLUDE_DIR}/subsurface.h
  ${INCLUDE_DIR}/testcase.h
  ${INCLUDE_DIR}/texture.h
  ${INCLUDE_DIR}/triaccel.h
  ${INCLUDE_DIR}/triaccel_sse.h
  ${INCLUDE_DIR}/trimesh.h
  ${INCLUDE_DIR}/util.h
  ${INCLUDE_DIR}/volume.h
  ${INCLUDE_DIR}/vpl.h
)

set(SRCS
  bsdf.cpp
  common.cpp
  emitter.cpp
  film.cpp
  gatherproc.cpp
  imageblock.cpp
  imageproc.cpp
  integrator.cpp
  intersection.cpp
  irrcache.cpp
  medium.cpp
  noise.cpp
  particleproc.cpp
  phase.cpp
  photon.cpp
  photonmap.cpp
  rectwu.cpp
  renderjob.cpp
  renderproc.cpp
  renderqueue.cpp
  sampler.cpp
  scene.cpp
  scenehandler.cpp
  sensor.cpp
  shader.cpp
  shape.cpp
  skdtree.cpp
  subsurface.cpp
  testcase.cpp
  texture.cpp
  trimesh.cpp
  util.cpp
  volume.cpp
  vpl.cpp
)

add_definitions(-DMTS_BUILD_MODULE=MTS_MODULE_RENDER)

include_directories(${ZLIB_INCLUDE_DIRS} ${XERCES_INCLUDE_DIRS})

add_mts_corelib(mitsuba-render ${HDRS} ${SRCS} LINK_LIBRARIES
  mitsuba-core ${XERCES_LIBRARIES})
# Make sure the schema is where expected by the scene handler
add_dependencies(mitsuba-render scene_schema)

mts_install_headers (${HDRS} FOLDER "mitsuba/render")
