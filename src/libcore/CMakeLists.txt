# mitsuba-core
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

set(MTS_MAIN_HEADER ${PROJECT_SOURCE_DIR}/include/mitsuba/mitsuba.h)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/mitsuba/core)
set(HDRS
  ${INCLUDE_DIR}/aabb.h
  ${INCLUDE_DIR}/aabb_sse.h
  ${INCLUDE_DIR}/appender.h
  ${INCLUDE_DIR}/atomic.h
  ${INCLUDE_DIR}/barray.h
  ${INCLUDE_DIR}/bitmap.h
  ${INCLUDE_DIR}/brent.h
  ${INCLUDE_DIR}/bsphere.h
  ${INCLUDE_DIR}/simplecache.h
  ${INCLUDE_DIR}/chisquare.h
  ${INCLUDE_DIR}/class.h
  ${INCLUDE_DIR}/cobject.h
  ${INCLUDE_DIR}/constants.h
  ${INCLUDE_DIR}/cstream.h
  ${INCLUDE_DIR}/formatter.h
  ${INCLUDE_DIR}/frame.h
  ${INCLUDE_DIR}/fresolver.h
  ${INCLUDE_DIR}/fstream.h
  ${INCLUDE_DIR}/fwd.h
  ${INCLUDE_DIR}/half.h
  ${INCLUDE_DIR}/kdtree.h
  ${INCLUDE_DIR}/lock.h
  ${INCLUDE_DIR}/logger.h
  ${INCLUDE_DIR}/lrucache.h
  ${INCLUDE_DIR}/math.h
  ${INCLUDE_DIR}/matrix.h
  ${INCLUDE_DIR}/matrix.inl
  ${INCLUDE_DIR}/mempool.h
  ${INCLUDE_DIR}/mmap.h
  ${INCLUDE_DIR}/mstream.h
  ${INCLUDE_DIR}/netobject.h
  ${INCLUDE_DIR}/normal.h
  ${INCLUDE_DIR}/object.h
  ${INCLUDE_DIR}/octree.h
  ${INCLUDE_DIR}/platform.h
  ${INCLUDE_DIR}/plugin.h
  ${INCLUDE_DIR}/pmf.h
  ${INCLUDE_DIR}/point.h
  ${INCLUDE_DIR}/properties.h
  ${INCLUDE_DIR}/qmc.h
  ${INCLUDE_DIR}/quad.h
  ${INCLUDE_DIR}/quat.h
  ${INCLUDE_DIR}/random.h
  ${INCLUDE_DIR}/ray.h
  ${INCLUDE_DIR}/ray_sse.h
  ${INCLUDE_DIR}/ref.h
  ${INCLUDE_DIR}/rfilter.h
  ${INCLUDE_DIR}/sched.h
  ${INCLUDE_DIR}/sched_remote.h
  ${INCLUDE_DIR}/serialization.h
  ${INCLUDE_DIR}/sfcurve.h
  ${INCLUDE_DIR}/shvector.h
  ${INCLUDE_DIR}/spectrum.h
  ${INCLUDE_DIR}/spline.h
  ${INCLUDE_DIR}/sse.h
  ${INCLUDE_DIR}/ssemath.h
  ${INCLUDE_DIR}/ssevector.h
  ${INCLUDE_DIR}/sshstream.h
  ${INCLUDE_DIR}/sstream.h
  ${INCLUDE_DIR}/statistics.h
  ${INCLUDE_DIR}/stream.h
  ${INCLUDE_DIR}/thread.h
  ${INCLUDE_DIR}/timer.h
  ${INCLUDE_DIR}/tls.h
  ${INCLUDE_DIR}/track.h
  ${INCLUDE_DIR}/transform.h
  ${INCLUDE_DIR}/triangle.h
  ${INCLUDE_DIR}/util.h
  ${INCLUDE_DIR}/vector.h
  ${INCLUDE_DIR}/version.h
  ${INCLUDE_DIR}/vmf.h
  ${INCLUDE_DIR}/warp.h
  ${INCLUDE_DIR}/zstream.h
)

set(SRCS
  aabb.cpp
  appender.cpp
  bitmap.cpp
  brent.cpp
  chisquare.cpp
  class.cpp
  cstream.cpp
  fmtconv.cpp
  formatter.cpp
  fresolver.cpp
  fstream.cpp
  half.cpp
  lock.cpp
  logger.cpp
  math.cpp
  mmap.cpp
  mstream.cpp
  object.cpp
  plugin.cpp
  properties.cpp
  qmc.cpp
  quad.cpp
  random.cpp
  rfilter.cpp
  sched.cpp
  sched_remote.cpp
  serialization.cpp
  shvector.cpp
  spectrum.cpp
  spline.cpp
  ssemath.cpp
  sshstream.cpp
  sstream.cpp
  statistics.cpp
  stream.cpp
  thread.cpp
  timer.cpp
  tls.cpp
  track.cpp
  transform.cpp
  triangle.cpp
  util.cpp
  vmf.cpp
  warp.cpp
  zstream.cpp
)

add_definitions(-DMTS_BUILD_MODULE=MTS_MODULE_CORE)

if (APPLE)
  list(APPEND SRCS platform_darwin.mm)
elseif (WIN32)
  list(APPEND HDRS ${INCLUDE_DIR}/getopt.h)
  list(APPEND SRCS getopt.c)
endif()

include_directories(
  ${ZLIB_INCLUDE_DIRS} ${PNG_INCLUDE_DIR} ${JPEG_INCLUDE_DIR}
  ${ILMBASE_INCLUDE_DIRS} ${OPENEXR_INCLUDE_DIRS})
  
set(LIBS ${ZLIB_LIBRARIES} ${PNG_LIBRARIES} ${JPEG_LIBRARIES}
  ${ILMBASE_LIBRARIES} ${OPENEXR_LIBRARIES}
  ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS})
  
if(MTS_FFTW)
  include_directories(${FFTW3_INCLUDE_DIRS})
  list(APPEND LIBS ${FFTW3_LIBRARIES})
endif()
  
if(WIN32)
  list(APPEND LIBS ws2_32 psapi)
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  list(APPEND LIBS rt)
elseif (APPLE)
  list(APPEND LIBS ${COCOA_LIBRARY})
endif()
add_mts_corelib(mitsuba-core ${MTS_MAIN_HEADER} ${HDRS} ${SRCS}
  LINK_LIBRARIES ${LIBS})
add_dependencies(mitsuba-core ior_database microfacet_data)

mts_install_headers (${MTS_MAIN_HEADER} FOLDER "mitsuba")
mts_install_headers (${HDRS} FOLDER "mitsuba/core")
