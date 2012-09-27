# mitsuba-bidir
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/mitsuba/bidir)
set(HDRS
  ${INCLUDE_DIR}/common.h
  ${INCLUDE_DIR}/edge.h
  ${INCLUDE_DIR}/geodist2.h
  ${INCLUDE_DIR}/manifold.h
  ${INCLUDE_DIR}/mempool.h
  ${INCLUDE_DIR}/mut_bidir.h
  ${INCLUDE_DIR}/mut_caustic.h
  ${INCLUDE_DIR}/mut_lens.h
  ${INCLUDE_DIR}/mut_manifold.h
  ${INCLUDE_DIR}/mut_mchain.h
  ${INCLUDE_DIR}/mutator.h
  ${INCLUDE_DIR}/path.h
  ${INCLUDE_DIR}/pathsampler.h
  ${INCLUDE_DIR}/rsampler.h
  ${INCLUDE_DIR}/util.h
  ${INCLUDE_DIR}/vertex.h
)

# Common sources
set(SRCS
  common.cpp
  edge.cpp
  manifold.cpp
  mut_bidir.cpp
  mut_caustic.cpp
  mut_lens.cpp
  mut_manifold.cpp
  mut_mchain.cpp
  path.cpp
  pathsampler.cpp
  rsampler.cpp
  util.cpp
  verification.cpp
  vertex.cpp
)

add_definitions(-DMTS_BUILD_MODULE=MTS_MODULE_BIDIR)
 
add_mts_corelib(mitsuba-bidir ${HDRS} ${SRCS}
  LINK_LIBRARIES mitsuba-core mitsuba-render)
  
mts_install_headers (${HDRS} FOLDER "mitsuba/bidir")
