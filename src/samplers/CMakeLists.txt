# Sample generators
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_sampler)
  add_mts_plugin(${ARGN} TYPE sampler)
endmacro()

add_sampler(independent independent.cpp)
add_sampler(stratified  stratified.cpp)
add_sampler(halton      halton.cpp faure.h faure.cpp)
add_sampler(hammersley  hammersley.cpp faure.h faure.cpp)
add_sampler(ldsampler   ldsampler.cpp)
add_sampler(sobol       sobol.cpp sobolseq.h sobolseq.cpp)
