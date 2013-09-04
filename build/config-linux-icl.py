import os, sys

BUILDDIR       = '#build/release'
DISTDIR        = '#dist'
CXX            = 'icpc'
CC             = 'icc'
CXXFLAGS       = ['-O3', '-Wall', '-g', '-pipe', '-O3', '-ipo', '-no-prec-div', '-xSSE3', '-fp-model', 'fast=2', '-openmp', '-mfpmath=sse', '-march=nocona', '-fno-math-errno', '-fomit-frame-pointer', '-DMTS_DEBUG', '-DSINGLE_PRECISION', '-DSPECTRUM_SAMPLES=3', '-DMTS_SSE', '-DMTS_HAS_COHERENT_RT', '-fopenmp', '-fvisibility=hidden', '-std=c++0x', '-wd2928', '-Qoption,cpp,--rvalue_ctor_is_not_copy_ctor']
LINKFLAGS      = []
SHLINKFLAGS    = ['-rdynamic', '-shared', '-fPIC', '-lstdc++']
BASEINCLUDE    = ['#include']
BASELIB        = ['dl', 'pthread', 'iomp5']
EIGENINCLUDE   = ['/usr/include/eigen3']
OEXRINCLUDE    = ['/usr/include/OpenEXR']
OEXRLIB        = ['Half', 'IlmImf', 'z']
PNGLIB         = ['png']
JPEGLIB        = ['jpeg']
XERCESINCLUDE  = []
XERCESLIB      = ['xerces-c']
GLLIB          = ['GL', 'GLU', 'GLEWmx', 'Xxf86vm', 'X11']
GLFLAGS        = ['-DGLEW_MX']
BOOSTLIB       = ['boost_system', 'boost_filesystem', 'boost_thread']
COLLADAINCLUDE = ['/usr/include/collada-dom', '/usr/include/collada-dom/1.4']
COLLADALIB     = ['collada14dom', 'xml2']

# The following assumes that the Mitsuba bindings should be built for the
# "default" Python version. It is also possible to build bindings for multiple
# versions at the same time by explicitly specifying e.g. PYTHON27INCLUDE,
# PYTHON27LIB, PYTHON27LIBDIR and PYTHON32INCLUDE, PYTHON32LIB, PYTHON32LIBDIR

pyver = os.popen("python --version 2>&1 | grep -oE '([[:digit:]].[[:digit:]])'").read().strip().replace('.', '')
env = locals()

env['PYTHON'+pyver+'INCLUDE']  = []
env['PYTHON'+pyver+'LIB']      = ['boost_python3' if pyver[0] == '3' else 'boost_python']

for entry in os.popen("python-config --cflags --libs").read().split():
	if entry[:2] == '-I':
		env['PYTHON'+pyver+'INCLUDE'] += [entry[2:]]
	if entry[:2] == '-l':
		env['PYTHON'+pyver+'LIB'] += [entry[2:]]
