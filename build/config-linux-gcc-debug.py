import os, sys

BUILDDIR       = '#build/debug'
DISTDIR        = '#dist'
CXX            = 'g++'
CC             = 'gcc'
CXXFLAGS       = ['-O0', '-Wall', '-g', '-pipe', '-march=nocona', '-msse2', '-ftree-vectorize', '-mfpmath=sse', '-funsafe-math-optimizations', '-fno-rounding-math', '-fno-signaling-nans', '-fno-math-errno', '-fno-omit-frame-pointer', '-DMTS_DEBUG', '-DSINGLE_PRECISION', '-DSPECTRUM_SAMPLES=3', '-DMTS_SSE', '-DMTS_HAS_COHERENT_RT', '-fopenmp', '-fvisibility=hidden']
LINKFLAGS      = []
SHLINKFLAGS    = ['-rdynamic', '-shared', '-fPIC', '-lstdc++']
BASEINCLUDE    = ['#include']
BASELIB        = ['dl', 'm', 'pthread', 'gomp']
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

