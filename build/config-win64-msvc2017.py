BUILDDIR        = '#build/release'
DISTDIR         = '#dist'
CXX             = 'cl'
CC              = 'cl'
# /O2=optimize for speed, global optimizations, intrinsic functions, favor fast code, frame pointer omission
# /EHsc=C++ exceptions, /fp:fast=Enable reasonable FP optimizations, /GS-=No buffer security checks, /GL=whole program optimizations
# To include debug information add '/Z7' to CXXFLAGS and '/DEBUG' to LINKFLAGS
# EDIT: MSVC2017 bugs Mitsuba when using /O2 because of /Ob2 - replaced it with the full set and /Ob2 for /Ob1
# EDIT: /Og is deprecated too
CXXFLAGS        = ['/nologo', #'/Og',
'/Oi', '/Ot', '/Oy', '/Ob1', #'/Ob2',
'/Gs', '/GF', '/Gy', #'/O2',
'/fp:fast', '/D', 'WIN32', '/D', 'WIN64', '/W3', '/EHsc', '/GS-', '/GL', '/MD', '/D', 'MTS_DEBUG', '/D', 'SINGLE_PRECISION', '/D', 'SPECTRUM_SAMPLES=3', '/D', 'MTS_SSE', '/D', 'MTS_HAS_COHERENT_RT', '/D', '_CONSOLE', '/D', 'NDEBUG', '/D', 'OPENEXR_DLL', '/openmp']
SHCXXFLAGS      = CXXFLAGS
TARGET_ARCH     = 'x86_64'
MSVC_VERSION    = '15.0'
LINKFLAGS       = ['/nologo', '/SUBSYSTEM:CONSOLE', '/MACHINE:X64', '/FIXED:NO', '/OPT:REF', '/OPT:ICF', '/LTCG', '/NODEFAULTLIB:LIBCMT', '/MANIFEST']
BASEINCLUDE     = ['#include', '#dependencies/include']
BASELIB         = ['msvcrt', 'ws2_32', 'Half', 'zlib']
BASELIBDIR      = ['#dependencies/lib']
OEXRINCLUDE     = ['#dependencies/include/openexr']
OEXRLIB         = ['IlmImf', 'IlmThread', 'Iex', 'zlib']
BOOSTLIB        = ['boost_system-vc141-mt-1_64', 'boost_filesystem-vc141-mt-1_64', 'boost_thread-vc141-mt-1_64']
# COLLADAINCLUDE  = ['#dependencies/include/collada-dom', '#dependencies/include/collada-dom/1.4']
# COLLADALIB      = ['libcollada14dom24']
XERCESLIB       = ['xerces-c_3']
PNGLIB          = ['libpng16']
JPEGLIB         = ['jpeg']
GLLIB           = ['opengl32', 'glu32', 'glew32mx', 'gdi32', 'user32']
GLFLAGS         = ['/D', 'GLEW_MX']
SHLIBPREFIX     = 'lib'
SHLIBSUFFIX     = '.dll'
LIBSUFFIX       = '.lib'
PROGSUFFIX      = '.exe'
PYTHON27LIB     = ['boost_python27-vc141-mt-1_64', 'python27']
PYTHON27INCLUDE = ['#dependencies/include/python27']
PYTHON35LIB     = ['boost_python35-vc141-mt-1_64', 'python35']
PYTHON35INCLUDE = ['#dependencies/include/python35']
PYTHON36LIB     = ['boost_python36-vc141-mt-1_64', 'python36']
PYTHON36INCLUDE = ['#dependencies/include/python36']
QTINCLUDE       = ['#dependencies/include']
QTDIR           = '#dependencies'
FFTWLIB         = ['libfftw3-3']
