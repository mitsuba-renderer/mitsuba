BUILDDIR        = '#build/release'
DISTDIR         = '#dist'
CXX             = 'cl'
CC              = 'cl'
# /O2=optimize for speed, global optimizations, intrinsic functions, favor fast code, frame pointer omission
# /EHsc=C++ exceptions, /fp:fast=Enable reasonable FP optimizations, /GS-=No buffer security checks, /GL=whole program optimizations
# To include debug information add '/Z7' to CXXFLAGS and '/DEBUG' to LINKFLAGS
CXXFLAGS        = ['/nologo', '/O2', '/fp:fast', '/D', 'WIN32', '/D', 'WIN64', '/W3', '/EHsc', '/GS-', '/GL', '/MD', '/D', 'MTS_DEBUG', '/D', 'SINGLE_PRECISION', '/D', 'SPECTRUM_SAMPLES=3', '/D', 'MTS_SSE', '/D', 'MTS_HAS_COHERENT_RT', '/D', '_CONSOLE', '/D', 'NDEBUG', '/D', 'OPENEXR_DLL', '/openmp']
SHCXXFLAGS      = CXXFLAGS
TARGET_ARCH     = 'x86_64'
MSVC_VERSION    = '10.0'
LINKFLAGS       = ['/nologo', '/SUBSYSTEM:CONSOLE', '/MACHINE:X64', '/FIXED:NO', '/OPT:REF', '/OPT:ICF', '/LTCG', '/NODEFAULTLIB:LIBCMT', '/MANIFEST']
BASEINCLUDE     = ['#include', '#dependencies/include']
BASELIB         = ['msvcrt', 'ws2_32', 'Half']
BASELIBDIR      = ['#dependencies/lib/x64']
OEXRINCLUDE     = ['#dependencies/include/openexr']
OEXRLIB         = ['IlmImf', 'IlmThread', 'Iex', 'zdll']
BOOSTLIB        = ['boost_system-vc100-mt-1_50', 'boost_filesystem-vc100-mt-1_50', 'boost_thread-vc100-mt-1_50']
COLLADAINCLUDE  = ['#dependencies/include/collada-dom', '#dependencies/include/collada-dom/1.4']
COLLADALIB      = ['libcollada14dom23']
XERCESLIB       = ['xerces-c_3']
PNGLIB          = ['libpng15']
JPEGLIB         = ['jpeg']
GLLIB           = ['opengl32', 'glu32', 'glew32mx', 'gdi32', 'user32']
GLFLAGS         = ['/D', 'GLEW_MX']
SHLIBPREFIX     = 'lib'
SHLIBSUFFIX     = '.dll'
LIBSUFFIX       = '.lib'
PROGSUFFIX      = '.exe'
PYTHON27LIB     = ['boost_python-vc100-mt-1_50', 'python27']
PYTHON27INCLUDE = ['#dependencies/include/python27']
PYTHON32LIB     = ['boost_python3-vc100-mt-1_50', 'python32']
PYTHON32INCLUDE = ['#dependencies/include/python32']
QTINCLUDE       = ['#dependencies/qt/include']
QTDIR           = '#dependencies/qt/x64'
