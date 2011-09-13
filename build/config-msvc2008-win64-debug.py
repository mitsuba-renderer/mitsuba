import sys, os

BUILDDIR       = '#build/debug'
DISTDIR        = '#dist'
CXX            = 'cl'
CC             = 'cl'
# /O2=optimize for speed, global optimizations, intrinsic functions, favor fast code, frame pointer omission
# /EHsc=C++ exceptions, /fp:fast=Enable reasonable FP optimizations, /GS-=No buffer security checks, /GL=whole program optimizations
# To include debug information add '/Z7' to CXXFLAGS and '/DEBUG' to LINKFLAGS
CXXFLAGS       = ['/nologo', '/Od', '/Z7', '/fp:fast', '/D', 'WIN32', '/D', 'WIN64', '/W3', '/EHsc', '/GS-', '/GL', '/MD', '/D', 'MTS_DEBUG', '/D', 'SINGLE_PRECISION', '/D', 'SPECTRUM_SAMPLES=3', '/D', 'MTS_SSE', '/D', 'MTS_HAS_COHERENT_RT', '/D', '_CONSOLE', '/D', 'NDEBUG', '/openmp']
SHCXXFLAGS     = CXXFLAGS
TARGET_ARCH    = 'x86_64'
MSVC_VERSION   = '9.0'
LINKFLAGS      = ['/nologo', '/SUBSYSTEM:CONSOLE', '/DEBUG', '/MACHINE:X64', '/FIXED:NO', '/OPT:REF', '/OPT:ICF', '/LTCG', '/NODEFAULTLIB:LIBCMT', '/MANIFEST']
BASEINCLUDE    = ['#include', '#dependencies/windows/include']
BASELIB        = ['pthreadVCE2', 'msvcrt', 'ws2_32']
OEXRINCLUDE    = ['#dependencies/windows/include/OpenEXR']
OEXRFLAGS      = ['/D', 'OPENEXR_DLL']
OEXRLIB        = ['IlmImf', 'IlmThread', 'zlib1', 'Half']
BOOSTINCLUDE   = ['#dependencies']
BOOSTLIB       = ['boost_system-vc90-mt-1_44', 'boost_filesystem-vc90-mt-1_44']
COLLADAINCLUDE = ['#dependencies/windows/include/colladadom', '#dependencies/windows/include/colladadom/1.4']
COLLADALIB     = ['libcollada14dom23']
XERCESLIB      = ['xerces-c_3']
PNGLIB         = ['libpng13']
JPEGLIB        = ['jpeg62']
GLLIB          = ['opengl32', 'glu32', 'glew32mx', 'gdi32', 'user32']
GLFLAGS        = ['/D', 'GLEW_MX']
BASELIBDIR     = ['#dependencies/windows/lib64', '#dependencies/windows/lib64/vc90']
PYTHONINCLUDE  = [os.path.join(os.path.split(sys.executable)[0], 'include')]
PYTHONLIBDIR   = [os.path.join(os.path.split(sys.executable)[0], 'libs')]
PYTHONLIB      = ['boost_python-vc90-mt-1_44', 'python26']
SHLIBPREFIX    = 'lib'
SHLIBSUFFIX    = '.dll'
LIBSUFFIX      = '.lib'
PROGSUFFIX     = '.exe'
