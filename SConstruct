import SCons
import sys
import glob
import os
import fnmatch
import multiprocessing

AddOption("--cfg", dest="cfg", type="string", nargs=1, action='store', help='Manually specify a configuration file')
configFile = GetOption('cfg')
if configFile == None:
	configFile = "config.py"

if not os.path.exists(configFile):
	print 'A configuration file must be selected! Have a look at \"README\"'
	Exit(1)

# Parse configuration options
vars = Variables(configFile);
vars.Add('CXX',           'C++ compiler')
vars.Add('CC',            'C compiler')
vars.Add('CXXFLAGS',      'C++ flags')
vars.Add('CCFLAGS',       'C compiler flags')
vars.Add('STRIP',         'Program for stripping away unused symbols')
vars.Add('SHCXXFLAGS',    'C++ flags (for shared libraries)')
vars.Add('LINK',          'Linker')
vars.Add('LINKFLAGS',     'Linker flags')
vars.Add('SHLINKFLAGS',   'Linker flags (dynamic libraries)')
vars.Add('BASEINCLUDE',   'Base include path')
vars.Add('BASELIB',       'Base libraries')
vars.Add('BASELIBDIR',    'Base library search path')
vars.Add('XERCESINCLUDE', 'Xerces-C include path')
vars.Add('XERCESLIB',     'Xerces-C libraries')
vars.Add('XERCESLIBDIR',  'Xerces-C library path')
vars.Add('OEXRINCLUDE',   'OpenEXR include path')
vars.Add('OEXRLIB',       'OpenEXR libraries')
vars.Add('OEXRFLAGS',     'OpenEXR-related compiler flags')
vars.Add('OEXRLIBDIR',    'OpenEXR library path')
vars.Add('PNGINCLUDE',    'libpng include path')
vars.Add('PNGLIB',        'libpng libraries')
vars.Add('PNGLIBDIR',     'libpng library path')
vars.Add('JPEGINCLUDE',   'libjpeg include path')
vars.Add('JPEGLIB',       'libjpeg libraries')
vars.Add('JPEGLIBDIR',    'libjpeg library path')
vars.Add('COLLADAINCLUDE','COLLADA DOM include path')
vars.Add('COLLADALIB',    'COLLADA DOM libraries')
vars.Add('COLLADALIBDIR', 'COLLADA DOM library path')
vars.Add('SHLIBPREFIX',   'Prefix for shared libraries')
vars.Add('SHLIBSUFFIX',   'Suffix for shared libraries')
vars.Add('PROGSUFFIX',    'Suffix for executables')
vars.Add('GLLIB',         'OpenGL+GLEW libraries')
vars.Add('GLINCLUDE',     'OpenGL+GLEW include path')
vars.Add('GLFLAGS',       'OpenGL+GLEW-related compiler flags')
vars.Add('GLLIBDIR',      'OpenGL+GLEW library path')
vars.Add('BOOSTINCLUDE',  'boost include path')
vars.Add('BOOSTLIB',      'boost libraries')
vars.Add('BOOSTLIBDIR',   'boost library path')
vars.Add('TARGET_ARCH',   'Target architecture')
vars.Add('MSVC_VERSION',  'MS Visual C++ compiler version')

try:
	env = Environment(options=vars, ENV = os.environ, tools=['default', 'qt4'], toolpath=['tools'])
	print 'Checking for Qt 4.x... yes'
	hasQt = True
except Exception:
	env = Environment(options=vars, ENV = os.environ, tools=['default'], toolpath=['tools'])
	print 'Unable to detect a Qt installation -- not building the GUI!'
	hasQt = False

hasCollada=True
env.Append(CPPPATH=env['BASEINCLUDE'])
env.Append(CPPFLAGS=[])
env.Append(LIBPATH=[])
env.Append(LIBS=env['BASELIB'])
if env.has_key('BOOSTINCLUDE'):
	env.Prepend(CPPPATH=env['BOOSTINCLUDE'])
if env.has_key('BOOSTLIBDIR'):
	env.Prepend(LIBPATH=env['BOOSTLIBDIR'])
if env.has_key('BOOSTLIB'):
	env.Prepend(LIBS=env['BOOSTLIB'])
if env.has_key('BASELIBDIR'):
	env.Append(LIBPATH=env['BASELIBDIR'])

env.Decider('MD5-timestamp')
env.SetOption('implicit_cache', 1)

#env.SetOption('num_jobs', multiprocessing.cpu_count())

AddOption("--dist", dest="dist", type="string", nargs=0, action='store', help='Make an official release')

# Check whether everything important is available
def CheckCXX(context):
	context.Message('Checking for ' + env['CXX'] + ' ...')
	ret = context.TryLink("#include <sstream>\n int main(int argc, char **argv) {\n std::ostringstream oss;\n return 0;\n }", '.cpp')
	context.Result(ret)
	return ret

# For running Uic & Moc (below)
def recursiveDirs(root) :
    return filter( (lambda a : a.rfind(".svn") == -1),  [ a[0] for a in os.walk(root)])

def unique(list) :
    return dict.fromkeys(list).keys()

def scanFiles(dir, accept=["*.cpp"], reject=[]) :
    sources = []
    paths = recursiveDirs(dir)
    for path in paths :
        for pattern in accept :
            sources+=glob.glob(path+"/"+pattern)
    for pattern in reject :
        sources = filter( (lambda a : a.rfind(pattern)==-1 ),  sources )
    return unique(sources)

conf = Configure(env, custom_tests = { 'CheckCXX' : CheckCXX })
cppPathPrevious = SCons.Util.semi_deepcopy(env['CPPPATH'])
libPathPrevious = SCons.Util.semi_deepcopy(env['LIBPATH'])
cppFlagsPrevious = SCons.Util.semi_deepcopy(env['CPPFLAGS'])
cxxFlagsPrevious = SCons.Util.semi_deepcopy(env['CXXFLAGS'])

if env.has_key('PNGINCLUDE'):
	env.Prepend(CPPPATH=env['PNGINCLUDE'])
if env.has_key('PNGLIBDIR'):
	env.Prepend(LIBPATH=env['PNGLIBDIR'])
if env.has_key('JPEGINCLUDE'):
	env.Prepend(CPPPATH=env['JPEGINCLUDE'])
if env.has_key('JPEGLIBDIR'):
	env.Prepend(LIBPATH=env['JPEGLIBDIR'])
if env.has_key('OEXRFLAGS'):
	env.Prepend(CPPFLAGS=env['OEXRFLAGS'])
if env.has_key('OEXRINCLUDE'):
	env.Prepend(CPPPATH=env['OEXRINCLUDE'])
if env.has_key('OEXRLIBDIR'):
	env.Prepend(LIBPATH=env['OEXRLIBDIR'])
if env.has_key('XERCESINCLUDE'):
	env.Prepend(CPPPATH=env['XERCESINCLUDE'])
if env.has_key('XERCESLIBDIR'):
	env.Prepend(LIBPATH=env['XERCESLIBDIR'])
if env.has_key('GLINCLUDE'):
	env.Prepend(CPPPATH=env['GLINCLUDE'])
if env.has_key('GLFLAGS'):
	env.Prepend(CPPFLAGS=env['GLFLAGS'])
if env.has_key('COLLADAINCLUDE'):
	env.Prepend(CPPPATH=env['COLLADAINCLUDE'])
if env.has_key('COLLADALIBDIR'):
	env.Prepend(LIBPATH=env['COLLADALIBDIR'])

if not conf.CheckCXX():
	print 'Could not compile a simple C++ fragment, verify that ' + env['CXX'] + ' is installed!'
	Exit(1)
if not conf.CheckCHeader(['png.h']):
	print 'libpng is missing (install libpng12-dev)'
	Exit(1)
if not conf.CheckCHeader(['stdio.h', 'jpeglib.h']):
	print 'libjpeg is missing (install libjpeg62-dev)'
	Exit(1)
if not conf.CheckCXXHeader('ImfRgba.h'):
	print 'OpenEXR is missing (install libopenexr-dev)'
	Exit(1)
if not conf.CheckCXXHeader('xercesc/dom/DOMLSParser.hpp'):
	print 'Xerces-C++ 3.x must be installed (install libxerces-c-dev)!'
	Exit(1)
if not conf.CheckCXXHeader('dae.h'):
	hasCollada = False
	print 'COLLADA DOM is missing: not building the COLLADA importer'
if not conf.CheckCXXHeader('boost/math/distributions/students_t.hpp'):
	print 'Boost is missing (install libboost1.40-dev and libboost-math1.40-dev)!'
	Exit(1)
if sys.platform == 'win32':
	if not (conf.CheckCHeader(['windows.h', 'GL/gl.h']) and conf.CheckCHeader(['windows.h', 'GL/glu.h']) and conf.CheckCHeader(['windows.h', 'GL/gl.h', 'GL/glext.h'])):
		print 'OpenGL headers are missing!'
		Exit(1)
	if not conf.CheckCHeader('GL/glew.h'):
		print 'GLEW headers are missing!'
		Exit(1)
elif sys.platform == 'linux2':
	if not (conf.CheckCHeader('GL/gl.h') and conf.CheckCHeader('GL/glu.h') and conf.CheckCHeader(['GL/gl.h', 'GL/glext.h'])):
		print 'OpenGL headers are missing!'
		Exit(1)
	if not conf.CheckCHeader('GL/glew.h'):
		print 'GLEW headers are missing (install libglewmx1.5-dev)!'
		Exit(1)
	if not conf.CheckType('GLEWContext', '#include <GL/glew.h>'):
		print 'GLEW-MX must be present!'
		Exit(1)
	if not conf.TryCompile("#include <GL/glew.h>\n int i = GL_VERTEX_ATTRIB_ARRAY_UNIFIED_NV;", '.cpp'):
		print 'Your version of GLEW-MX seems to be outdated!'
		Exit(1)
elif sys.platform == 'darwin':
	if not (conf.CheckCHeader('OpenGL/gl.h') and conf.CheckCHeader('OpenGL/glu.h') and conf.CheckCHeader(['OpenGL/gl.h', 'OpenGL/glext.h'])):
		print 'OpenGL headers are missing!'
		Exit(1)
	if not conf.CheckCHeader('OpenGL/glew.h'):
		print 'GLEW headers are missing!'
		Exit(1)
if sys.platform == 'linux2':
	if not (conf.CheckCHeader(['X11/Xlib.h', 'X11/extensions/xf86vmode.h'])):
		print 'X Video Mode selection library headers are missing! (Install libxxf86vm-dev)'
		Exit(1)
env.Replace(CPPPATH=cppPathPrevious)
env.Replace(LIBPATH=libPathPrevious)
env.Replace(CPPFLAGS=cppFlagsPrevious)
env.Replace(CXXFLAGS=cxxFlagsPrevious)
sys.stdout.write("Checking for Mitsuba version .. ")
file = open(env.GetBuildPath('include/mitsuba/mitsuba.h'), 'r')
MTS_VERSION=""
for line in file:
	if line.startswith("#define MTS_VERSION "):
		MTS_VERSION = line[21:len(line)-2]
if MTS_VERSION == "":
	print 'could not be determined!'
	Exit(1)
else:
	print MTS_VERSION
env = conf.Finish()
dist = GetOption('dist') != None

def stripinst_build_function(self, target, source, pkgname = None, use_own = None):
	inst = self.Install(target, source)
	self.AddPostAction(inst, env['STRIP'] + ' $TARGET')
	return inst

def osxlibinst_build_function(self, target, source, pkgname = None, use_own = None):
	inst = self.Install(target, source)
	prefix, name = os.path.split(source)
	self.AddPostAction(inst, 'install_name_tool -id @executable_path/../Frameworks/' + name + ' $TARGET')
	return inst

env.__class__.StripInst = stripinst_build_function
env.__class__.OSXLibInst = osxlibinst_build_function

if hasCollada:
	env.Append(CPPDEFINES = [['MTS_HAS_COLLADA', 1]] )

env.SConsignFile()

# MSVC: Embed the manifest
if sys.platform == 'win32':
	env['LINKCOM'] = [env['LINKCOM'], 'mt.exe -nologo -manifest ${TARGET}.manifest -outputresource:$TARGET;1']
	env['SHLINKCOM'] = [env['SHLINKCOM'], 'mt.exe -nologo -manifest ${TARGET}.manifest -outputresource:$TARGET;2']

try:
	os.mkdir('plugins')
except:
	pass

# Core library
coreEnv = env.Clone()
if coreEnv.has_key('OEXRLIBDIR'):
	coreEnv.Prepend(LIBPATH=env['OEXRLIBDIR'])
if coreEnv.has_key('OEXRINCLUDE'):
	coreEnv.Prepend(CPPPATH=env['OEXRINCLUDE'])
if coreEnv.has_key('OEXRFLAGS'):
	coreEnv.Prepend(CPPFLAGS=env['OEXRFLAGS'])
if coreEnv.has_key('OEXRLIB'):
	coreEnv.Prepend(LIBS=env['OEXRLIB'])
if coreEnv.has_key('PNGLIBDIR'):
	coreEnv.Prepend(LIBPATH=env['PNGLIBDIR'])
if coreEnv.has_key('PNGINCLUDE'):
	coreEnv.Prepend(CPPPATH=env['PNGINCLUDE'])
if coreEnv.has_key('PNGLIB'):
	coreEnv.Prepend(LIBS=env['PNGLIB'])
if coreEnv.has_key('JPEGLIBDIR'):
	coreEnv.Prepend(LIBPATH=env['JPEGLIBDIR'])
if coreEnv.has_key('JPEGINCLUDE'):
	coreEnv.Prepend(CPPPATH=env['JPEGINCLUDE'])
if coreEnv.has_key('JPEGLIB'):
	coreEnv.Prepend(LIBS=env['JPEGLIB'])

coreEnv.Prepend(CPPDEFINES = [['MTS_BUILD_MODULE', 'MTS_MODULE_CORE']])
libcore_objects = [
	'src/libcore/class.cpp', 'src/libcore/object.cpp', 
	'src/libcore/statistics.cpp', 'src/libcore/thread.cpp',
	'src/libcore/logger.cpp', 'src/libcore/appender.cpp',
	'src/libcore/formatter.cpp', 'src/libcore/lock.cpp', 
	'src/libcore/random.cpp', 'src/libcore/timer.cpp', 
	'src/libcore/util.cpp', 'src/libcore/properties.cpp', 
	'src/libcore/transform.cpp', 'src/libcore/spectrum.cpp', 
	'src/libcore/aabb.cpp', 'src/libcore/stream.cpp', 
	'src/libcore/fstream.cpp', 'src/libcore/plugin.cpp', 
	'src/libcore/triangle.cpp', 'src/libcore/bitmap.cpp',
	'src/libcore/serialization.cpp', 'src/libcore/sstream.cpp',
	'src/libcore/cstream.cpp', 'src/libcore/mstream.cpp', 
	'src/libcore/sched.cpp', 'src/libcore/sched_remote.cpp',
	'src/libcore/sshstream.cpp', 'src/libcore/wavelet.cpp',
	'src/libcore/zstream.cpp', 'src/libcore/shvector.cpp',
	'src/libcore/fresolver.cpp'
]

if sys.platform == 'darwin':
	coreEnv_osx = coreEnv.Clone();
	coreEnv_osx['CXXFLAGS'].remove('-fstrict-aliasing');
	coreEnv_osx['CXXFLAGS'].remove('-ftree-vectorize');
	coreEnv_osx['CXXFLAGS'].append('-fno-strict-aliasing');
	libcore_objects += coreEnv_osx.SharedObject('src/libcore/platform_darwin.mm')
elif sys.platform == 'win32':
	libcore_objects += coreEnv.SharedObject('src/libcore/getopt.c')
	libcore_objects += coreEnv.SharedObject('src/libcore/platform_win32.cpp')

libcore = coreEnv.SharedLibrary('src/libcore/mitsuba-core', libcore_objects);

if sys.platform == "darwin":
	coreEnv.AddPostAction(libcore, 'install_name_tool -id @executable_path/../Frameworks/libmitsuba-core.dylib $TARGET')

env = env.Clone()
env.Append(LIBS=['mitsuba-core'])
env.Append(LIBPATH=['src/libcore'])

# Rendering-specific library
renderEnv = env.Clone()
renderEnv.Append(CPPDEFINES = [['MTS_BUILD_MODULE', 'MTS_MODULE_RENDER']] )
if renderEnv.has_key('XERCESINCLUDE'):
	renderEnv.Prepend(CPPPATH=renderEnv['XERCESINCLUDE'])
if renderEnv.has_key('XERCESLIBDIR'):
	renderEnv.Prepend(LIBPATH=renderEnv['XERCESLIBDIR'])
if renderEnv.has_key('XERCESLIB'):
	renderEnv.Prepend(LIBS=renderEnv['XERCESLIB'])
librender = renderEnv.SharedLibrary('src/librender/mitsuba-render', [
	'src/librender/bsdf.cpp', 'src/librender/camera.cpp',
	'src/librender/film.cpp', 'src/librender/integrator.cpp',
	'src/librender/kdtree.cpp', 
	'src/librender/luminaire.cpp', 'src/librender/medium.cpp', 
	'src/librender/renderjob.cpp', 'src/librender/imageproc.cpp', 
	'src/librender/imageproc_wu.cpp', 'src/librender/renderproc.cpp', 
	'src/librender/renderproc_wr.cpp', 'src/librender/particleproc.cpp',
	'src/librender/renderqueue.cpp', 'src/librender/mipmap.cpp', 
	'src/librender/scene.cpp',  'src/librender/subsurface.cpp', 
	'src/librender/texture.cpp', 'src/librender/shape.cpp', 
	'src/librender/trimesh.cpp', 'src/librender/rfilter.cpp', 
	'src/librender/sampler.cpp', 'src/librender/util.cpp',
	'src/librender/irrcache.cpp', 'src/librender/testcase.cpp',
	'src/librender/preview.cpp', 'src/librender/photonmap.cpp',
	'src/librender/gatherproc.cpp', 'src/librender/mipmap3d.cpp',
	'src/librender/volume.cpp', 'src/librender/vpl.cpp',
	'src/librender/shader.cpp', 'src/librender/shandler.cpp',
	'src/librender/intersection.cpp'
])

if sys.platform == "darwin":
	renderEnv.AddPostAction(librender, 'install_name_tool -id @executable_path/../Frameworks/libmitsuba-render.dylib $TARGET')

env.Append(LIBS=['mitsuba-render'])
env.Append(LIBPATH=['src/librender'])

libhw_objects = ['src/libhw/session.cpp', 'src/libhw/device.cpp',
	 'src/libhw/gputexture.cpp', 'src/libhw/gpugeometry.cpp', 
	 'src/libhw/gpuprogram.cpp', 'src/libhw/renderer.cpp', 
	 'src/libhw/glrenderer.cpp', 'src/libhw/glprogram.cpp', 
	 'src/libhw/glgeometry.cpp', 'src/libhw/gltexture.cpp',
	 'src/libhw/gpusync.cpp', 'src/libhw/glsync.cpp',
	 'src/libhw/vpl.cpp', 'src/libhw/font.cpp',
	 'src/libhw/viewer.cpp']

if sys.platform == 'win32':
	libhw_objects += ['src/libhw/wglsession.cpp',
		'src/libhw/wgldevice.cpp',
		'src/libhw/wglrenderer.cpp']
elif sys.platform == 'linux2':
	libhw_objects += ['src/libhw/x11session.cpp',
		'src/libhw/x11device.cpp',
		'src/libhw/glxdevice.cpp', 
		'src/libhw/glxrenderer.cpp']

glEnv = env.Clone()
glEnv.Append(CPPDEFINES = [['MTS_BUILD_MODULE', 'MTS_MODULE_HW']] )
if glEnv.has_key('GLLIB'):
	glEnv.Prepend(LIBS=glEnv['GLLIB'])
if glEnv.has_key('GLLIBDIR'):
	glEnv.Prepend(LIBPATH=glEnv['GLLIBDIR'])
if glEnv.has_key('GLFLAGS'):
	glEnv.Prepend(CPPFLAGS=glEnv['GLFLAGS'])
if glEnv.has_key('GLINCLUDE'):
	glEnv.Prepend(CPPPATH=glEnv['GLINCLUDE'])

if sys.platform == 'darwin':
	glEnv_osx = glEnv.Clone();
	glEnv_osx['CXXFLAGS'].remove('-fstrict-aliasing');
	glEnv_osx['CXXFLAGS'].remove('-ftree-vectorize');
	glEnv_osx['CXXFLAGS'].append('-fno-strict-aliasing');
	libhw_objects += glEnv_osx.SharedObject(['src/libhw/nsglsession.mm',
		'src/libhw/nsgldevice.mm',
		'src/libhw/nsglrenderer.mm'])

libhw = glEnv.SharedLibrary('src/libhw/mitsuba-hw', libhw_objects)
if sys.platform == "darwin":
	glEnv.AddPostAction(libhw, 'install_name_tool -id @executable_path/../Frameworks/libmitsuba-hw.dylib $TARGET')

env = env.Clone()
env.Append(LIBS=['mitsuba-hw'])
env.Append(LIBPATH=['src/libhw'])
env['SHLIBPREFIX']=''

# Environment with Xerces + wxWidgets
mainEnv = env.Clone()
if mainEnv.has_key('XERCESINCLUDE'):
	mainEnv.Prepend(CPPPATH=mainEnv['XERCESINCLUDE'])
if mainEnv.has_key('XERCESLIBDIR'):
	mainEnv.Prepend(LIBPATH=mainEnv['XERCESLIBDIR'])
if mainEnv.has_key('XERCESLIB'):
	mainEnv.Prepend(LIBS=mainEnv['XERCESLIB'])
if mainEnv.has_key('GLLIB'):
	mainEnv.Prepend(LIBS=mainEnv['GLLIB'])
if mainEnv.has_key('GLLIBDIR'):
	mainEnv.Prepend(LIBPATH=mainEnv['GLLIBDIR'])
if mainEnv.has_key('GLFLAGS'):
	mainEnv.Prepend(CXXFLAGS=mainEnv['GLFLAGS'])
if mainEnv.has_key('GLINCLUDE'):
	mainEnv.Prepend(CPPPATH=mainEnv['GLINCLUDE'])

resources = []
darwinStub = []

if sys.platform == 'win32':
	resources += [env.RES('tools/windows/mitsuba_res.rc')]

# Build the command-line+GUI interface
mainEnv.Program('mtssrv', resources + ['src/mitsuba/mtssrv.cpp'])
mainEnv.Program('mitsuba', resources + ['src/mitsuba/mitsuba.cpp'])

if sys.platform == 'darwin':
	mainEnv_osx = mainEnv.Clone();
	mainEnv_osx['CXXFLAGS'].remove('-fstrict-aliasing');
	mainEnv_osx['CXXFLAGS'].remove('-ftree-vectorize');
	mainEnv_osx['CXXFLAGS'].append('-fno-strict-aliasing');
	darwinStub += [mainEnv_osx.StaticObject('src/mitsuba/darwin_stub.mm')]

mainEnv.Program('mtsutil', resources + darwinStub + ['src/mitsuba/mtsutil.cpp'])

env.Program('src/utils/joinrgb', ['src/utils/joinrgb.cpp'])
env.Program('src/utils/ttest', ['src/utils/ttest.cpp'])
env.Program('src/utils/createvol', ['src/utils/createvol.cpp'])

# COLLADA importer
if hasCollada:
	colladaEnv = mainEnv.Clone()
	temp = colladaEnv['CXXFLAGS']
	colladaEnv['CXXFLAGS'] = temp
	colladaEnv.Append(LIBS=['mitsuba-hw'])
	colladaEnv.Append(LIBPATH=['src/libhw'])
	if env.has_key('COLLADAINCLUDE'):
		colladaEnv.Prepend(CPPPATH=env['COLLADAINCLUDE'])
	if env.has_key('COLLADALIBDIR'):
		colladaEnv.Prepend(LIBPATH=env['COLLADALIBDIR'])
	if env.has_key('COLLADALIB'):
		colladaEnv.Prepend(LIBS=env['COLLADALIB'])
	converter_objects = [
		colladaEnv.StaticObject('src/converter/collada.cpp'),
		colladaEnv.StaticObject('src/converter/obj.cpp'),
		colladaEnv.StaticObject('src/converter/converter.cpp')
	]
	colladaEnv.Program('mtsimport', darwinStub + ['src/converter/mtsimport.cpp'] 
		+ resources + converter_objects)

if hasQt:
	qtEnv = mainEnv.Clone()
	qtEnv.Append(CPPPATH=['src/qtgui'])
	qtEnv.EnableQt4Modules(['QtGui', 'QtCore', 'QtOpenGL', 'QtXml', 'QtNetwork'])
	if sys.platform == 'win32':
		index = qtEnv['CXXFLAGS'].index('_CONSOLE')
		del qtEnv['CXXFLAGS'][index-1]
		del qtEnv['CXXFLAGS'][index-1]
		index = qtEnv['LINKFLAGS'].index('/SUBSYSTEM:CONSOLE')
		del qtEnv['LINKFLAGS'][index]
		qtEnv.Append(CXXFLAGS=['/D', '_WINDOWS'])
		qtEnv.Append(LINKFLAGS=['/SUBSYSTEM:WINDOWS'])
		qtEnv.Append(LIBS=['qtmain'])
	elif sys.platform == 'darwin':
		qtEnv.Append(LINKFLAGS=['-Ftools/darwin', '-framework', 'BWToolkitFramework'])

	qtInterfaces = [qtEnv.Uic4(uic) for uic in scanFiles('src/qtgui', ['*.ui'])]
	qtResources = [qtEnv.Qrc(qrc) for qrc in scanFiles('src/qtgui', ['*.qrc'])]

	qtgui_files = scanFiles('src/qtgui', ['*.cpp']) + qtResources + resources
	
	if hasCollada:
		qtgui_files += converter_objects
		if env.has_key('COLLADALIBDIR'):
			qtEnv.Prepend(LIBPATH=env['COLLADALIBDIR'])
		if env.has_key('COLLADALIB'):
			qtEnv.Prepend(LIBS=env['COLLADALIB'])

	if sys.platform == 'darwin':
		qtEnv_osx = qtEnv.Clone();
		# Objective C++ does not permit the following optimization flags
		qtEnv_osx['CXXFLAGS'].remove('-fstrict-aliasing');
		qtEnv_osx['CXXFLAGS'].remove('-ftree-vectorize');
		qtEnv_osx['CXXFLAGS'].append('-fno-strict-aliasing');
		qtEnv_osx['CXXFLAGS'].append(['-Ftools/darwin', '-framework', 'BWToolkitFramework'])
		qtgui_files += qtEnv_osx.StaticObject('src/qtgui/previewsettingsdlg_cocoa_impl.mm')
	else:
		qtgui_files = [x for x in qtgui_files if (not isinstance(x, str) or 'cocoa' not in x)]
	qtgui = qtEnv.Program('mtsgui', qtgui_files)
	if sys.platform == 'darwin':
		qtEnv.AddPostAction(qtgui, 'install_name_tool -change QtGui.framework/Versions/4/QtGui @executable_path/../Frameworks/QtGui $TARGET')
		qtEnv.AddPostAction(qtgui, 'install_name_tool -change QtCore.framework/Versions/4/QtCore @executable_path/../Frameworks/QtCore $TARGET')
		qtEnv.AddPostAction(qtgui, 'install_name_tool -change QtOpenGL.framework/Versions/4/QtOpenGL @executable_path/../Frameworks/QtOpenGL $TARGET')
		qtEnv.AddPostAction(qtgui, 'install_name_tool -change QtXml.framework/Versions/4/QtXml @executable_path/../Frameworks/QtXml $TARGET')
		qtEnv.AddPostAction(qtgui, 'install_name_tool -change QtNetwork.framework/Versions/4/QtNetwork @executable_path/../Frameworks/QtNetwork $TARGET')

plugins = []

# Build the plugins -- utilities
plugins += env.SharedLibrary('plugins/addimages', ['src/utils/addimages.cpp'])
plugins += env.SharedLibrary('plugins/cylclip', ['src/utils/cylclip.cpp'])
plugins += env.SharedLibrary('plugins/kdbench', ['src/utils/kdbench.cpp'])

# BSDFs
plugins += env.SharedLibrary('plugins/lambertian', ['src/bsdfs/lambertian.cpp'])
plugins += env.SharedLibrary('plugins/dielectric', ['src/bsdfs/dielectric.cpp'])
plugins += env.SharedLibrary('plugins/mirror', ['src/bsdfs/mirror.cpp'])
plugins += env.SharedLibrary('plugins/transparent', ['src/bsdfs/transparent.cpp'])
plugins += env.SharedLibrary('plugins/difftrans', ['src/bsdfs/difftrans.cpp'])
plugins += env.SharedLibrary('plugins/mask', ['src/bsdfs/mask.cpp'])
plugins += env.SharedLibrary('plugins/ward', ['src/bsdfs/ward.cpp'])
plugins += env.SharedLibrary('plugins/phong', ['src/bsdfs/phong.cpp'])
plugins += env.SharedLibrary('plugins/microfacet', ['src/bsdfs/microfacet.cpp'])
plugins += env.SharedLibrary('plugins/roughglass', ['src/bsdfs/roughglass.cpp'])
plugins += env.SharedLibrary('plugins/roughmetal', ['src/bsdfs/roughmetal.cpp'])
plugins += env.SharedLibrary('plugins/composite', ['src/bsdfs/composite.cpp'])

# Phase functions
plugins += env.SharedLibrary('plugins/isotropic', ['src/phase/isotropic.cpp'])
plugins += env.SharedLibrary('plugins/hg', ['src/phase/hg.cpp'])
plugins += env.SharedLibrary('plugins/kkay', ['src/phase/kkay.cpp'])

# Shapes and triangle mesh loaders
plugins += env.SharedLibrary('plugins/obj', ['src/shapes/obj.cpp'])
plugins += env.SharedLibrary('plugins/ply', ['src/shapes/ply/ply.cpp', 'src/shapes/ply/ply_parser.cpp'],
	CPPPATH = env['CPPPATH'] + ['src/shapes/ply'])
plugins += env.SharedLibrary('plugins/serialized', ['src/shapes/serialized.cpp'])
plugins += env.SharedLibrary('plugins/sphere', ['src/shapes/sphere.cpp'])
plugins += env.SharedLibrary('plugins/cylinder', ['src/shapes/cylinder.cpp'])
plugins += env.SharedLibrary('plugins/hair', ['src/shapes/hair.cpp'])
plugins += env.SharedLibrary('plugins/shapegroup', ['src/shapes/shapegroup.cpp'])
plugins += env.SharedLibrary('plugins/instance', ['src/shapes/instance.cpp'])

# Samplers
plugins += env.SharedLibrary('plugins/independent', ['src/samplers/independent.cpp'])
plugins += env.SharedLibrary('plugins/stratified', ['src/samplers/stratified.cpp'])
plugins += env.SharedLibrary('plugins/halton', ['src/samplers/halton.cpp'])
plugins += env.SharedLibrary('plugins/hammersley', ['src/samplers/hammersley.cpp'])
plugins += env.SharedLibrary('plugins/ldsampler', ['src/samplers/ldsampler.cpp'])

# Image reconstruction filters
plugins += env.SharedLibrary('plugins/box', ['src/rfilters/box.cpp'])
plugins += env.SharedLibrary('plugins/wsinc', ['src/rfilters/wsinc.cpp'])
plugins += env.SharedLibrary('plugins/mitchell', ['src/rfilters/mitchell.cpp'])
plugins += env.SharedLibrary('plugins/catmullrom', ['src/rfilters/catmullrom.cpp'])
plugins += env.SharedLibrary('plugins/gaussian', ['src/rfilters/gaussian.cpp'])

# Films
plugins += env.SharedLibrary('plugins/exrfilm', ['src/films/exrfilm.cpp'])
plugins += env.SharedLibrary('plugins/pngfilm', ['src/films/pngfilm.cpp'])
plugins += env.SharedLibrary('plugins/mfilm', ['src/films/mfilm.cpp'])

# Cameras
plugins += env.SharedLibrary('plugins/perspective', ['src/cameras/perspective.cpp'])
plugins += env.SharedLibrary('plugins/orthographic', ['src/cameras/orthographic.cpp'])

# Participating media
plugins += env.SharedLibrary('plugins/homogeneous', ['src/medium/homogeneous.cpp'])
plugins += env.SharedLibrary('plugins/heterogeneous', ['src/medium/heterogeneous.cpp'])
plugins += env.SharedLibrary('plugins/adaptive', ['src/medium/adaptive.cpp'])
plugins += env.SharedLibrary('plugins/flake', ['src/medium/flake.cpp'])
plugins += env.SharedLibrary('plugins/heterogeneous-stencil', ['src/medium/heterogeneous-stencil.cpp'])
plugins += env.SharedLibrary('plugins/heterogeneous-flake', ['src/medium/heterogeneous-flake.cpp'])

# Volumetric data sources
plugins += env.SharedLibrary('plugins/constvolume', ['src/volume/constvolume.cpp'])
plugins += env.SharedLibrary('plugins/gridvolume', ['src/volume/gridvolume.cpp'])
plugins += env.SharedLibrary('plugins/hgridvolume', ['src/volume/hgridvolume.cpp'])

# Sub-surface integrators
plugins += env.SharedLibrary('plugins/dipole', ['src/subsurface/dipole.cpp',
	'src/subsurface/irrproc.cpp', 'src/subsurface/irrtree.cpp'])
#plugins += env.SharedLibrary('plugins/adipole', ['src/subsurface/adipole.cpp',
#	'src/subsurface/irrproc.cpp', 'src/subsurface/irrtree.cpp'])

# Texture types
plugins += env.SharedLibrary('plugins/exrtexture', ['src/textures/exrtexture.cpp'])
plugins += env.SharedLibrary('plugins/ldrtexture', ['src/textures/ldrtexture.cpp'])
plugins += env.SharedLibrary('plugins/gridtexture', ['src/textures/gridtexture.cpp'])
plugins += env.SharedLibrary('plugins/checkerboard', ['src/textures/checkerboard.cpp'])
plugins += env.SharedLibrary('plugins/vertexcolors', ['src/textures/vertexcolors.cpp'])

# Light sources
plugins += env.SharedLibrary('plugins/area', ['src/luminaires/area.cpp'])
plugins += env.SharedLibrary('plugins/constant', ['src/luminaires/constant.cpp'])
plugins += env.SharedLibrary('plugins/envmap', ['src/luminaires/envmap.cpp'])
plugins += env.SharedLibrary('plugins/spot', ['src/luminaires/spot.cpp'])
plugins += env.SharedLibrary('plugins/point', ['src/luminaires/point.cpp'])
plugins += env.SharedLibrary('plugins/collimated', ['src/luminaires/collimated.cpp'])
plugins += env.SharedLibrary('plugins/directional', ['src/luminaires/directional.cpp'])

# Integrators
plugins += env.SharedLibrary('plugins/direct', ['src/integrators/direct/direct.cpp'])
plugins += env.SharedLibrary('plugins/errctrl', ['src/integrators/misc/errctrl.cpp'])
plugins += env.SharedLibrary('plugins/path', ['src/integrators/path/path.cpp'])
plugins += env.SharedLibrary('plugins/irrcache', ['src/integrators/misc/irrcache.cpp',
	'src/integrators/misc/irrcache_proc.cpp'])
plugins += env.SharedLibrary('plugins/volpath', ['src/integrators/path/volpath.cpp'])
plugins += env.SharedLibrary('plugins/volpath_simple', ['src/integrators/path/volpath_simple.cpp'])
plugins += env.SharedLibrary('plugins/ptracer', ['src/integrators/path/ptracer.cpp',
	'src/integrators/path/ptracer_proc.cpp'])
plugins += env.SharedLibrary('plugins/photonmapper', ['src/integrators/photonmapper/photonmapper.cpp'])
plugins += env.SharedLibrary('plugins/ppm', ['src/integrators/photonmapper/ppm.cpp'])
plugins += env.SharedLibrary('plugins/sppm', ['src/integrators/photonmapper/sppm.cpp'])
plugins += env.SharedLibrary('plugins/vpl', ['src/integrators/vpl/vpl.cpp'])
	
#camsampler = env.SharedObject('src/integrators/bidir/camsampler.cpp')
#pathvertex_bdpt = env.SharedObject('src/integrators/bidir/pathvertex_bdpt',
#	'src/integrators/bidir/pathvertex.cpp', CPPDEFINES = {'MTS_METHOD' : 'BDPT'});
#path_bdpt = env.SharedObject('src/integrators/bidir/path_bdpt',
#	'src/integrators/bidir/path.cpp', CPPDEFINES = {'MTS_METHOD' : 'BDPT'});
#pathvertex_mlt = env.SharedObject('src/integrators/bidir/pathvertex_mlt',
#	'src/integrators/bidir/pathvertex.cpp', CPPDEFINES = {'MTS_METHOD' : 'MLT'});
#path_mlt = env.SharedObject('src/integrators/bidir/path_mlt',
#	'src/integrators/bidir/path.cpp', CPPDEFINES = {'MTS_METHOD' : 'MLT'});

#plugins += env.SharedLibrary('plugins/bidir', [
#	camsampler, pathvertex_bdpt, path_bdpt,
#	'src/integrators/bidir/bidir.cpp',
#	'src/integrators/bidir/bidir_proc.cpp'], CPPDEFINES = {'MTS_METHOD' : 'BDPT'})

#plugins += env.SharedLibrary('plugins/kelemen', [
#	'src/integrators/bidir/kelemen.cpp',
#	'src/integrators/bidir/kelemen_sampler.cpp',
#	'src/integrators/bidir/kelemen_proc.cpp',
#	'src/integrators/bidir/mlt_sampler.cpp',
#	pathvertex_bdpt, path_bdpt], CPPDEFINES = {'MTS_METHOD' : 'BDPT'})

#plugins += env.SharedLibrary('plugins/mlt', [
#	camsampler,
#	'src/integrators/bidir/mlt.cpp',
#	'src/integrators/bidir/mlt_sampler.cpp',
#	'src/integrators/bidir/mlt_bidir.cpp',
#	'src/integrators/bidir/mlt_lens.cpp',
#	'src/integrators/bidir/mlt_caustic.cpp',
#	'src/integrators/bidir/mlt_mchain.cpp',
#	'src/integrators/bidir/mlt_proc.cpp',
#	pathvertex_mlt, path_mlt
#])

# Testcases
testEnv = env.Clone()
testEnv.Append(CPPDEFINES = [['MTS_TESTCASE', '1']])
for plugin in glob.glob('src/tests/test_*.cpp'):
	name = os.path.basename(plugin)
	plugins += testEnv.SharedLibrary('plugins/' + name[0:len(name)-4], plugin)

installTargets = []

# Windows build?
if sys.platform == 'win32':
	try:
		os.mkdir('dist')
		os.mkdir('dist/plugins')
		os.mkdir('dist/schema')
	except:
		pass
	for plugin in plugins:
		if '.dll' in plugin.__str__():
			installTargets += env.Install('dist/plugins', plugin)
	installTargets += env.Install('dist/schema', 'schema/scene.xsd')
	if 'WIN64' in env['CXXFLAGS']:
		dllprefix='tools/windows/lib64/'
	else:
		dllprefix='tools/windows/lib32/'
	installTargets += env.Install('dist', 'mitsuba.exe')
	installTargets += env.Install('dist', 'mtssrv.exe')
	installTargets += env.Install('dist', 'mtsutil.exe')
	installTargets += env.Install('dist', 'mtsimport.exe')
	installTargets += env.Install('dist', 'mtsgui.exe')
	installTargets += env.Install('dist', 'src/libcore/libmitsuba-core.dll')
	installTargets += env.Install('dist', 'src/libhw/libmitsuba-hw.dll')
	installTargets += env.Install('dist', 'src/librender/libmitsuba-render.dll')
	installTargets += env.Install('dist', dllprefix + 'Iex.dll')
	installTargets += env.Install('dist', dllprefix + 'Half.dll')
	installTargets += env.Install('dist', dllprefix + 'IlmThread.dll')
	installTargets += env.Install('dist', dllprefix + 'Imath.dll')
	installTargets += env.Install('dist', dllprefix + 'IlmImf.dll')
	installTargets += env.Install('dist', dllprefix + 'zlib1.dll')
	installTargets += env.Install('dist', dllprefix + 'libpng13.dll')
	installTargets += env.Install('dist', dllprefix + 'jpeg62.dll')
	installTargets += env.Install('dist', dllprefix + 'pthreadVCE2.dll')
	installTargets += env.Install('dist', dllprefix + 'xerces-c_3_0.dll')
	installTargets += env.Install('dist', dllprefix + 'glew32mx.dll')
	installTargets += env.Install('dist', dllprefix + 'libcollada14dom21.dll')
	if hasQt:
		installTargets += env.Install('dist', env['QT4_BINPATH']+'/QtCore4.dll')
		installTargets += env.Install('dist', env['QT4_BINPATH']+'/QtGui4.dll')
		installTargets += env.Install('dist', env['QT4_BINPATH']+'/QtXml4.dll')
		installTargets += env.Install('dist', env['QT4_BINPATH']+'/QtNetwork4.dll')
		installTargets += env.Install('dist', env['QT4_BINPATH']+'/QtOpenGL4.dll')
elif sys.platform == 'darwin':
	try:
		os.mkdir('Mitsuba.app')
		os.mkdir('Mitsuba.app/plugins')
		os.mkdir('Mitsuba.app/schema')
		os.mkdir('Mitsuba.app/Contents')
		os.mkdir('Mitsuba.app/Contents/Frameworks')
		os.mkdir('Mitsuba.app/Contents/MacOS')
		os.nkdir('Mitsuba.app/Contents/Resources/PreviewSettings.nib')
		os.mkdir('Mitsuba.app/Contents/Resources')
	except:
		pass
	for i in plugins:
		installTargets += env.Install('Mitsuba.app/plugins', i)
	installTargets += env.Install('Mitsuba.app/schema', 'schema/scene.xsd')
	installTargets += env.Install('Mitsuba.app/Contents/MacOS', 'mtssrv')
	installTargets += env.Install('Mitsuba.app/Contents/MacOS', 'mtsutil')
	installTargets += env.Install('Mitsuba.app/Contents/MacOS', 'mitsuba')
	installTargets += env.Install('Mitsuba.app/Contents/MacOS', 'mtsimport')
	plist = env.Install('Mitsuba.app/Contents', 'tools/darwin/Info.plist')
	installTargets += plist
	installTargets += env.AddPostAction(plist, 'perl -pi -e "s/MTS_VERSION/%s/" $TARGET' % MTS_VERSION)
	installTargets += env.Install('Mitsuba.app/Contents', 'tools/darwin/PkgInfo')
	installTargets += env.Install('Mitsuba.app/Contents/Resources', 'tools/darwin/Resources/mitsuba.icns')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'src/librender/libmitsuba-render.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'src/libcore/libmitsuba-core.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'src/libhw/libmitsuba-hw.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/GLEW.framework/Resources/libs/libGLEW.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/OpenEXR.framework/Resources/lib/libHalf.6.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/OpenEXR.framework/Resources/lib/libIex.6.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/OpenEXR.framework/Resources/lib/libImath.6.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/OpenEXR.framework/Resources/lib/libIlmThread.6.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/OpenEXR.framework/Resources/lib/libIlmImf.6.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/Xerces-C.framework/Resources/lib/libxerces-c-3.0.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/libpng.framework/Resources/lib/libpng.dylib')
	installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/libjpeg.framework/Resources/lib/libjpeg.dylib')
	if hasCollada:
		installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', 'tools/darwin/Collada14Dom.framework/Resources/lib/libCollada14Dom.dylib')
	if hasQt:
		installTargets += env.Install('Mitsuba.app/Contents/MacOS', 'mtsgui')
		installTargets += env.OSXLibInst('Mitsuba.app/Contents/Frameworks', '/Library/Frameworks/QtCore.framework/Versions/4/QtCore')
		opengl = env.OSXLibInst('Mitsuba.app/Contents/Frameworks', '/Library/Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL')
		xml = env.OSXLibInst('Mitsuba.app/Contents/Frameworks', '/Library/Frameworks/QtXml.framework/Versions/4/QtXml')
		network = env.OSXLibInst('Mitsuba.app/Contents/Frameworks', '/Library/Frameworks/QtNetwork.framework/Versions/4/QtNetwork')
		gui = env.OSXLibInst('Mitsuba.app/Contents/Frameworks', '/Library/Frameworks/QtGui.framework/Versions/4/QtGui')
		installTargets += env.AddPostAction(xml, 'install_name_tool -change QtCore.framework/Versions/4/QtCore @executable_path/../Frameworks/QtCore $TARGET')
		installTargets += env.AddPostAction(network, 'install_name_tool -change QtCore.framework/Versions/4/QtCore @executable_path/../Frameworks/QtCore $TARGET')
		installTargets += env.AddPostAction(gui, 'install_name_tool -change QtCore.framework/Versions/4/QtCore @executable_path/../Frameworks/QtCore $TARGET')
		installTargets += env.AddPostAction(opengl, 'install_name_tool -change QtCore.framework/Versions/4/QtCore @executable_path/../Frameworks/QtCore $TARGET')
		installTargets += env.AddPostAction(opengl, 'install_name_tool -change QtGui.framework/Versions/4/QtGui @executable_path/../Frameworks/QtGui $TARGET')
		installTargets += env.Install('Mitsuba.app/Contents/Resources', '/Library/Frameworks//QtGui.framework/Versions/4/Resources/qt_menu.nib')
		installTargets += env.Install('Mitsuba.app/Contents/Resources/PreviewSettings.nib', 'tools/darwin/PreviewSettings.nib/designable.nib')
		installTargets += env.Install('Mitsuba.app/Contents/Resources/PreviewSettings.nib', 'tools/darwin/PreviewSettings.nib/keyedobjects.nib')
		installTargets += env.Install('Mitsuba.app/Contents/Resources', 'tools/darwin/qt.conf')
		installTargets += env.Install('Mitsuba.app/Contents/Frameworks/BWToolkitFramework.framework/Versions/A', 'tools/darwin/BWToolkitFramework.framework/Versions/A/BWToolkitFramework')
		for file in os.listdir('tools/darwin/BWToolkitFramework.framework/Versions/A/Resources'):
			if fnmatch.fnmatch(file, '*.pdf') or fnmatch.fnmatch(file, '*.tiff') or fnmatch.fnmatch(file, '*.tif') or fnmatch.fnmatch(file, '*.png') or fnmatch.fnmatch(file, '*.rtf') or fnmatch.fnmatch(file, '*.plist'):
				installTargets += env.Install('Mitsuba.app/Contents/Frameworks/BWToolkitFramework.framework/Resources', 'tools/darwin/BWToolkitFramework.framework/Versions/A/Resources/' + file)

if dist:
	if sys.platform == 'win32':
		distTarget = env.Command("Mitsuba %s.zip" % MTS_VERSION, [], "tools\\windows\\build-dist.bat %s" % MTS_VERSION)
		Depends(distTarget, installTargets)
	elif sys.platform == 'darwin':
		distTarget = env.Command("Mitsuba %s.dmg" % MTS_VERSION, [], "tools/darwin/build-dmg.sh %s" % MTS_VERSION)
		Depends(distTarget, installTargets)
	elif sys.platform == 'linux2':
		env.Command("mitsuba-%s.tar.gz" % MTS_VERSION, [], "tools/linux/build-sourcedist.sh %s" % MTS_VERSION)
