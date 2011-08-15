import SCons
import sys
import glob
import os

resources = []
plugins = []
stubs = []

Export('SCons', 'sys', 'os', 'glob', 'resources', 
	'plugins', 'stubs')

# Configure the build framework
env = SConscript('build/SConscript.configure')

Export('env')

if sys.platform == 'win32':
	# Set an application icon on Windows
	resources += [ env.RES('data/windows/mitsuba_res.rc') ]

def build(scriptFile, exports = [], duplicate = 0):
	dirname = '/'.join(os.path.dirname(scriptFile).split('/')[1:])
	return SConscript(scriptFile, exports, 
		variant_dir=os.path.join(env['BUILDDIR'], dirname), duplicate=duplicate)

# ===== Build the support libraries ====

# Core support library
build('src/libcore/SConscript')
# Rendering-related APIs
build('src/librender/SConscript')
# Hardware acceleration
build('src/libhw/SConscript')
# Bidirectional support library
build('src/libbidir/SConscript')
# Python binding library
build('src/libpython/SConscript')

# ===== Build the applications =====
env = env.Clone()

# Build the command-line binaries
mainEnv = build('src/mitsuba/SConscript')

# Build the COLLADA converter
converter_objects = build('src/converter/SConscript', ['mainEnv'])

# Build the Qt-based GUI binaries
build('src/mtsgui/SConscript', ['mainEnv', 'converter_objects'], duplicate=True)

env['SHLIBPREFIX']=''

# ===== Build the plugins =====

Export('env')

# Utilities 
build('src/utils/SConscript')
# Surface scattering models
build('src/bsdfs/SConscript')
# Phase functions
build('src/phase/SConscript')
# Intersection shapes
build('src/shapes/SConscript')
# Sample generators
build('src/samplers/SConscript')
# Reconstruction filters
build('src/rfilters/SConscript')
# Film implementations
build('src/films/SConscript')
# Cameras
build('src/cameras/SConscript')
# Participating media
build('src/medium/SConscript')
# Volumetric data sources
build('src/volume/SConscript')
# Sub-surface integrators
build('src/subsurface/SConscript')
# Texture types
build('src/textures/SConscript')
# Light sources
build('src/luminaires/SConscript')
# Integrators
build('src/integrators/SConscript')
# Testcases
build('src/tests/SConscript')

# ===== Move everything to its proper place =====
SConscript('build/SConscript.install')
