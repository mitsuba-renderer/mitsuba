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
env = SConscript('config/SConscript.configure')

Export('env')

if sys.platform == 'win32':
	# Set an application icon on Windows
	resources += [ env.RES('data/windows/mitsuba_res.rc') ]

# ===== Build the support libraries ====

# Core support library
SConscript('src/libcore/SConscript')
# Rendering-related APIs
SConscript('src/librender/SConscript')
# Hardware acceleration
SConscript('src/libhw/SConscript')

# ===== Build the applications =====
env = env.Clone()

# Build the command-line binaries
mainEnv = SConscript('src/mitsuba/SConscript')

# Build the COLLADA converter
converter_objects = SConscript('src/converter/SConscript', ['mainEnv'])

# Build the Qt-based GUI binaries
SConscript('src/qtgui/SConscript', ['mainEnv', 'converter_objects'])

# ===== Build the plugins =====

env['SHLIBPREFIX']=''
Export('env')

# Utilities 
SConscript('src/utils/SConscript')
# Surface scattering models
SConscript('src/bsdfs/SConscript')
# Phase functions
SConscript('src/phase/SConscript')
# Intersection shapes
SConscript('src/shapes/SConscript')
# Sample generators
SConscript('src/samplers/SConscript')
# Reconstruction filters
SConscript('src/rfilters/SConscript')
# Film implementations
SConscript('src/films/SConscript')
# Cameras
SConscript('src/cameras/SConscript')
# Participating media
SConscript('src/medium/SConscript')
# Volumetric data sources
SConscript('src/volume/SConscript')
# Sub-surface integrators
SConscript('src/subsurface/SConscript')
# Texture types
SConscript('src/textures/SConscript')
# Light sources
SConscript('src/luminaires/SConscript')
# Integrators
SConscript('src/integrators/SConscript')
# Testcases
SConscript('src/tests/SConscript')

# ===== Move everything to its proper place =====
SConscript('config/SConscript.install')
