# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

import bpy, os, copy, subprocess, math, mathutils
from extensions_framework import util as efutil
from ..outputs import MtsLog

# From collada_internal.cpp

translate_start_name_map = list(map(chr, [
	95,  95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	65,  66,  67,  68,  69,  70,  71,  72,
	73,  74,  75,  76,  77,  78,  79,  80,
	81,  82,  83,  84,  85,  86,  87,  88,
	89,  90,  95,  95,  95,  95,  95,  95,
	97,  98,  99,  100,  101,  102,  103,  104,
	105,  106,  107,  108,  109,  110,  111,  112,
	113,  114,  115,  116,  117,  118,  119,  120,
	121,  122,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  192,
	193,  194,  195,  196,  197,  198,  199,  200,
	201,  202,  203,  204,  205,  206,  207,  208,
	209,  210,  211,  212,  213,  214,  95,  216,
	217,  218,  219,  220,  221,  222,  223,  224,
	225,  226,  227,  228,  229,  230,  231,  232,
	233,  234,  235,  236,  237,  238,  239,  240,
	241,  242,  243,  244,  245,  246,  95,  248,
	249,  250,  251,  252,  253,  254,  255]))

translate_name_map = list(map(chr, [
	95,  95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  45,  95,  95,  48,
	49,  50,  51,  52,  53,  54,  55,  56,
	57,  95,  95,  95,  95,  95,  95,  95,
	65,  66,  67,  68,  69,  70,  71,  72,
	73,  74,  75,  76,  77,  78,  79,  80,
	81,  82,  83,  84,  85,  86,  87,  88,
	89,  90,  95,  95,  95,  95,  95,  95,
	97,  98,  99,  100,  101,  102,  103,  104,
	105,  106,  107,  108,  109,  110,  111,  112,
	113,  114,  115,  116,  117,  118,  119,  120,
	121,  122,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  183,  95,
	95,  95,  95,  95,  95,  95,  95,  192,
	193,  194,  195,  196,  197,  198,  199,  200,
	201,  202,  203,  204,  205,  206,  207,  208,
	209,  210,  211,  212,  213,  214,  95,  216,
	217,  218,  219,  220,  221,  222,  223,  224,
	225,  226,  227,  228,  229,  230,  231,  232,
	233,  234,  235,  236,  237,  238,  239,  240,
	241,  242,  243,  244,  245,  246,  95,  248,
	249,  250,  251,  252,  253,  254,  255]))

def translate_id(name):
	# Doesn't handle duplicates at the moment
	result = ""
	if len(name) == 0:
		return name
	result += translate_start_name_map[ord(name[0])]
	for i in range(1, len(name)):
		result += translate_name_map[ord(name[i])]
	return result

class ParamSetItem(list):
	type		= None
	type_name	= None
	name		= None
	value		= None

	def __init__(self, *args):
		self.type, self.name, self.value = args
		self.type_name = "%s %s" % (self.type, self.name)
		self.append(self.type_name)
		self.append(self.value)
	
	def export(self, exporter):
		if self.type == "color":
			exporter.parameter('rgb', self.name,
				{ 'value' : "%s %s %s" % (self.value[0], self.value[1], self.value[2])})
		elif self.type == "point" or self.type == "vector":
			exporter.parameter(self.type, self.name,
				{ 'value' : "%s %s %s" % (self.value[0], self.value[1], self.value[2])})
		elif self.type == "integer" or self.type == "float" \
				or self.type ==	"string":
			exporter.parameter(self.type, self.name, { 'value' : "%s" % self.value })
	
	def export_ref(self, exporter):
		if self.type == "reference_texture" or self.type == "reference_material" or self.type == 'reference_medium':
			if self.name != "":
				exporter.element('ref', {'id' : translate_id(self.value), 'name' : self.name})
			else:
				exporter.element('ref', {'id' : translate_id(self.value)})

class ParamSet(list):
	names = []
	
	def update(self, other):
		for p in other:
			self.add(p.type, p.name, p.value)
		return self
	
	def add(self, type, name, value):
		if name in self.names:
			for p in self:
				if p.name == name:
					self.remove(p)
		
		self.append(
			ParamSetItem(type, name, value)
		)
		self.names.append(name)
		return self
	
	def add_float(self, name, value):
		self.add('float', name, value)
		return self
	
	def add_integer(self, name, value):
		self.add('integer', name, value)
		return self

	def add_reference(self, type, name, value):
		self.add('reference_%s' % type, name, value)
		return self

	def add_bool(self, name, value):
		self.add('bool', name, bool(value))
		return self

	def add_string(self, name, value):
		self.add('string', name, str(value))
		return self
	
	def add_vector(self, name, value):
		self.add('vector', name, [i for i in value])
		return self
	
	def add_point(self, name, value):
		self.add('point', name, [p for p in value])
		return self
	
	def add_color(self, name, value):
		self.add('color', name, [c for c in value])
		return self
	
	def export(self, exporter):
		for item in self:
			item.export(exporter)
		for item in self:
			item.export_ref(exporter)
	
def get_instance_materials(ob):
	obmats = []
	# Grab materials attached to object instances ...
	if hasattr(ob, 'material_slots'):
		for ms in ob.material_slots:
			obmats.append(ms.material)
	# ... and to the object's mesh data
	if hasattr(ob.data, 'materials'):
		for m in ob.data.materials:
			obmats.append(m)
	return obmats

def resolution(scene):
	'''
	scene		bpy.types.scene
	Calculate the output render resolution
	Returns		tuple(2) (floats)
	'''
	xr = scene.render.resolution_x * scene.render.resolution_percentage / 100.0
	yr = scene.render.resolution_y * scene.render.resolution_percentage / 100.0
	
	return xr, yr

def MtsLaunch(mts_path, commandline):
	env = copy.copy(os.environ)
	mts_render_libpath = os.path.join(mts_path, "src/librender")
	mts_core_libpath = os.path.join(mts_path, "src/libcore")
	mts_hw_libpath = os.path.join(mts_path, "src/libhw")
	mts_bidir_libpath = os.path.join(mts_path, "src/libbidir")
	env['LD_LIBRARY_PATH'] = mts_core_libpath + ":" + mts_render_libpath + ":" + mts_hw_libpath + ":" + mts_bidir_libpath
	commandline[0] = os.path.join(mts_path, commandline[0])
	return subprocess.Popen(commandline, env = env, cwd = mts_path)

class MtsExporter:
	'''
		Exports the scene using COLLADA and write additional information
		to an "adjustments" file. Thim mechanism is used to capture 
		any information that gets lost in translation when using the 
		Blender COLLADA exporter.
	'''

	def __init__(self, directory, filename, materials = None, textures = None):
		mts_basename = os.path.join(directory, filename)
		(path, ext) = os.path.splitext(mts_basename)
		if ext == '.xml':
			mts_basename = path
		self.dae_filename = mts_basename + ".dae"
		self.xml_filename = mts_basename + ".xml"
		self.adj_filename = mts_basename + "_adjustments.xml"
		self.meshes_dir = os.path.join(directory, "meshes")
		self.exported_materials = []
		self.exported_textures = []
		self.materials = materials if materials != None else bpy.data.materials
		self.textures = textures if textures != None else bpy.data.textures
		self.indent = 0
		self.stack = []

	def writeHeader(self):
		try:
			self.out = open(self.adj_filename, 'w')
		except IOError:
			MtsLog('Error: unable to write to file \"%s\"!' % self.adj_filename)
			return False
		self.out.write('<?xml version="1.0" encoding="utf-8"?>\n');
		self.openElement('scene')
		return True

	def writeFooter(self):
		self.closeElement()
		self.out.close()

	def openElement(self, name, attributes = {}):
		self.out.write('\t' * self.indent + '<%s' % name)
		for (k, v) in attributes.items():
			self.out.write(' %s=\"%s\"' % (k, v))
		self.out.write('>\n')
		self.indent = self.indent+1
		self.stack.append(name)

	def closeElement(self):
		self.indent = self.indent-1
		name = self.stack.pop()
		self.out.write('\t' * self.indent + '</%s>\n' % name)

	def element(self, name, attributes = {}):
		self.out.write('\t' * self.indent + '<%s' % name)
		for (k, v) in attributes.items():
			self.out.write(' %s=\"%s\"' % (k, v))
		self.out.write('/>\n')

	def parameter(self, paramType, paramName, attributes = {}):
		self.out.write('\t' * self.indent + '<%s name="%s"' % (paramType, paramName))
		for (k, v) in attributes.items():
			self.out.write(' %s=\"%s\"' % (k, v))
		self.out.write('/>\n')
	
	def exportWorldTrafo(self, trafo):
		self.openElement('transform', {'name' : 'toWorld'})
		value = ""
		for j in range(0,4):
			for i in range(0,4):
				value += "%f " % trafo[i][j]
		self.element('matrix', {'value' : value})
		self.closeElement()

	def exportLamp(self, lamp, idx):
		ltype = lamp.data.type
		name = translate_id(lamp.data.name)
		mult = lamp.data.mitsuba_lamp.intensity
		if ltype == 'POINT':
			self.openElement('luminaire', { 'type' : 'point', 'id' : '%s-light' % name })
			self.exportWorldTrafo(lamp.matrix_world)
			self.parameter('rgb', 'intensity', {'value' : 
				"%f %f %f" % (lamp.data.color.r*mult, lamp.data.color.g*mult,
					lamp.data.color.b*mult)})
			self.parameter('float', 'samplingWeight', {'value' : '%f' % lamp.data.mitsuba_lamp.samplingWeight})
			self.closeElement()
		elif ltype == 'AREA':
			self.element('remove', { 'id' : '%s-light' % name})
			self.openElement('shape', { 'type' : 'obj'} )
			(size_x, size_y) = (lamp.data.size, lamp.data.size)
			if lamp.data.shape == 'RECTANGLE':
				size_y = lamp.data.size_y
			mult = mult / (2 * size_x * size_y)
			filename = "area_luminaire_%d.obj" % idx
			try:
				os.mkdir(self.meshes_dir)
			except OSError:
				pass
			self.parameter('string', 'filename', { 'value' : 'meshes/%s' % filename})
			self.exportWorldTrafo(lamp.matrix_world)

			self.openElement('luminaire', { 'id' : '%s-arealight' % name, 'type' : 'area'})
			self.parameter('rgb', 'intensity', { 'value' : "%f %f %f"
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult)})
			self.closeElement()
			self.openElement('bsdf', { 'type' : 'lambertian'})
			self.parameter('spectrum', 'reflectance', {'value' : '0'})
			self.closeElement()
			self.closeElement()
			path = os.path.join(self.meshes_dir, filename)
			objFile = open(path, 'w')
			objFile.write('v %f %f 0\n' % (-size_x/2, -size_y/2))
			objFile.write('v %f %f 0\n' % ( size_x/2, -size_y/2))
			objFile.write('v %f %f 0\n' % ( size_x/2,  size_y/2))
			objFile.write('v %f %f 0\n' % (-size_x/2,  size_y/2))
			objFile.write('f 4 3 2 1\n')
			objFile.close()
		elif ltype == 'SUN':
			self.openElement('luminaire', { 'id' : '%s-light' % name, 'type' : 'directional'})
			scale = mathutils.Matrix.Scale(-1, 4, mathutils.Vector([0, 0, 1]))
			self.exportWorldTrafo(lamp.matrix_world * mathutils.Matrix.Scale(-1, 4, mathutils.Vector([0, 0, 1])))
			self.parameter('rgb', 'intensity', { 'value' : "%f %f %f"
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult)})
			self.parameter('float', 'samplingWeight', {'value' : '%f' % lamp.data.mitsuba_lamp.samplingWeight})
			self.closeElement()
		elif ltype == 'SPOT':
			self.openElement('luminaire', { 'id' : '%s-light' % name, 'type' : 'spot'})
			self.exportWorldTrafo(lamp.matrix_world * mathutils.Matrix.Scale(-1, 4, mathutils.Vector([0, 0, 1])))
			self.parameter('rgb', 'intensity', { 'value' : "%f %f %f"
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult)})
			self.parameter('float', 'cutoffAngle', {'value' : '%f' %  (lamp.data.spot_size * 180 / (math.pi * 2))})
			self.parameter('float', 'beamWidth', {'value' : '%f' % (lamp.data.spot_blend * lamp.data.spot_size * 180 / (math.pi * 2))})
			self.parameter('float', 'samplingWeight', {'value' : '%f' % lamp.data.mitsuba_lamp.samplingWeight})
			self.closeElement()
		elif ltype == 'HEMI':
			if lamp.data.mitsuba_lamp.envmap_type == 'constant':
				self.openElement('luminaire', { 'id' : '%s-light' % name, 'type' : 'constant'})
				self.parameter('float', 'samplingWeight', {'value' : '%f' % lamp.data.mitsuba_lamp.samplingWeight})
				self.parameter('rgb', 'intensity', { 'value' : "%f %f %f"
						% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult)})
				self.closeElement()
			elif lamp.data.mitsuba_lamp.envmap_type == 'envmap':
				self.openElement('luminaire', { 'id' : '%s-light' % name, 'type' : 'envmap'})
				self.parameter('string', 'filename', {'value' : efutil.filesystem_path(lamp.data.mitsuba_lamp.envmap_file)})
				self.exportWorldTrafo(lamp.matrix_world)
				self.parameter('float', 'intensityScale', {'value' : '%f' % lamp.data.mitsuba_lamp.intensity})
				self.parameter('float', 'samplingWeight', {'value' : '%f' % lamp.data.mitsuba_lamp.samplingWeight})
				self.closeElement()

	def exportIntegrator(self, integrator):
		self.openElement('integrator', { 'id' : 'integrator', 'type' : integrator.type})
		self.closeElement()

	def exportSampler(self, sampler):
		self.openElement('sampler', { 'id' : 'sampler', 'type' : sampler.type})
		self.parameter('integer', 'sampleCount', { 'value' : '%i' % sampler.sampleCount})
		self.closeElement()

	def findTexture(self, name):
		if name in self.textures:
			return self.textures[name]
		else:
			raise Exception('Failed to find texture "%s"' % name)
	
	
	def findMaterial(self, name):
		if name in self.materials:
			return self.materials[name]
		else:
			raise Exception('Failed to find material "%s" in "%s"' % (name,
				str(self.materials)))

	def exportTexture(self, tex):
		if tex.name in self.exported_textures:
			return
		self.exported_textures += [tex.name]
		params = tex.mitsuba_texture.get_params()

		for p in params:
			if p.type == 'reference_texture':
				self.exportTexture(self.findTexture(p.value))

		self.openElement('texture', {'id' : '%s' % translate_id(tex.name), 'type' : tex.mitsuba_texture.type})
		params.export(self)
		self.closeElement()

	def exportMaterial(self, mat):
		if not hasattr(mat, 'name') or mat.name in self.exported_materials:
			return
		self.exported_materials += [mat.name]
		mmat = mat.mitsuba_material
		params = mmat.get_params()
		twosided = False

		if mmat.twosided and mmat.type in ['lambertian', 'phong', 'ward', 
				'mirror', 'roughmetal', 'microfacet', 'composite']:
			twosided = True

		for p in params:
			if p.type == 'reference_material':
				self.exportMaterial(self.findMaterial(p.value))
			elif p.type == 'reference_texture':
				self.exportTexture(self.findTexture(p.value))

		if twosided:
			self.openElement('bsdf', {'id' : '%s-material' % translate_id(mat.name), 'type' : 'twosided'})
			self.openElement('bsdf', {'type' : mmat.type})
		else:
			self.openElement('bsdf', {'id' : '%s-material' % translate_id(mat.name), 'type' : mmat.type})

		params.export(self)
		self.closeElement()
		
		if twosided:
			self.closeElement()

	def exportEmission(self, obj):
			lamp = obj.data.materials[0].mitsuba_emission
			if obj.data.users > 1:
				MtsLog("Error: luminaires cannot be instantiated!")
				return
			mult = lamp.intensity
			name = translate_id(obj.data.name) + "-mesh_0"
			self.openElement('append', { 'id' : name})
			self.openElement('luminaire', { 'id' : '%s-emission' % name, 'type' : 'area'})
			self.parameter('float', 'samplingWeight', {'value' : '%f' % lamp.samplingWeight})
			self.parameter('rgb', 'intensity', { 'value' : "%f %f %f"
					% (lamp.color.r*mult, lamp.color.g*mult, lamp.color.b*mult)})
			self.closeElement()
			self.closeElement()

	def exportPreviewMesh(self, material):
		self.openElement('shape', {'id' : 'Exterior-mesh_0', 'type' : 'serialized'})
		self.parameter('string', 'filename', {'value' : 'matpreview.serialized'})
		self.parameter('integer', 'shapeIndex', {'value' : '1'})
		self.openElement('transform', {'name' : 'toWorld'})
		self.element('matrix', {'value' : '0.614046 0.614047 0 -1.78814e-07 -0.614047 0.614046 0 2.08616e-07 0 0 0.868393 1.02569 0 0 0 1'})
		self.closeElement()
		self.element('ref', {'name' : 'bsdf', 'id' : '%s-material' % translate_id(material.name)})
		lamp = material.mitsuba_emission
		if lamp and lamp.use_emission:
			mult = lamp.intensity
			self.openElement('luminaire', {'type' : 'area'})
			self.parameter('rgb', 'intensity', { 'value' : "%f %f %f"
					% (lamp.color.r*mult, lamp.color.g*mult, lamp.color.b*mult)})
			self.closeElement()
		self.closeElement()

	def exportCameraSettings(self, scene, camera):
		if scene.mitsuba_integrator.motionblur:
			frameTime = 1.0/scene.render.fps
			shuttertime = scene.mitsuba_integrator.shuttertime
			shutterOpen = (scene.frame_current - shuttertime/2) * frameTime
			shutterClose = (scene.frame_current + shuttertime/2) * frameTime
			self.openElement('prepend', {'id' : '%s-camera' % translate_id(camera.name)})
			self.parameter('float', 'shutterOpen', {'value' : str(shutterOpen)})
			self.parameter('float', 'shutterClose', {'value' : str(shutterClose)})
			self.closeElement()

	def exportMedium(self, medium):
		self.openElement('medium', {'id' : medium.name, 'type' : medium.type})
		if medium.g == 0:
			self.element('phase', {'type' : 'isotropic'})
		else:
			self.openElement('phase', {'type' : 'hg'})
			self.parameter('float', 'g', {'value' : str(medium.g)})
			self.closeElement()
		self.closeElement()

	def export(self, scene):
		if scene.mitsuba_engine.binary_path == '':
			MtsLog("Error: the Mitsuba binary path was not specified!")
			return False

		idx = 0
		# Force scene update; NB, scene.update() doesn't work
		scene.frame_set(scene.frame_current)
		efutil.export_path = self.xml_filename
	
		MtsLog('MtsBlend: Writing adjustments file to "%s"' % self.adj_filename)
		if not self.writeHeader():
			return False

		self.exportIntegrator(scene.mitsuba_integrator)
		self.exportSampler(scene.mitsuba_sampler)
		for medium in scene.mitsuba_media.media:
			self.exportMedium(medium)
		for obj in scene.objects:
			if obj.type == 'LAMP':
				self.exportLamp(obj, idx)
			elif obj.type == 'MESH':
				for mat in obj.data.materials:
					self.exportMaterial(mat)
				if len(obj.data.materials) > 0 and obj.data.materials[0] != None and obj.data.materials[0].mitsuba_emission.use_emission:
					self.exportEmission(obj)
			elif obj.type == 'CAMERA':
				self.exportCameraSettings(scene, obj)
			idx = idx+1
		self.writeFooter()
		(width, height) = resolution(scene)
		
		MtsLog('MtsBlend: Writing COLLADA file to "%s"' % self.dae_filename)
		scene.collada_export(self.dae_filename)

		MtsLog("MtsBlend: Launching mtsimport")
		command = ['mtsimport', '-r', '%dx%d' % (width, height),
			'-n', '-l', 'pngfilm', self.dae_filename, self.xml_filename, self.adj_filename]
		if scene.mitsuba_integrator.motionblur:
			command += ['-z']
		process = MtsLaunch(scene.mitsuba_engine.binary_path, command);
		if process.wait() != 0:
			return False
		return True
