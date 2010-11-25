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

import os, math, mathutils

import bpy

from extensions_framework import util as efutil
from mitsuba.export import translate_id

class MtsAdjustments:
	'''
		Writes scene information to an "adjustments" file. This
		mechanism is used to capture information which gets lost
		in translation when using the COLLADA exporter.
	'''

	def __init__(self, target_file, target_dir, materials = None, textures = None):
		self.target_file = target_file
		self.target_dir = target_dir
		self.exported_materials = []
		self.exported_textures = []
		self.materials = materials if materials != None else bpy.data.materials
		self.textures = textures if textures != None else bpy.data.textures

	def exportWorldtrafo(self, trafo):
		self.out.write('\t\t<transform name="toWorld">\n')
		self.out.write('\t\t\t<matrix value="')
		for j in range(0,4):
			for i in range(0,4):
				self.out.write("%f " % trafo[i][j])
		self.out.write('"/>\n\t\t</transform>\n')

	def exportLamp(self, lamp, idx):
		ltype = lamp.data.mitsuba_lamp.type
		name = translate_id(lamp.data.name)
		if ltype == 'POINT':
			self.out.write('\t<luminaire id="%s-light" type="point">\n' % name)
			mult = lamp.data.mitsuba_lamp.intensity
			self.exportWorldtrafo(lamp.matrix_world)
			self.out.write('\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			self.out.write('\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.samplingWeight)
			self.out.write('\t</luminaire>\n')
		elif ltype == 'AREA':
			self.out.write('\t<remove id="%s-light"/>\n' % name)
			self.out.write('\t<shape type="obj">\n')
			size_x = lamp.data.size
			size_y = lamp.data.size
			if lamp.data.shape == 'RECTANGLE':
				size_y = lamp.data.size_y
			mts_meshes_dir = os.path.join(self.target_dir, 'meshes')
			filename = "area_luminaire_%d.obj" % idx

			self.out.write('\t\t<string name="filename" value="meshes/%s"/>\n' % filename)
			self.exportWorldtrafo(lamp.matrix_world)
			self.out.write('\n\t\t<luminaire id="%s-arealight" type="area">\n' % name)
			mult = lamp.data.mitsuba_lamp.intensity / (2 * size_x * size_y)
			self.out.write('\t\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			self.out.write('\t\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.samplingWeight)
			self.out.write('\t\t</luminaire>\n')
			self.out.write('\t</shape>\n')
			
			try:
				os.mkdir(mts_meshes_dir)
			except OSError:
				pass
			path = os.path.join(mts_meshes_dir, filename)
			objFile = open(path, 'w')
			objFile.write('v %f %f 0\n' % (-size_x/2, -size_y/2))
			objFile.write('v %f %f 0\n' % ( size_x/2, -size_y/2))
			objFile.write('v %f %f 0\n' % ( size_x/2,  size_y/2))
			objFile.write('v %f %f 0\n' % (-size_x/2,  size_y/2))
			objFile.write('f 4 3 2 1\n')
			objFile.close()
		elif ltype == 'SUN':
			self.out.write('\t<luminaire id="%s-light" type="directional">\n' % name)
			mult = lamp.data.mitsuba_lamp.intensity
			scale = mathutils.Matrix.Scale(-1, 4, mathutils.Vector([0, 0, 1]))
			self.exportWorldtrafo(lamp.matrix_world * mathutils.Matrix.Scale(-1, 4, mathutils.Vector([0, 0, 1])))
			self.out.write('\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			self.out.write('\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.samplingWeight)
			self.out.write('\t</luminaire>\n')
		elif ltype == 'SPOT':
			self.out.write('\t<luminaire id="%s-light" type="spot">\n' % name)
			mult = lamp.data.mitsuba_lamp.intensity
			self.exportWorldtrafo(lamp.matrix_world * mathutils.Matrix.Scale(-1, 4, mathutils.Vector([0, 0, 1])))
			self.out.write('\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			self.out.write('\t\t<float name="cutoffAngle" value="%f"/>\n' % (lamp.data.spot_size * 180 / (math.pi * 2)))
			self.out.write('\t\t<float name="beamWidth" value="%f"/>\n' % (lamp.data.spot_blend * lamp.data.spot_size * 180 / (math.pi * 2)))
			self.out.write('\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.samplingWeight)
			self.out.write('\t</luminaire>\n')
		elif ltype == 'ENV':
			if lamp.data.mitsuba_lamp.envmap_type == 'constant':
				self.out.write('\t<luminaire id="%s-light" type="constant">\n' % name)
				mult = lamp.data.mitsuba_lamp.intensity
				self.out.write('\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
						% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
				self.out.write('\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.samplingWeight)
				self.out.write('\t</luminaire>\n')
			elif lamp.data.mitsuba_lamp.envmap_type == 'envmap':
				self.out.write('\t<luminaire id="%s-light" type="envmap">\n' % name)
				self.out.write('\t\t<string name="filename" value="%s"/>\n' % efutil.filesystem_path(lamp.data.mitsuba_lamp.envmap_file))
				self.exportWorldtrafo(lamp.matrix_world)
				self.out.write('\t\t<float name="intensityScale" value="%f"/>\n' % lamp.data.mitsuba_lamp.intensity)
				self.out.write('\t</luminaire>\n')

	def exportIntegrator(self, integrator):
		self.out.write('\t<integrator id="integrator" type="%s">\n' % integrator.type)
		self.out.write('\t</integrator>\n')

	def exportSampler(self, sampler):
		self.out.write('\t<sampler id="sampler" type="%s">\n' % sampler.type)
		self.out.write('\t\t<integer name="sampleCount" value="%i"/>\n' % sampler.sampleCount)
		self.out.write('\t</sampler>\n')

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

	def exportTexture(self, mat):
		if mat.name in self.exported_textures:
			return
		self.exported_textures += [mat.name]
		params = mat.mitsuba_texture.get_params()

		for p in params:
			if p.type == 'reference_texture':
				self.exportTexture(self.findTexture(p.value))

		self.out.write('\t<texture id="%s" type="%s">\n' % (translate_id(mat.name), mat.mitsuba_texture.type))
		self.out.write(params.to_string())
		self.out.write(params.to_string_ref())
		self.out.write('\t</texture>\n')

	def exportMaterial(self, mat):
		if mat.name in self.exported_materials:
			return
		self.exported_materials += [mat.name]
		params = mat.mitsuba_material.get_params()

		for p in params:
			if p.type == 'reference_material':
				self.exportMaterial(self.findMaterial(p.value))
			elif p.type == 'reference_texture':
				self.exportTexture(self.findTexture(p.value))

		self.out.write('\t<bsdf id="%s" type="%s">\n' % (translate_id(mat.name), mat.mitsuba_material.type))
		self.out.write(params.to_string())
		self.out.write(params.to_string_ref())
		self.out.write('\t</bsdf>\n')

	def exportEmission(self, obj):
			lamp = obj.data.materials[0].mitsuba_emission
			name = translate_id(obj.data.name)
			self.out.write('\t<append id="%s-mesh_0">\n' % name)
			self.out.write('\t\t<luminaire id="%s-emission" type="area">\n' % name)
			mult = lamp.intensity
			self.out.write('\t\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.color.r*mult, lamp.color.g*mult, lamp.color.b*mult))
			self.out.write('\t\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.samplingWeight)
			self.out.write('\t\t</luminaire>\n')
			self.out.write('\t</append>\n')

	def writeHeader(self):
		self.out = open(self.target_file, 'w')
		self.out.write('<scene>\n');

	def writeFooter(self):
		self.out.write('</scene>\n');
		self.out.close()

	def exportPreviewMesh(self, material):
		self.out.write('\t\t<shape id="Exterior-mesh_0" type="serialized">\n')
		self.out.write('\t\t\t<string name="filename" value="matpreview.serialized"/>\n')
		self.out.write('\t\t\t<integer name="shapeIndex" value="1"/>\n')
		self.out.write('\t\t\t<transform name="toWorld">\n')
		self.out.write('\t\t\t\t<matrix value="0.614046 0.614047 0 -1.78814e-07 -0.614047 0.614046 0 2.08616e-07 0 0 0.868393 1.02569 0 0 0 1"/>\n')
		self.out.write('\t\t\t</transform>\n')
		self.out.write('\t\t\t<ref id="%s" name="bsdf"/>\n' % translate_id(material.name))
		lamp = material.mitsuba_emission
		if lamp and lamp.use_emission:
			mult = lamp.intensity
			self.out.write('\t\t\t<luminaire type="area">\n')
			self.out.write('\t\t\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.color.r*mult, lamp.color.g*mult, lamp.color.b*mult))
			self.out.write('\t\t\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.samplingWeight)
			self.out.write('\t\t\t</luminaire>\n')
		self.out.write('\t\t</shape>\n')
		self.out.write('\n')

	def exportCameraSettings(self, scene, camera):
		if scene.mitsuba_integrator.motionblur:
			frameTime = 1.0/scene.render.fps
			shuttertime = scene.mitsuba_integrator.shuttertime
			shutterOpen = (scene.frame_current - shuttertime/2) * frameTime
			shutterClose = (scene.frame_current + shuttertime/2) * frameTime
			self.out.write('\t<prepend id="%s-camera">\n' % translate_id(camera.name))
			self.out.write('\t\t<float name="shutterOpen" value="%f"/>\n' % shutterOpen)
			self.out.write('\t\t<float name="shutterClose" value="%f"/>\n' % shutterClose)
			self.out.write('\t</prepend>\n')

	def export(self, scene):
		idx = 0
		self.writeHeader()
		self.exportIntegrator(scene.mitsuba_integrator)
		self.exportSampler(scene.mitsuba_sampler)
		for obj in scene.objects:
			if obj.type == 'LAMP':
				self.exportLamp(obj, idx)
			elif obj.type == 'MESH':
				for mat in obj.data.materials:
					self.exportMaterial(mat)
				if len(obj.data.materials) > 0 and obj.data.materials[0].mitsuba_emission.use_emission:
					self.exportEmission(obj)
			elif obj.type == 'CAMERA':
				self.exportCameraSettings(scene, obj)
			idx = idx+1
		self.writeFooter()

	
