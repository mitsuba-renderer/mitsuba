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

class MtsAdjustments:
	'''
		Writes scene information to an "adjustments" file. This
		mechanism is used to capture information which gets lost
		in translation when using the COLLADA exporter.
	'''

	def __init__(self, target_file, target_dir):
		self.target_file = target_file
		self.target_dir = target_dir
		self.exported_materials = []
		self.exported_textures = []

	def export_worldtrafo(self, adjfile, trafo):
		adjfile.write('\t\t<transform name="toWorld">\n')
		adjfile.write('\t\t\t<matrix value="')
		for j in range(0,4):
			for i in range(0,4):
				adjfile.write("%f " % trafo[i][j])
		adjfile.write('"/>\n\t\t</transform>\n')

	def export_lamp(self, adjfile, lamp, idx):
		ltype = lamp.data.mitsuba_lamp.type
		if ltype == 'POINT':
			adjfile.write('\t<luminaire id="%s-light" type="point">\n' % lamp.data.name)
			mult = lamp.data.mitsuba_lamp.intensity
			self.export_worldtrafo(adjfile, lamp.matrix_world)
			adjfile.write('\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			adjfile.write('\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.sampling_weight)
			adjfile.write('\t</luminaire>\n')
		elif ltype == 'AREA':
			adjfile.write('\t<shape type="obj">\n')
			size_x = lamp.data.size
			size_y = lamp.data.size
			if lamp.data.shape == 'RECTANGLE':
				size_y = lamp.data.size_y
			path = os.path.join(os.path.join(self.target_dir, 'meshes'), "_area_luminaire_%d.obj" % idx)

			adjfile.write('\t\t<string name="filename" value="%s"/>\n' % path)
			self.export_worldtrafo(adjfile, lamp.matrix_world)

			adjfile.write('\n\t\t<luminaire id="%s-light" type="area">\n' % lamp.data.name)
			mult = lamp.data.mitsuba_lamp.intensity / (2 * size_x * size_y)
			adjfile.write('\t\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			adjfile.write('\t\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.sampling_weight)
			adjfile.write('\t\t</luminaire>\n')
			adjfile.write('\t</shape>\n')
			objFile = open(path, 'w')
			objFile.write('v %f %f 0\n' % (-size_x/2, -size_y/2))
			objFile.write('v %f %f 0\n' % ( size_x/2, -size_y/2))
			objFile.write('v %f %f 0\n' % ( size_x/2,  size_y/2))
			objFile.write('v %f %f 0\n' % (-size_x/2,  size_y/2))
			objFile.write('f 4 3 2 1\n')
			objFile.close()
		elif ltype == 'SUN':
			adjfile.write('\t<luminaire id="%s-light" type="directional">\n' % lamp.data.name)
			mult = lamp.data.mitsuba_lamp.intensity
			scale = mathutils.Matrix.Scale(-1, 4, mathutils.Vector([0, 0, 1]))
			self.export_worldtrafo(adjfile, lamp.matrix_world * mathutils.Matrix.Scale(-1, 4, mathutils.Vector([0, 0, 1])))
			adjfile.write('\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			adjfile.write('\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.sampling_weight)
			adjfile.write('\t</luminaire>\n')
		elif ltype == 'SPOT':
			adjfile.write('\t<luminaire id="%s-light" type="spot">\n' % lamp.data.name)
			mult = lamp.data.mitsuba_lamp.intensity
			self.export_worldtrafo(adjfile, lamp.matrix_world * mathutils.Matrix.Scale(-1, 4, mathutils.Vector([0, 0, 1])))
			adjfile.write('\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			adjfile.write('\t\t<float name="cutoffAngle" value="%f"/>\n' % (lamp.data.spot_size * 180 / (math.pi * 2)))
			adjfile.write('\t\t<float name="beamWidth" value="%f"/>\n' % (lamp.data.spot_blend * lamp.data.spot_size * 180 / (math.pi * 2)))
			adjfile.write('\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.sampling_weight)
			adjfile.write('\t</luminaire>\n')
		elif ltype == 'ENV':
			if lamp.data.mitsuba_lamp.envmap_type == 'constant':
				adjfile.write('\t<luminaire id="%s-light" type="constant">\n' % lamp.data.name)
				mult = lamp.data.mitsuba_lamp.intensity
				adjfile.write('\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
						% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
				adjfile.write('\t\t<float name="samplingWeight" value="%f"/>\n' % lamp.data.mitsuba_lamp.sampling_weight)
				adjfile.write('\t</luminaire>\n')
			elif lamp.data.mitsuba_lamp.envmap_type == 'envmap':
				adjfile.write('\t<luminaire id="%s-light" type="envmap">\n' % lamp.data.name)
				adjfile.write('\t\t<string name="filename" value="%s"/>\n' % efutil.filesystem_path(lamp.data.mitsuba_lamp.envmap_file))
				self.export_worldtrafo(adjfile, lamp.matrix_world)
				adjfile.write('\t\t<float name="intensityScale" value="%f"/>\n' % lamp.data.mitsuba_lamp.intensity)
				adjfile.write('\t</luminaire>\n')

	def find_texture(self, name):
		if name in bpy.data.textures:
			return bpy.data.textures[name]
		else:
			raise Exception('Failed to find texture "%s"' % name)
	
	def find_material(self, name):
		if name in bpy.data.materials:
			return bpy.data.materials[name]
		else:
			raise Exception('Failed to find material "%s"' % name)

	def export_texture(self, adjfile, mat):
		if mat.name in self.exported_textures:
			return
		self.exported_textures += [mat.name]
		params = mat.mitsuba_texture.get_params()

		for p in params:
			if p.type == 'reference_texture':
				self.export_texture(adjfile, self.find_texture(p.value))

		adjfile.write('\t<texture id="%s" type="%s">\n' % (mat.name, mat.mitsuba_texture.type))
		adjfile.write(params.to_string())
		adjfile.write(params.to_string_ref())
		adjfile.write('\t</texture>\n')

	def export_material(self, adjfile, mat):
		if mat.name in self.exported_materials:
			return
		self.exported_materials += [mat.name]
		params = mat.mitsuba_material.get_params()

		for p in params:
			if p.type == 'reference_material':
				self.export_material(adjfile, self.find_material(p.value))
			elif p.type == 'reference_texture':
				self.export_texture(adjfile, self.find_texture(p.value))

		adjfile.write('\t<bsdf id="%s" type="%s">\n' % (mat.name, mat.mitsuba_material.type))
		adjfile.write(params.to_string())
		adjfile.write(params.to_string_ref())
		adjfile.write('\t</bsdf>\n')

	def export(self, scene):
		adjfile = open(self.target_file, 'w')
		adjfile.write('<adjustments>\n');
		idx = 0
		for obj in scene.objects:
			if obj.type == 'LAMP':
				self.export_lamp(adjfile, obj, idx)
			elif obj.type == 'MESH':
				for mat in obj.data.materials:
					self.export_material(adjfile, mat)
			idx = idx+1
		adjfile.write('</adjustments>\n');
		adjfile.close()


	
