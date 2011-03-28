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

# System Libs
import os, sys, subprocess, traceback, string, math

# Blender Libs
import bpy, bl_operators

# Extensions_Framework Libs
from extensions_framework import util as efutil

from .. import MitsubaAddon
from ..outputs import MtsLog
from ..export import MtsExporter

class MITSUBA_MT_base(bpy.types.Menu):
	preset_operator = "script.execute_preset"
	def draw(self, context):
		return self.draw_preset(context)

@MitsubaAddon.addon_register_class
class MITSUBA_MT_presets_engine(MITSUBA_MT_base):
	bl_label = "Mitsuba Engine Presets"
	preset_subdir = "mitsuba/engine"

@MitsubaAddon.addon_register_class
class MITSUBA_OT_preset_engine_add(bl_operators.presets.AddPresetBase, bpy.types.Operator):
	'''Save the current settings as a preset'''
	bl_idname = 'mitsuba.preset_engine_add'
	bl_label = 'Add Mitsuba Engine settings preset'
	preset_menu = 'MITSUBA_MT_presets_engine'
	preset_subdir = 'mitsuba/engine'
	
	def execute(self, context):
		self.preset_values = [
			'bpy.context.scene.mitsuba_engine.%s'%v['attr'] for v in bpy.types.mitsuba_engine.get_exportable_properties()
		]
		return super().execute(context)

@MitsubaAddon.addon_register_class
class MITSUBA_MT_presets_texture(MITSUBA_MT_base):
	bl_label = "Mitsuba Texture Presets"
	preset_subdir = "mitsuba/texture"

@MitsubaAddon.addon_register_class
class MITSUBA_OT_preset_texture_add(bl_operators.presets.AddPresetBase, bpy.types.Operator):
	'''Save the current settings as a preset'''
	bl_idname = 'mitsuba.preset_texture_add'
	bl_label = 'Add Mitsuba Texture settings preset'
	preset_menu = 'MITSUBA_MT_presets_texture'
	preset_values =  []
	preset_subdir = 'mitsuba/texture'
	
	def execute(self, context):
		pv = [
			'bpy.context.texture.mitsuba_texture.%s'%v['attr'] for v in bpy.types.mitsuba_texture.get_exportable_properties()
		]
		mts_type = context.texture.mitsuba_texture.type
		sub_type = getattr(bpy.types, 'mitsuba_tex_%s' % mts_type)

		pv.extend([
			'bpy.context.texture.mitsuba_texture.mitsuba_tex_%s.%s'%(mts_type, v['attr']) for v in sub_type.get_exportable_properties()
		])
		pv.extend([
			'bpy.context.texture.mitsuba_texture.mitsuba_tex_mapping.%s'%v['attr'] for v in bpy.types.mitsuba_tex_mapping.get_exportable_properties()
		])

		self.preset_values = pv
		return super().execute(context)

@MitsubaAddon.addon_register_class
class MITSUBA_MT_presets_material(MITSUBA_MT_base):
	bl_label = "Mitsuba Material Presets"
	preset_subdir = "mitsuba/material"

@MitsubaAddon.addon_register_class
class MITSUBA_OT_preset_material_add(bl_operators.presets.AddPresetBase, bpy.types.Operator):
	'''Save the current settings as a preset'''
	bl_idname = 'mitsuba.preset_material_add'
	bl_label = 'Add Mitsuba Material settings preset'
	preset_menu = 'MITSUBA_MT_presets_material'
	preset_values =  []
	preset_subdir = 'mitsuba/material'
	
	def execute(self, context):
		pv = [
			'bpy.context.material.mitsuba_material.%s'%v['attr'] for v in bpy.types.mitsuba_material.get_exportable_properties()
		] 

		# store only the sub-properties of the selected mitsuba material type
		mts_type = context.material.mitsuba_material.type
		sub_type = getattr(bpy.types, 'mitsuba_mat_%s' % mts_type)
		
		pv.extend([
			'bpy.context.material.mitsuba_material.mitsuba_mat_%s.%s'%(mts_type, v['attr']) for v in sub_type.get_exportable_properties()
		])
		pv.extend([
			'bpy.context.material.mitsuba_material.mitsuba_emission.%s'%v['attr'] for v in bpy.types.mitsuba_emission.get_exportable_properties()
		])

		self.preset_values = pv
		return super().execute(context)

@MitsubaAddon.addon_register_class
class EXPORT_OT_mitsuba(bpy.types.Operator):
	bl_idname = 'export.mitsuba'
	bl_label = 'Export Mitsuba Scene (.xml)'

	filename		= bpy.props.StringProperty(name='Target filename', subtype = 'FILE_PATH')
	directory		= bpy.props.StringProperty(name='Target directory')
	scene			= bpy.props.StringProperty(options={'HIDDEN'}, default='')

	def invoke(self, context, event):
		context.window_manager.fileselect_add(self)
		return {'RUNNING_MODAL'}

	def execute(self, context):
		try:
			if self.properties.scene == '':
				scene = context.scene
			else:
				scene = bpy.data.scenes[self.properties.scene]

			result = MtsExporter(
				directory = self.properties.directory,
				filename = self.properties.filename).export(scene)

			if not result:
				self.report({'ERROR'}, "Unsucessful export!");
				return {'CANCELLED'}

			return {'FINISHED'}
		except:
			typ, value, tb = sys.exc_info()
			elist = traceback.format_exception(typ, value, tb)
			MtsLog("Caught exception: %s" % ''.join(elist))
			self.report({'ERROR'}, "Unsucessful export!");
			return {'CANCELLED'}

def menu_func(self, context):
	default_path = os.path.splitext(os.path.basename(bpy.data.filepath))[0] + ".xml"
	self.layout.operator("export.mitsuba", text="Export Mitsuba scene...").filename = default_path
bpy.types.INFO_MT_file_export.append(menu_func)

@MitsubaAddon.addon_register_class
class MITSUBA_OT_material_slot_move(bpy.types.Operator):
	''' Rearrange the material slots '''
	bl_idname = 'mitsuba.material_slot_move'
	bl_label = 'Move a material entry up or down'
	type = bpy.props.StringProperty(name='type')

	def execute(self, context):
		obj = bpy.context.active_object
		index = obj.active_material_index
		new_index = index-1 if self.properties.type == 'UP' else index+1
		size = len(obj.material_slots)

		if new_index >= 0 and new_index < size:
			obj.active_material_index = 0
			# Can't write to material_slots, hence the kludge
			materials = []
			for i in range(0, size):
				materials += [obj.material_slots[i].name]
			for i in range(0, size):
				mat = obj.data.materials.pop(0)
				del(mat)
			temp = materials[index]
			materials[index] = materials[new_index]
			materials[new_index] = temp
			for i in range(0, size):
				obj.data.materials.append(bpy.data.materials[materials[i]])

			obj.active_material_index = new_index
		return {'FINISHED'}

@MitsubaAddon.addon_register_class
class MITSUBA_OT_material_add(bpy.types.Operator):
	''' Append a new material '''
	bl_idname = 'mitsuba.material_add'
	bl_label = 'Append a new material'
	type = bpy.props.StringProperty(name='type')

	def execute(self, context):
		obj = bpy.context.active_object
		index = obj.active_material_index
		if len(obj.material_slots) == 0:
			curName = 'Material'
		else:
			curName = obj.material_slots[index].name
		mat = bpy.data.materials.new(name=curName)
		obj.data.materials.append(mat)
		obj.active_material_index = len(obj.data.materials)-1
		return {'FINISHED'}

def material_converter(report, scene, blender_mat):
	try:
		mitsuba_mat = blender_mat.mitsuba_material

		mitsuba_mat.type = 'microfacet'
		mitsuba_mat.mitsuba_mat_microfacet.diffuseReflectance_color =  [blender_mat.diffuse_intensity*i for i in blender_mat.diffuse_color]
		mitsuba_mat.mitsuba_mat_microfacet.specularAmount = 1.0
		mitsuba_mat.mitsuba_mat_microfacet.diffuseAmount = 1.0
					
		logHardness = math.log(blender_mat.specular_hardness)
		specular_scale = 2.0 * max(0.0128415*logHardness**2 - 0.171266*logHardness + 0.575631, 0.0)
		mitsuba_mat.mitsuba_mat_microfacet.specularReflectance_color =  [specular_scale * blender_mat.specular_intensity*i for i in blender_mat.specular_color]
		mitsuba_mat.mitsuba_mat_microfacet.alphaB = min(max(0.757198 - 0.120395*logHardness, 0.0), 1.0)

		report({'INFO'}, 'Converted blender material "%s"' % blender_mat.name)
		return {'FINISHED'}
	except Exception as err:
		report({'ERROR'}, 'Cannot convert material: %s' % err)
		return {'CANCELLED'}

@MitsubaAddon.addon_register_class
class MITSUBA_OT_convert_all_materials(bpy.types.Operator):
	bl_idname = 'mitsuba.convert_all_materials'
	bl_label = 'Convert all Blender materials'
	
	def report_log(self, level, msg):
		MtsLog('Material conversion %s: %s' % (level, msg))
	
	def execute(self, context):
		for blender_mat in bpy.data.materials:
			# Don't convert materials from linked-in files
			if blender_mat.library == None:
				material_converter(self.report_log, context.scene, blender_mat)
		return {'FINISHED'}

@MitsubaAddon.addon_register_class
class MITSUBA_OT_convert_material(bpy.types.Operator):
	bl_idname = 'mitsuba.convert_material'
	bl_label = 'Convert selected Blender material'
	
	material_name = bpy.props.StringProperty(default='')
	
	def execute(self, context):
		if self.properties.material_name == '':
			blender_mat = context.material
		else:
			blender_mat = bpy.data.materials[self.properties.material_name]
		
		material_converter(self.report, context.scene, blender_mat)
		return {'FINISHED'}	
