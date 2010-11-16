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
import os, sys, copy, subprocess, traceback, string

# Blender Libs
import bpy
from presets import AddPresetBase

# Extensions_Framework Libs
from extensions_framework import util as efutil

from mitsuba.outputs import MtsLog
from mitsuba.export.adjustments import MtsAdjustments

def try_preset_path_create(preset_subdir):
	target_path = os.path.join(bpy.utils.preset_paths('')[0], preset_subdir)
	if not os.path.exists(target_path):
		os.makedirs(target_path)

class MITSUBA_MT_base(object):
	preset_operator = "script.execute_preset"
	def draw(self, context):
		try_preset_path_create(self.preset_subdir)
		return bpy.types.Menu.draw_preset(self, context)

class MITSUBA_OT_preset_base(AddPresetBase):
	def execute(self, context):
		try_preset_path_create(self.preset_subdir)
		return super().execute(context)

class MITSUBA_MT_presets_engine(MITSUBA_MT_base, bpy.types.Menu):
	bl_label = "Mitsuba Engine Presets"
	preset_subdir = "mitsuba/engine"

class MITSUBA_OT_preset_engine_add(MITSUBA_OT_preset_base, bpy.types.Operator):
	'''Save the current settings as a preset'''
	bl_idname = 'mitsuba.preset_engine_add'
	bl_label = 'Add Mitsuba Engine settings preset'
	preset_menu = 'MITSUBA_MT_presets_engine'
	preset_values = [
		'bpy.context.scene.mitsuba_engine.%s'%v['attr'] for v in bpy.types.mitsuba_engine.get_exportable_properties()
	]
	preset_subdir = 'mitsuba/engine'


class MITSUBA_MT_presets_texture(MITSUBA_MT_base, bpy.types.Menu):
	bl_label = "Mitsuba Texture Presets"
	preset_subdir = "mitsuba/texture"

class MITSUBA_OT_preset_texture_add(MITSUBA_OT_preset_base, bpy.types.Operator):
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


class MITSUBA_MT_presets_material(MITSUBA_MT_base, bpy.types.Menu):
	bl_label = "Mitsuba Material Presets"
	preset_subdir = "mitsuba/material"

class MITSUBA_OT_preset_material_add(MITSUBA_OT_preset_base, bpy.types.Operator):
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

class EXPORT_OT_mitsuba(bpy.types.Operator):
	bl_idname = 'export.mitsuba'
	bl_label = 'Export Mitsuba Scene (.xml)'

	filename		= bpy.props.StringProperty(name='Target filename')
	directory		= bpy.props.StringProperty(name='Target directory')
	scene			= bpy.props.StringProperty(options={'HIDDEN'}, default='')

	def invoke(self, context, event):
		context.window_manager.add_fileselect(self)
		return {'RUNNING_MODAL'}

	def execute(self, context):
		try:
			if self.properties.scene == '':
				scene = context.scene
			else:
				scene = bpy.data.scenes[self.properties.scene]

			if scene is None:
				self.report({'ERROR'}, 'Scene is not valid for export to %s' % self.properties.filename)
				return {'CANCELLED'}
			
			# Force scene update; NB, scene.update() doesn't work
			scene.frame_set(scene.frame_current)

			mts_basename = os.path.join(
				self.properties.directory,
				self.properties.filename)
			mts_dae_file = mts_basename + ".dae"
			mts_xml_file = mts_basename + ".xml"
			mts_adj_file = mts_basename + "_adjustments.xml"
			mts_meshes_dir = os.path.join(self.properties.directory, "meshes")

			efutil.export_path = mts_xml_file
			try:
				os.mkdir(mts_meshes_dir)
			except OSError:
				pass

			scene.collada_export(mts_dae_file)

			MtsLog('MtsBlend: Writing adjustments file to "%s"' % mts_adj_file)
			adj = MtsAdjustments(mts_adj_file, self.properties.directory)
			adj.export(scene)

			if scene.mitsuba_engine.binary_path == "":
				self.report({'ERROR'}, 'Mitsuba binary path must be specified!')
				return {'CANCELLED'}

			scene.mitsuba_engine.binary_path = efutil.filesystem_path(scene.mitsuba_engine.binary_path)
			efutil.write_config_value('mitsuba', 'defaults', 'binary_path', scene.mitsuba_engine.binary_path)

			(mts_path, tail) = os.path.split(bpy.path.abspath(scene.mitsuba_engine.binary_path))
			mtsimport_binary = os.path.join(mts_path, "mtsimport")
			env = copy.copy(os.environ)
			mts_render_libpath = os.path.join(mts_path, "src/librender")
			mts_core_libpath = os.path.join(mts_path, "src/libcore")
			mts_hw_libpath = os.path.join(mts_path, "src/libhw")
			env['LD_LIBRARY_PATH'] = mts_core_libpath + ":" + mts_render_libpath + ":" + mts_hw_libpath
			render = scene.render
			width = int(render.resolution_x * render.resolution_percentage * 0.01)
			height = int(render.resolution_y * render.resolution_percentage * 0.01)

			MtsLog("MtsBlend: Launching mtsimport")
			try:
				process = subprocess.Popen(
					[mtsimport_binary, '-r', '%dx%d' % (width, height),
						'-l', 'pngfilm', mts_dae_file, mts_xml_file, mts_adj_file],
					env = env,
					cwd = self.properties.directory
				)
				if process.wait() != 0:
					self.report({'ERROR'}, "mtsimport returned with a nonzero status!")
					return {'CANCELLED'}
			except OSError:
				self.report({'ERROR'}, "Could not execute '%s'" % mtsimport_binary)
				return {'CANCELLED'}

			return {'FINISHED'}
		except:
			typ, value, tb = sys.exc_info()
			elist = traceback.format_exception(typ, value, tb)
			MtsLog("Caught exception: %s" % ''.join(elist))
			return {'CANCELLED'}

menu_func = lambda self, context: self.layout.operator("export.mitsuba", text="Export Mitsuba scene...")
bpy.types.INFO_MT_file_export.append(menu_func)

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
