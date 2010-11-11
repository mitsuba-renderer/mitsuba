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
import os, copy, subprocess

# Blender Libs
import bpy
from presets import AddPresetBase

# Extensions_Framework Libs
from extensions_framework import util as efutil

from mitsuba.outputs import MtsLog
from mitsuba.export import MtsAdjustments

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

class MitsubaCheckOp(bpy.types.Operator):
	bl_idname = 'mts.check'
	bl_label = 'Check scene'

	def reportWarning(self, msg):
		self.report({'WARNING'}, msg)
		print("MtsBlend: %s" % msg)

	def _check_lamp(self, lamp):
		hasErrors = False
		if lamp.type == 'POINT' and lamp.falloff_type != 'INVERSE_SQUARE':
			self.reportWarning('Point light "%s" needs to have inverse square falloff' % lamp.name)
			hasErrors = True

		if hasErrors:
			self.reportWarning('Encountered one or more problems -- check the console')
		else:
			self.report({'INFO'}, "No problems found")

	def execute(self, context):
		scene = bpy.data.scenes[0]
		for obj in scene.objects:
			if obj.type == 'LAMP':
				self._check_lamp(obj.data)
		return {'FINISHED'}


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
		if self.properties.scene == '':
			scene = context.scene
		else:
			scene = bpy.data.scenes[self.properties.scene]

		if scene is None:
			self.report({'ERROR'}, 'Scene is not valid for export to %s' % self.properties.filename)
			return {'CANCELLED'}
		
		# Force scene update; NB, scene.update() doesn't work
		scene.frame_set(scene.frame_current)

		(self.properties.filename, _) = os.path.splitext(self.properties.filename)

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
		adj = MtsAdjustments(mts_adj_file)
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

menu_func = lambda self, context: self.layout.operator("export.mitsuba", text="Export Mitsuba scene...")
bpy.types.INFO_MT_file_export.append(menu_func)
