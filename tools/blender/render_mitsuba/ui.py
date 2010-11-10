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

import bpy

# Use some of the existing buttons.
import properties_render
properties_render.RENDER_PT_render.COMPAT_ENGINES.add('MITSUBA_RENDER')
properties_render.RENDER_PT_dimensions.COMPAT_ENGINES.add('MITSUBA_RENDER')
properties_render.RENDER_PT_antialiasing.COMPAT_ENGINES.add('MITSUBA_RENDER')
properties_render.RENDER_PT_output.COMPAT_ENGINES.add('MITSUBA_RENDER')
del properties_render

# Use only a subset of the world panels
import properties_world
properties_world.WORLD_PT_preview.COMPAT_ENGINES.add('MITSUBA_RENDER')
properties_world.WORLD_PT_context_world.COMPAT_ENGINES.add('MITSUBA_RENDER')
properties_world.WORLD_PT_world.COMPAT_ENGINES.add('MITSUBA_RENDER')
properties_world.WORLD_PT_mist.COMPAT_ENGINES.add('MITSUBA_RENDER')
del properties_world

# Example of wrapping every class 'as is'
import properties_material
for member in dir(properties_material):
	subclass = getattr(properties_material, member)
	try:
		subclass.COMPAT_ENGINES.add('MITSUBA_RENDER')
	except:
		pass
del properties_material

import properties_data_mesh
for member in dir(properties_data_mesh):
	subclass = getattr(properties_data_mesh, member)
	try:
		subclass.COMPAT_ENGINES.add('MITSUBA_RENDER')
	except:
		pass
del properties_data_mesh

import properties_texture
for member in dir(properties_texture):
	subclass = getattr(properties_texture, member)
	try:
		subclass.COMPAT_ENGINES.add('MITSUBA_RENDER')
	except:
		pass
del properties_texture

import properties_data_camera
for member in dir(properties_data_camera):
	subclass = getattr(properties_data_camera, member)
	try:
		subclass.COMPAT_ENGINES.add('MITSUBA_RENDER')
	except:
		pass
del properties_data_camera

import properties_data_lamp
for member in dir(properties_data_lamp):
	subclass = getattr(properties_data_lamp, member)
	try:
		subclass.COMPAT_ENGINES.add('MITSUBA_RENDER')
	except:
		pass
del properties_data_lamp

class RenderButtonsPanel():
	bl_space_type = 'PROPERTIES'
	bl_region_type = 'WINDOW'
	bl_context = "render"

	@classmethod
	def poll(cls, context):
		rd = context.scene.render
		return (rd.use_game_engine == False) and (rd.engine in cls.COMPAT_ENGINES)


class RENDER_PT_mitsuba_radiosity(RenderButtonsPanel, bpy.types.Panel):
	bl_label = "Mitsuba setup"
	COMPAT_ENGINES = {'MITSUBA_RENDER'}

	def draw(self, context):
		layout = self.layout
		scene = context.scene
		rd = scene.render

		split = layout.split()
		col = split.column();
		col.prop(scene, "mts_path", text="Executable")
		row = col.row();
		row.prop(scene, "mts_gui", text="In external GUI")
		row.operator("wm.save_homefile", text="Save", icon ='FILE_TICK')
		col.operator("mts.check", text="Check scene")
