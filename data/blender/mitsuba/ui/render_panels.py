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

import os, bpy, bl_ui

from .. import MitsubaAddon

from extensions_framework.ui import property_group_renderer
from extensions_framework import util as efutil

class render_panel(bl_ui.properties_render.RenderButtonsPanel, property_group_renderer):
	'''
	Base class for render engine settings panels
	'''
	
	COMPAT_ENGINES = { MitsubaAddon.BL_IDNAME }

@MitsubaAddon.addon_register_class
class setup_preset(render_panel, bpy.types.Panel):
	'''
	Engine settings presets UI Panel
	'''
	
	bl_label = 'Mitsuba Engine Presets'
	
	def draw(self, context):
		row = self.layout.row(align=True)
		row.menu("MITSUBA_MT_presets_engine", text=bpy.types.MITSUBA_MT_presets_engine.bl_label)
		row.operator("mitsuba.preset_engine_add", text="", icon="ZOOMIN")
		row.operator("mitsuba.preset_engine_add", text="", icon="ZOOMOUT").remove_active = True

		super().draw(context)

@MitsubaAddon.addon_register_class
class engine(render_panel, bpy.types.Panel):
	'''
	Engine settings UI Panel
	'''
	
	bl_label = 'Mitsuba Engine Settings'
	
	display_property_groups = [
		( ('scene',), 'mitsuba_engine' )
	]

@MitsubaAddon.addon_register_class
class integrator(render_panel, bpy.types.Panel):
	'''
	Integrator settings UI Panel
	'''
	
	bl_label = 'Mitsuba Integrator Settings'
	
	display_property_groups = [
		( ('scene',), 'mitsuba_integrator' )
	]

@MitsubaAddon.addon_register_class
class sampler(render_panel, bpy.types.Panel):
	'''
	Sampler settings UI Panel
	'''

	bl_label = 'Mitsuba Sampler Settings'
	
	display_property_groups = [
		( ('scene',), 'mitsuba_sampler' )
	]

