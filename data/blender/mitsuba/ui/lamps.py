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

import bpy, bl_ui

from .. import MitsubaAddon

from extensions_framework.ui import property_group_renderer

narrowui = 180

@MitsubaAddon.addon_register_class
class lamps(bl_ui.properties_data_lamp.DataButtonsPanel, property_group_renderer, bpy.types.Panel):
	bl_label = 'Mitsuba Lamps'
	COMPAT_ENGINES = {'mitsuba'}
	
	display_property_groups = [
		( ('lamp',), 'mitsuba_lamp' )
	]
	
	# Overridden to draw some of blender's lamp controls
	def draw(self, context):
		lamp = context.lamp
		if lamp is not None:
			layout = self.layout
			wide_ui = context.region.width > narrowui

			if wide_ui:
				layout.prop(lamp.mitsuba_lamp, "type", expand=True)
			else:
				layout.prop(lamp.mitsuba_lamp, "type", text="")
			
			if lamp.mitsuba_lamp.type == 'ENV':
				lamp.type = 'HEMI'
			else:
				lamp.type = lamp.mitsuba_lamp.type

			split = layout.split()
			
			col = split.column()

			if lamp.type == 'HEMI':
				layout.prop(lamp.mitsuba_lamp, "envmap_type", text="Type")

			if not(lamp.type == 'HEMI' and
					lamp.mitsuba_lamp.envmap_type == 'envmap'):
				layout.prop(lamp, "color", text="Color")
			else:
				layout.prop(lamp.mitsuba_lamp, "envmap_file", text="HDRI file")
			layout.prop(lamp.mitsuba_lamp, "intensity", text="Intensity")
			layout.prop(lamp.mitsuba_lamp, "samplingWeight", text = "Sampling weight")

			# SPOT LAMP: Blender Properties
			if lamp.type == 'SPOT':
				wide_ui = context.region.width > narrowui
				
				if wide_ui:
					#col = split.column()
					col=layout.row()
				else:
					col=layout.column()
				col.prop(lamp, "spot_size", text="Size")
				col.prop(lamp, "spot_blend", text="Blend", slider=True)
				col=layout.row()
				col.prop(lamp, "show_cone")

			# AREA LAMP: Blender Properties
			elif lamp.type == 'AREA':
				if wide_ui:
					col=layout.row()
				else:
					col=layout.column()
				col.row().prop(lamp, "shape", expand=True)
				sub = col.column(align=True)

				if (lamp.shape == 'SQUARE'):
					sub.prop(lamp, "size")
				elif (lamp.shape == 'RECTANGLE'):
					sub.prop(lamp, "size", text="Size X")
					sub.prop(lamp, "size_y", text="Size Y")
			elif wide_ui:
				col = split.column()
			
			if lamp.type == 'HEMI':
				layout.label('Note: these cover the whole sphere')
