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

class world_panel(bl_ui.properties_world.WorldButtonsPanel, property_group_renderer):
	COMPAT_ENGINES = { MitsubaAddon.BL_IDNAME }

@MitsubaAddon.addon_register_class
class media(world_panel):
	'''
	Participating Media Settings
	'''
	
	bl_label = 'Mitsuba Media'
	
	display_property_groups = [
		( ('scene',), 'mitsuba_media' )
	]

	def draw(self, context):
		super().draw(context)
		
		if context.world:
			row = self.layout.row(align=True)
			row.menu("MITSUBA_MT_presets_medium", text=bpy.types.MITSUBA_MT_presets_medium.bl_label)
			row.operator("mitsuba.preset_medium_add", text="", icon="ZOOMIN")
			row.operator("mitsuba.preset_medium_add", text="", icon="ZOOMOUT").remove_active = True

			if len(context.scene.mitsuba_media.media) > 0:
				current_vol_ind = context.scene.mitsuba_media.media_index
				current_vol = context.scene.mitsuba_media.media[current_vol_ind]
				# 'name' is not a member of current_vol.properties,
				# so we draw it explicitly
				self.layout.prop(
					current_vol, 'name'
				)
				# Here we draw the currently selected mitsuba_media_data property group
				for control in current_vol.controls:
					self.draw_column(
						control,
						self.layout,
						current_vol,
						context,
						property_group = current_vol
					)
		else:
			self.layout.label('No active World available!')
