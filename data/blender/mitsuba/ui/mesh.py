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

class world_panel(bl_ui.properties_data_mesh.MeshButtonsPanel, property_group_renderer):
	COMPAT_ENGINES = { MitsubaAddon.BL_IDNAME }

@MitsubaAddon.addon_register_class
class meshes(world_panel):
	'''
	Mesh Settings
	'''
	
	bl_label = 'Mitsuba Mesh Options'
	
	display_property_groups = [
		( ('mesh',), 'mitsuba_mesh' )
	]

	def draw(self, context):
		super().draw(context)
