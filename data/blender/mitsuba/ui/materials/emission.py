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
from ... import MitsubaAddon
from ...ui.materials import mitsuba_material_base

@MitsubaAddon.addon_register_class
class emission(mitsuba_material_base, bpy.types.Panel):
	'''
	Material Emission Settings
	'''
	
	bl_label = 'Mitsuba Material Emission'
	bl_options = {'DEFAULT_CLOSED'}
	MTS_PROPS = ['color', 'intensity', 'samplingWeight', 'use_emission']

	display_property_groups = [
		( ('material',), 'mitsuba_emission' )
	]

	def get_contents(self, mat):
		return mat.mitsuba_emission

	def draw_header(self, context):
		if hasattr(context, "material"):
			mat = bl_ui.properties_material.active_node_mat(context.material)
			self.layout.prop(mat.mitsuba_emission, "use_emission", text="")
			self.validate(context)
