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
from ... import MitsubaAddon
from ...ui.materials import mitsuba_material_base

@MitsubaAddon.addon_register_class
class main(mitsuba_material_base, bpy.types.Panel):
	'''
	Material Editor UI Panel
	'''

	bl_label	= 'Mitsuba Materials'

	display_property_groups = [
		( ('material',), 'mitsuba_material' )
	]
	
	def draw(self, context):
		row = self.layout.row(align=True)
		row.menu("MITSUBA_MT_presets_material", text=bpy.types.MITSUBA_MT_presets_material.bl_label)
		row.operator("mitsuba.preset_material_add", text="", icon="ZOOMIN")
		row.operator("mitsuba.preset_material_add", text="", icon="ZOOMOUT").remove_active = True
	
		row = self.layout.row(align=True)
		row.operator("mitsuba.convert_all_materials", icon='WORLD_DATA')
		row = self.layout.row(align=True)
		row.operator("mitsuba.convert_material", icon='MATERIAL_DATA')
		row = self.layout.row(align=True)

		row.menu('MATERIAL_MT_mitsuba_type', text=context.material.mitsuba_material.type_label)
		super().draw(context)
