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

from mitsuba.ui.materials import mitsuba_material_sub

class ui_material_composite(mitsuba_material_sub, bpy.types.Panel):
	bl_label = 'Mitsuba Composite Material'

	MTS_COMPAT = {'composite'}
	
	display_property_groups = [
		( ('material', 'mitsuba_material'), 'mitsuba_mat_composite' )
	]
		
	def draw(self, context):
		super().draw(context)

		mat = context.material.mitsuba_material.mitsuba_mat_composite
		weight = 0
		for i in range(1,mat.nElements+1):
			weight += getattr(mat, "mat%i_weight" % i)
		if weight > 1:
			row = self.layout.row()
			row.label("Warning: material weights sum to >1")
