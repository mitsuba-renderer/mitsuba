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
from properties_material import MaterialButtonsPanel

from extensions_framework.ui import property_group_renderer

class mitsuba_material_base(MaterialButtonsPanel, property_group_renderer):
	COMPAT_ENGINES	= {'mitsuba'}

class mitsuba_material_sub(MaterialButtonsPanel, property_group_renderer):
	COMPAT_ENGINES	= {'mitsuba'}
	MTS_COMPAT		= set()

	@classmethod
	def poll(cls, context):
		'''
		Only show Mitsuba panel if mitsuba_material.material in MTS_COMPAT
		'''
		
		return super().poll(context) and context.material.mitsuba_material.type in cls.MTS_COMPAT
