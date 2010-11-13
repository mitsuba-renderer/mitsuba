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

from properties_texture import TextureButtonsPanel

from extensions_framework.ui import property_group_renderer

class mitsuba_texture_base(TextureButtonsPanel, property_group_renderer):
	'''
	This is the base class for all Mitsuba texture sub-panels.
	'''
	
	COMPAT_ENGINES	= {'mitsuba'}
	MTS_COMPAT		= set()

	@classmethod
	def poll(cls, context):
		'''
		Only show panel if mitsuba_texture.type in MTS_COMPAT
		'''
		tex = context.texture
		return	tex and \
				(context.scene.render.engine in cls.COMPAT_ENGINES) and \
				context.texture.mitsuba_texture.type in cls.MTS_COMPAT
