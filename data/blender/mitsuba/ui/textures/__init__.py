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

from extensions_framework.ui import property_group_renderer

from ... import MitsubaAddon

@MitsubaAddon.addon_register_class
class TEXTURE_PT_context_texture_mts(bl_ui.properties_texture.TextureButtonsPanel, bpy.types.Panel):
	bl_label = ""
	bl_options = {'HIDE_HEADER'}
	COMPAT_ENGINES	= { MitsubaAddon.BL_IDNAME }

	@classmethod
	def poll(cls, context):
		engine = context.scene.render.engine
		if not hasattr(context, "texture_slot"):
			return False
		return ((context.material or context.world or context.lamp or context.brush or context.texture)
			and (engine in cls.COMPAT_ENGINES))

	def draw(self, context):
		layout = self.layout
		slot = context.texture_slot
		node = context.texture_node
		space = context.space_data
		tex = context.texture
		idblock = bl_ui.properties_texture.context_tex_datablock(context)
		tex_collection = space.pin_id is None and type(idblock) != bpy.types.Brush and not node

		if tex_collection:
			row = layout.row()

			row.template_list(idblock, "texture_slots", idblock, "active_texture_index", rows=2)

			col = row.column(align=True)
			col.operator("texture.slot_move", text="", icon='TRIA_UP').type = 'UP'
			col.operator("texture.slot_move", text="", icon='TRIA_DOWN').type = 'DOWN'
			col.menu("TEXTURE_MT_specials", icon='DOWNARROW_HLT', text="")

		split = layout.split(percentage=1)
		col = split.column()

		if tex_collection:
			col.template_ID(idblock, "active_texture", new="texture.new")
		elif node:
			col.template_ID(node, "texture", new="texture.new")
		elif idblock:
			col.template_ID(idblock, "texture", new="texture.new")

		if space.pin_id:
			col.template_ID(space, "pin_id")

		col = split.column()

class mitsuba_texture_base(bl_ui.properties_texture.TextureButtonsPanel, property_group_renderer):
	'''
	This is the base class for all Mitsuba texture sub-panels.
	'''
	
	COMPAT_ENGINES	= { MitsubaAddon.BL_IDNAME }
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

