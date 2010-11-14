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

class MATERIAL_PT_context_material_mts(MaterialButtonsPanel, bpy.types.Panel):
	bl_label = ""
	bl_options = {'HIDE_HEADER'}
	COMPAT_ENGINES = {'mitsuba'}

	@classmethod
	def poll(cls, context):
		# An exception, dont call the parent poll func because
		# this manages materials for all engine types

		engine = context.scene.render.engine
		return (context.material or context.object) and (engine in cls.COMPAT_ENGINES)

	def draw(self, context):
		layout = self.layout

		mat = context.material
		ob = context.object
		slot = context.material_slot
		space = context.space_data

		if ob:
			row = layout.row()

			row.template_list(ob, "material_slots", ob, "active_material_index", rows=2)

			col = row.column(align=True)
			col.operator("object.material_slot_add", icon='ZOOMIN', text="")
			col.operator("object.material_slot_remove", icon='ZOOMOUT', text="")

			col.menu("MATERIAL_MT_specials", icon='DOWNARROW_HLT', text="")

			if ob.mode == 'EDIT':
				row = layout.row(align=True)
				row.operator("object.material_slot_assign", text="Assign")
				row.operator("object.material_slot_select", text="Select")
				row.operator("object.material_slot_deselect", text="Deselect")

		split = layout.split(percentage=0.65)

		if ob:
			split.template_ID(ob, "active_material", new="material.new")
			row = split.row()

			if slot:
				row.prop(slot, "link", text="")
			else:
				row.label()
		elif mat:
			split.template_ID(space, "pin_id")
			split.separator()

