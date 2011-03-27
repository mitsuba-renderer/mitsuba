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
from ...outputs import MtsLog

from extensions_framework.ui import property_group_renderer
from extensions_framework import util as efutil

material_cache = {}
cached_spp = None
cached_depth = None

def copy(value):
	if value == None or isinstance(value, str) or isinstance(value, bool) \
		or isinstance(value, float) or isinstance(value, int):
		return value
	elif getattr(value, '__len__', False):
		return list(value)
	else:
		raise Exception("Copy: don't know how to handle '%s'" % str(vlaue))

class mitsuba_material_base(bl_ui.properties_material.MaterialButtonsPanel, property_group_renderer):
	COMPAT_ENGINES	= { MitsubaAddon.BL_IDNAME }
	MTS_PROPS	   = ['type']

	def validate(self, context):
		global material_cache
		mat = context.material
		if mat.name in material_cache:
			mat_cached = material_cache[mat.name]
		else:
			mat_cached = {}
			material_cache[mat.name] = mat_cached

		mat = self.get_contents(mat)
		repaint = False
		for prop in self.MTS_PROPS:
			prop_value = copy(getattr(mat, prop))
			prop_cache_value = mat_cached[prop] if prop in mat_cached else None
			if prop_cache_value != prop_value:
				mat_cached[prop] = prop_value
				repaint = True
		if repaint:
			# Cause a repaint
			MtsLog("Forcing a repaint")
			#context.material.preview_render_type = context.material.preview_render_type 

	def draw(self, context):
		if not hasattr(context, 'material'):
			return
		self.validate(context)
		return super().draw(context)

	def get_contents(self, mat):
		return mat.mitsuba_material

class mitsuba_material_sub(bl_ui.properties_material.MaterialButtonsPanel, property_group_renderer):
	COMPAT_ENGINES	= { MitsubaAddon.BL_IDNAME }
	MTS_COMPAT		= set()
	MTS_PROPS	   = []

	@classmethod
	def poll(cls, context):
		'''
		Only show Mitsuba panel if mitsuba_material.material in MTS_COMPAT
		'''
		if not hasattr(context, 'material'):
			return False

		return super().poll(context) and context.material.mitsuba_material.type in cls.MTS_COMPAT

	def draw(self, context):
		if not hasattr(context, 'material'):
			return
		mat = context.material
		sub_type = getattr(bpy.types, 'mitsuba_mat_%s' % mat.mitsuba_material.type)
		global material_cache
		if mat.name in material_cache:
			mat_cached = material_cache[mat.name]
		else:
			mat_cached = {}
			material_cache[mat.name] = mat_cached
	
		mat = getattr(mat.mitsuba_material, 'mitsuba_mat_%s' %
				mat.mitsuba_material.type)
		props = sub_type.get_exportable_properties()

		repaint = False
		for prop_entry in props:
			prop = prop_entry['attr']
			prop_value = copy(getattr(mat, prop) if hasattr(mat, prop) else None)
			prop_cache_value = mat_cached[prop] if prop in mat_cached else None
			if prop_cache_value != copy(prop_value):
				mat_cached[prop] = prop_value
				repaint = True
		if repaint:
			# Cause a repaint
			MtsLog("Forcing a repaint")
			#context.material.preview_render_type = context.material.preview_render_type 
		return super().draw(context)

@MitsubaAddon.addon_register_class
class MATERIAL_PT_preview_mts(bl_ui.properties_material.MaterialButtonsPanel, bpy.types.Panel):
	bl_label = "Preview"
	COMPAT_ENGINES	= { MitsubaAddon.BL_IDNAME }

	def draw(self, context):
		if not hasattr(context, 'material'):
			return
		self.layout.template_preview(context.material, show_buttons=False)
		engine = context.scene.mitsuba_engine
		row = self.layout.row(True)
		row.prop(engine, "preview_depth")
		row.prop(engine, "preview_spp")

		global cached_depth
		global cached_spp
		if engine.preview_depth != cached_depth or engine.preview_spp != cached_spp:
			actualChange = cached_depth != None
			cached_depth = engine.preview_depth
			cached_spp = engine.preview_spp
			if actualChange:
				MtsLog("Forcing a repaint")
				#context.material.preview_render_type = context.material.preview_render_type 
				efutil.write_config_value('mitsuba', 'defaults', 'preview_spp', str(cached_spp))
				efutil.write_config_value('mitsuba', 'defaults', 'preview_depth', str(cached_depth))

@MitsubaAddon.addon_register_class
class MATERIAL_PT_context_material_mts(bl_ui.properties_material.MaterialButtonsPanel, bpy.types.Panel):
	bl_label = ""
	bl_options = {'HIDE_HEADER'}
	COMPAT_ENGINES	= { MitsubaAddon.BL_IDNAME }

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

			row.template_list(ob, "material_slots", ob, "active_material_index", rows=4)

			col = row.column(align=True)
			col.operator("mitsuba.material_add", icon='ZOOMIN', text="")
			col.operator("object.material_slot_remove", icon='ZOOMOUT', text="")
			col.operator("mitsuba.material_slot_move", text="", icon='TRIA_UP').type = 'UP'
			col.operator("mitsuba.material_slot_move", text="", icon='TRIA_DOWN').type = 'DOWN'

			col.menu("MATERIAL_MT_specials", icon='DOWNARROW_HLT', text="")

			if ob.mode == 'EDIT':
				row = layout.row(align=True)
				row.operator("object.material_slot_assign", text="Assign")
				row.operator("object.material_slot_select", text="Select")
				row.operator("object.material_slot_deselect", text="Deselect")

		split = layout.split(percentage=0.75)

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

