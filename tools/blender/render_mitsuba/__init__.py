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

bl_addon_info = {
	"name": "Mitsuba",
	"author": "Wenzel Jakob",
	"version": (0, 1),
	"blender": (2, 5, 5),
	"api": 31667,
	"location": "Render > Engine > Mitsuba",
	"description": "Basic Mitsuba integration for Blender",
	"warning": "",
	"wiki_url": "http://wiki.blender.org/index.php/Extensions:2.5/Py/"\
		"Scripts/Render/Mitsuba",
	"tracker_url": "Unavailable",
	"category": "Render"}


if "bpy" in locals():
	reload(render)
	reload(ui)

else:
	import bpy
	from bpy.props import *
	from render_mitsuba import render
	from render_mitsuba import ui


def register():
	Scene = bpy.types.Scene

	# Not a real pov option, just to know if we should write
	Scene.mts_radio_enable = BoolProperty(
			name="Enable Radiosity",
			description="Enable mitsubas radiosity calculation",
			default=False)
	Scene.mts_radio_display_advanced = BoolProperty(
			name="Advanced Options",
			description="Show advanced options",
			default=False)

	# Real pov options
	Scene.mts_radio_adc_bailout = FloatProperty(
			name="ADC Bailout", description="The adc_bailout for radiosity rays. Use adc_bailout = 0.01 / brightest_ambient_object for good results",
			min=0.0, max=1000.0, soft_min=0.0, soft_max=1.0, default=0.01)

	Scene.mts_radio_always_sample = BoolProperty(
			name="Always Sample", description="Only use the data from the pretrace step and not gather any new samples during the final radiosity pass",
			default=True)

	Scene.mts_radio_brightness = FloatProperty(
			name="Brightness", description="Amount objects are brightened before being returned upwards to the rest of the system",
			min=0.0, max=1000.0, soft_min=0.0, soft_max=10.0, default=1.0)

	Scene.mts_radio_count = IntProperty(
			name="Ray Count", description="Number of rays that are sent out whenever a new radiosity value has to be calculated",
			min=1, max=1600, default=35)

	Scene.mts_radio_error_bound = FloatProperty(
			name="Error Bound", description="One of the two main speed/quality tuning values, lower values are more accurate",
			min=0.0, max=1000.0, soft_min=0.1, soft_max=10.0, default=1.8)

	Scene.mts_radio_gray_threshold = FloatProperty(
			name="Gray Threshold", description="One of the two main speed/quality tuning values, lower values are more accurate",
			min=0.0, max=1.0, soft_min=0, soft_max=1, default=0.0)

	Scene.mts_radio_low_error_factor = FloatProperty(
			name="Low Error Factor", description="If you calculate just enough samples, but no more, you will get an image which has slightly blotchy lighting",
			min=0.0, max=1.0, soft_min=0.0, soft_max=1.0, default=0.5)

	# max_sample - not available yet
	Scene.mts_radio_media = BoolProperty(
			name="Media", description="Radiosity estimation can be affected by media",
			default=False)

	Scene.mts_radio_minimum_reuse = FloatProperty(
			name="Minimum Reuse", description="Fraction of the screen width which sets the minimum radius of reuse for each sample point (At values higher than 2% expect errors)",
			min=0.0, max=1.0, soft_min=0.1, soft_max=0.1, default=0.015)

	Scene.mts_radio_nearest_count = IntProperty(
			name="Nearest Count", description="Number of old ambient values blended together to create a new interpolated value",
			min=1, max=20, default=5)

	Scene.mts_radio_normal = BoolProperty(
			name="Normals", description="Radiosity estimation can be affected by normals",
			default=False)

	Scene.mts_radio_recursion_limit = IntProperty(
			name="Recursion Limit", description="how many recursion levels are used to calculate the diffuse inter-reflection",
			min=1, max=20, default=3)


def unregister():
	import bpy
	Scene = bpy.types.Scene

	del Scene.mts_radio_enable
	del Scene.mts_radio_display_advanced
	del Scene.mts_radio_adc_bailout
	del Scene.mts_radio_always_sample
	del Scene.mts_radio_brightness
	del Scene.mts_radio_count
	del Scene.mts_radio_error_bound
	del Scene.mts_radio_gray_threshold
	del Scene.mts_radio_low_error_factor
	del Scene.mts_radio_media
	del Scene.mts_radio_minimum_reuse
	del Scene.mts_radio_nearest_count
	del Scene.mts_radio_normal
	del Scene.mts_radio_recursion_limit

if __name__ == "__main__":
	register()
