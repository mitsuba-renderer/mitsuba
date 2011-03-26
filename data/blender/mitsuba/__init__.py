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

import os

bl_info = {
	"name": "Mitsuba",
	"author": "Wenzel Jakob",
	"version": (0, 2, 1),
	"blender": (2, 5, 6),
	"api": 35669,
	"category": "Render",
	"location": "Render > Engine > Mitsuba",
	"description": "Basic Mitsuba integration for Blender",
	"warning": "",
	"wiki_url": "http://wiki.blender.org/index.php/Extensions:2.5/Py/"\
		"Scripts/Render/Mitsuba",
	"tracker_url": "https://www.mitsuba-renderer.org/bugtracker/projects/mitsuba"}

def plugin_path():
	return os.path.dirname(os.path.realpath(__file__))

if 'core' in locals():
	import imp
	imp.reload(core)
else:
	import bpy
	
	from extensions_framework import Addon
	MitsubaAddon = Addon(bl_info)
	register, unregister = MitsubaAddon.init_functions()

	# Importing the core package causes extensions_framework managed
	# RNA class registration via @MitsubaAddon.addon_register_class
	from . import core
