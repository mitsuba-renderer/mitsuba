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

from .. import MitsubaAddon

from extensions_framework import declarative_property_group
from extensions_framework import util as efutil

@MitsubaAddon.addon_register_class
class mitsuba_mesh(declarative_property_group):
	ef_attach_to = ['Mesh', 'SurfaceCurve', 'TextCurve', 'Curve']

	controls = [
		'normals'
	]

	visibility = {
	}

	properties = [
		{
			'type': 'enum',
			'attr': 'normals',
			'name': 'Normal mode',
			'description': 'Specifies how Mitsuba obtains normal information',
			'items' : [
				('dihedralangle','Smooth vertex normals (using angle constraint)', 'dihedralangle'),
				('vertexnormals','Smooth vertex normals (using connectivity)', 'vertexnormals'),
				('facenormals','Flat face normals', 'facenormals'),
				('default','Use normals from Blender', 'default')
			],
			'default': 'default'
		}
	]

