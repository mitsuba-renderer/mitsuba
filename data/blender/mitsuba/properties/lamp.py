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

@MitsubaAddon.addon_register_class
class mitsuba_lamp(declarative_property_group):
	ef_attach_to = ['Lamp']

	controls = [
		'samplingWeight',
		'envmap_type',
		'envmap_file'
	]

	visibility = {
		'envmap_type': { 'type': 'ENV' },
		'envmap_file': { 'type': 'ENV', 'envmap_type' : 'envmap' }
	}

	properties = [
		{
			'type': 'float',
			'attr': 'samplingWeight',
			'name': 'Sampling weight',
			'description': 'Relative amount of samples to place on this light source (e.g. the "importance")',
			'default': 1.0,
			'min': 1e-3,
			'soft_min': 1e-3,
			'max': 1e3,
			'soft_max': 1e3,
			'save_in_preset': True
		},
		{
			'type': 'float',
			'attr': 'intensity',
			'name': 'Intensity',
			'description': 'Specifies the intensity of the light source',
			'default': 10.0,
			'min': 1e-3,
			'soft_min': 1e-3,
			'max': 1e5,
			'soft_max': 1e5,
			'save_in_preset': True
		},
		{
			'type': 'enum',
			'attr': 'envmap_type',
			'name': 'Environment map type',
			'description': 'Environment map type',
			'default': 'constant',
			'items': [
				('constant', 'Constant background source', 'constant'),
				('envmap', 'HDRI environment map', 'envmap')
			],
			'save_in_preset': True
		},
		{
			'type': 'string',
			'subtype': 'FILE_PATH',
			'attr': 'envmap_file',
			'name': 'HDRI Map',
			'description': 'EXR image to use for lighting (in latitude-longitude format)',
			'default': '',
			'save_in_preset': True
		}
	]

