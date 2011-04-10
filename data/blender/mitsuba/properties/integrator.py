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
class mitsuba_integrator(declarative_property_group):
	'''
	Storage class for Mitsuba Integrator settings.
	This class will be instantiated within a Blender scene
	object.
	'''

	ef_attach_to = ['Scene']

	controls = [
		'type',
		'maxdepth',
		['motionblur',
		'shuttertime']
	]
	
	visibility = {
		'shuttertime':		{ 'motionblur': True }
	}

	properties = [
		{
			'type': 'enum',
			'attr': 'type',
			'name': 'Type',
			'description': 'Specifies the type of integrator to use',
			'default': 'direct',
			'items': [
				('volpath', 'Volumetric path tracer', 'volpath'),
				('path', 'Path tracer', 'path'),
				('direct', 'Direct Illumination', 'direct'),
				('ptracer', 'Adjoint Particle Tracer', 'ptracer')
			],
			'save_in_preset': True
		},
		{
			'type': 'bool',
			'attr': 'motionblur',
			'name': 'Motion Blur',
			'description': 'Should motion blur be enabled?',
			'default' : False,
			'save_in_preset': True
		},
		{
			'type': 'float',
			'attr': 'shuttertime',
			'name': 'Shutter time',
			'description': 'Amount of time, for which the shutter remains open (measured in frames)',
			'save_in_preset': True,
			'min': 0,
			'max': 100,
			'default': 1
		},
		{
			'type': 'int',
			'attr': 'maxdepth',
			'name': 'Max. path depth',
			'description': 'Maximum path depth to be rendered. 2 corresponds to direct illumination, 3 is 1-bounce indirect illumination, etc.',
			'save_in_preset': True,
			'min': 2,
			'max': 100,
			'default': 4
		}
	]

