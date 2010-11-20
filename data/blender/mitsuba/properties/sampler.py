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

from extensions_framework import declarative_property_group
from extensions_framework import util as efutil

class mitsuba_sampler(declarative_property_group):
	'''
	Storage class for Mitsuba Sampler settings.
	This class will be instantiated within a Blender scene
	object.
	'''

	controls = [
		'type',
		'sampleCount'
	]

	properties = [
		{
			'type': 'enum',
			'attr': 'type',
			'name': 'Type',
			'description': 'Specifies the type of sampler to use',
			'default': 'ldsampler',
			'items': [
				('independent', 'Independent', 'independent'),
				('stratified', 'Stratified', 'stratified'),
				('ldsampler', 'Low discrepancy', 'ldsampler')
			],
			'save_in_preset': True
		},
		{
			'type': 'int',
			'attr': 'sampleCount',
			'name': 'Pixel samples',
			'description': 'Number of samples to use for estimating the illumination at each pixel',
			'default': 8,
			'min': 1,
			'max': 10240,
			'save_in_preset': True
		}
	]

