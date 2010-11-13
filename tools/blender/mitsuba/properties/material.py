import math
from copy import deepcopy

import bpy

from extensions_framework import declarative_property_group
from extensions_framework import util as efutil
from mitsuba.properties.texture import TextureParameter
from mitsuba.export import ParamSet

param_reflectance = TextureParameter('reflectance', 'Reflectance', \
		'Diffuse reflectance value', default=(0.5, 0.5, 0.5))

def dict_merge(*args):
	vis = {}
	for vis_dict in args:
		vis.update(deepcopy(vis_dict))
	return vis


class mitsuba_material(declarative_property_group):
	'''
	Storage class for Mitsuba Material settings.
	This class will be instantiated within a Blender Material
	object.
	'''
	
	controls = [
		'type',
	] 

	properties = [
		# Material Type Select
		{
			'type': 'enum',
			'attr': 'type',
			'name': 'Type',
			'description': 'Mitsuba material type',
			'default': 'matte',
			'items': [
				('lambertian', 'Lambertian', 'Lambertian (i.e. ideally diffuse) material')
			],
			'save_in_preset': True
		}
	]

	def get_params(self):
		sub_type = getattr(self, 'mitsuba_mat_%s' % self.type)
		return sub_type.get_params()

class mitsuba_mat_lambertian(declarative_property_group):
	controls = [
	] + param_reflectance.controls
	
	properties = [
	] + param_reflectance.properties
	
	visibility = dict_merge(
		param_reflectance.visibility
	)

	def get_params(self):
		params = ParamSet()
		params.update(param_reflectance.get_params(self))
		return params

