import math
from copy import deepcopy

import bpy

from extensions_framework import declarative_property_group
from extensions_framework import util as efutil
from mitsuba.properties.texture import TextureParameter
from mitsuba.export import ParamSet

param_reflectance = TextureParameter('reflectance', 'Reflectance', \
		'Diffuse reflectance value', default=(0.5, 0.5, 0.5))
param_diffuseReflectance = TextureParameter('diffuseReflectance', 'Diffuse reflectance', \
		'Diffuse reflectance value', default=(0.5, 0.5, 0.5))
param_specularReflectance = TextureParameter('specularReflectance', 'Specular reflectance', \
		'Specular reflectance value', default=(1.0, 1.0, 1.0))

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
				('lambertian', 'Lambertian', 'Lambertian (i.e. ideally diffuse) material'),
				('phong', 'Phong', 'Modified Phong BRDF'),
				('ward', 'Anisotropic Ward', 'Anisotropic Ward BRDF'),
				('dielectric', 'Ideal dielectric', 'Ideal dielectric material (e.g. glass)'),
				('mirror', 'Ideal mirror', 'Ideal mirror material'),
				('roughglass', 'Rough glass', 'Rough dielectric material (e.g. sand-blasted glass)'),
				('roughmetal', 'Rough metal', 'Rough conductor (e.g. sand-blasted metal)'),
				('difftrans', 'Diffuse transmitter', 'Material with an ideally diffuse transmittance'),
				('microfacet', 'Microfacet', 'Microfacet material (like the rough glass material, but without transmittance)')
			],
			'save_in_preset': True
		}
	]

	def get_params(self):
		sub_type = getattr(self, 'mitsuba_mat_%s' % self.type)
		return sub_type.get_params()

class mitsuba_mat_lambertian(declarative_property_group):
	controls = param_reflectance.controls
	
	properties = param_reflectance.properties
	
	visibility = dict_merge(param_reflectance.visibility)

	def get_params(self):
		params = ParamSet()
		params.update(param_reflectance.get_params(self))
		return params

class mitsuba_mat_phong(declarative_property_group):
	controls = [
		'diffuseAmount',
		'specularAmount',
		'exponent'
	] + param_diffuseReflectance.controls \
	  + param_specularReflectance.controls

	properties = [
		{
			'attr': 'diffuseAmount',
			'type': 'float',
			'description' : 'Diffuse reflection lobe multiplier',
			'name' : 'Diffuse amount',
			'default' : 1.0,
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'specularAmount',
			'type': 'float',
			'description' : 'Specular reflection lobe multiplier',
			'name' : 'Specular amount',
			'default' : 1.0,
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'exponent',
			'type': 'float',
			'description' : 'Phong exponent',
			'name' : 'Exponent',
			'default' : 10.0,
			'min': 0.001,
			'max': 10000.0,
			'save_in_preset': True
		}
	] + param_diffuseReflectance.properties \
	  + param_specularReflectance.properties
	
	visibility = dict_merge(
		param_diffuseReflectance.visibility, 
		param_specularReflectance.visibility
	)

	def get_params(self):
		params = ParamSet()
		params.update(param_diffuseReflectance.get_params(self))
		params.update(param_specularReflectance.get_params(self))
		params.add_float('diffuseAmount', self.diffuseAmount)
		params.add_float('specularAmount', self.specularAmount)
		params.add_float('exponent', self.exponent)
		return params

class mitsuba_mat_ward(declarative_property_group):
	controls = [
		'diffuseAmount',
		'specularAmount',
		['alphaX', 'alphaY']
	] + param_diffuseReflectance.controls \
	  + param_specularReflectance.controls

	properties = [
		{
			'attr': 'diffuseAmount',
			'type': 'float',
			'description' : 'Diffuse reflection lobe multiplier',
			'name' : 'Diffuse amount',
			'default' : 1.0,
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'specularAmount',
			'type': 'float',
			'description' : 'Specular reflection lobe multiplier',
			'name' : 'Specular amount',
			'default' : 1.0,
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'alphaX',
			'type': 'float',
			'description' : 'Roughness value along U (0.3=coarse, 0.001=very fine)',
			'name' : 'U Roughness',
			'default' : 0.1,
			'min': 0.001,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'alphaY',
			'type': 'float',
			'description' : 'Roughness value along V (0.3=coarse, 0.001=very fine)',
			'name' : 'V Roughness',
			'default' : 0.1,
			'min': 0.001,
			'max': 1.0,
			'save_in_preset': True
		}
	] + param_diffuseReflectance.properties \
	  + param_specularReflectance.properties
	
	visibility = dict_merge(
		param_diffuseReflectance.visibility,
		param_specularReflectance.visibility
	)

	def get_params(self):
		params = ParamSet()
		params.update(param_diffuseReflectance.get_params(self))
		params.update(param_specularReflectance.get_params(self))
		params.add_float('diffuseAmount', self.diffuseAmount)
		params.add_float('specularAmount', self.specularAmount)
		params.add_float('alphaX', self.alphaX)
		params.add_float('alphaY', self.alphaY)
		return params

class mitsuba_mat_microfacet(declarative_property_group):
	controls = [
		'diffuseAmount',
		'specularAmount',
		'alphaB',
		['extIOR', 'intIOR']
	] + param_diffuseReflectance.controls \
	  + param_specularReflectance.controls

	properties = [
		{
			'attr': 'diffuseAmount',
			'type': 'float',
			'description' : 'Diffuse reflection lobe multiplier',
			'name' : 'Diffuse amount',
			'default' : 0.0,
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'specularAmount',
			'type': 'float',
			'description' : 'Specular reflection lobe multiplier',
			'name' : 'Specular amount',
			'default' : 1.0,
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'alphaB',
			'type': 'float',
			'name': 'Roughness',
			'description' : 'Roughness value (0.3=coarse, 0.001=very fine)',
			'default' : 0.1,
			'min': 0.001,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'extIOR',
			'type': 'float',
			'name' : 'Ext. IOR',
			'description' : 'Exterior index of refraction (e.g. air=1, glass=1.5 approximately)',
			'default' : 1,
			'min': 1.0,
			'max': 10.0,
			'save_in_preset': True
		},
		{
			'attr': 'intIOR',
			'type': 'float',
			'name' : 'Int. IOR',
			'description' : 'Interior index of refraction (e.g. air=1, glass=1.5 approximately)',
			'default' : 1.5,
			'min': 1.0,
			'max': 10.0,
			'save_in_preset': True
		}
	] + param_diffuseReflectance.properties \
	  + param_specularReflectance.properties
	
	visibility = dict_merge(
		param_diffuseReflectance.visibility,
		param_specularReflectance.visibility
	)

	def get_params(self):
		params = ParamSet()
		params.update(param_diffuseReflectance.get_params(self))
		params.update(param_specularReflectance.get_params(self))
		params.add_float('diffuseAmount', self.diffuseAmount)
		params.add_float('specularAmount', self.specularAmount)
		params.add_float('alphaB', self.alphaB)
		params.add_float('extIOR', self.extIOR)
		params.add_float('intIOR', self.intIOR)
		return params

class mitsuba_mat_roughglass(declarative_property_group):
	controls = [
		'specularReflectance',
		'specularTransmittance',
		'alphaB',
		['extIOR', 'intIOR']
	]

	properties = [
		{
			'attr': 'specularReflectance',
			'type': 'float_vector',
			'subtype': 'COLOR',
			'description' : 'Weight of the specular reflectance',
			'name' : 'Specular reflectance',
			'default' : (1.0, 1.0, 1.0),
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'specularTransmittance',
			'type': 'float_vector',
			'subtype': 'COLOR',
			'description' : 'Weight of the specular transmittance',
			'name' : 'Specular transmittance',
			'default' : (1.0, 1.0, 1.0),
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'alphaB',
			'type': 'float',
			'name' : 'Roughness',
			'description' : 'Roughness value (0.3=coarse, 0.001=very fine)',
			'default' : 0.1,
			'min': 0.001,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'extIOR',
			'type': 'float',
			'name' : 'Ext. IOR',
			'description' : 'Exterior index of refraction (e.g. air=1, glass=1.5 approximately)',
			'default' : 1,
			'min': 1.0,
			'max': 10.0,
			'save_in_preset': True
		},
		{
			'attr': 'intIOR',
			'type': 'float',
			'name' : 'Int. IOR',
			'description' : 'Interior index of refraction (e.g. air=1, glass=1.5 approximately)',
			'default' : 1.5,
			'min': 1.0,
			'max': 10.0,
			'save_in_preset': True
		}
	]

	def get_params(self):
		params = ParamSet()
		params.add_color('specularReflectance', self.specularReflectance)
		params.add_color('specularTransmittance', self.specularTransmittance)
		params.add_float('alphaB', self.alphaB)
		params.add_float('extIOR', self.extIOR)
		params.add_float('intIOR', self.intIOR)
		return params

class mitsuba_mat_roughmetal(declarative_property_group):
	controls = [
		'alphaB',
		'ior', 'k'
	] + param_specularReflectance.controls

	properties = [
		{
			'attr': 'alphaB',
			'type': 'float',
			'name' : 'Roughness',
			'description' : 'Roughness value (0.3=coarse, 0.001=very fine)',
			'default' : 0.1,
			'min': 0.001,
			'max': 1.0,
			'expand' : False,
			'save_in_preset': True
		},
		{
			'attr': 'ior',
			'type': 'float_vector',
			'name' : 'IOR',
			'description' : 'Per-channel index of refraction of the conductor',
			'default' : (0.370, 0.370, 0.370),
			'min': 1.0,
			'max': 10.0,
			'expand' : False,
			'save_in_preset': True
		},
		{
			'attr': 'k',
			'type': 'float_vector',
			'name' : 'Absorption coefficient',
			'description' : 'Per-channel absorption coefficient of the conductor',
			'default' : (2.820, 2.820, 2.820),
			'min': 1.0,
			'max': 10.0,
			'save_in_preset': True
		}
	] + param_specularReflectance.properties

	visibility = dict_merge(param_specularReflectance.visibility)

	def get_params(self):
		params = ParamSet()
		params.update(param_specularReflectance.get_params(self))
		params.add_float('alphaB', self.alphaB)
		params.add_color('ior', self.ior)
		params.add_color('k', self.k)
		return params

class mitsuba_mat_dielectric(declarative_property_group):
	controls = [
		'specularReflectance',
		'specularTransmittance',
		['extIOR', 'intIOR']
	]

	properties = [
		{
			'attr': 'specularReflectance',
			'type': 'float_vector',
			'subtype': 'COLOR',
			'description' : 'Weight of the specular reflectance',
			'name' : 'Specular reflectance',
			'default' : (1.0, 1.0, 1.0),
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'specularTransmittance',
			'type': 'float_vector',
			'subtype': 'COLOR',
			'description' : 'Weight of the specular transmittance',
			'name' : 'Specular transmittance',
			'default' : (1.0, 1.0, 1.0),
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'extIOR',
			'type': 'float',
			'name' : 'Ext. IOR',
			'description' : 'Exterior index of refraction (e.g. air=1, glass=1.5 approximately)',
			'default' : 1,
			'min': 1.0,
			'max': 10.0,
			'save_in_preset': True
		},
		{
			'attr': 'intIOR',
			'type': 'float',
			'name' : 'Int. IOR',
			'description' : 'Interior index of refraction (e.g. air=1, glass=1.5 approximately)',
			'default' : 1.5,
			'min': 1.0,
			'max': 10.0,
			'save_in_preset': True
		}
	]

	def get_params(self):
		params = ParamSet()
		params.add_color('specularReflectance', self.specularReflectance)
		params.add_color('specularTransmittance', self.specularTransmittance)
		params.add_float('extIOR', self.extIOR)
		params.add_float('intIOR', self.intIOR)
		return params

class mitsuba_mat_mirror(declarative_property_group):
	controls = [
		'specularReflectance',
	]

	properties = [
		{
			'attr': 'specularReflectance',
			'type': 'float_vector',
			'subtype': 'COLOR',
			'description' : 'Weight of the specular reflectance',
			'name' : 'Specular reflectance',
			'default' : (1.0, 1.0, 1.0),
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		}
	]

	def get_params(self):
		params = ParamSet()
		params.add_color('specularReflectance', self.specularReflectance)
		return params

class mitsuba_mat_difftrans(declarative_property_group):
	controls = [
		'transmittance',
	]

	properties = [
		{
			'attr': 'transmittance',
			'type': 'float_vector',
			'subtype': 'COLOR',
			'description' : 'Diffuse transmittance value',
			'name' : 'Diffuse transmittance',
			'default' : (0.5, 0.5, 0.5),
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		}
	]

	def get_params(self):
		params = ParamSet()
		params.add_color('transmittance', self.transmittance)
		return params
