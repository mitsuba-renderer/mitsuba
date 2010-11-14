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
from extensions_framework.validate import Logic_OR as O

from mitsuba.export import ParamSet

#------------------------------------------------------------------------------ 
# Texture property group construction helpers
#------------------------------------------------------------------------------ 

class TextureParameterBase(object):
	attr				= None
	name				= None
	default				= (0.8, 0.8, 0.8)
	min					= 0.0
	max					= 1.0
	
	texture_collection	= 'texture_slots'
	
	controls			= None
	visibility			= None
	properties			= None
	
	def __init__(self, attr, name, description, default=None, min=None, max=None):
		self.attr = attr
		self.name = name
		self.description = description
		if default is not None:
			self.default = default
		if min is not None:
			self.min = min
		if max is not None:
			self.max = max

		self.controls = self.get_controls()
		self.visibility = self.get_visibility()
		self.properties = self.get_properties()
	
	def texture_collection_finder(self):
		return lambda s,c: s.object.material_slots[s.object.active_material_index].material

	def texture_slot_set_attr(self):
		def set_attr(s,c):
			if type(c).__name__ == 'mitsuba_material':
				return getattr(c, 'mitsuba_mat_%s'%c.type)
			else:
				return getattr(c, 'mitsuba_tex_%s'%c.type)
		return set_attr
	
	def get_controls(self):
		'''
		Subclasses can override this for their own needs
		'''	
		return []
	
	def get_visibility(self):
		'''
		Subclasses can override this for their own needs
		'''	
		return {}
	
	def get_properties(self):
		'''
		Subclasses can override this for their own needs
		'''	
		return []
	
	def get_extra_controls(self):
		'''
		Subclasses can override this for their own needs
		'''	
		return []
	
	def get_extra_visibility(self):
		'''
		Subclasses can override this for their own needs
		'''	
		return {}
	
	def get_extra_properties(self):
		'''
		Subclasses can override this for their own needs
		'''	
		return []
	
	def get_params(self, context):
		'''
		Return a Mitsuba ParamSet of the properties
		defined in this Texture, getting parameters
		from property group 'context'
		'''
		return ParamSet()

class TextureParameter(TextureParameterBase):
	def get_controls(self):
		return [
			[ 0.9, [0.375,'%s_colorlabel' % self.attr, '%s_color' % self.attr], '%s_usetexture' % self.attr ],
			'%s_texture' % self.attr
		] + self.get_extra_controls()
	
	def get_visibility(self):
		vis = {
			'%s_texture' % self.attr: { '%s_usetexture' % self.attr: True },
		}
		vis.update(self.get_extra_visibility())
		return vis
	
	# colour for each material type. If the property name is
	# not set, then the color won't be changed.
	master_color_map = {
		'lambertian': 'reflectance'
	}

	def set_master_colour(self, s, c):
		'''
		This neat little hack will set the blender material colour to the value
		given in the material panel via the property's 'draw' lambda function.
		'''

		if c.type in self.master_color_map.keys() and self.attr == self.master_color_map[c.type]:
			submat = getattr(c, 'mitsuba_mat_%s'%c.type)
			submat_col = getattr(submat, self.attr+'_color')
			if s.material.diffuse_color != submat_col:
				s.material.diffuse_color = submat_col
	
	def get_properties(self):
		return [
			{
				'attr': self.attr,
				'type': 'string',
				'default': 'mts_color_texture'
			},
			{
				'attr': '%s_usetexture' % self.attr,
				'type': 'bool',
				'name': 'T',
				'description': 'Textured %s' % self.name,
				'default': False,
				'toggle': True,
				'save_in_preset': True
			},
			{
				'type': 'text',
				'attr': '%s_colorlabel' % self.attr,
				'name': self.name
			},
			{
				'type': 'float_vector',
				'attr': '%s_color' % self.attr,
				'name': '', #self.name,
				'description': self.description,
				'default': self.default,
				'min': self.min,
				'soft_min': self.min,
				'max': self.max,
				'soft_max': self.max,
				'subtype': 'COLOR',
				'draw': lambda s,c: self.set_master_colour(s, c),
				'save_in_preset': True
			},
			{
				'attr': '%s_texturename' % self.attr,
				'type': 'string',
				'name': '%s_texturename' % self.attr,
				'description': '%s Texture' % self.name,
				'save_in_preset': True
			},
			{
				'type': 'prop_search',
				'attr': '%s_texture' % self.attr,
				'src': self.texture_collection_finder(),
				'src_attr': self.texture_collection,
				'trg': self.texture_slot_set_attr(),
				'trg_attr': '%s_texturename' % self.attr,
				'name': self.name
			}
		] + self.get_extra_properties()

	def get_params(self, context):
		params = ParamSet()
		if hasattr(context, '%s_usetexture' % self.attr) \
			and getattr(context, '%s_usetexture' % self.attr):
			params.add_reference('texture', self.attr, getattr(context, '%s_texturename' % self.attr))
		else:
			params.add_color(
				self.attr,
				getattr(context, '%s_color' % self.attr)
			)
		return params

class mitsuba_texture(declarative_property_group):
	'''
	Storage class for Mitsuba Texture settings.
	This class will be instantiated within a Blender Texture
	object.
	'''

	controls = [
		'type'
	]

	properties = [
		{
			'attr': 'type',
			'name': 'Texture type',
			'type': 'enum',
			'items': [
				('ldrtexture', 'Bitmap', 'ldrtexture'),
				('checkerboard', 'Checkerboard', 'checkerboard')
			],
			'default' : 'ldrtexture',
			'save_in_preset': True
		},
	]

	def get_params(self):
		if hasattr(self, 'mitsuba_tex_%s' % self.type):
			mts_texture = getattr(self, 'mitsuba_tex_%s' % self.type) 
			params = mts_texture.get_params()
			params.update(self.mitsuba_tex_mapping.get_params())
			return params
		else:
			return ParamSet()

class mitsuba_tex_mapping(declarative_property_group):
	controls = [
		['uscale', 'vscale'],
		['uoffset', 'voffset']
	]

	properties = [
		{
			'attr': 'uscale',
			'type': 'float',
			'name': 'U Scale',
			'default': 1.0,
			'min': -100.0,
			'soft_min': -100.0,
			'max': 100.0,
			'soft_max': 100.0,
			'save_in_preset': True
		},
		{
			'attr': 'vscale',
			'type': 'float',
			'name': 'V Scale',
			'default': 1.0,
			'min': -100.0,
			'soft_min': -100.0,
			'max': 100.0,
			'soft_max': 100.0,
			'save_in_preset': True
		},
		{
			'attr': 'uoffset',
			'type': 'float',
			'name': 'U Offset',
			'default': 0.0,
			'min': -100.0,
			'soft_min': -100.0,
			'max': 100.0,
			'soft_max': 100.0,
			'save_in_preset': True
		},
		{
			'attr': 'voffset',
			'type': 'float',
			'name': 'V Offset',
			'default': 0.0,
			'min': -100.0,
			'soft_min': -100.0,
			'max': 100.0,
			'soft_max': 100.0,
			'save_in_preset': True
		}
	]

	def get_params(self):
		mapping_params = ParamSet()
		mapping_params.add_float('uscale', self.uscale)
		mapping_params.add_float('vscale', self.vscale)
		mapping_params.add_float('uoffset', self.uoffset)
		mapping_params.add_float('voffset', self.voffset)
		return mapping_params

class mitsuba_tex_ldrtexture(declarative_property_group):
	controls = [
		'filename',
		'wrapMode',
		'filterType', 
		'maxAnisotropy',
		['srgb', 'gamma']
	]
	
	visibility = {
		'maxAnisotropy': { 'filterType': 'ewa' },
		'gamma': { 'srgb': False }
	}

	properties = [
		{
			'type': 'string',
			'subtype': 'FILE_PATH',
			'attr': 'filename',
			'description' : 'Path to a PNG/JPG/TGA/BMP file',
			'name': 'File Name',
			'save_in_preset': True
		},
		{
			'type': 'enum',
			'attr': 'filterType',
			'name': 'Filter type',
			'description' : 'Specifies the type of texture filtering',
			'items': [
				('isotropic', 'Isotropic (Trilinear MipMap)', 'isotropic'),
				('ewa', 'Anisotropic (EWA Filtering)', 'ewa')
			],
			'default' : 'ewa',
			'save_in_preset': True
		},
		{
			'type': 'bool',
			'attr': 'srgb',
			'name': 'sRGB',
			'description' : 'Is the texture stored in sRGB color space?',
			'default': True,
			'save_in_preset': True
		},
		{
			'type': 'float',
			'attr': 'gamma',
			'name': 'Gamma',
			'description' : 'Specifies the texture gamma value',
			'default': 2.2,
			'min': 0.01,
			'soft_min': 0.01,
			'max': 6.0,
			'soft_max': 6.0,
			'save_in_preset': True
		},
		{
			'type': 'float',
			'description' :  'Maximum allowed anisotropy when using the EWA filter',
			'attr': 'maxAnisotropy',
			'name': 'Max. Anisotropy',
			'default': 8.0,
			'save_in_preset': True
		},
		{
			'type': 'enum',
			'attr': 'wrapMode',
			'name': 'Wrapping',
			'description' : 'What should be done when encountering UV coordinates outside of the range [0,1]',
			'items': [
				('repeat', 'Repeat', 'repeat'),
				('black', 'Black outside of [0, 1]', 'black'),
				('white', 'White outside of [0, 1]', 'white'),
				('clamp', 'Clamp to edges', 'clamp')
			],
			'save_in_preset': True
		}
	]

	def get_params(self):
		params = ParamSet()

		params.add_string('filename', efutil.path_relative_to_export(self.filename) ) \
			  .add_string('filterType', self.filterType) \
			  .add_float('maxAnisotropy', self.maxAnisotropy) \
			  .add_string('wrapMode', self.wrapMode) \
			  .add_float('gamma', -1 if self.srgb else self.gamma)

		return params

class mitsuba_tex_checkerboard(declarative_property_group):
	controls = [
		'darkColor',
		'brightColor'
	]

	properties = [
		{
			'attr': 'darkColor',
			'type': 'float_vector',
			'subtype': 'COLOR',
			'name' : 'Dark color',
			'description' : 'Color of the dark patches',
			'default' : (0.2, 0.2, 0.2),
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		},
		{
			'attr': 'brightColor',
			'type': 'float_vector',
			'subtype': 'COLOR',
			'description' : 'Color of the bright patches',
			'name' : 'Bright color',
			'default' : (0.4, 0.4, 0.4),
			'min': 0.0,
			'max': 1.0,
			'save_in_preset': True
		}
	]

	def get_params(self):
		params = ParamSet()
		
		params.add_color('darkColor', self.darkColor) 
		params.add_color('brightColor', self.brightColor) 

		return params

