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

# From collada_internal.cpp

translate_start_name_map = list(map(chr, [
	95,  95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	65,  66,  67,  68,  69,  70,  71,  72,
	73,  74,  75,  76,  77,  78,  79,  80,
	81,  82,  83,  84,  85,  86,  87,  88,
	89,  90,  95,  95,  95,  95,  95,  95,
	97,  98,  99,  100,  101,  102,  103,  104,
	105,  106,  107,  108,  109,  110,  111,  112,
	113,  114,  115,  116,  117,  118,  119,  120,
	121,  122,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  192,
	193,  194,  195,  196,  197,  198,  199,  200,
	201,  202,  203,  204,  205,  206,  207,  208,
	209,  210,  211,  212,  213,  214,  95,  216,
	217,  218,  219,  220,  221,  222,  223,  224,
	225,  226,  227,  228,  229,  230,  231,  232,
	233,  234,  235,  236,  237,  238,  239,  240,
	241,  242,  243,  244,  245,  246,  95,  248,
	249,  250,  251,  252,  253,  254,  255]))

translate_name_map = list(map(chr, [
	95,  95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  45,  95,  95,  48,
	49,  50,  51,  52,  53,  54,  55,  56,
	57,  95,  95,  95,  95,  95,  95,  95,
	65,  66,  67,  68,  69,  70,  71,  72,
	73,  74,  75,  76,  77,  78,  79,  80,
	81,  82,  83,  84,  85,  86,  87,  88,
	89,  90,  95,  95,  95,  95,  95,  95,
	97,  98,  99,  100,  101,  102,  103,  104,
	105,  106,  107,  108,  109,  110,  111,  112,
	113,  114,  115,  116,  117,  118,  119,  120,
	121,  122,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  95,  95,
	95,  95,  95,  95,  95,  95,  183,  95,
	95,  95,  95,  95,  95,  95,  95,  192,
	193,  194,  195,  196,  197,  198,  199,  200,
	201,  202,  203,  204,  205,  206,  207,  208,
	209,  210,  211,  212,  213,  214,  95,  216,
	217,  218,  219,  220,  221,  222,  223,  224,
	225,  226,  227,  228,  229,  230,  231,  232,
	233,  234,  235,  236,  237,  238,  239,  240,
	241,  242,  243,  244,  245,  246,  95,  248,
	249,  250,  251,  252,  253,  254,  255]))

def translate_id(name):
	# Doesn't handle duplicates at the moment
	result = ""
	if len(name) == 0:
		return name
	result += translate_start_name_map[ord(name[0])]
	for i in range(1, len(name)):
		result += translate_name_map[ord(name[i])]
	return result

class ParamSetItem(list):
	type		= None
	type_name	= None
	name		= None
	value		= None

	def __init__(self, *args):
		self.type, self.name, self.value = args
		self.type_name = "%s %s" % (self.type, self.name)
		self.append(self.type_name)
		self.append(self.value)
	
	def to_string(self):
		if self.type == "color":
			return '\t\t<rgb name="%s" value="%s %s %s"/>\n' % (self.name,
					self.value[0], self.value[1], self.value[2])
		elif self.type == "point" or self.type == "vector":
			return '\t\t<%s name="%s" value="%s %s %s"/>\n' % (self.type,
					self.name, self.value[0], self.value[1], self.value[2])
		elif self.type == "integer" or self.type == "float" \
				or self.type ==	"string":
			return '\t\t<%s name="%s" value="%s"/>\n' % (self.type, self.name, self.value)
		else:
			return ""
	
	def to_string_ref(self):
		if self.type == "reference_texture" or self.type == "reference_material":
			if self.name == "":
				return '\t\t<ref id="%s"/>\n' % translate_id(self.value)
			else:
				return '\t\t<ref name="%s" id="%s"/>\n' % (self.name, translate_id(self.value))
		else:
			return ""

class ParamSet(list):
	names = []
	
	def update(self, other):
		for p in other:
			self.add(p.type, p.name, p.value)
		return self
	
	def add(self, type, name, value):
		if name in self.names:
			for p in self:
				if p.name == name:
					self.remove(p)
		
		self.append(
			ParamSetItem(type, name, value)
		)
		self.names.append(name)
		return self
	
	def add_float(self, name, value):
		self.add('float', name, value)
		return self
	
	def add_integer(self, name, value):
		self.add('integer', name, value)
		return self

	def add_reference(self, type, name, value):
		self.add('reference_%s' % type, name, value)
		return self

	def add_bool(self, name, value):
		self.add('bool', name, bool(value))
		return self

	def add_string(self, name, value):
		self.add('string', name, str(value))
		return self
	
	def add_vector(self, name, value):
		self.add('vector', name, [i for i in value])
		return self
	
	def add_point(self, name, value):
		self.add('point', name, [p for p in value])
		return self
	
	def add_color(self, name, value):
		self.add('color', name, [c for c in value])
		return self
	
	def to_string(self):
		return ''.join(item.to_string() for item in self)

	def to_string_ref(self):
		return ''.join(item.to_string_ref() for item in self)

def get_instance_materials(ob):
	obmats = []
	# Grab materials attached to object instances ...
	if hasattr(ob, 'material_slots'):
		for ms in ob.material_slots:
			obmats.append(ms.material)
	# ... and to the object's mesh data
	if hasattr(ob.data, 'materials'):
		for m in ob.data.materials:
			obmats.append(m)
	return obmats
