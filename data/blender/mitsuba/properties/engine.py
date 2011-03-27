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
class mitsuba_engine(declarative_property_group):
	ef_attach_to = ['Scene']

	controls = [
		'export_mode',
		'render_mode',
		'binary_path',
		'refresh_interval'
	]

	visibility = {
		'render_mode':		{ 'export_mode': 'render' },
		'refresh_interval':	{ 'export_mode': 'render', 'render_mode' : 'cli' }
	}

	properties = [
		{
			'type': 'enum',
			'attr': 'export_mode',
			'name': 'Export mode',
			'description': 'Specifies whether or not to launch the renderer after exporting the scene',
			'default': 'render',
			'items': [
				('render', 'Export + Render', 'render'),
				('exportonly', 'Only export', 'exportonly')
			],
			'save_in_preset': True
		},
		{
			'type': 'enum',
			'attr': 'render_mode',
			'name': 'Rendering mode',
			'description': 'Launch the external GUI or use the command-line renderer?',
			'default': 'cli',
			'items': [
				('cli', 'Mitsuba CLI', 'cli'),
				('gui', 'Mitsuba GUI', 'gui')
			],
			'save_in_preset': True
		},
		{
			'type': 'string',
			'subtype': 'DIR_PATH',
			'attr': 'binary_path',
			'name': 'Executable path',
			'description': 'Path to the Mitsuba install',
			'default': efutil.find_config_value('mitsuba', 'defaults', 'binary_path', '')
		},
		{
			'type': 'int',
			'attr': 'refresh_interval',
			'name': 'Refresh interval',
			'description': 'Period for updating rendering on screen (in seconds)',
			'default': 5,
			'min': 1,
			'soft_min': 1,
			'save_in_preset': True
		},
		{
			'type': 'int',
			'attr': 'preview_depth',
			'name': 'Depth',
			'description': 'Max. path depth used when generating the preview (2: direct illumination, 3: one-bounce indirect, etc.)',
			'default': int(efutil.find_config_value('mitsuba', 'defaults', 'preview_depth', '2')),
			'min': 2,
			'max': 10,
			'save_in_preset': True
		},
		{
			'type': 'int',
			'attr': 'preview_spp',
			'name': 'SPP',
			'description': 'Samples per pixel used to generate the preview',
			'default': int(efutil.find_config_value('mitsuba', 'defaults', 'preview_spp', '16')),
			'min': 1,
			'max': 128,
			'save_in_preset': True
		},
	]

