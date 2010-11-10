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

import bpy
import subprocess
import os
import sys
import math
import copy
import shutil
import time

# Basic Mitsuba integration based on the POV-Ray add-on
# Piggybacks on the COLLADA exporter to get most things done
class MitsubaRender(bpy.types.RenderEngine):
	bl_idname = 'MITSUBA_RENDER'
	bl_label = "Mitsuba"

	def _export(self, scene):
		import tempfile
		self._temp_dir = tempfile.mkdtemp(prefix='mitsuba_');
		self._temp_dae = os.path.join(self._temp_dir, 'scene.dae')
		self._temp_xml = os.path.join(self._temp_dir, 'scene.xml')
		self._temp_adj = os.path.join(self._temp_dir, 'scene_adjustments.xml')
		self._temp_out = os.path.join(self._temp_dir, 'scene.png')
		print("MtsBlend: Writing COLLADA file")
		while True:
			try:
				bpy.ops.wm.collada_export(filepath=self._temp_dae, check_existing=False)
				break
			except SystemError:
				# Weird SystemError (Operator bpy.ops.wm.collada_export.poll() 
				# failed, context is incorrect) -> try again
				print("MtsBlend: Retrying")
				time.sleep(0.1)
		print("MtsBlend: Writing adjustments file")
		adjfile = open(self._temp_adj, 'w')
		adjfile.write('<adjustments>\n');
		adjfile.write('</adjustments>\n');
		adjfile.close()

	def _render(self):
		mts_path = '/home/wenzel/mitsuba'
		mtsimport_binary = os.path.join(mts_path, "mtsimport")
		mitsuba_binary = os.path.join(mts_path, "mitsuba")
		mts_render_libpath = os.path.join(mts_path, "src/librender")
		mts_core_libpath = os.path.join(mts_path, "src/libcore")
		mts_hw_libpath = os.path.join(mts_path, "src/libhw")
		env = copy.copy(os.environ)
		env['LD_LIBRARY_PATH'] = mts_core_libpath + ":" + mts_render_libpath + ":" + mts_hw_libpath
		scene = bpy.data.scenes[0]
		render = scene.render
		width = int(render.resolution_x * render.resolution_percentage * 0.01)
		height = int(render.resolution_y * render.resolution_percentage * 0.01)

		try:
			print("MtsBlend: Launching mtsimport")
			process = subprocess.Popen(
				[mtsimport_binary, '-s', '-r', '%dx%d' % (width, height),
					'-l', 'pngfilm', 
					self._temp_dae, self._temp_xml,
					self._temp_adj],
				env = env,
				cwd = self._temp_dir
			)
			if process.wait() != 0:
				print("MtsBlend: mtsimport returned with a nonzero status")
				return False
			
			self._process = subprocess.Popen(
				[mitsuba_binary, self._temp_xml, '-o', self._temp_out],
				env = env,
				cwd = self._temp_dir
			)
		except OSError:
			print("MtsBlend: Could not execute '%s', possibly Mitsuba isn't installed" % mtsimport_binary)
			return False

		return True

	def _cleanup(self):
	 	#shutil.rmtree(self._temp_dir)
		print("Not cleaning up")

	def render(self, scene):
		self._export(scene)
		if not self._render():
			self.update_stats("", "MtsBlend: Unable to render (please check the console)")
			return

		r = scene.render
		x = int(r.resolution_x * r.resolution_percentage * 0.01)
		y = int(r.resolution_y * r.resolution_percentage * 0.01)
		DELAY = 0.02

		# Wait for the file to be created
		while not os.path.exists(self._temp_out):
			if self.test_break():
				try:
					self._process.terminate()
				except:
					pass
				break

			if self._process.poll() != None:
				self.update_stats("", "MtsBlend: Failed to render (check console)")
				break

			time.sleep(DELAY)

		if os.path.exists(self._temp_out):
			self.update_stats("", "MtsBlend: Rendering")
			prev_size = -1

			def update_image():
				result = self.begin_result(0, 0, x, y)
				layer = result.layers[0]
				try:
					print("Loading %s" % self._temp_out)
					layer.load_from_file(self._temp_out)
				except:
					pass
				self.end_result(result)

			while True:
				if self._process.poll() is not None:
					update_image()
					break

				if self.test_break():
					try:
						self._process.terminate()
					except:
						pass

					break

				new_size = os.path.getsize(self._temp_out)
				if new_size != prev_size:
					update_image()
					prev_size = new_size

				time.sleep(DELAY)

		self._cleanup()
