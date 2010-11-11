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

class MitsubaCheckOp(bpy.types.Operator):
	bl_idname = 'mts.check'
	bl_label = 'Check scene'

	def reportWarning(self, msg):
		self.report({'WARNING'}, msg)
		print("MtsBlend: %s" % msg)

	def _check_lamp(self, lamp):
		hasErrors = False
		if lamp.type == 'POINT' and lamp.falloff_type != 'INVERSE_SQUARE':
			self.reportWarning('Point light "%s" needs to have inverse square falloff' % lamp.name)
			hasErrors = True

		if hasErrors:
			self.reportWarning('Encountered one or more problems -- check the console')
		else:
			self.report({'INFO'}, "No problems found")

	def execute(self, context):
		scene = bpy.data.scenes[0]
		for obj in scene.objects:
			if obj.type == 'LAMP':
				self._check_lamp(obj.data)
		return {'FINISHED'}

# Basic Mitsuba integration based on the POV-Ray add-on
# Piggybacks on the COLLADA exporter to get most things done
class MitsubaRender(bpy.types.RenderEngine):
	bl_idname = 'MITSUBA_RENDER'
	bl_label = "Mitsuba"

	def _export_worldtrafo(self, adjfile, trafo):
		adjfile.write('\t\t<transform name="toWorld">\n')
		adjfile.write('\t\t\t<matrix value="')
		for j in range(0,4):
			for i in range(0,4):
				adjfile.write("%f " % trafo[i][j])
		adjfile.write('"/>\n\t\t</transform>\n')


	def _export_lamp(self, adjfile, lamp, idx):
		if lamp.data.type == 'POINT':
			adjfile.write('\t<luminaire id="%s-light" type="point">\n' % lamp.data.name)
			mult = lamp.data.energy * lamp.data.distance * lamp.data.distance / 2
			self._export_worldtrafo(adjfile, lamp.matrix_world)
			adjfile.write('\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			adjfile.write('\t</luminaire>\n')
		elif lamp.data.type == 'AREA':
			adjfile.write('\t<shape type="obj">\n')
			size_x = lamp.data.size
			size_y = lamp.data.size
			if lamp.data.shape == 'RECTANGLE':
				size_y = lamp.data.size_y
			path = os.path.join(os.path.join(self._temp_dir, 'meshes'), "_area_luminaire_%d.obj" % idx)

			adjfile.write('\t\t<string name="filename" value="%s"/>\n' % path)
			self._export_worldtrafo(adjfile, lamp.matrix_world)

			adjfile.write('\n\t\t<luminaire id="%s-light" type="area">\n' % lamp.data.name)
			mult = lamp.data.energy * lamp.data.distance * lamp.data.distance / (size_x * size_y)
			adjfile.write('\t\t\t<rgb name="intensity" value="%f %f %f"/>\n' 
					% (lamp.data.color.r*mult, lamp.data.color.g*mult, lamp.data.color.b*mult))
			adjfile.write('\t\t</luminaire>\n')
			adjfile.write('\t</shape>\n')
		
			objFile = open(path, 'w')
			objFile.write('v %f %f 0\n' % (-size_x/2, -size_y/2))
			objFile.write('v %f %f 0\n' % ( size_x/2, -size_y/2))
			objFile.write('v %f %f 0\n' % ( size_x/2,  size_y/2))
			objFile.write('v %f %f 0\n' % (-size_x/2,  size_y/2))
			objFile.write('f 4 3 2 1\n')
			objFile.close()

	def _export(self, scene):
		import tempfile
		self._temp_dir = tempfile.mkdtemp(prefix='mitsuba_');
		self._temp_dae = os.path.join(self._temp_dir, 'scene.dae')
		self._temp_xml = os.path.join(self._temp_dir, 'scene.xml')
		self._temp_adj = os.path.join(self._temp_dir, 'scene_adjustments.xml')
		self._temp_out = os.path.join(self._temp_dir, 'scene.png')
		os.mkdir(os.path.join(self._temp_dir, 'meshes'))
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
		idx = 0
		for obj in scene.objects:
			if obj.type == 'LAMP':
				self._export_lamp(adjfile, obj, idx)
			idx = idx+1
		adjfile.write('</adjustments>\n');
		adjfile.close()

	def _render(self):
		scene = bpy.data.scenes[0]
		(mts_path, tail) = os.path.split(bpy.path.abspath(scene.mts_path))
		mtsimport_binary = os.path.join(mts_path, "mtsimport")
		mitsuba_binary = os.path.join(mts_path, "mitsuba")
		mts_render_libpath = os.path.join(mts_path, "src/librender")
		mts_core_libpath = os.path.join(mts_path, "src/libcore")
		mts_hw_libpath = os.path.join(mts_path, "src/libhw")
		env = copy.copy(os.environ)
		env['LD_LIBRARY_PATH'] = mts_core_libpath + ":" + mts_render_libpath + ":" + mts_hw_libpath
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
		print("Not cleaning up")
	 	#shutil.rmtree(self._temp_dir)

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
