class MtsAdjustments:
	def __init__(self, target_file):
		self.target_file = target_file

	def export_worldtrafo(self, adjfile, trafo):
		adjfile.write('\t\t<transform name="toWorld">\n')
		adjfile.write('\t\t\t<matrix value="')
		for j in range(0,4):
			for i in range(0,4):
				adjfile.write("%f " % trafo[i][j])
		adjfile.write('"/>\n\t\t</transform>\n')

	def export_lamp(self, adjfile, lamp, idx):
		if lamp.data.type == 'POINT':
			adjfile.write('\t<luminaire id="%s-light" type="point">\n' % lamp.data.name)
			mult = lamp.data.energy * lamp.data.distance * lamp.data.distance / 2
			self.export_worldtrafo(adjfile, lamp.matrix_world)
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
			self.export_worldtrafo(adjfile, lamp.matrix_world)

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

	def export(self, scene):
		adjfile = open(self.target_file, 'w')
		adjfile.write('<adjustments>\n');
		idx = 0
		for obj in scene.objects:
			if obj.type == 'LAMP':
				self.export_lamp(adjfile, obj, idx)
			idx = idx+1
		adjfile.write('</adjustments>\n');
		adjfile.close()
