def resolution(scene):
	'''
	scene		bpy.types.scene
	Calculate the output render resolution
	Returns		tuple(2) (floats)
	'''
	xr = scene.render.resolution_x * scene.render.resolution_percentage / 100.0
	yr = scene.render.resolution_y * scene.render.resolution_percentage / 100.0
	
	return xr, yr
