import os
import bpy
from extensions_framework import log
from extensions_framework.util import TimerThread

def MtsLog(*args, popup=False):
	'''
	Send string to AF log, marked as belonging to Mitsuba module.
	Accepts variable args 
	'''
	if len(args) > 0:
		log(' '.join(['%s'%a for a in args]), module_name='Mitsuba', popup=popup)

class MtsFilmDisplay(TimerThread):
	'''
	Periodically update render result with Mituba's framebuffer
	'''

	STARTUP_DELAY = 1

	def begin(self, renderer, output_file, resolution):
		(xres, yres) = resolution
		self.result = renderer.begin_result(0, 0, int(xres), int(yres))
		self.renderer = renderer
		self.output_file = output_file
		self.resolution = resolution
		self.start()

	def shutdown(self):
		self.renderer.end_result(self.result)

	def kick(self, render_end=False):
		if not bpy.app.background or render_end:
			if render_end:
				MtsLog('Final render result %ix%i' % self.resolution)
			else:
				MtsLog('Updating render result %ix%i' % self.resolution)

			if os.path.exists(self.output_file):
				try:
					self.result.layers[0].load_from_file(self.output_file)
					self.renderer.update_result(self.result)
				except:
					pass
			else:
				err_msg = 'ERROR: Could not load render result from %s' % self.output_file
				MtsLog(err_msg)
