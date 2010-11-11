import os
import bpy
from extensions_framework import log
from extensions_framework.util import TimerThread

def MtsLog(*args, popup=False):
	'''
	Send string to AF log, marked as belonging to Mitsuba module.
	Accepts variable args (can be used as pylux.errorHandler)
	'''
	if len(args) > 0:
		log(' '.join(['%s'%a for a in args]), module_name='Mitsuba', popup=popup)

class MtsFilmDisplay(TimerThread):
	'''
	Periodically update render result with Mituba's framebuffer
	'''

	STARTUP_DELAY = 1

	def kick(self, render_end=False):
		xres, yres = self.LocalStorage['resolution']

		if render_end:
			MtsLog('Final render result %ix%i' % (xres,yres))
		else:
			MtsLog('Updating render result %ix%i' % (xres,yres))

		result = self.LocalStorage['RE'].begin_result(0, 0, int(xres), int(yres))
		if os.path.exists(self.LocalStorage['RE'].output_file):
			bpy.ops.ef.msg(msg_text='Updating RenderResult')
			lay = result.layers[0]
			lay.load_from_file(self.LocalStorage['RE'].output_file)
		else:
			err_msg = 'ERROR: Could not load render result from %s' % self.LocalStorage['RE'].output_file
			MtsLog(err_msg)
			bpy.ops.ef.msg(msg_type='ERROR', msg_text=err_msg)
		self.LocalStorage['RE'].end_result(result)
