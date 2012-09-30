#!/usr/bin/env python
#
#  This file is part of Mitsuba, a physically based rendering system.
#
#  Copyright (c) 2007-2012 by Wenzel Jakob and others.
#
#  Mitsuba is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License Version 3
#  as published by the Free Software Foundation.
#
#  Mitsuba is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program. If not, see <http://www.gnu.org/licenses/>.

from PySide.QtCore import QFile, QTextStream
from xml.etree import ElementTree as et
from lepl import Regexp, Space, Literal, Separator, List
import numpy as np
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import sys

Real = lambda: Regexp(r'[\+\-]?(?:[0-9]*\.[0-9]+|[0-9]+\.' +
		r'|[0-9]+)(?:[eE][\+\-]?[0-9]+)?')

def make_path_grammar():
	sep = ~(Space() | Literal(','))[:]
	with Separator(sep):
		num = Real() >> float
		# Moveto
		M = ((Literal('M') | Literal('m')) & num[2][:]) 
		# Horizontal straight lines
		H = (Literal('H') | Literal('h')) & num[:]
		# Vertical straight lines
		V = (Literal('V') | Literal('v')) & num[:]
		# General straight lines
		L = (Literal('L') | Literal('l')) & num[2][:] 
		# Cubic bezier curves (curveto)
		C = (Literal('C') | Literal('c')) & num[6][:]
		# Cubic bezier curves (smooth curveto)
		S = (Literal('S') | Literal('s')) & num[4][:]
		# Close the path
		z = Literal('z') | Literal('Z')
		grammar = sep & ((M|H|V|L|C|S|z) > List)[:] & sep 
		grammar.config.no_compile_to_regexp()
		return grammar

def make_polygon_grammar():
	sep = ~(Space() | Literal(','))[:]
	with Separator(sep):
		num = Real() >> float
		grammar = sep & num[2][:] & sep 
		grammar.config.no_compile_to_regexp()
		return grammar

class BezierSpline(object):
	def __init__(self, *args):
		if len(args) == 4:
			self.start = args[0]
			self.cp1 = args[1]
			self.cp2 = args[2]
			self.end = args[3]
		elif len(args) == 2:
			d = args[1] - args[0]
			self.start = args[0]
			self.cp1 = args[0] + 1.0/3.0 * d
			self.cp2 = args[0] + 2.0/3.0 * d
			self.end = args[1]
		else:
			raise Exception("Invalid constructor call")

	def _eval(self, t):
		tmp = 1 - t
		tmp2, t2 = tmp*tmp, t*t
		return self.start * (tmp*tmp2) + self.cp1 * (3*tmp2*t) + \
			self.cp2 * (3*tmp*t2) + self.end * (t*t2) 

	def drawGL(self):
		steps = 10
		p = self.start
		glVertex2f(p[0], p[1])
		for i in range(1, steps):
			p = self._eval(float(i)/(steps-1))
			glVertex2f(p[0], p[1])

	def drawGL_tess(self):
		steps = 10
		p = self.start
		gluTessVertex(tobj, [p[0], p[1], 0], [p[0], p[1], 0])
		for i in range(1, steps):
			p = self._eval(float(i)/(steps-1))
			gluTessVertex(tobj, [p[0], p[1], 0], [p[0], p[1], 0])

class AABB(object):
	def __init__(self):
		inf = float("inf")
		self.min = np.array([inf,  inf])
		self.max = np.array([-inf, -inf])
	
	def expand_by(self, p):
		self.min[0] = min(self.min[0], p[0])
		self.min[1] = min(self.min[1], p[1])
		self.max[0] = max(self.max[0], p[0])
		self.max[1] = max(self.max[1], p[1])

	def expand_by_aabb(self, aabb):
		self.min[0] = min(self.min[0], aabb.min[0])
		self.min[1] = min(self.min[1], aabb.min[1])
		self.max[0] = max(self.max[0], aabb.max[0])
		self.max[1] = max(self.max[1], aabb.max[1])

	def size(self):
		return self.max - self.min

	def __repr__(self):
		return "AABB[min=%s, max=%s]" % \
			(repr(self.min), repr(self.max))

class Path(object):
	PathGrammar = make_path_grammar()
	PolygonGrammar = make_polygon_grammar()

	def __init__(self, node):
		self.pos = None
		self.start = None
		self.cp2 = None
		self.splines = []
		self.aabb = AABB()

		def getflt(key):
			value = node.get(key)
			return float(value) if value != None else 0

		if 'path' in node.tag:
			instructions = Path.PathGrammar.parse(node.get("d"))
		elif 'polygon' in node.tag:
			points = Path.PolygonGrammar.parse(node.get("points"))
			instructions = [['M'] + points, ['z']]
		elif 'rect' in node.tag:
			x, y = getflt("x"), getflt("y")
			width, height = getflt("width"), getflt("height")
			instructions = [['M', x, y ], ['h', width], ['v', height], ['h', -width], ['z']]
		elif 'line' in node.tag:
			x1, x2 = getflt("x1"), getflt("x2")
			y1, y2 = getflt("y1"), getflt("y2")
			instructions = [['M', x1, y1, x2, y2 ]]
		else:
			raise Exception("Unknown tag!")

		self.id = node.get("id")

		self.stroke = self._color(node.get('stroke'))
		self.fill = self._color(node.get('fill'))

		commandList = {
			'm' : Path._moveto,
			'h' : Path._hlineto,
			'v' : Path._vlineto,
			'l' : Path._lineto,
			'c' : Path._curveto,
			's' : Path._scurveto,
			'z' : Path._close
		}

		self.index = 1
		for item in instructions:
			cmd, args = item[0], item[1:]
			commandList[cmd.lower()](self, cmd, args)
			self.lastcmd = cmd
			self.index += 1

		for spline in self.splines:
			if spline:
				self.aabb.expand_by(spline.start)
				self.aabb.expand_by(spline.cp1)
				self.aabb.expand_by(spline.cp2)
				self.aabb.expand_by(spline.end)

	def _moveto(self, cmd, args):
		if cmd == 'M':
			self.pos = np.array(args[0:2])
		else:
			self.pos =+ np.array(args[0:2])
		if self.start is None or self.lastcmd in ['z', 'Z']:
			self.start = self.pos
		if len(args) > 2:
			# Implicit lineto
			cmd = 'L' if cmd == 'M' else 'l'
			self._lineto(cmd, args[2:])

	def _lineto(self, cmd, args):
		if cmd == 'L':
			end = np.array(args[0:2])
		else:
			end = self.pos + np.array(args[0:2])
		d = end - self.pos
		self.splines.append(BezierSpline(self.pos, end))
		self.pos = end
		if len(args) > 2:
			self._lineto(cmd, args[2:])

	def _hlineto(self, cmd, args):
		if cmd == 'H':
			end = np.array([args[0], self.pos[1]])
		else:
			end = np.array([args[0] + self.pos[0], self.pos[1]])
		self.splines.append(BezierSpline(self.pos, end))
		self.pos = end
		if len(args) > 1:
			self._hlineto(cmd, args[1:])

	def _vlineto(self, cmd, args):
		if cmd == 'V':
			end = np.array([self.pos[0], args[0]])
		else:
			end = np.array([self.pos[0], args[0] + self.pos[1]])
		self.splines.append(BezierSpline(self.pos, end))
		self.pos = end
		if len(args) > 1:
			self._vlineto(cmd, args[1:])

	def _curveto(self, cmd, args):
		start = self.pos
		cp1 = np.array(args[0:2])
		cp2 = np.array(args[2:4])
		end = np.array(args[4:6])

		if cmd == 'c':
			cp1 += start
			cp2 += start
			end += start
		
		self.splines.append(BezierSpline(start, cp1, cp2, end))
		self.pos = end
		self.cp2 = cp2

		if len(args) > 6:
			self._curveto(cmd, args[6:])

	def _scurveto(self, cmd, args):
		start = self.pos
		cp2 = np.array(args[0:2])
		end = np.array(args[2:4])

		if cmd == 's':
			cp2 += start
			end += start

		if self.cp2 != None and self.lastcmd in ['C', 'c', 'S', 's']:
			cp1 = 2*start - self.cp2
		else:
			cp1 = start
		
		self.splines.append(BezierSpline(start, cp1, cp2, end))
		self.pos = end
		self.cp2 = cp2

		if len(args) > 4:
			self._scurveto(cmd, args[4:])

	def _close(self, cmd, args):
		self._lineto('L', self.start)
		self.splines.append(None)

	def _color(self, value):
		if value == None or value.lower() == 'none':
			return None
		r = int(value[1:3], 16) / 255.0
		g = int(value[3:5], 16) / 255.0
		b = int(value[5:7], 16) / 255.0
		return [r, g, b]

	def drawGL(self):
		if self.fill:
			glColor4f(self.fill[0], self.fill[1], self.fill[2], 1.0)
			gluTessBeginPolygon(tobj, None)
			gluTessBeginContour(tobj) 
			for spline in self.splines:
				if spline:
					spline.drawGL_tess()
				else:
					gluTessEndContour(tobj)
					gluTessBeginContour(tobj) 
			gluTessEndContour(tobj)
			gluTessEndPolygon(tobj)

		if self.stroke:
			glColor4f(self.stroke[0], self.stroke[1], self.stroke[2], 1.0)
			glBegin(GL_LINE_STRIP)
			for spline in self.splines:
				if spline:
					spline.drawGL()
				else:
					glEnd()
					glBegin(GL_LINE_LOOP)
			glEnd()

	def export(self, out):
		if self.id:
			out.write("obj %s\n" % self.id)
		else:
			out.write("obj\n")

		if self.stroke:
			out.write("stroke %f %f %f\n" % (self.stroke[0], self.stroke[1], self.stroke[2]))
		if self.fill:
			out.write("fill %f %f %f\n" % (self.fill[0], self.fill[1], self.fill[2]))
		for spline in self.splines:
			if spline is None:
				out.write("skip\n")
			else:
				out.write("spline %f %f %f %f %f %f %f %f\n" %
					(spline.start[0], spline.start[1], spline.cp1[0], spline.cp1[1],
					 spline.cp2[0], spline.cp2[1], spline.end[0], spline.end[1]))

		out.write("\n")

class Scene(object):
	def __init__(self, filename):
		ns = '{http://www.w3.org/2000/svg}'
		self.paths = []
		self.aabb = AABB()
		def register(node):
			if node.tag == ns + 'path' or node.tag == ns + 'polygon' \
				or node.tag == ns + 'line' or node.tag == ns + 'rect':
				path = Path(node)
				self.aabb.expand_by_aabb(path.aabb)
				self.paths.append(path)
			else:
				for child in node:
					register(child)
				return
		register(et.parse(filename).getroot())

	def drawGL(self):
		glClearColor(.3, .3, .3, 1.0)
		glClear(GL_COLOR_BUFFER_BIT)
		size = self.aabb.size().max()
		xs = self.aabb.min[0] - size * 0.1
		ys = self.aabb.min[1] - size * 0.1
		size = size * 1.2
		glOrtho(xs, size, size*6.0/8.0, ys, 0, 1)
		for path in self.paths:
			path.drawGL()
		glutSwapBuffers()

	def export(self, target):
		with open(target, "w") as f:
			for path in self.paths:
				path.export(f)

def keyboard(key, x, y):
	if key == 'q' or key == "\x1b":
		sys.exit(0)

def initGL():
	global tobj
	tobj = gluNewTess()
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
	glEnable(GL_LINE_SMOOTH)
	glEnable(GL_POLYGON_SMOOTH)
	glEnable(GL_BLEND)
	glLineWidth(2.0)
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
	gluTessCallback(tobj, GLU_TESS_VERTEX, glVertex3fv)
	gluTessCallback(tobj, GLU_TESS_BEGIN, lambda x: glBegin(x))
	gluTessCallback(tobj, GLU_TESS_END, lambda: glEnd())
	gluTessCallback(tobj, GLU_TESS_COMBINE, lambda pos, data, weights: pos)

if __name__ == '__main__':
	if len(sys.argv) != 3:
		print('Syntax: importsvg.py <input SVG file> <output SC2 file>')
	scene = Scene(sys.argv[1])
	scene.export(sys.argv[2])
	glutInit(sys.argv)
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB)
	glutInitWindowSize(800, 600)
	glutCreateWindow("Bezier spline importer")
	glutDisplayFunc(lambda: scene.drawGL())
	glutKeyboardFunc(keyboard)
	initGL()
	glutMainLoop()
