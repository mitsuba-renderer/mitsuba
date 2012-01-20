#! /usr/bin/python
# 
# This script walks through all plugin files and 
# extracts documentation that should go into the
# reference manual

import os, platform, re

def findOrderID(filename):
	f = open(filename)
	for line in f.readlines():
		match = re.match(r'.*\\order{([^}]*)}.*', line)
		if match != None:
			return int(match.group(1))
	return 1000

def process(target, filename):
	f = open(filename)
	inheader = False
	for line in f.readlines():
		match = re.match(r'^/\*! ?(.*)$', line)
		if match != None:
			print("Processing %s" % filename)
			line = match.group(1).replace('%', '\%')
			target.write(line + '\n')
			inheader = True
			continue
		if not inheader:
			continue
		if re.search(r'^[\s\*]*\*/$', line):
			inheader = False
			continue
		match = re.match(r'^\s*\** ?(.*)$', line)
		if match != None:
			line = match.group(1).replace('%', '\%')
			target.write(line + '\n')
	f.close()

# Traverse all source directories and find any plugin code
def traverse(target, dirname, files):
	suffix = os.path.split(dirname)[1]
	if 'lib' in suffix or suffix == 'tests' \
		or suffix == 'mitsuba' or suffix == 'utils' \
		or suffix == 'converter' or suffix == 'qtgui':
		return

	ordering = []
	for filename in files:
		if '.cpp' == os.path.splitext(filename)[1]:
			fname = os.path.join(dirname, filename)
			ordering = ordering + [(findOrderID(fname), fname)]
	ordering = sorted(ordering, key = lambda entry: entry[0])

	for entry in ordering:
		process(target, entry[1])

# Wrap the walk function to make this work in python 2 and 3.
pyVer = int(platform.python_version_tuple()[0])
def walk(path, visit, arg):
	if pyVer >= 3:
		os.walk(path, visit, arg)
	else:
		os.path.walk(path, visit, arg)

os.chdir(os.path.dirname(__file__))
f = open('plugins_generated.tex', 'w')
f.write('\section{Plugin reference}\n')
f.write('\input{section_shapes}\n')
walk('../src/shapes', traverse, f)
f.write('\input{section_bsdf}\n')
walk('../src/bsdfs', traverse, f)
f.write('\input{section_textures}\n')
walk('../src/textures', traverse, f)
f.write('\input{section_subsurface}\n')
walk('../src/subsurface', traverse, f)
f.write('\input{section_media}\n')
walk('../src/medium', traverse, f)
f.write('\input{section_phase}\n')
walk('../src/phase', traverse, f)
f.write('\input{section_volumes}\n')
walk('../src/volume', traverse, f)
f.write('\input{section_luminaires}\n')
walk('../src/luminaires', traverse, f)
f.write('\input{section_integrators}\n')
walk('../src/integrators', traverse, f)
f.write('\input{section_films}\n')
walk('../src/films', traverse, f)
f.close()
os.system('bibtex main.aux')
os.system('pdflatex main.tex | grep -i warning')
