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

def extract(target, filename):
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

pyVer = int(platform.python_version_tuple()[0])

# Traverse source directories and process any found plugin code
def process(path, target):
	def capture(fileList, dirname, files):
		suffix = os.path.split(dirname)[1]
		if 'lib' in suffix or suffix == 'tests' \
			or suffix == 'mitsuba' or suffix == 'utils' \
			or suffix == 'converter' or suffix == 'mtsgui':
			return
		for filename in files:
			if '.cpp' == os.path.splitext(filename)[1]:
				fname = os.path.join(dirname, filename)
				fileList += [fname]

	fileList = []
	# Wrap the walk function to make this work in python 2 and 3.
	if pyVer >= 3:
		os.walk(path, capture, fileList)
	else:
		os.path.walk(path, capture, fileList)

	ordering = [(findOrderID(fname), fname) for fname in fileList]
	ordering = sorted(ordering, key = lambda entry: entry[0])

	for entry in ordering:
		extract(target, entry[1])

os.chdir(os.path.dirname(__file__))
f = open('plugins_generated.tex', 'w')
f.write('\input{section_shapes}\n')
process('../src/shapes', f)
f.write('\input{section_bsdf}\n')
process('../src/bsdfs', f)
f.write('\input{section_textures}\n')
process('../src/textures', f)
f.write('\input{section_subsurface}\n')
process('../src/subsurface', f)
f.write('\input{section_media}\n')
process('../src/medium', f)
f.write('\input{section_phase}\n')
process('../src/phase', f)
f.write('\input{section_volumes}\n')
process('../src/volume', f)
f.write('\input{section_emitters}\n')
process('../src/emitters', f)
f.write('\input{section_sensors}\n')
process('../src/sensors', f)
f.write('\input{section_integrators}\n')
process('../src/integrators', f)
f.write('\input{section_samplers}\n')
process('../src/samplers', f)
f.write('\input{section_films}\n')
process('../src/films', f)
f.write('\input{section_rfilters}\n')
f.close()
os.system('bibtex main.aux')
os.system('pdflatex main.tex')
#os.system('pdflatex main.tex | grep -i warning | grep -v "Package \(typearea\|hyperref\)"')
