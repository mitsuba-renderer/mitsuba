#! /usr/bin/python
# 
# This script walks through all plugin files and 
# extracts documentation that should go into the
# reference manual

import os, re

def findOrderID(filename):
	f = open(filename)
	for line in f.readlines():
		match = re.match(r'.*\\order{([^}])}.*', line)
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

os.chdir(os.path.dirname(__file__))
f = open('plugins_generated.tex', 'w')
f.write('\section{Plugin reference}\n')
f.write('\input{section_bsdf}\n')
os.path.walk('../src/bsdfs', traverse, f)
f.write('\input{section_subsurface}\n')
os.path.walk('../src/subsurface', traverse, f)
f.write('\input{section_media}\n')
os.path.walk('../src/media', traverse, f)
f.write('\input{section_integrators}\n')
os.path.walk('../src/integrators', traverse, f)
f.close()
os.system('bibtex main.aux')
os.system('pdflatex main.tex')
