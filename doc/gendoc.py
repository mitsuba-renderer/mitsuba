#!/usr/bin/env python
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
	for (dirname, subdirs, files) in os.walk(path):
	    capture(fileList, dirname, files)

	ordering = [(findOrderID(fname), fname) for fname in fileList]
	ordering = sorted(ordering, key = lambda entry: entry[0])

	for entry in ordering:
		extract(target, entry[1])

def process_src(target, src_subdir, section=None):
	if section is None:
		section = "section_" + src_subdir
	target.write('\input{{{0}}}\n'.format(section))
	process('../src/{0}'.format(src_subdir), target)

def texify(texfile):
	from subprocess import Popen, PIPE, check_call
	version = Popen(["pdflatex", "-version"], stdout=PIPE).communicate()[0]
	# Call decode() to convert from bytes to string, required in Python 3
	if re.match('.*MiKTeX.*', version.decode()):
		# MiKTeX's "texify" calls latex/bibtex in tandem automatically
		print("Running texify on {0}...".format(texfile))
		check_call(['texify', '-pq', texfile])
	else:
		check_call(['pdflatex', texfile])
		check_call(['bibtex',   texfile.replace('.tex', '.aux')])
		check_call(['pdflatex', texfile])
		check_call(['pdflatex', texfile])

os.chdir(os.path.dirname(os.path.abspath(__file__)))
with open('plugins_generated.tex', 'w') as f:
	process_src(f, 'shapes')
	process_src(f, 'bsdfs', 'section_bsdf')
	process_src(f, 'textures')
	process_src(f, 'subsurface')
	process_src(f, 'medium', 'section_media')
	process_src(f, 'phase')
	process_src(f, 'volume', 'section_volumes')
	process_src(f, 'emitters')
	process_src(f, 'sensors')
	process_src(f, 'integrators')
	process_src(f, 'samplers')
	process_src(f, 'films')
	process_src(f, 'rfilters')

texify('main.tex')
