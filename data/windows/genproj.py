from lxml import etree
import os, shutil

doc = etree.parse('data/windows/mitsuba.vcproj.template')

headers = etree.ETXPath('/VisualStudioProject/Files/Filter[@Name="Header Files"]')(doc)[0]
sources = etree.ETXPath('/VisualStudioProject/Files/Filter[@Name="Source Files"]')(doc)[0]


def traverse(dirname, base):
	for file in [file for file in os.listdir(dirname) if not file in ['.', '..']]:
		filename = os.path.join(dirname, file)
		if os.path.isdir(filename):
			if filename == '.\\include\\mitsuba':
				traverse(filename, base)
			else:
				node = etree.SubElement(base, 'Filter')
				node.set('Name', os.path.split(filename)[1])
				node.tail = '\n\t\t'
				node.text = '\n\t\t'
				traverse(filename, node)
		else:
			ext = os.path.splitext(filename)[1]
			if ext == '.cpp' or ext == '.c' or ext == '.h':
				node = etree.SubElement(base, 'File')
				node.set('RelativePath', filename)
				node.tail = '\n\t\t'

traverse('.\\src', sources)
traverse('.\\include', headers)
of = open('mitsuba.vcproj', 'w')
of.write(etree.tostring(doc, pretty_print=True))
of.close()
shutil.copyfile('data/windows/mitsuba.sln.template', 'mitsuba.sln')
