from lxml import etree
import os, string, shutil, uuid

doc2017 = etree.parse('data/windows/mitsuba-msvc2017.vcxproj.template')
doc2017_filters = etree.parse('data/windows/mitsuba-msvc2017.vcxproj.filters.template')

ns = {'n' : 'http://schemas.microsoft.com/developer/msbuild/2003'}
headers2017 = etree.XPath('/n:Project/n:ItemGroup[@Label="Header Files"]', namespaces = ns)(doc2017)[0]
sources2017 = etree.XPath('/n:Project/n:ItemGroup[@Label="Source Files"]', namespaces = ns)(doc2017)[0]
headers2017_filters = etree.XPath('/n:Project/n:ItemGroup[@Label="Header Files"]', namespaces = ns)(doc2017_filters)[0]
sources2017_filters = etree.XPath('/n:Project/n:ItemGroup[@Label="Source Files"]', namespaces = ns)(doc2017_filters)[0]
filters2017 = etree.XPath('/n:Project/n:ItemGroup[@Label="Filters"]', namespaces = ns)(doc2017_filters)[0]
def traverse(dirname, prefix):
	for file in [file for file in os.listdir(dirname) if not file in ['.', '..']]:
		filename = os.path.join(dirname, file)
		if os.path.isdir(filename):
			lastname = os.path.split(filename)[1]
			# Visual Studio 2017 nodes
			subprefix = os.path.join(prefix, lastname)
			node2017 = etree.SubElement(filters2017, 'Filter')
			node2017.set('Include', subprefix.replace('/', '\\'))
			ui = etree.SubElement(node2017, 'UniqueIdentifier')
			ui.text = '{' + str(uuid.uuid4()) + '}'
			ui.tail = '\n\t\t'
			node2017.tail = '\n\t\t'
			node2017.text = '\n\t\t\t'

			traverse(filename, subprefix)
		else:
			ext = os.path.splitext(filename)[1]

			filename = '..\\' + filename.replace('/', '\\')
			prefix = prefix.replace('/', '\\')

			# Visual Studio 2017 nodes
			if ext == '.cpp' or ext == '.c':
				node = etree.SubElement(sources2017, 'ClCompile')
				node.set('Include', filename)
				node.tail = '\n\t\t'
				node.text = '\n\t\t\t'
				node = etree.SubElement(sources2017_filters, 'ClCompile')
				node.set('Include', filename)
				node.tail = '\n\t\t'
				node.text = '\n\t\t\t'
				filter = etree.SubElement(node, 'Filter')
				filter.text = prefix
				filter.tail = '\n\t\t'
			elif ext == '.h' or ext == '.inl':
				node = etree.SubElement(headers2017, 'ClInclude')
				node.set('Include', filename)
				node.tail = '\n\t\t'
				node.text = '\n\t\t\t'
				node = etree.SubElement(headers2017_filters, 'ClInclude')
				node.set('Include', filename)
				node.tail = '\n\t\t'
				node.text = '\n\t\t\t'
				filter = etree.SubElement(node, 'Filter')
				filter.text = prefix
				filter.tail = '\n\t\t'


traverse('src', 'Source Files')
traverse('include', 'Header Files')

of = open('build/mitsuba-msvc2017.vcxproj', 'w')
of.write(etree.tostring(doc2017, pretty_print=True))
of.close()

of = open('build/mitsuba-msvc2017.vcxproj.filters', 'w')
of.write(etree.tostring(doc2017_filters, pretty_print=True))
of.close()
