from lxml import etree
import os, string, shutil, uuid

doc2010 = etree.parse('data/windows/mitsuba-msvc2010.vcxproj.template')
doc2010_filters = etree.parse('data/windows/mitsuba-msvc2010.vcxproj.filters.template')

ns = {'n' : 'http://schemas.microsoft.com/developer/msbuild/2003'}
headers2010 = etree.XPath('/n:Project/n:ItemGroup[@Label="Header Files"]', namespaces = ns)(doc2010)[0]
sources2010 = etree.XPath('/n:Project/n:ItemGroup[@Label="Source Files"]', namespaces = ns)(doc2010)[0]
headers2010_filters = etree.XPath('/n:Project/n:ItemGroup[@Label="Header Files"]', namespaces = ns)(doc2010_filters)[0]
sources2010_filters = etree.XPath('/n:Project/n:ItemGroup[@Label="Source Files"]', namespaces = ns)(doc2010_filters)[0]
filters2010 = etree.XPath('/n:Project/n:ItemGroup[@Label="Filters"]', namespaces = ns)(doc2010_filters)[0]
def traverse(dirname, prefix):
	for file in [file for file in os.listdir(dirname) if not file in ['.', '..']]:
		filename = os.path.join(dirname, file)
		if os.path.isdir(filename):
			lastname = os.path.split(filename)[1]
			# Visual Studio 2010 nodes
			subprefix = os.path.join(prefix, lastname)
			node2010 = etree.SubElement(filters2010, 'Filter')
			node2010.set('Include', subprefix.replace('/', '\\'))
			ui = etree.SubElement(node2010, 'UniqueIdentifier')
			ui.text = '{' + str(uuid.uuid4()) + '}'
			ui.tail = '\n\t\t'
			node2010.tail = '\n\t\t'
			node2010.text = '\n\t\t\t'

			traverse(filename, subprefix)
		else:
			ext = os.path.splitext(filename)[1]

			filename = '..\\' + filename.replace('/', '\\')
			prefix = prefix.replace('/', '\\')

			# Visual Studio 2010 nodes
			if ext == '.cpp' or ext == '.c':
				node = etree.SubElement(sources2010, 'ClCompile')
				node.set('Include', filename)
				node.tail = '\n\t\t'
				node.text = '\n\t\t\t'
				node = etree.SubElement(sources2010_filters, 'ClCompile')
				node.set('Include', filename)
				node.tail = '\n\t\t'
				node.text = '\n\t\t\t'
				filter = etree.SubElement(node, 'Filter')
				filter.text = prefix
				filter.tail = '\n\t\t'
			elif ext == '.h' or ext == '.inl':
				node = etree.SubElement(headers2010, 'ClInclude')
				node.set('Include', filename)
				node.tail = '\n\t\t'
				node.text = '\n\t\t\t'
				node = etree.SubElement(headers2010_filters, 'ClInclude')
				node.set('Include', filename)
				node.tail = '\n\t\t'
				node.text = '\n\t\t\t'
				filter = etree.SubElement(node, 'Filter')
				filter.text = prefix
				filter.tail = '\n\t\t'


traverse('src', 'Source Files')
traverse('include', 'Header Files')

of = open('build/mitsuba-msvc2010.vcxproj', 'w')
of.write(etree.tostring(doc2010, pretty_print=True))
of.close()

of = open('build/mitsuba-msvc2010.vcxproj.filters', 'w')
of.write(etree.tostring(doc2010_filters, pretty_print=True))
of.close()
