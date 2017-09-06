# Detect present Python & Boost-python libraries on Linux
import os, struct

class PkgConfig(dict):
    _paths = [
        '/usr/lib/pkgconfig',
        '/usr/lib/%s-linux-gnu/pkgconfig' % (os.uname()[4]),
        '/usr/lib%i/pkgconfig' % (struct.calcsize('P')*8)
    ]

    def __init__(self, name):
        for path in self._paths:
            fn = os.path.join(path, '%s.pc' % name)
            if os.path.exists(fn):
                self._parse(fn)
                break

    def _parse(self, filename):
        from string import Template

        lines = open(filename).readlines()
        localVariables = {}

        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            elif ':' in line:
                name, val = line.split(':')
                val = val.strip()
                if '$' in val:
                    val = Template(val).substitute(localVariables)
                self[name] = val
            elif '=' in line:
                name, val = line.split('=')
                val = val.strip()
                if '$' in val:
                    val = Template(val).substitute(localVariables)
                localVariables[name] = val

def find_boost_python(version):
    libnames = [
        'boost_python-mt-py%s' % version,
        'boost_python-py%s' % version,
        'boost_python' + ('3' if version.startswith('3') else '')
    ]
    basepaths = [
        '/usr/lib',
        '/usr/lib/%s-linux-gnu' % (os.uname()[4]),
        '/usr/lib%i' % (struct.calcsize('P')*8)
    ]

    for basepath in basepaths:
        for libname in libnames:
            if os.path.isfile(os.path.join(basepath, "lib" + libname + ".so")):
                return libname
    return None

def detect_python():
    pyver = ['2.6', '2.7', '3.0', '3.1', '3.2', '3.3', '3.4', '3.5', '3.6']
    pyenv = {}

    for version in pyver:
        pkgconfig = PkgConfig('python-%s' % version)
        version = version.replace('.', '')
        flags = []
        if 'Cflags' in pkgconfig:
            flags += pkgconfig['Cflags'].split()
        if 'Libs' in pkgconfig:
            flags += pkgconfig['Libs'].split()
        if len(flags) == 0:
            continue
        boost_libname = find_boost_python(version)
        if boost_libname == None:
            continue
        pyenv['PYTHON' + version + 'INCLUDE'] = []
        pyenv['PYTHON' + version + 'LIBDIR'] = []
        pyenv['PYTHON' + version + 'LIB'] = [ boost_libname ]
        for flag in flags:
            if flag.startswith('-I'):
                pyenv['PYTHON' + version + 'INCLUDE'] += [flag[2:]]
            elif flag.startswith('-L'):
                pyenv['PYTHON' + version + 'LIBDIR'] += [flag[2:]]
            elif flag.startswith('-l'):
                pyenv['PYTHON' + version + 'LIB'] += [flag[2:]]
    return pyenv

if __name__ == '__main__':
    import pprint
    pprint.pprint(detect_python())
