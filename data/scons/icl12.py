import os, sys, subprocess, copy, re

def get_output(script, args = None, shellenv = None):
	if sys.platform == 'win32':
		cmdLine = '"%s" %s & set' % (script, (args if args else ''))
		shell = False
	elif sys.platform.startswith('linux'):
		cmdLine = 'source "%s" %s ; set' % (script, (args if args else ''))
		shell = True
	else:
		raise Exception("Unsuported OS type: " + sys.platform)

	popen = subprocess.Popen(cmdLine, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=shellenv)

	# Use the .stdout and .stderr attributes directly because the
	# .communicate() method uses the threading module on Windows
	# and won't work under Pythons not built with threading.
	stdout = popen.stdout.read()
	if popen.wait() != 0:
		raise IOError(popen.stderr.read())

	output = stdout
	return output


def parse_output(output, keep = None):
	ret={} #this is the data we will return

	## parse everything
	reg=re.compile('(\\w*)=(.*)', re.I)
	for line in output.splitlines():
		m=reg.match(line)
		if m:
			if keep is not None:
				#see if we need to filter out data
				k=m.group(1)
				if k in keep:
					ret[k]=m.group(2)#.split(os.pathsep)
			else:
				# take everything
				ret[m.group(1)]=m.group(2)#.split(os.pathsep)

	#see if we need to filter out data
	if keep is not None:
		pass

	return ret

def normalize_env(shellenv, keys):
	"""Given a dictionary representing a shell environment, add the variables
	from os.environ needed for the processing of .bat files; the keys are
	controlled by the keys argument.

	It also makes sure the environment values are correctly encoded.

	Note: the environment is copied"""
	normenv = {}
	if shellenv:
		if sys.platform=='win32':
			for k in shellenv.keys():
				normenv[k] = copy.deepcopy(shellenv[k]).encode('mbcs')

		for k in keys:
			if os.environ.has_key(k):
				normenv[k] = os.environ[k]

	return normenv

def get_script_env(env,script,args=None,vars=None):
	'''
	this function returns a dictionary of all the data we want to merge
	or process in some other way.
	'''
	if sys.platform=='win32':
		nenv = normalize_env(env['ENV'], ['COMSPEC'])
	else:
		nenv = normalize_env(env['ENV'], [])
	output = get_output(script,args,nenv)
	vars = parse_output(output, vars)

	return vars


def merge_script_vars(env,script,args=None,vars=None):
	'''
	This merges the data retieved from the script in to the Enviroment
	by prepending it.
	script is the name of the script, args is optional arguments to pass
	vars are var we want to retrieve, if None it will retieve everything found
	'''
	shell_env=get_script_env(env,script,args,vars)
	for k, v in shell_env.iteritems():
		env.PrependENVPath(k, v, delete_existing=1)

def generate(env):
	if 'INTEL_COMPILER' not in env or env['INTEL_COMPILER'] != True:
		return
	if env['TARGET_ARCH'] == 'x86':
		arch = 'ia32'
		arch_redist = 'ia32'
	elif env['TARGET_ARCH'] == 'x86_64' or env['TARGET_ARCH'] == 'amd64':
		arch = 'ia32_intel64'
		arch_redist = 'intel64'
	else:
		raise Exception('Unknown architecture ' + env['TARGET_ARCH'])

	if env['MSVC_VERSION'] == '9.0':
		vsrelease = 'vs2008'
	elif env['MSVC_VERSION'] == '10.0':
		vsrelease = 'vs2010'
	else:
		raise Exception('Unknown version of visual studio!')

	if 'ICPP_COMPOSER2014' in os.environ:
		icpp_path = os.environ.get('ICPP_COMPOSER2014')
	elif 'ICPP_COMPILER14' in os.environ:
		icpp_path = os.environ.get('ICPP_COMPILER14')
	elif 'ICPP_COMPOSER2013' in os.environ:
		icpp_path = os.environ.get('ICPP_COMPOSER2013')
	elif 'ICPP_COMPILER13' in os.environ:
		icpp_path = os.environ.get('ICPP_COMPILER13')
	elif 'ICPP_COMPOSER2011' in os.environ:
		icpp_path = os.environ.get('ICPP_COMPOSER2011')
	elif 'ICPP_COMPILER12' in os.environ:
		icpp_path = os.environ.get('ICPP_COMPILER12')
	else:
		raise Exception('Could not find any of the ICCPP_* environment variables!')

	merge_script_vars(env, os.path.join(icpp_path, 'bin/iclvars.bat'), arch + ' ' + vsrelease)
	env['REDIST_PATH'] = os.path.join(os.path.join(os.path.join(icpp_path, 'redist'), arch_redist), 'compiler')

def exists(env):
	if 'INTEL_COMPILER' not in env or env['INTEL_COMPILER'] != True:
		return False
	return 'ICPP_COMPOSER2011' in os.environ

