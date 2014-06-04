#
# This script adds Mitsuba to the current path.
# It works with both Bash and Zsh.
#
# NOTE: this script must be sourced and not run, i.e.
#    . setpath.sh        for Bash
#    source setpath.sh   for Zsh or Bash
#

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "The setpath.sh script must be sourced, not executed. In other words, run\n"
	echo "$ source setpath.sh\n"
	echo "If you wish to use the Mitsuba Python bindings, you should also specify"
	echo "your Python version /before/ sourcing setpath.sh, e.g.\n"
	echo "$ export MITSUBA_PYVER=3.3"
	echo "$ source setpath.sh"
	exit 0
fi

if [ "$BASH_VERSION" ]; then
	MITSUBA_DIR=$(dirname "$BASH_SOURCE")
	export MITSUBA_DIR=$(builtin cd "$MITSUBA_DIR"; builtin pwd)
elif [ "$ZSH_VERSION" ]; then
	export MITSUBA_DIR=$(dirname "$0:A")
fi

if [ "$MITSUBA_PYVER" ]; then
	pyver=$MITSUBA_PYVER
else
	pyver=`python --version 2>&1 | grep -oE '([[:digit:]].[[:digit:]])'`
fi

if [[ "$(uname)" == 'Darwin' ]]; then
	export PYTHONPATH="$MITSUBA_DIR/Mitsuba.app/python/$pyver:$PYTHONPATH"
	mitsuba_plugin_dir="$MITSUBA_DIR/Mitsuba.app/plugins"
else
	export PYTHONPATH="$MITSUBA_DIR/dist/python:$MITSUBA_DIR/dist/python/$pyver:$PYTHONPATH"
	mitsuba_plugin_dir="$MITSUBA_DIR/dist/plugins"
fi
unset pyver

if [[ ! -z "$ZSH_VERSION" && -d "$mitsuba_plugin_dir" ]]; then
	# Zsh autocomplete for mitsuba, mtsutil, and mtssrv
	mitsuba_plugins=$(ls -1 "$mitsuba_plugin_dir" | grep -oE '\w+(\.so|\.dylib)' | sed 's/.so$//;s/.dylib$//')
	compdef "_arguments '-c[connect to host(s)]:host:_hosts' '-s[connect to list of hosts in a file]:hostfile:_files' '-o[output file]:out:_files' '*:scene:_files -g \*.\(xml\|XML\)'" mitsuba
	compdef "_arguments '-c[connect to host(s)]:host:_hosts' '-s[connect to list of hosts in a file]:hostfile:_files' '1:plugins:($mitsuba_plugins)' '*:utilargs:_files'" mtsutil
	compdef "_arguments '-c[connect to host(s)]:host:_hosts' '-s[connect to list of hosts in a file]:hostfile:_files'" mtssrv
	unset mitsuba_plugins
fi

unset mitsuba_plugin_dir

if [[ "$(uname)" == 'Darwin' ]]; then
	export PATH="$MITSUBA_DIR/Mitsuba.app/Contents/MacOS:$PATH"
else
	export LD_LIBRARY_PATH="$MITSUBA_DIR/dist:$LD_LIBRARY_PATH"
	export PATH="$MITSUBA_DIR/dist:$PATH"
fi
