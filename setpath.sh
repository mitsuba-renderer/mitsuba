#
# This script adds Mitsuba to the current path.
# It works with both Bash and Zsh.
#
# NOTE: this script must be sourced and not run, i.e.
#    . setpath.sh        for Bash
#    source setpath.sh   for Zsh or Bash
#

if [ "$BASH_VERSION" ]; then
	MITSUBA_DIR=$(dirname "$BASH_SOURCE")
	export MITSUBA_DIR=$(builtin cd "$MITSUBA_DIR"; builtin pwd)
elif [ "$ZSH_VERSION" ]; then
	export MITSUBA_DIR=$(dirname "$0")

	# Zsh autocomplete for mitsuba, mtsutil, and mtssrv
	mitsuba_plugins=$(ls -1 "$MITSUBA_DIR/dist/plugins" | grep -oE '\w+.so' | sed 's/.so$//')
	compdef "_arguments '-c[connect to host(s)]:host:_hosts' '-s[connect to list of hosts in a file]:hostfile:_files' '-o[output file]:out:_files' '*:scene:_files -g \*.\(xml\|XML\)'" mitsuba
	compdef "_arguments '-c[connect to host(s)]:host:_hosts' '-s[connect to list of hosts in a file]:hostfile:_files' '1:plugins:($mitsuba_plugins)' '*:utilargs:_files'" mtsutil
	compdef "_arguments '-c[connect to host(s)]:host:_hosts' '-s[connect to list of hosts in a file]:hostfile:_files'" mtssrv
	unset mitsuba_plugins
fi

if [[ "$(uname)" == 'Darwin' ]]; then
	export PATH="$MITSUBA_DIR/Mitsuba.app/Contents/MacOS:$PATH"
	export PYTHONPATH="$MITSUBA_DIR/Mitsuba.app/python/2.7:$PYTHONPATH"
else
	export LD_LIBRARY_PATH="$MITSUBA_DIR/dist:$LD_LIBRARY_PATH"
	export PATH="$MITSUBA_DIR/dist:$PATH"

	# Add Mitsuba to PYTHONPATH if there is only one version of Python
	mitsuba_python=$(ls -1 "$MITSUBA_DIR/dist/python")
	if [[ "$(echo $mitsuba_python | wc -l)" == "1" ]]; then
		export PYTHONPATH="$MITSUBA_DIR/dist:$MITSUBA_DIR/dist/python/$mitsuba_python:$PYTHONPATH"
	fi
	unset mitsuba_python

	# Generate core dumps if something goes wrong
	ulimit -c 1000000000
fi
