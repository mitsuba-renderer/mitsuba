This directory contains a very basic plugin for Mitsuba <-> Blender
integration. It is based on the excellent LuxBlend 2.5 code.

MtsBlend piggybacks on Blender's internal COLLADA exporter to transfer most
of the scene information into Mitsuba. Whatever gets "lost in translation"
is added back in using Mitsuba's adjustment mechanism.
