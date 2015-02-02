# TODO: Add plugin name in lower case.
# Plugin must have a function called Pluginname_Init (plugin name part
# capitalized).
PLUGINNAME = filters

# TODO: Subdirectory name. Usually there is no need to change this.
DIRNAME    = ${PLUGINNAME}

# TODO: Change this if there is more than just one source file.
OBJECTS    = ${PLUGINNAME}.${OBJSUFFIX}

# Always include this.
# TODO: If you want to compile outside gOpenMol directories
#       uncomment the first line, comment the second line out and
#       set environment variable GOM_ROOT to /path/to/gOpenMol/directory.
#include ${GOM_ROOT}/src/plugins/config.mk
include ../config.mk

# The rest are only needed if you don't want to write your own make rules.

# TODO: Select the source code language (used during linkage only).
#LANG       = CC
LANG       = CXX

# TODO: Set compiler flags.
PLUG_LIBS  = ${TCL_LIBS} ${LIBS}

include ${gomsrcdir}/plugins/rules.mk
