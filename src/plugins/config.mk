# -*- mode: makefile -*-
# Directory entries.
prefix           = /usr/local
exec_prefix      = ${prefix}
datadir          = ${prefix}/share
includedir       = ${prefix}/include
libdir           = ${exec_prefix}/lib
gomsubdir        = gOpenMol-3.00
gomroot          = ${libdir}/${gomsubdir}
gomdataroot      = ${datadir}/${gomsubdir}
gombindir        = ${gomroot}/bin
gomincdir        = ${gomsrcdir}/include/gomlib
#INS gomincdir        = ${includedir}/${gomsubdir}
gomlibdir        = ${gomroot}/bin
gomsrcdir        = /Users/shkljo14/src/gOpenMol-3.00-osx/src
gombltdir        = /Users/shkljo14/src/gOpenMol-3.00-osx/src
#INS gomsrcdir        = ${gomroot}/src
#INS gombltdir        = ${gomroot}/src
gompluginsdir    = ${gomroot}/plugins
gompluginstcldir = ${gomdataroot}/data/pkgtcl/plugins

# File suffixies.
OBJSUFFIX        = lo
ARSUFFIX         = la
# Cached user values.
CPPFLAGS         = 
CFLAGS           = -g -O2
CXXFLAGS         = -g -O2
FFLAGS           = 
LDFLAGS          = 
LTLDFLAGS        = 
DEFS             = 
LIBS             = 
# Some useful C/C++ compiler entries found by gOpenMol configure.
# C library entries.
TCL_LIBS         =  -ltcl8.6
X11_LIBS         = -lXmu -lXi  -lSM -lICE -lX11 
TK_LIBS          =  -ltk8.6
GL_LIBS          = -lGL 
GLU_LIBS         = -lGLU 
GLUT_LIBS        = -lglut 
JPEG_LIBS        = -ljpeg 
PY_LIBS          = -lpython2.7 
# C compiler entries for libraries.
TCL_CPPFLAGS     = -I/usr/local/include 
TCL_LDFLAGS      =  -L/usr/local/lib
TCL_LTLDFLAGS    =  -R/usr/local/lib
X11_CPPFLAGS     = -I/opt/X11/include/ 
X11_CFLAGS       = 
X11_LDFLAGS      =  -L/opt/X11/lib/
TK_CPPFLAGS      = -I/usr/local/include 
TK_LDFLAGS       =  -L/usr/local/lib
TK_LTLDFLAGS     =  -R/usr/local/lib
PY_CPPFLAGS      = -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
PY_LDFLAGS       = -L/System/Library/Frameworks/Python.framework/Versions/2.7/lib -L/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config
PY_LTLDFLAGS     =  -R/System/Library/Frameworks/Python.framework/Versions/2.7/lib -R/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config
COMMON_CPPFLAGS  = 
COMMON_LDFLAGS   = 
COMMON_LTLDFLAGS = 
# Values used by gOpenMol.
# Not quite... Don't show private gOpenMol headers to plugins.
GOM_CPPFLAGS     =  -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7     -I/usr/local/include  -I/opt/X11/include/  -I/usr/local/include     -I${gomincdir} -I${gombltdir}/include/gomlib
GOM_CFLAGS       =            -g -O2
GOM_FFLAGS       =            
GOM_LDFLAGS      =  -L/System/Library/Frameworks/Python.framework/Versions/2.7/lib -L/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config      -L/usr/local/lib  -L/opt/X11/lib/  -L/usr/local/lib   
GOM_LTLDFLAGS    =   -R/System/Library/Frameworks/Python.framework/Versions/2.7/lib -R/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config      -R/usr/local/lib   -R/usr/local/lib   
GOM_DEFS         = -DHAVE_CONFIG_H -DOS_NAME='"Darwin"'
GOM_LIBS         =  -lpython2.7  -ljpeg  -lglut  -lGLU  -lGL   -ltk8.6 -lXmu -lXi  -lSM -lICE -lX11   -ltcl8.6  -lutil -lm -lc  
# Default values for plugin.
PLUG_CPPFLAGS    = ${GOM_CPPFLAGS}
PLUG_CFLAGS      = ${CFLAGS}
PLUG_CXXFLAGS    = ${CXXFLAGS}
PLUG_LDFLAGS     = ${GOM_LDFLAGS}
PLUG_LTLDFLAGS   = ${GOM_LTLDFLAGS}
PLUG_DEFS        = ${DEFS}
PLUG_LIBS        = ${LIBS}
# C/C++ compiler entries.
CPP              = gcc -E
CC               = gcc
CXX              = g++
CCDEPMODE        = depmode=gcc3
CCDEPCOMP        = $(CCDEPMODE) $(DEPCOMP)
CXXDEPMODE       = depmode=gcc3
CXXDEPCOMP       = $(CXXDEPMODE) $(DEPCOMP)
GOM_COMPILE_CC   = ${LIBTOOL} --tag=CC  --tag=disable-static --mode=compile ${CC}
GOM_COMPILE_CXX  = ${LIBTOOL} --tag=CXX --tag=disable-static --mode=compile ${CXX}
GOM_LINK_CC      = ${LIBTOOL} --tag=CC  --tag=disable-static --mode=link ${CC}
GOM_LINK_CXX     = ${LIBTOOL} --tag=CXX --tag=disable-static --mode=link ${CXX}
# Tools.
SHELL            = /bin/sh
LN_S             = ln -s
DEPCOMP          = ${SHELL} ${gomsrcdir}/unix/depcomp
LIBTOOL          = ${SHELL} ${gombltdir}/libtool
MKDIR_P          = ${SHELL} ${gomsrcdir}/unix/install-sh -d -m 755
RM               = rm -f
RM_R             = ${RM} -r
GOM_REMOVE       = ${LIBTOOL} --mode=clean ${RM}
GOM_COMPILEFLAGS =
GOM_LINKFLAGS    = -avoid-version -rpath ${gompluginsdir}
# Install library using following commands:
#  ${GOM_PREINSTALL} ${gompluginsdir}
#  ${GOM_INSTALL} libexample.${ARSUFFIX} ${gompluginsdir}
#  ${GOM_POSTINSTALL} ${gompluginsdir}
GOM_PREINSTALL   = :
GOM_INSTALL      = ${LIBTOOL} --mode=install ${SHELL} ${gomsrcdir}/unix/install-sh -c -m 555
GOM_POSTINSTALL  = ${LIBTOOL} --mode=finish
# To install data files such as Tcl scripts.
GOM_INSTALL_DATA = ${SHELL} ${gomsrcdir}/unix/install-sh -c -m 644
