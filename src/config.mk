# -*- mode: makefile -*-
# Directories.
SUBDIRS          = general graphics parser readwrite
# Source directories.
incdir           = $(top_srcdir)/include
rootdir          = $(top_srcdir)/..
makedir          = $(top_srcdir)/unix/make
make_config      = non-gmake
# Target directories.
prefix           = /usr/local
exec_prefix      = ${prefix}
bindir           = ${exec_prefix}/bin
datadir          = ${prefix}/share
includedir       = ${prefix}/include
libdir           = ${exec_prefix}/lib
gomsubdir        = gOpenMol-3.00
gomroot          = $(libdir)/$(gomsubdir)
gomdataroot      = $(datadir)/$(gomsubdir)
gombindir        = $(gomroot)/bin
gomincdir        = $(includedir)/$(gomsubdir)
gomlibdir        = $(gomroot)/bin
# Tcl/Tk entries.
TCL_LIBRARY      = /Library/Frameworks/Tcl.framework/Versions/8.6/Resources/
TK_LIBRARY       = /Library/Frameworks/Tcl.framework/Versions/8.6/Resources/
# Cached user values.
CPPFLAGS         = 
CFLAGS           = -g -O2
CXXFLAGS         = -g -O2
FFLAGS           = 
LDFLAGS          = 
LTLDFLAGS        = 
DEFS             = 
LIBS             = 
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
# C compiler entries.
CPP              = gcc -E
CC               = gcc
CXX              = g++
F77              = 
CCTAGS           = --tag=disable-shared
CCDEPMODE        = depmode=gcc3
CCDEPCOMP        = $(CCDEPMODE) $(DEPCOMP)
GOM_CPPFLAGS     =  -I/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7     -I/usr/local/include  -I/opt/X11/include/  -I/usr/local/include     -I$(incdir) -Iinclude -I.
GOM_CFLAGS       = ${COMMON_CFLAGS} ${PY_CFLAGS} ${JPEG_CFLAGS} ${GLUT_CFLAGS} ${GLU_CFLAGS} ${GL_CFLAGS} ${TK_CFLAGS} ${X11_CFLAGS} ${TCL_CFLAGS} ${COMMON_CFLAGS} ${STD_CFLAGS}  ${CFLAGS}
GOM_FFLAGS       = ${COMMON_FFLAGS} ${PY_FFLAGS} ${JPEG_FFLAGS} ${GLUT_FFLAGS} ${GLU_FFLAGS} ${GL_FFLAGS} ${TK_FFLAGS} ${X11_FFLAGS} ${TCL_FFLAGS} ${COMMON_FFLAGS} ${STD_FFLAGS}  ${FFLAGS}
GOM_LDFLAGS      = ${COMMON_LDFLAGS} ${PY_LDFLAGS} ${JPEG_LDFLAGS} ${GLUT_LDFLAGS} ${GLU_LDFLAGS} ${GL_LDFLAGS} ${TK_LDFLAGS} ${X11_LDFLAGS} ${TCL_LDFLAGS} ${COMMON_LDFLAGS} ${STD_LDFLAGS}  ${LDFLAGS}
GOM_LTLDFLAGS    = ${COMMON_LTLDFLAGS} ${PY_LTLDFLAGS} ${JPEG_LTLDFLAGS} ${GLUT_LTLDFLAGS} ${GLU_LTLDFLAGS} ${GL_LTLDFLAGS} ${TK_LTLDFLAGS} ${X11_LTLDFLAGS} ${TCL_LTLDFLAGS} ${COMMON_LTLDFLAGS} ${STD_LTLDFLAGS}  ${LTLDFLAGS}
GOM_DEFS         = ${COMMON_DEFS} ${PY_DEFS} ${JPEG_DEFS} ${GLUT_DEFS} ${GLU_DEFS} ${GL_DEFS} ${TK_DEFS} ${X11_DEFS} ${TCL_DEFS} ${COMMON_DEFS} ${STD_DEFS}  -DHAVE_CONFIG_H ${DEFS} -DOS_NAME='"Darwin"'
GOM_LIBS         = ${COMMON_LIBS} ${PY_LIBS} ${JPEG_LIBS} ${GLUT_LIBS} ${GLU_LIBS} ${GL_LIBS} ${TK_LIBS} ${X11_LIBS} ${TCL_LIBS} ${COMMON_LIBS} ${STD_LIBS}  ${LIBS}
STD_LIBS         = -lutil -lm -lc 
# Allow static linkage.
AR               = ar
ARFLAGS          = cru
RANLIB           = 
# File suffixies.
ARSUFFIX         = a
EXEEXT           = 
PATH_SEPARATOR   = :
# Tools.
set_show         = set x $${MAKEFLAGS}; show='echo'; \
		case " $$* " in " -s ") show=':' ;; esac; \
		case $$2 in *=*) ;; *s*) show=':' ;; esac
SHELL            = /bin/sh
LN_S             = ln -s

MAKEOLD          = touch -t 199001010000
MKDIR_P          = ${SHELL} ${top_srcdir}/unix/install-sh -d -m 755
RM               = rm -f
RM_R             = $(RM) -r
TCLSH            = tclsh8.6
DEPCOMP          = $(SHELL) $(top_srcdir)/unix/depcomp
LIBTOOL          = $(SHELL) $(top_builddir)/libtool
REMOVE           = $(LIBTOOL) --mode=clean $(RM)
COMPILE          = $(LIBTOOL) $(CCTAGS) --mode=compile $(CC) $(GOM_CPPFLAGS) $(GOM_CFLAGS) $(GOM_DEFS) -c -o
COMPILEFLAGS     =
LINK             = $(LIBTOOL) --mode=link $(CC) $(GOM_CFLAGS) -o
LINK_SO          = $(LINK)
LINK_A           = $(LINK)
LINKFLAGS        = $(GOM_LDFLAGS) $(GOM_LIBS) $(GOM_LTLDFLAGS) -rpath $(gomlibdir)
LINKFLAGS_SO     = $(LINKFLAGS) -avoid-version
LINKFLAGS_A      =
LTPREINSTALL     = :
LTINSTALL        = $(LIBTOOL) --mode=install $(INSTALL_PROGRAM)
LTPOSTINSTALL    = $(LIBTOOL) --mode=finish
INSTALL          = /usr/bin/install -c
INSTALL_PROGRAM  = ${INSTALL}
INSTALL_SCRIPT   = ${INSTALL}
INSTALL_DATA     = ${INSTALL} -m 644
# Build target:
#     dynamic     Build a single binary which exports global symbols.
#     shared      Build a shareable library if possible.
#     noshared    Wont't use libraries at all.
#     static      Build static libraries using $(AR) or libtool.
BUILDTARGET      = dynamic
