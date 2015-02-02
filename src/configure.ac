#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_COPYRIGHT([
############################################################################
#                        Copyright (c) 2002 - 2005 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
############################################################################
])dnl

AC_PREREQ(2.57)
AC_INIT([gOpenMol],[GOM_PACKAGE_VERSION])
AC_CONFIG_AUX_DIR([unix])
AC_CONFIG_SRCDIR([include/gomlib/gopenmolext.h])
AC_CONFIG_HEADER([config.h])

############################################################################
#
# Configure Posix build environment ...
#
############################################################################

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_F77
GOM_PROG_DEPCOMP([CXX])
GOM_PROG_DEPCOMP([CC])
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Keep flags specified by user or AC_PROG_CC separately.
GOM_SAVE_USER_OPTIONS

# Checks for libraries.
for gom_lib in c m util
do
    AC_CHECK_LIB(
	[$gom_lib],[main],[STD_LIBS="-l$gom_lib $STD_LIBS"],[],[$STD_LIBS])
done
GOM_PREPEND_PACKAGE_OPTIONS([STD])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([sys/procfs.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_C_VOLATILE
GOM_CHECK_SIZEOF([intmax_t])
GOM_CHECK_SIZEOF([intptr_t])
GOM_CHECK_SIZEOF([ptrdiff_t])
GOM_CHECK_SIZEOF([size_t])
GOM_CHECK_SIZEOF([void *])
GOM_CHECK_SIZEOF([int])
GOM_CHECK_SIZEOF([long int])
GOM_CHECK_SIZEOF([long long int])

GOM_DEFINE_OS_NAME

# Checks for library functions.
AC_TYPE_SIGNAL
GOM_FUNC_C99_PRINTF
GOM_FUNC_C99_SCANF
# Following functions are in ANSI C althought autoscan warns about them:
#     floor, malloc, memmove, memset, pow, realloc, sqrt, strchr, strerror,
#     strstr
AC_CHECK_FUNCS([asprintf vasprintf getcwd getrusage nearbyint])

AC_CACHE_SAVE

############################################################################
#
# Configure gOpenMol ...
#
############################################################################

dnl Help strings.
AC_ARG_ENABLE([plugins],
    [AC_HELP_STRING([--disable-plugins],
	[Disable plugin support])])
AC_ARG_ENABLE([debug],
    [AC_HELP_STRING([--enable-debug],
	[Enable some debug printing and assertions lines.
        Disabled by default.])],
    [AC_DEFINE([DEBUG], [],[Define if you want to enable debugging code.])],
    [AC_DEFINE([NDEBUG],[],[Define if you want to disable debugging code.])])
AC_ARG_WITH([mesa],
    [AC_HELP_STRING([--with-mesa],
	[First try Mesa. If not found then try OpenGL.])],
    [GL_LIB='MesaGL GL'; GLU_LIB='MesaGLU GLU'])
AC_ARG_WITH([gmake],
    [AC_HELP_STRING([--with-gmake],
	[Use GNU Make makefile extensions.])])
AC_ARG_WITH([dirs],
    [AC_HELP_STRING([--with-dirs=dirs],
	[List of existing directories. Separator is the same as for
	PATH (usually a colon).
	For each directory <dir> add -I<dir>/include to CPPFLAGS,
	-L<dir>/lib to LDFLAGS, pass -R<dir>/lib to libtool and
	search for programs in <dir>/bin.])])
AC_ARG_WITH([subdir],
    [AC_HELP_STRING([--with-subdir=DIR],
	[Install all files to PREFIX/share/DIR and EPREFIX/lib/DIR
	@<:@default=gOpenMol-]AC_PACKAGE_VERSION[@:>@])])

# Set subdirectory name.
case ${with_subdir-yes} in
yes) gomsubdir=[gOpenMol-]AC_PACKAGE_VERSION ;;
*)   gomsubdir="$with_subdir" ;;
esac
AC_SUBST([gomsubdir])

case $with_gmake in
yes) make_config=gmake     ;;
*)   make_config=non-gmake ;;
esac
AC_SUBST([make_config])

gom_enable_shared=$enable_shared
if test -f libtool; then
    AC_MSG_CHECKING([for libtool])
    AC_SUBST([LIBTOOL],['$(SHELL) $(top_builddir)/libtool'])
    AC_MSG_RESULT([./libtool])
    if $CONFIG_SHELL ./libtool --features |
	    grep 'enable shared libraries' >/dev/null; then
	enable_shared=yes
    else
	enable_shared=no
    fi
else
    AC_PROG_LIBTOOL
fi

# Variable gom_enable_shared is `yes' if user gave --enable-shared.
# Variable enable_shared is `yes' if libtool can build shared libraries.
case $gom_enable_shared-$enable_shared in
yes-yes) ARSUFFIX=la; BUILDTARGET=shared;  CCTAGS= ;;
*)       ARSUFFIX=a;  BUILDTARGET=dynamic; CCTAGS="--tag=disable-shared" ;;
esac
case $enable_shared in
no) enable_plugins=no ;;
esac
case $enable_plugins in
no) BUILDTARGET=static ;;
*)  enable_plugins=yes
    AC_DEFINE([ENABLE_EXTENSIONS],[],
	[Define if you want to enable API functions for extensions.])
    ;;
esac
AC_SUBST([ARSUFFIX])
AC_SUBST([BUILDTARGET])
AC_SUBST([CCTAGS])

AC_CACHE_SAVE

# Prepend common library directories.
gom_PATH="$PATH"
GOM_PREPEND_ITEMS([COMMON_CPPFLAGS],
    [-I],[$with_dirs],[/include], [],                [$PATH_SEPARATOR])
GOM_PREPEND_ITEMS([COMMON_LDFLAGS],
    [-L],[$with_dirs],[/lib],     [],                [$PATH_SEPARATOR])
GOM_PREPEND_ITEMS([COMMON_LTLDFLAGS],
    [-R],[$with_dirs],[/lib],     [],                [$PATH_SEPARATOR])
GOM_PREPEND_ITEMS([gom_PATH],
    [],  [$with_dirs],[/bin],     [$PATH_SEPARATOR], [$PATH_SEPARATOR])
GOM_PREPEND_PACKAGE_OPTIONS([COMMON])

GOM_PROG_TCLSH([$gom_PATH])
case $TCLSH in
:) AC_MSG_ERROR([gOpenMol cannot be compiled without Tcl.]) ;;
esac

# Checks for packages.
GOM_PATH_TCL(
    [AC_MSG_ERROR([gOpenMol cannot be compiled without Tcl.])],
    [Tcl_StringCaseMatch],[$gom_PATH])
GOM_PATH_X11
gom_ALL_LIBS="$ALL_LIBS"
gom_LIBS="$LIBS"
# Search for other graphics libraries.
if test "$no_x" = yes; then
    gom_use_graphics=no
else
    gom_use_graphics=yes
    # Search for Tk library.
    GOM_PATH_TK([gom_use_graphics=no],[Tk_Init])
    # Search for OpenGL and Mesa libraries.
    GOM_CHECK_PACKAGE([GL],[glClear],[GL/gl.h],[GL],,[gom_use_graphics=no])
    GOM_CHECK_PACKAGE([GLU],[gluDisk],[GL/glu.h],[GLU],,[gom_use_graphics=no])
    GOM_CHECK_PACKAGE([glut],[glutInit],[GL/glut.h],[GLUT],,[gom_use_graphics=no])
fi
if test "$gom_use_graphics" = yes; then
    AC_DEFINE([ENABLE_GRAPHICS],[],[Define to enable graphics sub system.])
else
    ALL_LIBS="${gom_ALL_LIBS}"
    LIBS="$gom_LIBS"
    AC_MSG_WARN([Will not build graphical user interface.])
fi

GOM_CHECK_PACKAGE([jpeg],[jpeg_start_compress],[jpeglib.h],[JPEG])
GOM_PATH_PYTHON([],[$gom_PATH])

AC_CACHE_SAVE

# Prefer directories specified by user with --with-dirs. Prepend again.
GOM_PREPEND_PACKAGE_OPTIONS([COMMON])

AC_DEFINE([GOMLIB],[],
    [Produce a library (will produce a gom_main entry point).])

AC_CONFIG_FILES(
    config.mk:config.mk.ac
    Makefile:unix/make/top.mk.ac
    init.mk:unix/make/init.mk.ac:config.mk.ac:unix/make/compile.mk:unix/make/gomapi.mk
    general/Makefile:unix/make/subdir.mk.ac
    graphics/Makefile:unix/make/subdir.mk.ac
    parser/Makefile:unix/make/subdir.mk.ac
    readwrite/Makefile:unix/make/subdir.mk.ac
    plugins/config.mk:plugins/config.mk.ac
    utility/Makefile:unix/make/utility.mk.ac
)

AC_CONFIG_COMMANDS(
    [graphics/mk],
    [# Initialize dependency files.
    ${MAKE-make} -f init.mk >/dev/null && rm -f init.mk])

AC_CONFIG_COMMANDS(
    [plugins],
    [# Create plugin directories.
    for gom_dir in `cd "$srcdir" && echo \
	plugins \
	plugins/* \
	plugins/*/* \
	plugins/*/*/*`
    do
	test -d "$srcdir/$gom_dir" && test ! -d "$gom_dir" && mkdir "$gom_dir"
    done])

# Create Makefile links.
echo dummy >conf$$
if test ! -f "$srcdir/conf$$"
then
    for gom_mk in `cd "\$srcdir" && echo \
	plugins/Makefile \
	plugins/*/Makefile`
    do
	AC_CONFIG_LINKS($gom_mk:$gom_mk)
    done
fi
rm -f conf$$

AC_OUTPUT

echo "
Summary
=======
tclsh:    $TCLSH
Tcl:      $TCL_LIBS
Tk:       $TK_LIBS
OpenGL:   $GL_LIBS
GLU:      $GLU_LIBS
glut:     $GLUT_LIBS
JPEG:     $JPEG_LIBS
Python:   $PY_LIBS
Graphics: $gom_use_graphics
Plugins:  $enable_plugins"
