dnl -*- mode: autoconf -*-
AC_COPYRIGHT([
############################################################################
#                       Copyright (c) 2003 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
############################################################################
])dnl

dnl # _GOM_HAVE_PYTHON_SYSCONFIG()
dnl # ----------------------------
dnl # Sets gom4_python_have_sysconfig to `yes' or `no' depending on
dnl # whether python has module distutils.sysconfig.
AC_DEFUN([_GOM_HAVE_PYTHON_SYSCONFIG],
[AC_REQUIRE([GOM_PROG_PYTHON])dnl
AC_CACHE_CHECK([for python module distutils.sysconfig],
    [gom4_cv_python_have_sysconfig],
    [test : != "$PYTHON" &&
	$PYTHON -c 'import distutils.sysconfig' 2>/dev/null &&
	gom4_cv_python_have_sysconfig=yes ||
	gom4_cv_python_have_sysconfig=no])])

dnl # _GOM_PYTHON_SYSCONFIG(VARIABLE, KEYLIST)
dnl # ----------------------------------------
dnl # Set VARIABLE to the values of Python configure variables in KEYLIST.
AC_DEFUN([_GOM_PYTHON_SYSCONFIG],
[AC_REQUIRE([_GOM_HAVE_PYTHON_SYSCONFIG])dnl
case $gom4_cv_python_have_sysconfig in yes)
    AC_CACHE_VAL([gom4_cv_$1],
	[if test "${$1+set}" = set; then
	    gom4_cv_$1=$$1
	else
	    # Generade summed list of <get_config_var>s.
	    gom4_$0_gets=`echo $2 | \
		sed 's/[[^ ]][[^ ]]*/get_config_var("&")/g;s/ /+" "+/g'`
	    # Get the value.
	    gom4_cv_$1=`$PYTHON -c \
"from distutils.sysconfig import get_config_var
print $gom4_$0_gets"`
	fi])
    $1=$gom4_cv_$1 ;;
esac])

dnl # GOM_PROG_PYTHON([PATH])
dnl # --------------------------------------
dnl # Search for working Python interpreter.
dnl # Sets: PYTHON
AC_DEFUN([GOM_PROG_PYTHON],
[_GOM_SET_PATH([PYTHON],[$2])dnl
AC_ARG_VAR([PYTHON],[Python interpreter command])dnl
# Get a working Python interpreter.
_GOM_CHECK_WORKING_PROG_VER([PYTHON],[python],
    [],[.py],[print "works"],[works],[:],[_GOM_PYTHON_PATH])])

dnl # GOM_PATH_PYTHON([ACTION-IF-NOT-FOUND],[PATH])
dnl # -----------------
dnl # Sets PY_CPPFLAGS, PY_LIBS, PY_LDFLAGS, PY_LTLDFLAGS,
dnl #      CPPFLAGS, LIBS, LDFLAGS and LTLDFLAGS.
dnl # Defines HAVE_LIBPYTHON.
AC_DEFUN([GOM_PATH_PYTHON],
[_GOM_SET_PATH([PYTHON],[$2])dnl
AC_REQUIRE([GOM_PROG_PYTHON])dnl
case $PYTHON in
:) $1 ;;
*)
    # Get flags.
    AC_MSG_CHECKING([for Python])
    _GOM_PYTHON_SYSCONFIG([PY_INC_DIR],    [INCLUDEPY])
    _GOM_PYTHON_SYSCONFIG([PY_LIB_DIR_SH], [LIBDIR])
    _GOM_PYTHON_SYSCONFIG([PY_LIB_DIR_ST], [LIBPL])
    _GOM_PYTHON_SYSCONFIG([PY_LIBS],       [LIBS SYSLIBS])
    _GOM_PYTHON_SYSCONFIG([PY_VERSION],    [VERSION])
    _GOM_EXCLUDE_FLAGS([PY_LIBS],[$USER_LIBS])
    GOM_PREPEND_ITEMS([PY_CPPFLAGS],[-I],[$PY_INC_DIR])
    GOM_PREPEND_ITEMS([PY_LDFLAGS],[-L],[$PY_LIB_DIR_SH $PY_LIB_DIR_ST])
    AC_MSG_RESULT([version $PY_VERSION, libraries $PY_LDFLAGS $PY_LIBS, headers $PY_CPPFLAGS])
    # Check package.
    GOM_CHECK_PACKAGE(
	[python${PY_VERSION}],[Py_Initialize],[Python.h],[PY],
	[GOM_HAVE_PACKAGE([PYTHON],[python])],[$1])
    ;;
esac])
