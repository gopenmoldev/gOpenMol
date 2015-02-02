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

dnl # _GOM_INIT_EXEC_PREFIX_LIST()
dnl # ----------------------------
dnl # Sets gom4_exec_prefix_list to list of installation prefixes.
AC_DEFUN([_GOM_INIT_EXEC_PREFIX_LIST],
    [gom4_$0_IFS="$IFS"
    # Get parent directories of library directories in LDFLAGS.
    for gom4_$0_flag in $LDFLAGS
    do
	for gom4_$0_dir in `expr x"$gom4_$0_flag" : x'-L\(.*\)[[\\\\/]]lib$'`
	do
	    gom4_exec_prefix_list="$gom4_exec_prefix_list${gom4_exec_prefix_list+$PATH_SEPARATOR}$gom4_$0_dir"
	done
    done
    # Get parent directories of bin directories in PATH.
    IFS="$PATH_SEPARATOR"
    for gom4_$0_dir in _GOM_TCL_PATH
    do
	for gom4_$0_dir in `expr x"$gom4_$0_dir" : x'\(.*\)[[\\\\/]]bin$'`
	do
	    gom4_exec_prefix_list="$gom4_exec_prefix_list${gom4_exec_prefix_list+$PATH_SEPARATOR}$gom4_$0_dir"
	done
    done
    IFS="$gom4_$0_IFS"])

dnl # _GOM_TCLTK_FILE_LIST({Tcl | Tk}, DIR, File)
dnl # -------------------------------------------
dnl # List of Tcl/Tk directories. Highest version number comes last.
AC_DEFUN([_GOM_TCLTK_FILE_LIST],
    [$2/m4_tolower($1)[[0-9]]*$3 \
    $2/m4_tolower($1)[[0-9]][[0-9]]*$3 \
    $2/m4_tolower($1)${m4_toupper($1)_VERSION-[[0-9]][[0-9]]*}$3])

dnl # _GOM_TCLTK_FIND_FILE(VARIABLE, {Tcl | Tk}, SUBDIRS, TEST, FILE)
dnl # --------------------------------------------------------------
AC_DEFUN([_GOM_TCLTK_FIND_FILE],
    [AC_REQUIRE([_GOM_INIT_EXEC_PREFIX_LIST])dnl
    gom4_$0_IFS="$IFS"
    IFS="$PATH_SEPARATOR"
    # Get the last one from the first possible directory.
    # That way we might get one with the highest version number, and
    # user installed packages take precedence over system packages.
    for gom4_$0_dir in \
	$TCL_EXEC_PREFIX $TCL_PREFIX $TK_EXEC_PREFIX  $TK_PREFIX  \
	$gom4_exec_prefix_list
    do
	for gom4_$0_file in \
	    $gom4_$0_dir$3$5 _GOM_TCLTK_FILE_LIST([$2],[$gom4_$0_dir$3],[$5])
	do
	    test $4 "$gom4_$0_file" || continue
	    $1="$gom4_$0_file"
	    echo "$as_me:$LINENO: found $gom4_$0_file" >&AS_MESSAGE_LOG_FD
	done
	test -n "$$1" && break
    done
    IFS="$gom4_$0_IFS"])

dnl # _GOM_TCLTK_LOAD_CONFIG({Tcl | Tk})
dnl # --------------------------------------------------------
AC_DEFUN([_GOM_TCLTK_LOAD_CONFIG],
    [gom4_$0_IFS="$IFS"

    AC_ARG_WITH(m4_tolower($1),
	[AC_HELP_STRING([--with-]m4_tolower($1)[=DIR],
	    [Directory containing $1 configuration (]m4_tolower($1)[Config.sh)])],
	[m4_toupper($1)_BIN_DIR="$withval"])

    # Find $1Config.sh
    AC_CACHE_CHECK([for $1 configuration file],[gom4_cv_$1_config],
	[if test -n "[$]m4_toupper($1)_BIN_DIR"; then
	    gom4_cv_$1_config="[$]m4_toupper($1)_BIN_DIR/m4_tolower($1)Config.sh"
	else
	    _GOM_TCLTK_FIND_FILE([gom4_cv_$1_config],[$1],[/lib],[-f],
		[/m4_tolower($1)Config.sh])
	fi])

    # Load $1Config.sh
    if test -f "$gom4_cv_$1_config"; then
	m4_toupper($1)_BIN_DIR=`_GOM_DIRNAME([$gom4_cv_$1_config])`
	AC_MSG_NOTICE([reading $1 configurarions])
	. "$gom4_cv_$1_config"
    else
	AC_MSG_WARN([No configuration file found for $1.])
	AC_MSG_WARN([Consider using --with-]m4_tolower($1)[=DIR.])
	m4_toupper($1)_BIN_DIR=
    fi
])

AC_DEFUN([_GOM_TCLTK_LOAD_TCL_CONFIG],[_GOM_TCLTK_LOAD_CONFIG([Tcl])])
AC_DEFUN([_GOM_TCLTK_LOAD_TK_CONFIG], [_GOM_TCLTK_LOAD_CONFIG([Tk])])

dnl # _GOM_TCLTK_CHECK_PACKAGE({Tcl | Tk}, FUNCTION)
dnl # ----------------------------------------------
AC_DEFUN([_GOM_TCLTK_CHECK_PACKAGE],
    [AC_REQUIRE([GOM_PROG_TCLSH])

    # Find the library (.tcl) files of $1.
    AC_ARG_VAR(m4_toupper($1)[_LIBRARY],[$1 library directory])dnl
    AC_MSG_CHECKING([for $1 library directory])
    if test -z "[$]m4_toupper($1)_LIBRARY"; then
	gom4_$0_dir=m4_tolower($1)${m4_toupper($1)_VERSION}
	if test -f "[$]m4_toupper($1)_PREFIX/lib/gom4_$0_dir/tclIndex"; then
	    # This is the default location.
	    m4_toupper($1)_LIBRARY=[$]m4_toupper($1)_PREFIX/lib/$gom4_$0_dir
	elif test -f "[$]m4_toupper($1)_PREFIX/share/gom4_$0_dir/tclIndex"; then
	    m4_toupper($1)_LIBRARY=[$]m4_toupper($1)_PREFIX/share/$gom4_$0_dir
	elif echo ['puts [info tclversion]'] | $TCLSH |
		grep "${TCL_VERSION-dummy}" >/dev/null; then
	    # If we retrieve TCL_LIBRARY from tclsh, versions MUST match.
	    m4_toupper($1)_LIBRARY=`echo ['puts [info library]'] |
		$TCLSH | sed \
		-e ['s,\([\\\\/][tcltk0123456789.]*\)[^\\\\/]*$,\1,'] \
		m4_if([$1],Tk,[-e ['s,\([\\\\/]\)tcl\([^\\\\/]*\)$,\1tk\2,']])`
	    if test -z "[$]m4_toupper($1)_LIBRARY"; then
		AC_MSG_RESULT([not found])
	    fi
	else
	    AC_MSG_RESULT([not found])
	    AC_MSG_WARN([The $1 configuration file and the interpreter are not the same version.])
	    AC_MSG_WARN([Consider setting ]m4_toupper($1)[_LIBRARY and/or TCLSH.])
	fi
    fi
    if test -n "[$]m4_toupper($1)_LIBRARY"; then
       AC_MSG_RESULT([[$]m4_toupper($1)_LIBRARY])

	# Find headers.
	AC_ARG_WITH(m4_tolower($1)[include],
	    [AC_HELP_STRING([--with-]m4_tolower($1)[include=DIR],
		[Directory containing $1 header (]m4_tolower($1)[.h)])],
	    [m4_toupper($1)_INCLUDE_SPEC="-I$withval"])

	AC_CACHE_CHECK([for $1 header directive],[gom4_cv_$1_CPPFLAGS],
	    [if test -n "[$]m4_toupper($1)_INCLUDE_SPEC"; then
		# We have a Tcl >= 8.4.4 or --with-tclinclude.
		gom4_cv_$1_CPPFLAGS=[$]m4_toupper($1)_INCLUDE_SPEC
	    else
		# We have a Tcl < 8.4.
		gom4_$0_header=
		_GOM_TCLTK_FIND_FILE(
		    [gom4_$0_header],[$1],[/include],[-f],[/m4_tolower($1).h])
		test -n "$gom4_$0_header" &&
		    gom4_cv_$1_CPPFLAGS="-I"`_GOM_DIRNAME([$gom4_$0_header])`
	    fi])
	m4_toupper($1)_CPPFLAGS="$gom4_cv_$1_CPPFLAGS [$]m4_toupper($1)_CPPFLAGS"

	# Find linkage options.
	AC_MSG_CHECKING([for $1 linkage options])
	# TCL_LIB_SPEC contains ${TCL_DBGX}. We have to eval.
	eval "gom4_$0_specs=\"[$]m4_toupper($1)_LIB_SPEC \[$]m4_toupper($1)_LIBS \[$]m4_toupper($1)_LDFLAGS\""
	m4_toupper($1)_LIBS=
	m4_toupper($1)_LDFLAGS=
	for gom4_$0_flag in $gom4_$0_specs
	do
	    case $gom4_$0_flag in
	    -l* | -Wl,-bI:*)
		m4_toupper($1)_LIBS="[$]m4_toupper($1)_LIBS $gom4_$0_flag" ;;
	    *)  m4_toupper($1)_LDFLAGS="[$]m4_toupper($1)_LDFLAGS $gom4_$0_flag" ;;
	    esac
	done
	AC_MSG_RESULT([[$]m4_toupper($1)_LIBS [$]m4_toupper($1)_LDFLAGS])

	# Check that we have correct values.
	GOM_CHECK_PACKAGE(
	    [],[$2],[m4_tolower($1).h],m4_toupper($1),[],
	    [GOM_INIT_PACKAGE_OPTIONS(m4_toupper($1))
	    m4_toupper($1)_LIBRARY=])
    else
	GOM_INIT_PACKAGE_OPTIONS(m4_toupper($1))
    fi
])

dnl # GOM_PATH_TCL([ACTION-IF-NOT-FOUND],[FUNCTION],[PATH])
dnl # -----------------------------------------------------
dnl # Sets TCL_CPPFLAGS, TCL_LIBS, TCL_LDFLAGS, TCL_LTLDFLAGS, TCL_LIBRARY,
dnl #      CPPFLAGS, LIBS, LDFLAGS and LTLDFLAGS.
dnl # Defines HAVE_LIBTCL.
AC_DEFUN([GOM_PATH_TCL],
    [_GOM_SET_PATH([TCL],[$3])dnl
    AC_REQUIRE([_GOM_TCLTK_LOAD_TCL_CONFIG])dnl
    if test -n "$TCL_BIN_DIR"; then
	# Exclude common flags from TCL_LIBS.
	_GOM_EXCLUDE_FLAGS([TCL_LIBS],[$USER_LIBS])
	_GOM_TCLTK_CHECK_PACKAGE(Tcl,[$2])
    fi
    if test -z "$TCL_LIBRARY" || test -z "$TCL_BIN_DIR"; then
       :;$1
    fi])

dnl # GOM_PATH_TK([ACTION-IF-NOT-FOUND],[FUNCTION],[PATH])
dnl # -----------------------------------------------------
dnl # Sets TK_CPPFLAGS, TK_LIBS, TK_LDFLAGS, TK_LTLDFLAGS, TK_LIBRARY,
dnl #      CPPFLAGS, LIBS, LDFLAGS and LTLDFLAGS.
dnl # Defines HAVE_LIBTCL.
AC_DEFUN([GOM_PATH_TK],
    [_GOM_SET_PATH([TCL],[$3])dnl
    AC_REQUIRE([GOM_PATH_TCL])dnl
    AC_REQUIRE([GOM_PATH_X11])dnl
    AC_REQUIRE([_GOM_TCLTK_LOAD_TK_CONFIG])dnl
    if test -n "$TK_BIN_DIR"; then
	# Exclude common flags from TK_LIBS and TK_LDFLAGS (not splitted yet).
	_GOM_EXCLUDE_FLAGS(
	    [TK_LIBS],
	    [$USER_LIBS $TCL_LIBS $X11_LIBS $TCL_LDFLAGS $X11_LDFLAGS])
	_GOM_TCLTK_CHECK_PACKAGE(Tk,[$2])
    fi
    if test -z "$TK_LIBRARY" || test -z "$TK_BIN_DIR"; then
       :;$1
    fi])

dnl # GOM_PROG_TCLSH([PATH])
dnl # ----------------------
dnl # Search for working Tcl interpreter.
dnl # Sets: TCLSH
AC_DEFUN([GOM_PROG_TCLSH],
    [_GOM_SET_PATH([TCL],[$1])dnl
    AC_REQUIRE([_GOM_TCLTK_LOAD_TCL_CONFIG])dnl
    AC_ARG_VAR([TCLSH],[Tcl interpreter command])dnl
    # Get working Tcl interpreter.
    # Try to get from the same Tcl installation.
    _GOM_CHECK_WORKING_PROG_VER([TCLSH],[tclsh],[],[.tcl],
	[puts works],[works],[:],
	[$TCL_EXEC_PREFIX/bin${PATH_SEPARATOR}$TCL_BIN_DIR/../bin${PATH_SEPARATOR}_GOM_TCL_PATH],
	[${TCL_VERSION-[[0-9][0-9]]}])])
