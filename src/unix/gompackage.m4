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

m4_define(_GOM_OPTION_VARIABLES,
[CPPFLAGS CFLAGS CXXFLAGS FFLAGS DEFS LIBS LDFLAGS LTLDFLAGS])

dnl # GOM_INIT_PACKAGE_OPTIONS(PREFIX)
dnl # ---------------------------------
dnl # Set package options variable to empty values.
AC_DEFUN([GOM_INIT_PACKAGE_OPTIONS],
[m4_bpatsubst(_GOM_OPTION_VARIABLES,[\w+],[$1_\&=;])])

dnl # GOM_PREPEND_PACKAGE_OPTIONS(PREFIX)
dnl # -------------------------------------------------------------
dnl # Prepend the value of package variable to a global variable.
dnl # _GOM_UNEXPAND_VARIABLES will later unexpand added options.
dnl # In Makefile.in use:
dnl # GL_LIBS = @GL_LIBS@  # could be `-lGLU -lGL'
dnl # LIBS    = @ALL_LIBS@ # could be `${GL_LIBS} -lm -lc'
AC_DEFUN([GOM_PREPEND_PACKAGE_OPTIONS],
[AC_REQUIRE([_GOM_CHECK_PACKAGE_OPTIONS_DELAY])dnl
m4_bpatsubst(_GOM_OPTION_VARIABLES,[\w+],
    [\&="$$1_\& $\&"
     ALL_\&='${$1_\&}'" $ALL_\&"
     AC_SUBST([\&])dnl
     AC_SUBST([$1_\&])dnl
     AC_SUBST([ALL_\&])])])dnl

dnl # GOM_SAVE_USER_OPTIONS()
dnl # -----------------------
dnl # Sets the value of USER_<VARIABLE> to the value of <VARIABLE>.
dnl # In Makefile.in use:
dnl # CFLAGS    = @USER_CFLAGS@  # could be `-g -O2'
dnl # MY_CFLAGS = @ALL_CPPFLAGS@ # could be `-fPIC ${CPPFLAGS}'
AC_DEFUN([GOM_SAVE_USER_OPTIONS],
[AC_REQUIRE([_GOM_CHECK_PACKAGE_OPTIONS_DELAY])dnl
m4_bpatsubst(_GOM_OPTION_VARIABLES,[\w+],
[USER_\&=$\&
AC_SUBST([USER_\&])])])dnl

dnl # _GOM_CHECK_PACKAGE_OPTIONS_DELAY()
dnl # -------------------------------
AC_DEFUN([_GOM_CHECK_PACKAGE_OPTIONS_DELAY],
[AC_CONFIG_COMMANDS_PRE([_GOM_CHECK_PACKAGE_OPTIONS])])

dnl # _GOM_CHECK_PACKAGE_OPTIONS()
dnl # ----------------------------
# At the moment, ALL_ variables contain references to package variables.
# Append options which are in global variable but not in package
# variables literally to ALL_ variables.
# Finally, append ${USER_<VAR>} to ALL_<VAR> variables.
AC_DEFUN([_GOM_CHECK_PACKAGE_OPTIONS],
[for gom4_$0_var in _GOM_OPTION_VARIABLES; do
    eval "gom4_$0_glb=\${${gom4_$0_var}}"
    eval "gom4_$0_usr=\${USER_${gom4_$0_var}}"
    eval "gom4_$0_all=\${ALL_${gom4_$0_var}}"
    eval "gom4_$0_pkgs=\"${gom4_$0_all}\""
    for gom4_$0_opt in $gom4_$0_glb; do
	case " $gom4_$0_pkgs $gom4_$0_usr " in
	*" $gom4_$0_opt "*) ;;
	*)
	    gom4_$0_all="$gom4_$0_all $gom4_$0_opt"
	    gom4_$0_pkgs="$gom4_$0_pkgs $gom4_$0_opt"
	    ;;
	esac
    done
    eval "ALL_${gom4_$0_var}=\"\${gom4_$0_all} \\\${${gom4_$0_var}}\""
done])

dnl # GOM_PREPEND_ITEMS(VARIABLE, PREFIX, LIST, SUFFIX, [SEPARATOR], [IFS])
dnl # ---------------------------------------------------------------------
AC_DEFUN([GOM_PREPEND_ITEMS],
[gom4_$0_IFS="$IFS"
m4_ifval([$6],[IFS="$6"])
set x
for gom4_$0_item in $3
do
    case $gom4_$0_item in "") continue ;; esac
    set "[$]@" "$2${gom4_$0_item}$4"
done
set "[$]@" $$1
shift # x
IFS="m4_default([$5],[ ])"
$1="[$]*"
IFS="$gom4_$0_IFS"])

dnl # _GOM_EXCLUDE_OPTIONS(VARIABLE, OPTIONS)
dnl # ---------------------------------------
dnl # Remove options OPTIONS from the value of VARIABLE.
dnl # NOTE: two parameter options are not supported.
AC_DEFUN([_GOM_EXCLUDE_FLAGS],
[set x $$1; shift
for gom4_$0_opt
do
    case " $2 " in
    *" $gom4_$0_flag "*) ;;
    *) set x "[$]@" "$gom4_$0_opt"; shift ;;
    esac
    shift # $gom4_$0_opt
done
$1="[$]*"])

dnl # GOM_CHECK_PACKAGE([LIBRARY-LIST], FUNCTION, HEADERS, VARIABLE-PREFIX,
dnl #                   [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl # ---------------------------------------------------------------------
dnl # If one the libraries in <LIBRARY-LIST> and all the headers are
dnl # found, set <VARIABLE-PREFIX>_LIBS to
dnl # `-l<library> $<VARIABLE_PREFIX>_LIBS' and prepend that value to
dnl # LIBS and execute ACTION-IF-FOUND (default is GOM_HAVE_PACKAGE).
dnl # Otherwise execute ACTION-IF-NOT-FOUND.
AC_DEFUN([GOM_CHECK_PACKAGE],
[AC_REQUIRE([_GOM_CHECK_PACKAGE_SOURCE])dnl
if { set x "m4_ifval([$1],[${$4_LIB-$1}])" '$2' '$3' '$4'; shift
     . cftest.gom/chkpkg
}; then
    GOM_PREPEND_PACKAGE_OPTIONS([$4])
    m4_default([$5],[GOM_HAVE_PACKAGE([$4],m4_tolower([$4]))])
else
    GOM_INIT_PACKAGE_OPTIONS([$4])
    $6
fi])
AC_DEFUN([_GOM_CHECK_PACKAGE_SOURCE],
[AC_REQUIRE([_GOM_TMPDIR])dnl
# This is a shell function emulation:
# set LIBRARY-LIST FUNCTION HEADERS VARIABLE-PREFIX
# . cftest.gom/chkpkg
test -d cftest.gom || mkdir cftest.gom
cat >cftest.gom/chkpkg <<'_GOM_EOF'
gom4_$0_liblist=[$]1
gom4_$0_func=[$]2
gom4_$0_headers=[$]3
gom4_$0_prefix=[$]4
eval "
    gom4_$0_CPPFLAGS=\${${gom4_$0_prefix}_CPPFLAGS}
    gom4_$0_LIBS=\${${gom4_$0_prefix}_LIBS}
    gom4_$0_LDFLAGS=\${${gom4_$0_prefix}_LDFLAGS}
    gom4_$0_LTLDFLAGS=\${${gom4_$0_prefix}_LTLDFLAGS}
"
# Check for headers.
gom4_$0_found=yes
gom4_$0_CPPFLAGS_save=$CPPFLAGS
CPPFLAGS="$gom4_$0_CPPFLAGS $CPPFLAGS"
for gom4_$0_header in $gom4_$0_headers
do
    AC_CHECK_HEADER([$gom4_$0_header],[:],[gom4_$0_found=no])
done
CPPFLAGS=$gom4_$0_CPPFLAGS_save
if test "$gom4_$0_found" = yes; then
    gom4_$0_found=no
    # Check for a library.
    gom4_$0_LIBS_save="$LIBS"
    gom4_$0_LDFLAGS_save=$LDFLAGS
    LIBS="$gom4_$0_LIBS $LIBS"
    LDFLAGS="$gom4_$0_LDFLAGS $LDFLAGS"
    if test -n "$gom4_$0_liblist"; then
	for gom4_$0_lib in $gom4_$0_liblist; do
	    AC_CHECK_LIB([$gom4_$0_lib],[$gom4_$0_func],
		[eval "${gom4_$0_prefix}_LIBS=\"-l\$gom4_$0_lib \$${gom4_$0_prefix}_LIBS\""
		 gom4_$0_found=yes;break])
	done
    else
	AC_CACHE_CHECK(
	    [for $gom4_$0_func in $gom4_$0_LIBS],
	    [gom4_cv_lib_$4_ok],
	    [AC_TRY_LINK_FUNC([$gom4_$0_func],
		[gom4_cv_lib_$4_ok=yes],[gom4_cv_lib_$4_ok=no])])
	gom4_$0_found=$gom4_cv_lib_$4_ok
    fi
    LIBS=$gom4_$0_LIBS_save
    LDFLAGS=$gom4_$0_LDFLAGS_save
fi
if test "$gom4_$0_found" = yes; then
    # Convert library directories to -R options for libtool.
    for gom4_$0_flag in $gom4_$0_LDFLAGS $gom4_$0_LIBS
    do
	case $gom4_$0_flag in
	-L*) gom4_$0_LTLDFLAGS="$gom4_$0_LTLDFLAGS -R"`expr x"$gom4_$0_flag" : 'x-L\(.*\)'` ;;
	esac
    done
    eval "${gom4_$0_prefix}_LTLDFLAGS=\$gom4_$0_LTLDFLAGS"
fi
test "$gom4_$0_found" = yes
_GOM_EOF
])dnl

dnl # GOM_HAVE_PACKAGE(LIBRARY-MACRO, LIBRARY-NAME)
dnl # -------------------------------------------
AC_DEFUN([GOM_HAVE_PACKAGE],
[AC_DEFINE([HAVE_LIB$1], [1],
    [Define to 1 if you have the `$2' library (-l$2) and
    corresponding header files.])])dnl `

dnl # GOM_PATH_X11()
dnl # --------------
dnl # Set variables X11_CPPFLAGS, X11_CFLAGS, X11_LIBS and X11_LDFLAGS.
dnl # Prepends these values to CPPFLAGS, CFLAGS, LIBS and LDFLAGS.
AC_DEFUN([GOM_PATH_X11],
[AC_REQUIRE([AC_PATH_XTRA])dnl
if test "$no_x" != yes; then
    # Variable names used by AC_PATH_XTRA are quite odd.
    # X_CFLAGS: C preprocessor and C flags.
    for gom4_$0_flag in $X_CFLAGS
    do
        case "$gom4_$0_flag" in
        -I* | -D*) X11_CPPFLAGS="$gom4_$0_flag $X11_CPPFLAGS" ;;
	*)         X11_CFLAGS="$gom4_$0_flag $X11_CFLAGS" ;;
	esac
    done
    # X_PRE_LIBS, X_EXTRA_LIBS: libraries
    # X_LIBS: linker flags
    X11_LIBS="$X_PRE_LIBS -lX11 $X_EXTRA_LIBS"
    X11_LDFLAGS="$X_LIBS"
    for gom4_$0_lib in Xi Xmu
    do
	AC_CHECK_LIB(
	    [$gom4_$0_lib],[main],[X11_LIBS="-l$gom4_$0_lib $X11_LIBS"],
            [],[$X11_LDFLAGS $X11_LIBS])
    done
    GOM_PREPEND_PACKAGE_OPTIONS([X11])
fi])
