dnl -*- mode: autoconf -*-
AC_COPYRIGHT([
############################################################################
#                           Copyright (c) 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
############################################################################
])dnl

dnl # GOM_CHECK_SIZEOF(TYPE)
dnl # ----------------------
dnl # Check size of a type using shell function emulation.
AC_DEFUN([GOM_CHECK_SIZEOF],
[AC_REQUIRE([_GOM_CHECK_SIZEOF_SOURCE])dnl
AH_TEMPLATE(
	[SIZEOF_]_GOM_CPP_VAR([$1]),
	[The size of a `$1', as computed by sizeof.])dnl
sed \
    -e 's/gom dummy sizeof \*/$1/g' \
    -e 's/gom_dummy_sizeof_p/_GOM_CPP_VAR([$1])/g' \
    -e 's/GOM_DUMMY_SIZEOF_P/_GOM_CPP_VAR([$1])/g' \
    cftest.gom/sizeof >cftest.gom/sizeof.$$
. cftest.gom/sizeof.$$])
AC_DEFUN([_GOM_CHECK_SIZEOF_SOURCE],
[AC_REQUIRE([_GOM_TMPDIR])dnl
cat >cftest.gom/sizeof <<'_GOM_EOF'
AC_CHECK_SIZEOF([gom dummy sizeof *])
_GOM_EOF
])

dnl # _GOM_TMPDIR()
dnl # -------------
AC_DEFUN([_GOM_TMPDIR],
[test -d cftest.gom || mkdir cftest.gom
AC_CONFIG_COMMANDS_POST([_GOM_REMOVE_TMPDIR])])
AC_DEFUN([_GOM_REMOVE_TMPDIR],
[rm -fr cftest.gom])

dnl # _GOM_CPP_VAR(STRING)
dnl # ---------------------
AC_DEFUN([_GOM_CPP_VAR],
[m4_translit(
    [$1],
    [ *abcdefghijklmnopqrstuvwxyz],
    [_PABCDEFGHIJKLMNOPQRSTUVWXYZ])])

dnl # GOM_HAVE_VA_COPY()
dnl # ------------------
dnl # C99 does not enforce va_copy to be a macro.
dnl # On the other hand, it can be a macro, thus AC_CHECK_FUNCS will
dnl # not work.
AC_DEFUN([GOM_HAVE_VA_COPY],
[AC_CACHE_CHECK([for va_copy],[gom4_cv_$0],
    [gom4_cv_$0=no
     for gom4_$0_copy in va_copy __va_copy; do
	AC_TRY_LINK(
	    [],[va_arg ap1;va_arg ap2;$gom4_$0_copy(ap1,ap2);],
	    [gom4_cv_$0=$gom4_$0_copy; break])
     done])
if test "$gom4_cv_$0" != no; then
    AC_DEFINE([HAVE_VA_COPY],[],[Define to 1 if you have `va_copy'.])
    if test x"$gom4_cv_$0" != x"va_copy"; then
	AC_DEFINE([va_copy],[$gom4_cv_$0],
	    [Define if you have `va_copy' with a different name.])
    fi
fi])

dnl # GOM_DEFINE_OS_NAME()
dnl # --------------------
dnl # Set OS_DEFS to -DOS_NAME='"<os name>"'
AC_DEFUN([GOM_DEFINE_OS_NAME],
[AC_CACHE_CHECK([for OS name],[gom4_cv_$0],
    [gom4_cv_$0=`uname -o 2>/dev/null || uname -s 2>/dev/null`
    test -n "$gom4_cv_$0" || gom4_cv_$0=='Unknown OS'])
OS_DEFS="-DOS_NAME='\"$gom4_cv_$0\"'"
AC_SUBST([OS_DEFS])])dnl
