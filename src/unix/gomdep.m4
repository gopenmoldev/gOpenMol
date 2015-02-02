dnl -*- mode: autoconf -*-
AC_COPYRIGHT([
############################################################################
#                           Copyright (c) 2003 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
############################################################################
])dnl

dnl # _GOM_DEP_TRACKING()
dnl # -------------------
AC_DEFUN([_GOM_DEP_TRACKING],
[AC_ARG_ENABLE([deps],
[AC_HELP_STRING([--disable-deps],
	[Do not determine dependencies of object files.
	Will speed up one time build.])
AC_HELP_STRING([--enable-deps],
	[Do determine dependencies of object files even if
	significantly more work is required.])],
    [gom4_dep_tracking="$enableval"])])

dnl # GOM_PROG_DEPCOMP(COMPILER)
dnl # --------------------------
dnl # Sets: <COMPILER>DEPMODE to `depmode=<depmode>'
dnl #       <COMPILER>DEPCOMP to `$(<COMPILER>DEPMODE) $(DEPCOMP)' or to `'
AC_DEFUN([GOM_PROG_DEPCOMP],
[AC_REQUIRE([_GOM_DEP_TRACKING])
AC_CACHE_CHECK([dependency mode for $$1],[gom4_cv_$1_depmode],
    [case $gom4_dep_tracking in
    no) gom4_cv_$1_depmode=none ;; # Disable totally.
    *)
	# ac_aux_dir may be relative.
	gom4_$0_depcomp=`cd "${ac_aux_dir-.}" && pwd`/depcomp
	# Create a sub directory for testing.
	rm -fr conf$$
	mkdir conf$$
	cd conf$$
	for gom4_cv_$1_depmode in `\
	    sed -n ['s/^\#*\([^ *()]*\)).*/\\1/p'] <"\$gom4_$0_depcomp"`
	do
	    case $gom4_cv_$1_depmode in
	    nosideeffect)
		case $gom4_dep_tracking in
		yes) ;; # Significantly more work is needed. Still continue.
		*) gom4_cv_$1_depmode=none; break ;; # Stop here.
		esac
		;;
	    esac
	    # Recreate test files.
	    echo "int test;" >test.h
	    echo \#include \"test.h\" >test.c
	    echo \# dummy >test.d
	    # Try to compile.
	    depmode="$gom4_cv_$1_depmode" \
	    source=test.c \
	    object=test.o \
	    depfile=test.d \
	    $SHELL $gom4_$0_depcomp $$1 -c -o test.o test.c >/dev/null 2>&1 &&
		grep 'test\.h' test.d >/dev/null 2>&1 &&
		break
	done
	# Go back.
	cd ..
	# Clean up.
	rm -fr conf$$
	;;
    esac])
$1DEPMODE="depmode=$gom4_cv_$1_depmode"
case $gom4_cv_$1_depmode in
none) $1DEPCOMP= ;;
*)    $1DEPCOMP='$($1DEPMODE) $(DEPCOMP)' ;;
esac
AC_SUBST([$1DEPCOMP])dnl
AC_SUBST([$1DEPMODE])])dnl
