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

dnl # _GOM_SET_PATH(CATEGORY,[PATH])
dnl # -------------------------
dnl # Defines M4 macro _GOM_<CATEGORY>_PATH.
dnl # May be called multiple times. Warns if called with different PATHs.
AC_DEFUN([_GOM_SET_PATH],
[m4_case(_GOM_$1_PATH,
    [_GOM_$1_PATH],[m4_define([_GOM_$1_PATH],m4_default([$2],[$PATH]))],
    m4_default([$2],[$PATH]),[],
    [m4_ifval([$2],[m4_warning(
        [multiple values for PATH for $1 given: `]_GOM_$1_PATH[' != `$2']
)])])])

dnl # _GOM_BASENAME_FOR_PATH_PROG(VARIABLE)
dnl # -------------------------------------
AC_DEFUN([_GOM_BASENAME_FOR_PATH_PROG],
[gom4_$0_dir=`expr "\$$1" : '\(.*\)[[\\\\/]][[^\\\\/]]*\$'`
gom4_$0_prog=`expr "\$$1" : '.*[[\\\\/]]\([[^\\\\/]]*\)\$'`
case "$PATH_SEPARATOR$PATH$PATH_SEPARATOR" in
*"$PATH_SEPARATOR$gom4_$0_dir$PATH_SEPARATOR"*) $1=$gom4_$0_prog ;; # In PATH.
esac])

dnl # _GOM_DIRNAME(DIR)
dnl # -----------------
AC_DEFUN([_GOM_DIRNAME],[echo x"$1" | sed 's,^x\(.*\)[[\\\\/]].*,\1,'])

dnl # GOM_CHECK_PROG_VER(VARIABLE, PROG-TO-CHECK-FOR, [VALUE-IF-NOT-FOUND],
dnl #                    [PATH], [PREFERRED-VERSION], [REJECT])
dnl # -----------------------------------------------
dnl # Search from the PATH for program called PROG-TO-CHECK-FOR.
dnl # Program may have a version number suffix.
dnl # Subversions are separated by dashes and periods.
dnl # Also, the version suffix may begin with the dash or period.
dnl # Otherwise: similar to AC_CHECK_PROG.
AC_DEFUN([GOM_CHECK_PROG_VER],
[if test -z "$$1"
then
    # Get the first word of the program (that is: omit arguments).
    set x $2; gom4_$0_prog=$[2]
    # Search for a possible program.
    gom4_$0_result=
    gom4_$0_IFS="$IFS"
    IFS="$PATH_SEPARATOR"
    # Use directories in the PATH (may be overwritten by programmer).
    for gom4_$0_dir in m4_default([$4],[$PATH])
    do
	IFS="$gom4_$0_IFS"
	# We want an absolute directory.
	case $gom4_$0_dir in "") gom4_$0_dir=. ;; esac
	gom4_$0_dir=`cd "$gom4_$0_dir" 2>/dev/null && pwd`
	case $gom4_$0_dir in "") continue ;; esac
	for gom4_$0_file in \
	    "$gom4_$0_dir/$gom4_$0_prog"* \
	    "$gom4_$0_dir/$gom4_$0_prog"$5
	do
	    test -f "$gom4_$0_file" || continue
	    # Do not accept programs like <dir>/<proc>-info.
	    # Only dashes and periods allowed between program name and
	    # version number. Version part is optional.
	    case `expr x"$gom4_$0_file" : '.*[[\\\\/]]\(.*\)'` in
	    "$gom4_$0_prog") ;;
	    "$gom4_$0_prog$EXEEXT") ;;
	    "$gom4_$0_prog"[[0-9]]*) ;;
	    "$gom4_$0_prog"[[-.]][[0-9]]*) ;;
	    *) continue ;;
	    esac
	    # Skip rejected (probably non-working) programs.
	    m4_ifval([$6],
		[case "$PATH_SEPARATOR$6$PATH_SEPARATOR" in
		*"$PATH_SEPARATOR$gom4_$0_file$PATH_SEPARATOR"*) continue ;;
		esac])
	    echo "$as_me:$LINENO: found $gom4_$0_file" >&AS_MESSAGE_LOG_FD
	    gom4_$0_result="$gom4_$0_file"
	    # 8.5 comes after 8.2 in sort order, so continue.
	done
	test -z "$gom4_$0_result" || break
    done
    IFS="$gom4_$0_IFS"
    # Let's see what we got.
    case $gom4_$0_result in
    "") m4_ifval([$3],[$1="$3"]) ;; # Nothing found.
    *)
	# Got the command name.
	$1=$gom4_$0_result
	# Directory may be omitted if program is the first one in the PATH.
	m4_ifval([$6],
	    [case "$6" in "") _GOM_BASENAME_FOR_PATH_PROG([$1]) ;; esac],
	    [_GOM_BASENAME_FOR_PATH_PROG([$1])])
	;;
    esac
fi
AC_SUBST([$1])])dnl

dnl # _GOM_CHECK_WORKING_PROG_VER(VARIABLE, PROG-TO-CHECK-FOR,
dnl #                             TEST-ARGS, FILE-SUFFIX, TEST-INPUT, OUTPUT,
dnl #                             [VALUE-IF-NOT-FOUND],
dnl #                             [PATH], [PREFERRED-VERSION])
dnl # -----------------------------------------------------------------------
AC_DEFUN([_GOM_CHECK_WORKING_PROG_VER],
[# Get the first word of the program (that is: omit arguments).
set x $2; gom4_$0_prog=$[2]
AC_CACHE_CHECK([for $gom4_$0_prog],[gom4_cv_prog_$1],
    [if test "${$1+set}" = set
    then
	# User knows best. Use the value supplied by user.
	gom4_cv_prog_$1="$$1"
    else
	# Accept only working interpreters.
	gom4_$0_reject=
	while test "${gom4_cv_prog_$1+set}" != set
	do
	    # Get the newest interpreter.
	    gom4_$0_prog=
	    GOM_CHECK_PROG_VER(
		[gom4_$0_prog],[$2],[:],[$8],[$9],
		[$gom4_$0_reject$PATH_SEPARATOR])
	    case $gom4_$0_prog in :) break ;; esac
	    # Test if it works.
	    cat >conftest$4 <<_GOM4_EOF
$5
_GOM4_EOF
	    set x $2; shift; shift
	    AC_TRY_COMMAND(
		[$gom4_$0_prog ${1+"[$]@"} $3 conftest$4 >conftest.out])
	    case $?-`cat conftest.out` in
	    0-$6)
		# It works!
		gom4_cv_prog_$1="$gom4_$0_prog"
		case $gom4_$0_reject in
		"") _GOM_BASENAME_FOR_PATH_PROG([gom4_cv_prog_$1]) ;;
		esac
		;;
	    *)
		echo "$as_me: failed script was:" >&AS_MESSAGE_LOG_FD
		sed 's/^/| /' < conftest$4 >&AS_MESSAGE_LOG_FD
		# It doesn't work. Get the next one.
		gom4_$0_reject="$gom4_$0_reject$PATH_SEPARATOR$gom4_$0_prog"
		;;
	    esac
	    rm -f conftest$4
	done
    fi])
if test "${gom4_cv_prog_$1+set}" = "set"
then
    $1="$gom4_cv_prog_$1"
m4_ifval([$7],
[else
    $1="$7"])
fi])
