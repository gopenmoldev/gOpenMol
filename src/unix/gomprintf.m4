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

dnl # GOM_FUNC_C99_PRINTF()
dnl # ---------------------
dnl # Check if we have snprintf and vsnprintf which conform C99.
AC_DEFUN([GOM_FUNC_C99_PRINTF],[AC_REQUIRE([_GOM_FUNC_C99_PRINTF])])
AC_DEFUN([_GOM_FUNC_C99_PRINTF],
[AC_REQUIRE([GOM_FUNC_C99_SCANF])dnl
if test -z "$gom4_cv_$0"; then
    # If sscanf supports all type modifiers, snprintf probably does, too.
    gom4_$0_ok=$gom4_cv_scanf_mod
    if test "$gom4_$0_ok" = yes; then
	# Check a return value and that snprintf honours the count parameter.
	AC_MSG_CHECKING([for snprintf])
	AC_RUN_IFELSE(
	    [AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT([])],
	    [[char b[5]="";
	      int n;
	      n=snprintf(b,1,"123");
	      return ( n==3 && !*b ) ? 0 : 1;]])],[],[gom4_$0_ok=no])
	AC_MSG_RESULT([$gom4_$0_ok])
	if test "$gom4_$0_ok" = yes; then
	    # It is probably enough to check only an existence of vsnprint.
	    AC_CHECK_FUNC([vsnprintf],[],[gom4_$0_ok=no])
	fi
    fi
    AC_CACHE_CHECK([for C99 vsnprintf],[gom4_cv_$0],[gom4_cv_$0=$gom4_$0_ok])
fi
gom4_GOM_FUNC_C99_PRINTF_ok=$gom4_cv_$0
if test "$gom4_GOM_FUNC_C99_PRINTF_ok" = yes; then
    AC_DEFINE(
	[HAVE_C99_XPRINTF],[],
	[Define to 1 if you have C99 [v][s[n]|f]printf.])
    $1
fi])

dnl # GOM_FUNC_C99_SCANF()
dnl # --------------------
dnl # Check if we have sscanf and vsscanf which supports all C99 type
dnl # modifiers.
AC_DEFUN([GOM_FUNC_C99_SCANF],[AC_REQUIRE([_GOM_FUNC_C99_SCANF])])
AC_DEFUN([_GOM_FUNC_C99_SCANF],
[AC_REQUIRE([_GOM_CHECK_SCANF_MOD])dnl
gom4_GOM_FUNC_C99_SCANF_ok=$gom4_cv_scanf_mod
if test "$gom4_GOM_FUNC_C99_SCANF_ok" = yes; then
    AC_CHECK_FUNC([vsscanf],[:],[gom4_GOM_FUNC_C99_SCANF_ok=no])
fi
if test "$gom4_GOM_FUNC_C99_SCANF_ok" = yes; then
    AC_DEFINE(
	[HAVE_C99_XSCANF],[],
	[Define to 1 if you have C99 [v][s|f]scanf.])
else
    GOM_CHECK_SIZEOF([int])
    GOM_CHECK_SIZEOF([long int])
    GOM_CHECK_SIZEOF([long long int])
    GOM_CHECK_SIZEOF([size_t])
    GOM_CHECK_SIZEOF([intmax_t])
    GOM_CHECK_SIZEOF([ptrdiff_t])
fi])

dnl # GOM_FUNC_C99_SCANF_MOD()
dnl # ------------------------
dnl # Check if we have sscanf supports all C99 type modifiers.
AC_DEFUN([_GOM_CHECK_SCANF_MOD],
[if test -z "$gom4_cv_scanf_mod"; then
    gom4_$0_ok=yes
    for gom4_$0_spec in \
	'hhd hhn char' \
	'hd  hn  short int' \
	'd   n   int' \
	'ld  ln  long int' \
	'lld lln long long int' \
	'jd  jn  intmax_t' \
	'td  tn  ptrdiff_t' \
	'zd  zn  size_t' \
	'f   n   float' \
	'lf  n   double' \
	'Lf  n   long double'
    do
	set $gom4_$0_spec
	gom4_$0_fmt=%[$]1%[$]2
	AC_MSG_CHECKING([for sscanf, %[$]1 and %[$]2])
	shift 2
	case $[1] in
	(*f) gom4_$0_code="int  n=0;" ;;
	(*)  gom4_$0_code="[$]* n=0;" ;;
	esac
	gom4_$0_code="
	    $gom4_$0_code
	    $[*] x;
	    int  nf;
	    nf = sscanf(\"0_\",\"$gom4_$0_fmt\", &x, &n);
	    return ( nf == 1 && n > 0 ) ? 0 : 1;"
	AC_RUN_IFELSE(
	    [AC_LANG_PROGRAM([AC_INCLUDES_DEFAULT([])],[$gom4_$0_code])],
	    [AC_MSG_RESULT([yes])],
	    [AC_MSG_RESULT([no])
	     gom4_$0_ok=no])
    done
    AC_CACHE_CHECK(
	[for C99 sscanf],[gom4_cv_scanf_mod],[gom4_cv_scanf_mod=$gom4_$0_ok])
fi
])dnl
