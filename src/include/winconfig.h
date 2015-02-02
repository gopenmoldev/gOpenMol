/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define if you have the asprintf function.  */
/* #undef HAVE_ASPRINTF */

/* Define to 1 if you have C99 [v][s[n]|f]printf. */
#undef HAVE_C99_XPRINTF

/* Define to 1 if you have C99 [v][s|f]scanf. */
#undef HAVE_C99_XSCANF

/* Define to 1 if you have the `asprintf' function. */
#undef HAVE_ASPRINTF

/* Define to 1 if you have the `vasprintf' function. */
#undef HAVE_VASPRINTF

/* Use GLUT in graphics sub system. */
#define GLUT 

/* Define to 1 if you have the `jpeg' library (-ljpeg) and corresponding
   header files. */
#define HAVE_LIBJPEG 1

/* Define to 1 if you have the `pythonX.X' library (-lpythonX.X) and
   corresponding header files. */
#define HAVE_LIBPYTHON 1

/* Define if you want to enable API functions for extensions. */
#define ENABLE_EXTENSIONS 

/* Produce a library (will produce a gom_main entry point). */
#define GOPENMOL

/* Define to enable graphics sub system. */
#define ENABLE_GRAPHICS 

/* MSVC++ accepts ANSI C code. */
#define __STDC__ 1

/* Define to operation system name. */
#define OS_NAME "WINDOWS"

/* Define as `__inline' if that's what the C compiler calls it, or to nothing
   if it is not supported. */
#define inline __inline

#define popen  _popen
#define pclose _pclose
