#ifndef INC_GOPENMOL_GOMSTDIO_H
#define INC_GOPENMOL_GOMSTDIO_H

/*
Copyright (c) 2004 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of
Leif Laaksonen
All rights reserved

Coded 2004 - 2005 by:
Eero Häkkinen
*/

#include <stdarg.h>
#include <stdio.h>

#include "gomcext.h"

#if defined(HAVE_ASPRINTF) && defined(HAVE_C99_XPRINTF)
#ifndef asprintf
extern int asprintf(char **, const char *, ...)
    CHECK_FORMAT_PRINTF_1_EXTRA;
#endif
#else
extern int gomp_asprintf(char **, const char *, ...)
    CHECK_FORMAT_PRINTF_1_EXTRA;
#undef  asprintf
#define asprintf gomp_asprintf
#endif

#if defined(HAVE_VASPRINTF) && defined(HAVE_C99_XPRINTF)
#ifndef vasprintf
extern int vasprintf(char **, const char *, va_list)
    CHECK_FORMAT_VPRINTF_1_EXTRA;
#endif
#else
extern int gomp_vasprintf(char **, const char *, va_list)
    CHECK_FORMAT_VPRINTF_1_EXTRA;
#undef  vasprintf
#define vasprintf gomp_vasprintf
#endif

/* We want functions conforming C99. */

#ifndef HAVE_C99_XPRINTF
extern int gomp_printf(const char *, ...)
    CHECK_FORMAT_PRINTF;
extern int gomp_fprintf(FILE *, const char *, ...)
    CHECK_FORMAT_PRINTF_1_EXTRA;
extern int gomp_sprintf(char *, const char *, ...)
    CHECK_FORMAT_PRINTF_1_EXTRA;
extern int gomp_snprintf(char *, size_t, const char *, ...)
    CHECK_FORMAT_PRINTF_2_EXTRA;
extern int gomp_vprintf(const char *, va_list)
    CHECK_FORMAT_VPRINTF;
extern int gomp_vfprintf(FILE *, const char *, va_list)
    CHECK_FORMAT_VPRINTF_1_EXTRA;
extern int gomp_vsprintf(char *, const char *, va_list)
    CHECK_FORMAT_VPRINTF_1_EXTRA;
extern int gomp_vsnprintf(char *, size_t, const char *, va_list)
    CHECK_FORMAT_VPRINTF_2_EXTRA;

#ifndef KEEP_STDIO_MACROS
#undef  printf
#undef  fprintf
#undef  sprintf
#undef  snprintf
#undef  asprintf

#undef  vprintf
#undef  vfprintf
#undef  vsprintf
#undef  vsnprintf
#undef  vasprintf

#define printf    gomp_printf
#define fprintf   gomp_fprintf
#define sprintf   gomp_sprintf
#define snprintf  gomp_snprintf
#define asprintf  gomp_asprintf

#define vprintf   gomp_vprintf
#define vfprintf  gomp_vfprintf
#define vsprintf  gomp_vsprintf
#define vsnprintf gomp_vsnprintf
#define vasprintf gomp_vasprintf
#endif /* ! KEEP_STDIO_MACROS */

#endif /* ! HAVE_C99_XPRINTF */

#ifndef HAVE_C99_XSCANF
extern int gomp_scanf(const char *, ...)
    CHECK_FORMAT_SCANF;
extern int gomp_fscanf(FILE *, const char *, ...)
    CHECK_FORMAT_SCANF_1_EXTRA;
extern int gomp_sscanf(const char *, const char *, ...)
    CHECK_FORMAT_SCANF_1_EXTRA;
extern int gomp_vscanf(const char *, va_list)
    CHECK_FORMAT_VSCANF;
extern int gomp_vfscanf(FILE *, const char *, va_list)
    CHECK_FORMAT_VSCANF_1_EXTRA;
extern int gomp_vsscanf(const char *, const char *, va_list)
    CHECK_FORMAT_VSCANF_1_EXTRA;

#ifndef KEEP_STDIO_MACROS
#undef  scanf
#undef  fscanf
#undef  sscanf
#undef  vscanf
#undef  vfscanf
#undef  vsscanf

#define scanf    gomp_scanf
#define fscanf   gomp_fscanf
#define sscanf   gomp_sscanf
#define vscanf   gomp_vscanf
#define vfscanf  gomp_vfscanf
#define vsscanf  gomp_vsscanf
#endif /* ! KEEP_STDIO_MACROS */

#endif /* ! HAVE_C99_XSCANF */

extern FILE* gomp_tmpfile(void);

#ifndef KEEP_STDIO_MACROS
#undef  tmpfile
#define tmpfile gomp_tmpfile
#endif

#endif /* ! INC_GOPENMOL_GOMSTDIO_H */
