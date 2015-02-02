/*  

Copyright (c) 2004 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Coded by:
Eero Häkkinen
*/

#include "maindefs.h"

#include <stdarg.h>
#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>

#include "stdafx.h"

#include "gomfmt.h"

#if !defined(HAVE_VASPRINTF) || !defined(HAVE_C99_XPRINTF)
#ifndef HAVE_VA_COPY
static int do_gomp_vasprintf(
    const char *string, size_t len, vxprintf_data *data)
{
    vsnprintf_data *sn = (vsnprintf_data*)data;
    size_t min_count   = data->index + len + 1;
    if ( sn->count < min_count ) {
        char *buffer;
        size_t count = 2 * sn->count + 1;
        if ( count < min_count )
            count = min_count;
        buffer = realloc(sn->dst, count);
        if ( ! buffer ) {
            buffer = realloc(sn->dst, min_count);
            if ( ! buffer ) {
                /** Once failed, there is no clean way to continue. */
                free(sn->dst);
                sn->dst = NULL;
                return 0;
            }
        }
    }
    strncpy(sn->dst + data->index, string, len);
    return 1;
}
#endif /* ! HAVE_VA_COPY */
int gomp_vasprintf(char **dst, const char *format, va_list args)
{
#ifdef HAVE_VA_COPY
    int len;
    va_list args2;

    va_copy(args,args2);
    len = vsnprintf(NULL, 0, format, args2);
    va_end(args2);
    *dst = malloc(len + 1);
    if ( ! *dst )
        return -1;
    return vsprintf(dst, format, args);
#else /* ! HAVE_VA_COPY */
    vsnprintf_data data = { { NULL, 0 } };
    data.data.format    = format;
    data.data.args      = args;
    data.count          = 63;
    data.dst            = malloc(data.count);
    if ( data.dst )
        gomp_vxprintf(&data.data);
    *dst = data.dst;
    if ( ! *dst )
        return -1;
    (*dst)[data.data.index] = '\0';
    return data.data.index;
#endif /* ! HAVE_VA_COPY */
}
#else /* HAVE_VASPRINTF && HAVE_C99_XPRINTF */
#define gomp_vasprintf vasprintf
#endif /* HAVE_VASPRINTF && HAVE_C99_XPRINTF */

#if !defined(HAVE_VASPRINTF) || !defined(HAVE_C99_XPRINTF)
DELEGATE_TO_VFMT(
    asprint,
    (char **dst, char const *format, ...),
    (dst, format, args))
#endif
