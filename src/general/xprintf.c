/*  

Copyright (c) 2004 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Coded by:
Eero Häkkinen
*/

#include "maindefs.h"

#if STDC_HEADERS
#  include <stddef.h>
#endif
#if HAVE_STDINT_H
#  include <stdint.h>
#endif
#if HAVE_INTTYPES_H
#  include <inttypes.h>
#endif

#define KEEP_STDIO_MACROS

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>

#include "stdafx.h"

#include "gomfmt.h"

#ifndef HAVE_C99_XPRINTF

/**
 * Function:
 *   - Mimic following C99 functions:
 *         int printf(const char *format, ...);
 *         int fprintf(FILE *stream, const char *format, ...);
 *         int sprintf(char *dst, const char *format, ...);
 *         int snprintf(char *dst, size_t count, const char *format, ...);
 *
 *         int vprintf(const char *format, va_list args);
 *         int vfprintf(FILE *stream, const char *format, va_list args);
 *         int vsprintf(char *dst, const char *format, va_list args);
 *         int vsnprintf(
 *             char *dst, size_t count, const char *format, va_list args);
 *   - Implement following BSD and GNU functions:
 *         int asprintf(char **dst, const char *format, ...);
 *         int vasprintf(char **dst, const char *format, va_list args);
 * Return:
 *   - Number of characters which would have been written (excluding the
 *     terminating nul character) if there would have been enough
 *     space.
 * Recognized format specification fields:
 *     Fields are in form of
 *     %[<flags>][<width>][.<precision>][<size>]<type>
 *     Asterisk (*) may be used as a <width> and/or <precision> in
 *     which case the width and/or precision is taken from the next
 *     argument which must be of a type int.
 *
 *     Integer conversions (d, i, o, u, x, X):
 *       - System sprintf must recognized these.
 *       - All flags recognized by a system sprintf are recognized.
 *
 *       - Allowed <size> modifiers are:
 *           * hh: signed char / unsigned char
 *           * h:  short int / unsigned short int
 *           * l   long int / unsigned long int
 *           * ll: long long int / unsigned long long int /
 *                 __int64 / unsigned __int64
 *           * j:  intmax_t / uintmax_t
 *           * z:  size_t
 *           * t:  ptrdiff_t
 *       - intmax_t, uintmax_t and size_t are passed to system sprintf
 *         as long long int or unsigned long long int.
 *       - ptrdiff_t is passed to sprintf as
 *         int, long int or long long int.
 *     Floating point conversions (e, E, f, F, g, G, a, A), character
 *     conversion (c) and a pointer conversion (p):
 *       - System sprintf must recognized these.
 *       - All flags recognized by a system sprintf are recognized.
 *       - <size> modifier 'L' (long double) is supported but the value is
 *         converted to double before passing to system sprintf.
 *     String conversion (s):
 *       - The only flag which is recognized is '-' (left justified).
 *     Wide character conversions (lc, ls):
 *       - system sprintf must recognized these.
 *       - All flags recognized by a system sprintf are recognized.
 *     Position conversion (n):
 *       - An argument must be a pointer to an integer.
 *       - The same <size> modifier as for integer conversions are allowed.
 *
 *     If a resulting string from a single conversion is very long or
 *     hard to estimate, a temporary file and fprintf is used instead
 *     of sprintf.
 */

int gomp_vfprintf(FILE *file, char const *format, va_list args)
{
    vxprintf_data data = { NULL, 0 };
    data.format = format;
    data.args   = args;
    data.file   = file;
    gomp_vxprintf(&data);
    return ferror(file) ? -1 : (int)data.index;
}

int gomp_vprintf(char const *format, va_list args)
{
    return gomp_vfprintf(stdout, format, args);
}

/**
 * Copy from string to ((vsnprintf_data*)data)->dst.
 * There is no limit test in sprintf.
 */
static int do_gomp_vsprintf(
    const char *string, size_t len, vxprintf_data *data)
{
    /* Do range test for sprintf.
     */
    strncpy(((vsnprintf_data*)data)->dst + data->index, string, len);
    return 1;
}
int gomp_vsprintf(char *dst, char const *format, va_list args)
{
    vsnprintf_data data = { { do_gomp_vsprintf, 0 } };
    data.data.format = format;
    data.data.args   = args;
    data.dst         = dst;
    gomp_vxprintf(&data.data);
    dst[data.data.index] = '\0';
    return data.data.index;
}

/**
 * Copy from string to ((vsnprintf_data*)data)->dst if there is space.
 * Ignore the rest. Never fail because the total length is needed.
 */
static int do_gomp_vsnprintf(
    const char *string, size_t len, vxprintf_data *data)
{
    vsnprintf_data *sn = (vsnprintf_data*)data;
    size_t i = data->index;
    while ( len > 0 && i < sn->count ) {
        sn->dst[i++] = *string++;
        len--;
    }
    return 1;
}
int gomp_vsnprintf(
    char *dst, size_t count, char const *format, va_list args)
{
    vsnprintf_data data = { { do_gomp_vsprintf, 0 } };
    data.data.format = format;
    data.data.args   = args;
    data.dst         = dst;
    data.count       = count;
    gomp_vxprintf(&data.data);
    if ( dst && count > 0 )
        dst[data.data.index < count ? data.data.index : ( count - 1 )] = '\0';
    return data.data.index;
}

DELEGATE_TO_VFMT(print, (char const *format,  ...), (format, args))
DELEGATE_TO_VFMT(
    fprint, (FILE *file, char const *format, ...), (file, format,args))
DELEGATE_TO_VFMT(
    sprint, (char *dst, char const *format, ...), (dst, format, args))
DELEGATE_TO_VFMT(
    snprint,
    (char *dst, size_t count, char const *format, ...),
    (dst, count, format, args))

#else /* HAVE_C99_XPRINTF */
#define gomp_vsnprintf vsnprintf
#define gomp_vsprintf  vsprintf
#endif /* HAVE_C99_XPRINTF */

#ifdef NEED_GOMP_VXPRINTF

static int vxprintf_arg(vxprintf_data *);
static int vxprintf_strarg(vxprintf_data *, const char *, const char *);

/**
 * Format arguments from args according to format to *pDst.
 * Arguments pDst and count are arguments passed to append_text.
 * Function append_text may use them in any way. This function doesn't use
 * them directly.
 * */
int gomp_vxprintf(vxprintf_data* data)
{
    size_t length;
    
    for ( ;; ) {
        switch ( *data->format ) {
        case '\0':
            return 1;
        case '%':
            switch ( *(data->format + 1) ) {
            case '%':
                data->format++;
                /* fall through */
            case '\0':
                /* Add a '%' character. */
                if ( data->file )
                    fputc('%', data->file);
                else if ( ! data->append_text("%", 1, data) )
                    return 0;
                data->format++;
                data->index++;
                break;
            default:
                /* Process the next argument. */
                if ( ! vxprintf_arg(data) ) {
                    /* Unknown conversion. We can't continue safely. */
                    /* Append the format directly.                   */
                    length = strlen(data->format);
                    if ( data->file )
                        fwrite(data->format, 1, length, data->file);
                    else if ( ! data->append_text(data->format, length, data) )
                        return 0;
                    data->format += length;
                    data->index  += length;
                    return 0;
                }
            }
            break;
        default:
            /* Append normal text. */
            length = 0;
            while ( data->format[length] && data->format[length] != '%' )
                length++;
            if ( data->file )
                fwrite(data->format, 1, length, data->file);
            else if ( ! data->append_text(data->format, length, data) )
                return 0;
            data->format += length;
            data->index  += length;
        }
    }
}

/**
 * Format one argument from *pArgs according to *pFormat to *pDst.
 * Arguments pDst and pCount are arguments passed to append_text.
 * Function append_text may use them in any way. This function doesn't
 * use them directly.
 */
int vxprintf_arg(vxprintf_data* data)
{
    char temp[PRINTF_BUFF_LEN];
    char subformat[FORMAT_FIELD_LEN];
    int width, size, count;
    int min_width = 0, i= 0;
    FILE *file = data->file;
    const char *format = data->format;

    /**
     * Format is "%[<flags>][<width>][.<precision>][<size>][<type>]".
     * Copy a "%<flags><width>.<precision>" into the subformat buffer.
     */
    subformat[i++] = *format++;
    while ( *format && !isalpha((unsigned char)*format) ) {
        if ( *format == '*' ) {
            /* Get width or precision from the arg list. */
            width = va_arg(data->args,int);
            format++;
        }
        else if ( isdigit((unsigned char)*format) && *format != '0' ) {
            /* Get width or precision from the format string. */
            width = atoi(format);
            while ( isdigit((unsigned char)*format) )
                format++;
        }
        else {
            subformat[i++] = *format++;
            continue;
        }
        if ( min_width < width )
            min_width = width;
        sprintf(subformat + i, "%d", width);
        i += strlen(subformat + i);
    }

    /** Parse "<size>" and "<type>". */
    if ( *format == 's' ) {
        const char *string = va_arg(data->args, const char *);
        subformat[i++] = *format++;
        subformat[i  ] = '\0';
        if ( data->file )
            fprintf(data->file, subformat, string);
        else if ( ! vxprintf_strarg(data, subformat, string) )
            /* Unknown format. */
            return 0;
        goto conversion_done;
    }

    if ( min_width >= MAX_PRINTF_WIDTH ) {
        /** Cannot be done safely by using sprintf. Let's use fprintf. */
        if ( ! file && ! ( file = gomp_tmpfile() ) )
            return 0;
    }

    size = gomp_parse_format_size(
        subformat, &i, &format, &data->args, data->index);

#define PRINT_ARG2(type, as_type) \
    if ( file ) \
        fprintf(file, subformat, (as_type)va_arg(data->args, type), &count); \
    else \
        sprintf(temp, subformat, (as_type)va_arg(data->args, type), &count)
#define PRINT_ARG(type) PRINT_ARG2(type, type)

    subformat[i] = *format++;
    strcpy(subformat + i + 1, "%n");
    count = -1;
    switch ( subformat[i] ) {
    case 'n':
        /** Handled by gomp_parse_format_size(...). */
        goto conversion_done;
    case INT_TYPE:
        switch ( size ) {
        case pfsz_hh: PRINT_ARG2(int,       char     ); break;
        case pfsz_h:  PRINT_ARG2(int,       short int); break;
        case pfsz:    PRINT_ARG(int      ); break;
        case pfsz_l:  PRINT_ARG(long int ); break;
        case pfsz_ll: PRINT_ARG(LLONG    ); break;
        case pfsz_j:
            switch ( PFSZ_DELEGATE_AS(intmax_t) ) {
            case pfsz:    PRINT_ARG2(intmax_t, int      ); break;
            case pfsz_l:  PRINT_ARG2(intmax_t, long int ); break;
            case pfsz_ll: PRINT_ARG2(intmax_t, LLONG    ); break;
            default: assert(0);
            }
            break;
        case pfsz_t:
            switch ( PFSZ_DELEGATE_AS(ptrdiff_t) ) {
            case pfsz:    PRINT_ARG2(ptrdiff_t, int      ); break;
            case pfsz_l:  PRINT_ARG2(ptrdiff_t, long int ); break;
            case pfsz_ll: PRINT_ARG2(ptrdiff_t, LLONG    ); break;
            default: assert(0);
            }
            break;
        case pfsz_z:
            switch ( PFSZ_DELEGATE_AS(size_t) ) {
            case pfsz:    PRINT_ARG2(size_t, int      ); break;
            case pfsz_l:  PRINT_ARG2(size_t, long int ); break;
            case pfsz_ll: PRINT_ARG2(size_t, LLONG    ); break;
            default: assert(0);
            }
            break;
        default: goto error;
        }
        break;
    case UINT_TYPE:
        switch ( size ) {
        case pfsz_hh: PRINT_ARG2(int,       unsigned char     ); break;
        case pfsz_h:  PRINT_ARG2(int,       unsigned short int); break;
        case pfsz:    PRINT_ARG(unsigned int      ); break;
        case pfsz_l:  PRINT_ARG(unsigned long int ); break;
        case pfsz_ll: PRINT_ARG(ULLONG            ); break;
        case pfsz_j:
            switch ( PFSZ_DELEGATE_AS(intmax_t) ) {
            case pfsz:    PRINT_ARG2(uintmax_t, unsigned int      ); break;
            case pfsz_l:  PRINT_ARG2(uintmax_t, unsigned long int ); break;
            case pfsz_ll: PRINT_ARG2(uintmax_t, ULLONG            ); break;
            default: assert(0);
            }
            break;
        case pfsz_t:
            switch ( PFSZ_DELEGATE_AS(ptrdiff_t) ) {
            case pfsz:    PRINT_ARG2(ptrdiff_t, unsigned int      ); break;
            case pfsz_l:  PRINT_ARG2(ptrdiff_t, unsigned long int ); break;
            case pfsz_ll: PRINT_ARG2(ptrdiff_t, ULLONG            ); break;
            default: assert(0);
            }
            break;
        case pfsz_z:
            switch ( PFSZ_DELEGATE_AS(size_t) ) {
            case pfsz:    PRINT_ARG2(size_t, unsigned int      ); break;
            case pfsz_l:  PRINT_ARG2(size_t, unsigned long int ); break;
            case pfsz_ll: PRINT_ARG2(size_t, ULLONG            ); break;
            default: assert(0);
            }
            break;
        default: goto error;
        }
        break;
    case FLOAT_TYPE:
        {
            double val;
            switch ( size ) {
            case pfsz:   val = va_arg(data->args, double); break;
            case pfsz_L: val = va_arg(data->args, long double); break;
            default: goto error;
            }
            switch ( subformat[i] ) {
            case 'f':
            case 'F':
                /**
                 * Format in style of "%.17f" may result a very long
                 * string.
                 */
                if ( ! file &&
                     fabs(val) > MAX_PRINTF_F_VALUE &&
                     ! ( file = tmpfile() ) )
                    return 0;
                break;
            }
            if ( file )
                fprintf(file, subformat, val, &count);
            else
                sprintf(temp, subformat, val, &count);
        }
        break;
    case 'p':
        switch ( size ) {
        case pfsz: PRINT_ARG(void *); break;
        default: goto error;
        }
        break;
    case 'c':
        switch ( size ) {
        case pfsz:   PRINT_ARG(int);     break;
        case pfsz_l: PRINT_ARG(wchar_t); break;
        default: goto error;
        }
        break;
    case 's':
        if ( ! file && ! ( file = tmpfile() ) )
            return 0;
        switch ( size ) {
        case pfsz_l: PRINT_ARG(wchar_t *); break;
        default: goto error;
        }
        break;
    default: goto error;
    }

    if ( count < 0 )
        goto error;
    if ( file ) {
        if ( file != data->file ) {
            int ch;
            if ( ferror(file) )
                goto error;
            rewind(file);
            while ( ( ch = getc(file) ) != EOF ) {
                char c = ch;
                if ( ! data->append_text(&c, 1, data) )
                    goto error;
            }
        }
    }
    else if ( ! data->append_text(temp, strlen(temp), data) )
        goto error;
    data->index += count;

  conversion_done:
    if ( file && file != data->file )
        fclose(file);
    data->format = format;
    return 1;

  error:
    if ( file && file != data->file )
        fclose(file);
    return 0;
}

/**
 * Format one string argument from *pArgs according to format to *pDst.
 * Format must be in form of "%[-][width][.precision]s".
 * Arguments pDst and pCount are arguments passed to append_text.
 * Function append_text may use them in any way. This function doesn't use
 * them directly.
 * */
int vxprintf_strarg(
    vxprintf_data *data, const char *format, const char *string)
{
    int iWidth = 0, iLength = 0;
    int left_justified = 0;
    size_t length, pad = 0;

    /* Format of the format must be */
    /* "%[-][<width>][.<precision>]s"    */
    format++; /* Skip the % character. */
    if ( *format == '-' ) {
        left_justified = 1;
        format++;
    }
    if ( strcmp(format, "s") == 0 )
        /* Easy case: "%s" or "%-s". */
        length = strlen(string);
    else {
        /* Parse width and length. */
        if ( sscanf(format, "%d.%ds", &iWidth, &iLength) == 2 ) {
            length = iLength < 0 ? 0 : iLength;
            pad    = iWidth  < (int)length ? 0 : ( iWidth - length );
        }
        else if ( sscanf(format, ".%ds", &iLength) == 1 ) {
            length = iLength < 0 ? 0 : iLength;
        }
        else if ( sscanf(format, "%ds", &iWidth) == 1 ) {
            length = strlen(string);
            pad    = iWidth  < (int)length ? 0 : ( iWidth - length );
        }
        else
            /* Unknown format. */
            return 0;
    }
    if ( left_justified ) {
        /* Append the string. */
        if ( ! data->append_text(string, length, data) )
            return 0;
        data->index += length;
        /* Add padding. */
        while ( pad > 0 ) {
            if ( ! data->append_text(" ", 1, data) )
                return 0;
            data->index++;
            pad--;
        }
    }
    else {
        /* Add padding. */
        while ( pad > 0 ) {
            if ( ! data->append_text(" ", 1, data) )
                return 0;
            data->index++;
            pad--;
        }
        /* Append the string. */
        if ( ! data->append_text(string, length, data) )
            return 0;
        data->index += length;
    }
    return 1;
}

#endif /* NEED_GOMP_VXPRINTF */

#if defined(NEED_GOMP_VXPRINTF) || !defined(HAVE_C99_XSCANF)

int gomp_parse_format_size(
    char *subformat, int *pIndex, const char **pFormat,
    va_list *pArgs, size_t index)
{
    int size;
    switch ( **pFormat ) {
    case 'h':
        /** Do not pass "h"s to a system [fs]print nor [fs]scanf. */
        ++*pFormat;
        switch ( **pFormat ) {
        case 'h':
            ++*pFormat;
            size = pfsz_hh;
            break;
        default:
            size = pfsz_h;
            break;
        }
        break;
    case 'l':
        ++*pFormat;
        switch ( **pFormat ) {
        case 'l':
            ++*pFormat;
            size = pfsz_ll;
fmt_lld:
            strcpy(&subformat[*pIndex], PRINTF_FMT_LL);
            *pIndex += strlen(&subformat[*pIndex]);
            break;
        default:
            size = pfsz_l;
fmt_ld:
            subformat[(*pIndex)++] = 'l';
            break;
        }
        break;
    case 'j':
        ++*pFormat;
        size = pfsz_j;
        switch ( PFSZ_DELEGATE_AS(intmax_t) ) {
        case pfsz:    goto fmt_d;
        case pfsz_l:  goto fmt_ld;
        case pfsz_ll: goto fmt_lld;
        default: assert(0);
        }
        break;
    case 't':
        ++*pFormat;
        size = pfsz_t;
        switch ( PFSZ_DELEGATE_AS(ptrdiff_t) ) {
        case pfsz:    goto fmt_d;
        case pfsz_l:  goto fmt_ld;
        case pfsz_ll: goto fmt_lld;
        default: assert(0);
        }
        break;
    case 'z':
        ++*pFormat;
        size = pfsz_z;
        switch ( PFSZ_DELEGATE_AS(size_t) ) {
        case pfsz:    goto fmt_d;
        case pfsz_l:  goto fmt_ld;
        case pfsz_ll: goto fmt_lld;
        default: assert(0);
        }
        break;
    case 'L':
        size = pfsz_L;
        subformat[(*pIndex)++] = 'l';
        break;
    default:
        size = pfsz;
    }
fmt_d:
    if ( **pFormat == 'n' ) {
        switch ( size ) {
        case pfsz_hh: *va_arg((*pArgs), char *)      = index; break;
        case pfsz_h:  *va_arg((*pArgs), short int *) = index; break;
        case pfsz:    *va_arg((*pArgs), int *)       = index; break;
        case pfsz_l:  *va_arg((*pArgs), long int *)  = index; break;
        case pfsz_ll: *va_arg((*pArgs), LLONG *)     = index; break;
        case pfsz_j:  *va_arg((*pArgs), intmax_t *)  = index; break;
        case pfsz_t:  *va_arg((*pArgs), ptrdiff_t *) = index; break;
        case pfsz_z:  *va_arg((*pArgs), size_t *)    = index; break;
        default: *pFormat = " "; /** Error */
        }
    }
    return size;
}

#endif
