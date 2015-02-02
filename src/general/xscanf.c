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
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "gomstdio.h"

#include "stdafx.h"

#include "gomfmt.h"

#ifndef HAVE_C99_XSCANF

static int vxscanf(FILE *, const char *, const char *, va_list);

/**
 * Function:
 *   - Mimic following C99 functions:
 *         int scanf(const char *format, ...);
 *         int fscanf(FILE *stream, const char *format, ...);
 *         int sscanf(const char *dst, const char *format, ...);
 *
 *         int vscanf(const char *format, va_list args);
 *         int vfscanf(FILE *stream, const char *format, va_list args);
 *         int vsscanf(char *dst, const char *format, va_list args);
 * Return:
 *   - Number of successfully converted fields or EOF.
 * Recognized format specification fields:
 *     Fields are in form of
 *     %[<flags>][<width>][<size>]<type>
 *     The only accepted <flag> is an asterisk (*) which suppresses
 *     assignment.
 *
 *     Integer conversions (d, i, o, u, x, X):
 *       - System sscanf must recognized these.
 *       - Allowed <size> modifiers are:
 *           * hh: signed char / unsigned char
 *           * h:  short int / unsigned short int
 *           * l   long int / unsigned long int
 *           * ll: long long int / unsigned long long int /
 *                 __int64 / unsigned __int64
 *           * j:  intmax_t / uintmax_t
 *           * z:  size_t
 *           * t:  ptrdiff_t
 *       - Internal temporary variables are used to circumvent limitations
 *         of a system sscanf.
 *       - intmax_t, uintmax_t and size_t are passed to system sscanf
 *         as a pointer to long long int or unsigned long long int.
 *       - ptrdiff_t is passed to sscanf as pointer to
 *         int, long int or long long int.
 *       - size_t is passed to sscanf as a pointer to
 *         unsigned int, unsigned long int or unsigned long long int.
 *     Floating point conversions (e, E, f, F, g, G, a, A) and a
 *     pointer conversion (p):
 *       - System sscanf must recognized these.
 *       - <size> modifier 'L' (long double) is supported but not passed to
 *         a system sscanf.
 *     Character conversion (c, s, [seq], [^seq]):
 *       - Fully supported
 *     Wide character convertsions (lc, ls):
 *       - system sscanf must recognized these.
 *     Wide character sequence convertsions (l[seq], l[^seq]):
 *       - NOT SOPPORTED at all.
 *     Position conversion (n):
 *       - An argument must be a pointer to an integer.
 *       - The same <size> modifier as for integer conversions are allowed.
 */

int gomp_vfscanf(FILE *file, char const *format, va_list args)
{
    return vxscanf(file, NULL, format, args);
}

int gomp_vscanf(char const *format, va_list args)
{
    return vxscanf(stdout, NULL, format, args);
}

int gomp_vsscanf(const char *src, char const *format, va_list args)
{
    return vxscanf(NULL, src, format, args);
}

DELEGATE_TO_VFMT(scan, (char const *format,  ...), (format, args))
DELEGATE_TO_VFMT(
    fscan, (FILE *file, char const *format, ...), (file, format,args))
DELEGATE_TO_VFMT(
    sscan, (const char *src, char const *format, ...), (src, format, args))

static int vxscanf_literal(FILE *, const char **, char);
static int vxscanf_scanset(
    char *, FILE *, const char **, const char *, const char **);

int vxscanf(FILE *file, const char *src, const char *format, va_list args)
{
    int  index = 0, fields = 0;
    int  count, ignore, i, size;
    char subformat[FORMAT_FIELD_LEN];
    char ch;

    while ( *format ) {
        switch ( *format ) {
        case '%':
            if ( format[1] ) {
                if ( format[1] == '%' )
                    format++;
                else
                    break;
            }
            /** Fall through. */
        default:
            count = vxscanf_literal(file, &src, *format++);
            if ( count < 0 )
                goto eof;
            index += count;
            continue;
        }
        /**
         * Format is "%[<flags>][<width>][<size>][<type>]".
         * Copy "%<flags><width>" into the subformat buffer.
         */
        i      = 0;
        ignore = 0;
        subformat[i++] = *format++;
        if ( *format == '*') {
            ignore = 1;
            subformat[i++] = *format++;
        }
        while ( isdigit((unsigned char)*format) )
            subformat[i++] = *format++;
        /** Parse "<size>". */
        size = gomp_parse_format_size(subformat, &i, &format, &args, index);

        /** Parse "<type>". */
        subformat[i] = *format++;
        strcpy(subformat + i + 1, "%n");
        count = -1;
        if ( ignore ) {
            switch ( subformat[i] ) {
            case '[':
                switch ( size ) {
                case pfsz:
                    count = vxscanf_scanset(
                        NULL, file, &src, subformat, &format);
                    break;
                default: goto error;
                }
                break;
            default:
                if ( file )
                    fscanf(file, subformat, &count);
                else {
                    sscanf(src,  subformat, &count);
                    src += count;
                }
            }
        }
        else {
#define GET_ARG(type) va_arg(args, type)
#define SCAN_ARG(type) { \
        if ( file ) \
            fscanf(file, subformat, GET_ARG(type *), &count); \
        else { \
            sscanf(src,  subformat, GET_ARG(type *), &count); \
            src += count; \
        } \
    }
#define SCAN_ARG2(type, as_type) { \
        as_type v; \
        if ( file ) \
            fscanf(file, subformat, &v, &count); \
        else { \
            sscanf(src,  subformat, &v, &count); \
            src += count; \
        } \
        if ( count >= 0 ) \
            *GET_ARG(type *) = v; \
    }
            switch ( subformat[i] ) {
            case 'n':
                /** Handled by gomp_parse_format_size(...). */
                continue;
            case INT_TYPE:
                switch ( size ) {
                case pfsz_hh: SCAN_ARG(char     ); break;
                case pfsz_h:  SCAN_ARG(short int); break;
                case pfsz:    SCAN_ARG(int      ); break;
                case pfsz_l:  SCAN_ARG(long int ); break;
                case pfsz_ll: SCAN_ARG(LLONG    ); break;
                case pfsz_j:
                    switch ( PFSZ_DELEGATE_AS(intmax_t) ) {
                    case pfsz:    SCAN_ARG2(intmax_t, int      ); break;
                    case pfsz_l:  SCAN_ARG2(intmax_t, long int ); break;
                    case pfsz_ll: SCAN_ARG2(intmax_t, LLONG    ); break;
                    default: assert(0);
                    }
                    break;
                case pfsz_t:
                    switch ( PFSZ_DELEGATE_AS(ptrdiff_t) ) {
                    case pfsz:    SCAN_ARG2(ptrdiff_t, int      ); break;
                    case pfsz_l:  SCAN_ARG2(ptrdiff_t, long int ); break;
                    case pfsz_ll: SCAN_ARG2(ptrdiff_t, LLONG    ); break;
                    default: assert(0);
                    }
                    break;
                case pfsz_z:
                    switch ( PFSZ_DELEGATE_AS(size_t) ) {
                    case pfsz:    SCAN_ARG2(size_t, int      ); break;
                    case pfsz_l:  SCAN_ARG2(size_t, long int ); break;
                    case pfsz_ll: SCAN_ARG2(size_t, LLONG    ); break;
                    default: assert(0);
                    }
                    break;
                default: goto error;
                }
                break;
            case UINT_TYPE:
                switch ( size ) {
                case pfsz_hh: SCAN_ARG(unsigned char     ); break;
                case pfsz_h:  SCAN_ARG(unsigned short int); break;
                case pfsz:    SCAN_ARG(unsigned int      ); break;
                case pfsz_l:  SCAN_ARG(unsigned long int ); break;
                case pfsz_ll: SCAN_ARG(ULLONG            ); break;
                case pfsz_j:
                    switch ( PFSZ_DELEGATE_AS(intmax_t) ) {
                    case pfsz:    SCAN_ARG2(intmax_t, unsigned int      ); break;
                    case pfsz_l:  SCAN_ARG2(intmax_t, unsigned long int ); break;
                    case pfsz_ll: SCAN_ARG2(intmax_t, ULLONG            ); break;
                    default: assert(0);
                    }
                    break;
                case pfsz_t:
                    switch ( PFSZ_DELEGATE_AS(ptrdiff_t) ) {
                    case pfsz:    SCAN_ARG2(ptrdiff_t, unsigned int      ); break;
                    case pfsz_l:  SCAN_ARG2(ptrdiff_t, unsigned long int ); break;
                    case pfsz_ll: SCAN_ARG2(ptrdiff_t, ULLONG            ); break;
                    default: assert(0);
                    }
                    break;
                case pfsz_z:
                    switch ( PFSZ_DELEGATE_AS(size_t) ) {
                    case pfsz:    SCAN_ARG2(size_t, unsigned int      ); break;
                    case pfsz_l:  SCAN_ARG2(size_t, unsigned long int ); break;
                    case pfsz_ll: SCAN_ARG2(size_t, ULLONG            ); break;
                    default: assert(0);
                    }
                    break;
                default: goto error;
                }
                break;
            case FLOAT_TYPE:
                switch ( size ) {
                case pfsz:   SCAN_ARG(float ); break;
                case pfsz_l: SCAN_ARG(double); break;
                case pfsz_L: SCAN_ARG2(long double, double); break;
                default: goto error;
                }
                break;
            case 'p':
                switch ( size ) {
                case pfsz: SCAN_ARG(void *); break;
                default: goto error;
                }
                break;
            case 'c':
            case 's':
                switch ( size ) {
                case pfsz:   SCAN_ARG(char);    break;
                case pfsz_l: SCAN_ARG(wchar_t); break;
                default: goto error;
                }
                break;
            case '[':
                switch ( size ) {
                case pfsz:
                    count = vxscanf_scanset(
                        va_arg(args, char *), file, &src, subformat, &format);
                    break;
                default: goto error;
                }
                break;
            }
        }
        if ( count < 0 )
            goto eof;
        else
            index += count;
        fields++;
    }
  eof:
    return fields;
  error:
    return EOF;
}

int vxscanf_literal(FILE *file, const char **pSrc, char fc)
{
    int ch;

    if ( isspace((unsigned char)fc) ) {
        size_t count = 0;
        if ( file ) {
            while ( ( ch = getc(file) ) != EOF ) {
                if ( ! isspace((unsigned char)ch) ) {
                    ungetc(ch, file);
                    break;
                }
                count++;
            }
        }
        else {
            while ( **pSrc && isspace((unsigned char)**pSrc) ) {
                ++*pSrc;
                count++;
            }
        }
        return count;
    }
    else {
        if ( file ) {
            if ( ( ch = getc(file) ) == EOF )
                return -1;
            else if ( ch != fc ) {
                ungetc(ch, file);
                return -1;
            }
        }
        else {
            if ( **pSrc == fc )
                ++*pSrc;
            else
                return -1;
        }
        return 1;
    }
}

int vxscanf_scanset(
    char *dst,
    FILE *file, const char **pSrc, const char *subformat, const char **pFormat)
{
    int count, ch;
    int not    =  0;
    int width  = -1;
    const char *end, *p;

    if ( **pFormat == '^' ) {
        not = 1;
        ++*pFormat;
    }
    if ( ! **pFormat )
        return -1;
    end = strchr(*pFormat + 1, ']');

    if ( *++subformat == '*' )
        subformat++;
    if ( isdigit((unsigned char)*subformat) )
        width = atoi(subformat);

    for ( count = 0 ; count != width ; count++ ) {
        if ( file ) {
            if ( ( ch = getc(file) ) == EOF )
                break;
        }
        else {
            if ( ! ( ch = *(*pSrc)++ ) )
                break;
        }
        p = strchr(*pFormat, ch);
        if ( ( p && p < end ) ? not : ! not ) {
            if ( file )
                ungetc(ch, file);
            else
                --*pSrc;
            break;
        }
        if ( dst )
            *dst++ = ch;
    }
    if ( ! count )
        return -1;
    if ( dst )
        *dst = '\0';
    return count;
}

#endif /* ! HAVE_C99_XSCANF */
