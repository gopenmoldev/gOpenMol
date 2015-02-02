#include "maindefs.h"

#include <stdarg.h>
#include <stdlib.h>
#include "gomstdio.h"
#include <tcl.h>

#include "parser.h"
#include "printmsg.h"
/*#include "tclutils.h"*/

#include "stdafx.h"

#define IMPLEMENT_FORMAT_MSG(type) \
int gomp_Format##type(const char *format, ...) \
{ \
    int result; \
    va_list args1,args2; \
    va_start(args1,format); \
    va_start(args2,format); \
    result = DoPrintMessage(gomp_Print##type, format, args1, args2); \
    va_end(args1); \
    va_end(args2); \
    return result; \
}

static int DoPrintMessage(int (*)(const char *),const char *, va_list, va_list)
    CHECK_FORMAT_VPRINTF_1_EXTRA;

int DoPrintMessage(int (*PrintMsg)(const char *),
                   const char *format, va_list args1,va_list args2)
{
    char  message[BUFF_LEN];
    char *msg = NULL;
    int   result;

    result = vsnprintf(message, BUFF_LEN, format, args1);
    if ( result >= 0 && result < BUFF_LEN )
        /* No need to dynamically allocate memory. */
        return PrintMsg(message);
    /* We need a longer buffer. */
    /* It is not portable to pass the same va_list twice for v...printf. */
    /* So we have two va_list arguments which points to the same set of  */
    /* arguments.                                                        */
    if ( vasprintf(&msg, format, args2) < 0 )
        /* Out of memory. */
        return PrintMsg(format);
    result = PrintMsg(msg);
    free(msg);
    return result;
}

IMPLEMENT_FORMAT_MSG(Message)
IMPLEMENT_FORMAT_MSG(WARNING)
IMPLEMENT_FORMAT_MSG(ERROR)
IMPLEMENT_FORMAT_MSG(EXIT)

static int FormatTclReturn(Tcl_Interp *, const char *, va_list, va_list)
    CHECK_FORMAT_VPRINTF_1_EXTRA;

int FormatTclReturn(Tcl_Interp *interp, const char *format,
                    va_list args1, va_list args2)
{
    char  message[BUFF_LEN];
    char *msg = NULL;
    int   count;

    count = vsnprintf(message, BUFF_LEN, format, args1);
    if ( count < BUFF_LEN ) {
        /* No need to dynamically allocate memory (by us). */
        Tcl_SetResult(interp, message, TCL_VOLATILE);
        return 0;
    }
    
    /**
     * We need a longer buffer.  It is not portable to pass the same
     * va_list twice for v...printf.  So we have two va_list arguments
     * which points to the same set of arguments.
     */
    /* Allocate the buffer. */
    msg = Tcl_Alloc( count + 1 );
    if ( ! msg )
        /* Out of memory. */
        return(-1);
    vsprintf(msg, format, args2);
    Tcl_SetResult(interp, msg, TCL_DYNAMIC);
    return 0;
}

int gomp_ParserFormatReturn(Tcl_Interp *interp, const char *format, ...)
{
    int result;
    va_list args1,args2;
    va_start(args1,format);
    va_start(args2,format);
    result = FormatTclReturn(interp, format, args1, args2);
    va_end(args1);
    va_end(args2);
    return result;
}
