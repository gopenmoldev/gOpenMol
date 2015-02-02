/*
   Print utility functions 

   Leif Laaksonen, CSC , 1995 - 2004
   Eero Häkkinen 2003 - 2004
*/

#ifndef INC_GOPENMOL_PRINTMSG
#define INC_GOPENMOL_PRINTMSG

#include <stdarg.h>

#include "gomcext.h"

/* message levels            */
extern int  gomp_PrintMessage(const char *);     /* print message */
extern int  gomp_PrintWARNING(const char *);     /* warning message */
extern int  gomp_PrintERROR(const char *);       /* error message */
extern int  gomp_PrintEXIT(const char *);        /* program has to halt */

extern int  gomp_FormatMessage(const char *, ...) /* print message */
                CHECK_FORMAT_PRINTF;
extern int  gomp_FormatWARNING(const char *, ...) /* warning message */
                CHECK_FORMAT_PRINTF;
extern int  gomp_FormatERROR(const char *, ...)   /* error message */
                CHECK_FORMAT_PRINTF;
extern int  gomp_FormatEXIT(const char *, ...)    /* program has to halt */
                CHECK_FORMAT_PRINTF;

extern int  gomp_PrintBANNER(void);
extern int  gomp_PrintUSAGE(void);
extern int  gomp_SuppressErrors(int);
extern int  gomp_GetDebugLevel(void);

#endif /* INC_GOPENMOL_PRINTMSG */
