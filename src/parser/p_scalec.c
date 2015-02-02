/*

Copyright (c) 1994 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <sys/types.h>

#include "gomstring.h"
#include "gomtcl.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#include "stdafx.h"

int gomp_ScaleCommand(ClientData , Tcl_Interp *, int , const char **);

/*********************************************************************/
int gomp_ScaleCommand(ClientData clientdata, Tcl_Interp *interp,
                    int argc, const char **argv)
/*********************************************************************/
{
    static char   Text[BUFF_LEN];
#ifdef ENABLE_GRAPHICS
    static float  Scale;
#endif


/* #1   scale  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "gsc$ale")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

#ifdef ENABLE_GRAPHICS
/* display */
        if(gomp_StringMatch(Text , "disp$lay")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] != (char)NULL) {
                Scale = atof(Text);
                if(!gomp_ScaleDisplay(Scale , Scale , Scale))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
        }
        else 
#endif /* ENABLE_GRAPHICS */
        {
            gomp_PrintERROR("command 'gscale' not recognized");
            return(TCL_ERROR);
        }

    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'gscale' command not recognized");
    return(TCL_ERROR);

}

