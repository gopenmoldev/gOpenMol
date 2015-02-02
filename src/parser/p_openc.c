/*

Copyright (c) 1994 - 2005 by:
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
#include "model_file.h"
#include "plot.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_MOpenCommand(ClientData clientdata, Tcl_Interp *interp,
                    int argc, const char **argv)
/*********************************************************************/
{
    static char Text[BUFF_LEN];
    static int  Wstr;
    static int  Level;

/* mopen */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "mope$n")) {

        Wstr = 0;

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

        if(Text[0] == (char)NULL) {
            gomp_PrintERROR("model file name missing");
            return(TCL_ERROR);
        }
        Level = gomp_ReadOldModel(Text);

        if(!Level)
            return(TCL_OK);
        else
            return(TCL_ERROR);
    }

/*  E R R O R command not recognized         */
    gomp_PrintERROR("'mopen' command not recognized");
    return(TCL_ERROR);
}

