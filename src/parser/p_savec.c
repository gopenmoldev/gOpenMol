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
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_MSaveCommand(ClientData clientdata, Tcl_Interp *interp,
                    int argc, const char **argv)
/*********************************************************************/
{
    static char Text[BUFF_LEN];
    static int  Wstr;

/* mopen */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "msav$e")) {

        Wstr = 0;

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

        if(strlen(Text) == 0) {
            gomp_PrintERROR("model file name missing");
            return(TCL_ERROR);
        }
        if(!gomp_WriteOldModel(Text))
            return(TCL_OK);
        else
            return(TCL_ERROR);
    }

/*  E R R O R command not recognized         */
    gomp_PrintERROR("'msave' command not recognized");
    return(TCL_ERROR);

}

