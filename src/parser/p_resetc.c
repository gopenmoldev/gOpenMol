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

#include "gommain.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "molecule.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_ResetCommand(ClientData clientdata, Tcl_Interp *interp,
                    int argc, const char **argv)
/*********************************************************************/
{
    static char   Text[BUFF_LEN];

/* #1   reset  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);

    if(gomp_StringMatch(Text , "rese$t")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

/* view */
        if(gomp_StringMatch(Text , "view")) {
            if(!gomp_ResetView())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "gope$nmol")) {
            if(!gomp_ResetgOpenMol())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "atomc$olours") ||
                gomp_StringMatch(Text , "atomc$olors")) {
            if(!gomp_IdentifyAtomColoursAllStructures())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "conn$ectivity") ||
                gomp_StringMatch(Text , "conn$ections")) {
            if(!gomp_ResetAtmConn(1))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else {
            gomp_PrintERROR("command 'reset' not recognized");
            return(TCL_ERROR);
        }

    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'reset' command not recognized");
    return(TCL_ERROR);

}

