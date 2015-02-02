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

#include "gommonitor.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "plot.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_FillCommand(ClientData clientdata, Tcl_Interp *interp,
                   int argc, const char **argv)
/*********************************************************************/
{
    static int  ITemp1,ITemp2;
    static char Text[BUFF_LEN];

/* #1   fill  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "fill")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
/* distance */
        if(gomp_StringMatch(Text , "dist$ance")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text , "arra$y") ||
               gomp_StringMatch(Text , "list")) {
                if(!gomp_FillDistanceSeries())
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "-dis$tance")) {
            if(!gomp_DeleteDistanceSeries())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* angle */
        else if(gomp_StringMatch(Text , "angl$e")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text , "arra$y") ||
               gomp_StringMatch(Text , "list")) {
                if(!gomp_FillAngleSeries())
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "-ang$le")) {
            if(!gomp_DeleteAngleSeries())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* torsion */
        else if(gomp_StringMatch(Text , "tors$ion")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            if(gomp_StringMatch(Text , "arra$y") ||
               gomp_StringMatch(Text , "list")) {
                if(!gomp_FillTorsionSeries())
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "-tor$sion")) {
            if(!gomp_DeleteTorsionSeries())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* structure */
        else if(gomp_StringMatch(Text , "stru$cture")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] == (char)NULL) {
                ITemp1 = 1;
            }
            else 
                ITemp1 = atoi(Text);
            if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("structure index out of range");
                return(TCL_ERROR);
            }
            ITemp1--;
            if(!gomp_IdentifyAtoms(ITemp1))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* atom */
        else if(gomp_StringMatch(Text , "atom")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] == (char)NULL) {
                gomp_PrintERROR("atom index is missing");
                return(TCL_ERROR);
            }
            else
                ITemp2 = atoi(Text);
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] == (char)NULL) {
                ITemp1 = 1;
            }
            else
                ITemp1 = atoi(Text);
            if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("structure number out of range");
                return(TCL_ERROR);
            }
            ITemp1--;
            if(ITemp2 < 1 || ITemp2 > gomp_GetNumAtomsInMolecStruct(ITemp1)) {
                gomp_PrintERROR("atom number is out of range");
                return(TCL_ERROR);
            }
            ITemp2--;
            if(!gomp_IdentifyAtom(ITemp1,ITemp2,NULL))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else {
            gomp_PrintERROR("'fill' command not recognized");
            return(TCL_ERROR);
        }

    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'fill' command not recognized");
    return(TCL_ERROR);

}

