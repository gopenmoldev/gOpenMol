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

#include "coord_man.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "molecule.h"
#include "molecstruct.h"
#include "plot.h"
#include "printmsg.h"
#include "selection.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_SelectCommand(ClientData clientdata, Tcl_Interp *interp,
                     int argc, const char **argv)
/*********************************************************************/
{
    static int    i;
    static char   Text[BUFF_LEN];
    static char   Text1[BUFF_LEN];
    static char   Text2[BUFF_LEN];
    static char   Text3[BUFF_LEN];
    static char   Text4[BUFF_LEN];
    static int    ITemp;
    static int    Loop;
    static const char *Value;

/* #1   select  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "sele$ct")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

/* atoms */
        if(gomp_StringMatch(Text , "atom$s")) {
            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no structures defined");
                return(TCL_ERROR);
            }

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
        
            if(strlen(Text4) == 0) {
                ITemp = 1;
            }
            else {
                ITemp = atoi(Text4);
                if(ITemp < 1 || ITemp > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure index out of range");
                    return(TCL_ERROR);
                }
            }
            if(!gomp_ParseSelectionList((ITemp - 1), 
                                      Text1, 
                                      Text2, 
                                      Text3))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-ato$ms")) {
            if(!gomp_DeleteSelectionList())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* structures */
        if(gomp_StringMatch(Text , "stru$ctures")) {

            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no structures defined");
                return(TCL_ERROR);
            }

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            if(gomp_StringMatch(Text1 , "all")) {
                for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {
                    (void)gomp_ActivateSelectedStructure(i , 
                                                       STRUCTURE_SELECTION_ON);

                    sprintf(Text,"gomSelectStruct(%d)",(i+1));
                    Value = Tcl_SetVar(gomp_GetTclInterp(),Text,"1",TCL_GLOBAL_ONLY);
                    if(!Value) {
                        gomp_PrintERROR("can't set tcl variable 'gomSelectStruct'");
                        return(TCL_ERROR);
                    }
                }
                return(TCL_OK);
            }

            Loop = 0;

            while(strlen(Text1) != 0) {
                ITemp      = atoi(Text1);
                if(ITemp < 1 || ITemp > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("index is out of range");
                    return(TCL_ERROR);
                }
                (void)gomp_ActivateSelectedStructure(ITemp - 1 ,
                                                   STRUCTURE_SELECTION_ON);
                sprintf(Text,"gomSelectStruct(%d)",ITemp);
                Value = Tcl_SetVar(gomp_GetTclInterp(),Text,"1",TCL_GLOBAL_ONLY);
                if(!Value) {
                    gomp_PrintERROR("can't set tcl variable 'gomSelectStruct'");
                    return(TCL_ERROR);
                }
                Loop++;
                strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
            }
            return(TCL_OK);
        }
/* structures */
        if(gomp_StringMatch(Text , "-str$uctures")) {

            if(!gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("no structures defined");
                return(TCL_ERROR);
            }

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            if(gomp_StringMatch(Text1 , "all")) {
                for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {
                    (void)gomp_ActivateSelectedStructure(i , 
                                                       STRUCTURE_SELECTION_OFF);

                    sprintf(Text,"gomSelectStruct(%d)",(i+1));
                    Value = Tcl_SetVar(gomp_GetTclInterp(),Text,"0",TCL_GLOBAL_ONLY);
                    if(!Value) {
                        gomp_PrintERROR("can't set tcl variable 'gomSelectStruct'");
                        return(TCL_ERROR);
                    }
                }
                return(TCL_OK);
            }

            Loop = 0;

            while(strlen(Text1) != 0) {
                ITemp      = atoi(Text1);
                if(ITemp < 1 || ITemp > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("index is out of range");
                    return(TCL_ERROR);
                }
                (void)gomp_ActivateSelectedStructure(ITemp - 1,
                                                   STRUCTURE_SELECTION_OFF);
                sprintf(Text,"gomSelectStruct(%d)",ITemp);
                Value = Tcl_SetVar(gomp_GetTclInterp(),Text,"0",TCL_GLOBAL_ONLY);
                if(!Value) {
                    gomp_PrintERROR("can't set tcl variable 'gomSelectStruct'");
                    return(TCL_ERROR);
                }
                Loop++;
                strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
            }
            return(TCL_OK);
        }
        else {
            gomp_PrintERROR("command 'select' not recognized");
            return(TCL_ERROR);
        }

    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'select' command not recognized");
    return(TCL_ERROR);

}

