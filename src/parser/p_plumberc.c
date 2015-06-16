/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2001 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <sys/types.h>

#include "colouring.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "plot.h"
#include "molecstruct.h"
#include "plumber.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_PlumberCommand(ClientData clientdata, Tcl_Interp *interp,
                      int argc, const char **argv)
/*********************************************************************/
{
    static int    Type,ITemp,ITemp1,Glue;
    static char   Text[BUFF_LEN];
    static char   Text1[BUFF_LEN];
    static char   Text2[BUFF_LEN];
    static char   Text3[BUFF_LEN];
    static char   Text4[BUFF_LEN];
    static char   Text5[BUFF_LEN];
    static char   Text6[BUFF_LEN];
    static char   Text9[BUFF_LEN];
    static char   Text10[BUFF_LEN];
    static char   Text11[BUFF_LEN];
    static float  RedC;
    static float  GreenC;
    static float  BlueC;
    static float  Rad;
    static float  Width;
    static float  Thickness;

/* #1   plumber  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);

    if(gomp_StringMatch(Text , "plum$ber")) {

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

/* atoms */
        if(gomp_StringMatch(Text , "atom$s")) {

            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            if(strlen(Text4) == 0) {
                gomp_PrintERROR("radius value missing");
                return(TCL_ERROR);
            }
            Rad = atof(Text4);

            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            if(strlen(Text5) == 0) {
                gomp_PrintERROR("colour name (value) missing");
                return(TCL_ERROR);
            }
            if(gomp_ColourName2RGB(Text5 , &RedC , &GreenC , &BlueC)) {
                sprintf(Text,"can't resolve colour name '%s'",Text5);
                gomp_PrintERROR(Text5);
                return(TCL_ERROR);
            }
 
            gomp_CopyString(Text6,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            Type = CYLINDER_TYPE;
            if(gomp_StringMatch(Text6 , "cyli$nder") ||
               gomp_StringMatch(Text6 , "tube"))
                Type = CYLINDER_TYPE;
            else if(gomp_StringMatch(Text6 , "ribb$on")) 
                Type = RIBBON_TYPE;
            else if(gomp_StringMatch(Text6 , "shel$ix"))
                Type = SOLID_HELIX_TYPE;
            else if(gomp_StringMatch(Text6 , "fhel$ix"))
                Type = FLAT_HELIX_TYPE;
            else if(gomp_StringMatch(Text6 , "arro$w"))
                Type = ARROW_TYPE;
            else if(gomp_StringMatch(Text6 , "stra$nd"))
                Type = STRAND_TYPE;
            else if(gomp_StringMatch(Text6 , "trac$e"))
                Type = TRACE_TYPE;
            else {
                gomp_PrintERROR("unrecognized plumber type");
                return(TCL_ERROR);
            }

            if( Type == SOLID_HELIX_TYPE ||
                Type == FLAT_HELIX_TYPE ||
                Type == ARROW_TYPE ||
                Type == STRAND_TYPE ) {
                gomp_CopyString(Text9,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text10,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text9) == 0) {
                    sprintf(Text9,"2.0");
                    gomp_PrintMessage("using default width value 2.0");
                }
                if(strlen(Text10) == 0) {
                    sprintf(Text10,"0.5");
                    gomp_PrintMessage("using default thickness value 0.5");
                }
                Width = atof(Text9);
                Thickness = atof(Text10);
            }

            Glue = 0;
            do {
                gomp_CopyString(Text11,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_StringMatch(Text11 , "glue"))
                    Glue = 1;
            } while( strlen(Text11) != 0);

            (void)gomp_ParsePlumberList( Text1,
                                       Text2,
                                       Text3,
                                       RedC,
                                       GreenC,
                                       BlueC,
                                       Rad  ,
                                       Type ,
                                       Width,
                                       Thickness,
                                       Glue);

            return(TCL_OK);
        }
/* -atoms */
        else if(gomp_StringMatch(Text , "-ato$ms")) {
            gomp_SetPlumberDisplay(0);
            if(!gomp_DeletePlumbers())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* display */
        else if(gomp_StringMatch(Text , "disp$lay")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "on")) {
                if(!gomp_GetPlumberSets()) {
                    gomp_PrintERROR("no plumber sets are defined");
                    return(TCL_ERROR);
                }
                if(!gomp_SetPlumberDisplay(1))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "off")) {
                if(!gomp_SetPlumberDisplay(0))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("unrecognized display parameter");
                return(TCL_ERROR);
            }
        }
/* -display */
        else if(gomp_StringMatch(Text , "-dis$play")) {
            if(!gomp_SetPlumberDisplay(0))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "colo$r") ||
                gomp_StringMatch(Text , "colo$ur")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
        
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("colour name (value) missing");
                return(TCL_ERROR);
            }
            if(gomp_ColourName2RGB(Text1 , &RedC , &GreenC , &BlueC)) {
                sprintf(Text,"can't resolve colour name '%s'",Text1);
                gomp_PrintERROR(Text1);
                return(TCL_ERROR);
            }

            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text2) == 0) {
                ITemp = 1;
            }
            else
                ITemp = atoi(Text2);
            if(ITemp < 1 || ITemp > gomp_GetPlumberSets()) {
                gomp_PrintERROR("plumber number out of range");
                return(TCL_ERROR);
            } 
            ITemp--;
            gomp_SetPlumberRed  (ITemp, RedC  );
            gomp_SetPlumberGreen(ITemp, GreenC);
            gomp_SetPlumberBlue (ITemp, BlueC );
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "stru$cture")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
        
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("structure number missing");
                return(TCL_ERROR);
            }
            else
                ITemp1 = atoi(Text1);
            if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("structure number out of range");
                return(TCL_ERROR);
            }
            ITemp1--;

            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text2) == 0) {
                ITemp = 1;
            }
            else
                ITemp = atoi(Text2);
            if(ITemp < 1 || ITemp > gomp_GetPlumberSets()) {
                gomp_PrintERROR("plumber number out of range");
                return(TCL_ERROR);
            } 
            ITemp--;
            gomp_SetPlumberStructure(ITemp, ITemp1);
            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text , "dest$roy")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                ITemp = 1;
            }
            else
                ITemp = atoi(Text1);
            if(ITemp < 1 || ITemp > gomp_GetPlumberSets()) {
                gomp_PrintERROR("plumber number out of range");
                return(TCL_ERROR);
            } 
            ITemp--;
            if(!gomp_DeletePlumber(ITemp))
                return(TCL_OK);
            else
                return (TCL_ERROR);
        }
        else {
            gomp_PrintERROR("command 'plumber' not recognized");
            return(TCL_ERROR);
        }

    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'plumber' command not recognized");

    return(TCL_ERROR);

}

