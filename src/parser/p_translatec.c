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
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_TranslateCommand(ClientData clientdata, Tcl_Interp *interp,
                        int argc, const char **argv)
/*********************************************************************/
{
    static char  Text[BUFF_LEN];
    static float a;
    static float b;
    static float c;

/* #1   translate ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "tran$slate")) {

        if(!gomp_GetNumMolecStructs()) {
            gomp_PrintERROR("no atom structures defined");
            return(TCL_ERROR);
        }

        a = 0.0;
        b = 0.0;
        c = 0.0;

/* #1.1 primary */
        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        if(gomp_StringMatch(Text , "prim$ary") || 
           gomp_StringMatch(Text , "disp$lay")) {

#ifdef ENABLE_GRAPHICS
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] != (char)NULL) {
                a = atof(Text);
            }

            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] != (char)NULL) {
                b = atof(Text);
            }

            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] != (char)NULL) {
                c = atof(Text);
            }

            gomp_Translate(a , b , c);
#endif /* ENABLE_GRAPHICS */
       
            return(TCL_OK);
           }

/* #1.1 selection */
        else if(gomp_StringMatch(Text , "sele$ction")) {

            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] != (char)NULL) {
                a = atof(Text);
            }

            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] != (char)NULL) {
                b = atof(Text);
            }

            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(Text[0] != (char)NULL) {
                c = atof(Text);
            }
       
            gomp_TranslateCoordinatesX( a , b , c);
            return(TCL_OK);
        }
    }

/*  E R R O R command not recognized         */
    gomp_PrintERROR("'translate' command not recognized");

    return(TCL_ERROR);
}

