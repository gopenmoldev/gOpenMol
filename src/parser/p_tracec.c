/*

Copyright (c) 1996 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2004 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <tcl.h>

#include "gomstring.h"
#include "gomtcl.h"
#include "printmsg.h"
#include "tclutils.h"
#include "trace.h"

#include "stdafx.h"

#define TRACE_ON  1
#define TRACE_OFF 0

/*********************************************************************/
int gomp_TraceCommand(ClientData clientdata, Tcl_Interp *interp,
                    int argc, const char **argv)
/*********************************************************************/
{
    static char Text[BUFF_LEN];
    static char Text1[BUFF_LEN];
    static char Text2[BUFF_LEN];
    static char Text3[BUFF_LEN];
    static char Text4[BUFF_LEN];
    static char Text5[BUFF_LEN];
    static int iappend;
    static FILE *Trace_p;

/* fill command */

    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);

    if(gomp_StringMatch(Text , "ptra$ce")) {

        gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

        if(gomp_StringMatch(Text1 , "-atom$s")) {
            (void)gomp_SetDisplayTraceAtoms(TRACE_OFF);
            if(!gomp_DeleteTrace())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        if(gomp_StringMatch(Text1 , "writ$e")) {
            strncpy(Text2,
                    gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN-1);
            if(gomp_StringMatch(Text2 , "prob$esurf") ||
               gomp_StringMatch(Text2 , "prob$surf")) {
                strncpy(Text3,
                        gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN-1);
                if(Text3[0] == '\0') {
                    gomp_PrintERROR("file name is missing");
                    return(TCL_ERROR);
                }

                Trace_p = fopen(Text3 , "w");
                if(Trace_p == (FILE *)NULL) {
                    gomp_PrintERROR("can't open trace file");
                    return(TCL_ERROR);
                }

                if(!gomp_WriteAtomTrace(PROBESURF_INPUT , Trace_p)) {
                    fclose(Trace_p);
                    return(TCL_OK);
                }
                else {
                    fclose(Trace_p);
                    return(TCL_ERROR);
                }
            }
            else {
                gomp_PrintERROR("unrecognized 'ptrace write' command");
                return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text1 , "atom$s")) {

            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

            iappend = 0;
            if(gomp_StringMatch(Text5 , "appe$nd")) 
                iappend = 1;

            if(gomp_TraceAtoms(Text2 , Text3 , Text4 , iappend))
                return(TCL_ERROR);

            return(TCL_OK);
        }
        else if(gomp_StringMatch(Text1 , "disp$lay")) {
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
            if(gomp_StringMatch(Text2 , "on")) {
                (void)gomp_SetDisplayTraceAtoms(TRACE_ON);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text2 , "off")) {
                (void)gomp_SetDisplayTraceAtoms(TRACE_OFF);
                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("unrecognized 'ptrace display' option");
                return(TCL_ERROR);
            }
        }
        else {
            gomp_PrintERROR("unrecognized 'ptrace atoms' command");
            return(TCL_ERROR);
        }
    }
    else {
        gomp_PrintERROR("unrecognized 'ptrace' command");
        return(TCL_ERROR);
    }

    return(TCL_OK);
}


