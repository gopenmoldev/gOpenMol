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

#ifdef ENABLE_GRAPHICS
#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <sys/types.h>

#if defined(WIN32)
#include <windows.h>
#endif

#include "gomclipbrd.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "gomwindow.h"
#include "printmsg.h"
#include "tclutils.h"

#endif /* ENABLE_GRAPHICS */

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS

#define DIST_TYPE  0
#define ANG_TYPE   1
#define TORS_TYPE  2


/*********************************************************************/
int gomp_CopyCommand(ClientData clientdata, Tcl_Interp *interp,
                   int argc, const char **argv)
/*********************************************************************/
{
    static char Text[BUFF_LEN];
    static char Text1[BUFF_LEN];
    static char Text2[BUFF_LEN];
    static int  Type;
    static int  SaveWindowID;

/* #1   copy ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "copy")) {

        gomp_CopyString(Text ,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

        if(gomp_StringMatch(Text , "bitm$ap"))  {

            strncpy(Text ,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);

            SaveWindowID = gomp_GetWindowID();

/* default window is window == 1 */
            if(strlen(Text) == 0) {
                if(gomp_SetWindow(gomp_GetWindowIDFromStack(0))) {
                    gomp_SetWindow(SaveWindowID);
                    return(1);
                }
            } else {
                Type = atoi(Text) - 1;
                if(gomp_SetWindow(gomp_GetWindowIDFromStack(Type))) {
                    gomp_SetWindow(SaveWindowID);
                    return(1);
                }
            }

            Type = gomp_CopyBitmap2Clipboard();

            gomp_SetWindow(SaveWindowID);

            if(!Type)
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-bit$map")) {
            if(!gomp_DeleteClipboardData())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "atom$s")) {
            if(!gomp_CopyText2Clipboard("Dummy text"))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-ato$ms")) {
            if(!gomp_DeleteClipboardData())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "corr$elation")) {
            if(!gomp_CopyCorrelationArray2Clipboard())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-corr$elation")) {
            if(!gomp_DeleteClipboardData())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "rdf")) {
            if(!gomp_CopyRDFarray2Clipboard())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-rdf")) {
            if(!gomp_DeleteClipboardData())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "times$eries")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("time series type is missing");
                return(TCL_ERROR);
            }
            Type = -1;
            if(gomp_StringMatch(Text1 , "dist$ance")) {
                Type = DIST_TYPE;
            }
            else if(gomp_StringMatch(Text1 , "angl$e")) {
                Type = ANG_TYPE;
            }
            else if(gomp_StringMatch(Text1 , "tors$ion")) {
                Type = TORS_TYPE;
            }

            if(Type < 0) {
                gomp_PrintERROR("unrecognized type of list (distance, angle, torsion)");
                return(TCL_ERROR);
            }
            strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            if(strlen(Text2) == 0) {
                gomp_PrintERROR("time series index is missing");
                return(TCL_ERROR);
            }
            if(!gomp_CopyTimeseries2Clipboard(Type , atoi(Text2)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-times$eries")) {
            if(!gomp_DeleteClipboardData())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "text")) {
            if(!gomp_CopyText2Clipboard(argv[2]))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "-text")) {
            if(!gomp_DeleteClipboardData())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "msdi$splacement")) {
            if(!gomp_CopyMSDarray2Clipboard())
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else {
            gomp_PrintERROR("'copy' command not recognized");
            return(TCL_ERROR);
        }
    }

/*  E R R O R command not recognized         */
    return(TCL_ERROR);
}
#else /* ENABLE_GRAPHICS */
extern int i;
#endif

