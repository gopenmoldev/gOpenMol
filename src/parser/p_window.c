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

#include "gomstring.h"
#include "gomtcl.h"
#include "gomwindow.h"
#include "printmsg.h"
#include "tclutils.h"
#endif /* ENABLE_GRAPHICS */

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS

/*********************************************************************/
int gomp_WindowCommand(ClientData clientdata, Tcl_Interp *interp,
                     int argc, const char **argv)
/*********************************************************************/
{
    static int  WinID;
    static char Text[BUFF_LEN];
    static char Text1[BUFF_LEN];
    static char Text2[BUFF_LEN];
    static char Text3[BUFF_LEN];

/* #1   window ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "window")) {

        gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        if(Text3[0] == (char)NULL) {
            gomp_PrintERROR("window ID (number) missing");
            return(TCL_OK);
        }

        WinID = atoi(Text3);

        gomp_CopyString(Text ,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

        if(gomp_StringMatch(Text , "resi$ze"))  {

            if(gomp_SetWindow(WinID)) return(TCL_ERROR);

            strncpy(Text1 ,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
            strncpy(Text2 ,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);

            if(Text1[0] == (char)NULL || Text2[0] == (char)NULL) {
                gomp_PrintERROR("new window size parameter(s) missing");
                return(TCL_ERROR);
            }
            if(!gomp_ResizeWindow(atoi(Text1) , atoi(Text2)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* move window */
        else if(gomp_StringMatch(Text , "move")) {

            if(gomp_SetWindow(WinID)) return(TCL_ERROR);

            strncpy(Text1 ,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);
            strncpy(Text2 ,gomp_GetNextFromParserStack(argc,(const char **)NULL),
                    BUFF_LEN-1);

            if(Text1[0] == (char)NULL || Text2[0] == (char)NULL) {
                gomp_PrintERROR("new window size parameter(s) missing");
                return(TCL_ERROR);
            }
            if(!gomp_MoveWindow(atoi(Text1) , atoi(Text2)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* iconify window */
        else if(gomp_StringMatch(Text , "icon$ify")) {

            if(WinID == 0) {
/* iconify the Tcl control widget ... */
                if(TCL_OK != gomp_SendCommand2Parser("wm iconify .")) {
                    gomp_PrintERROR("can't iconify main control widget");
                    return(TCL_ERROR);
                }
            } else {
                if(gomp_SetWindow(WinID)) return(TCL_ERROR);
                gomp_IconifyWindow();
            }

            return(TCL_OK);
        }
/* deiconify window */
        else if(gomp_StringMatch(Text , "deic$onify")) {

            if(WinID == 0) {
/* iconify the Tcl control widget ... */
                if(TCL_OK != gomp_SendCommand2Parser("wm deiconify .")) {
                    gomp_PrintERROR("can't deiconify main control widget");
                    return(TCL_ERROR);
                }
            } else {
                if(gomp_SetWindow(WinID)) return(TCL_ERROR);
                gomp_DeIconifyWindow();
            }

            return(TCL_OK);
        }
/* fullscreen window */
        else if(gomp_StringMatch(Text , "full$screen")) {

            if(WinID == 0) {
/* can't fullscreen the Tcl control widget ... */
                gomp_PrintERROR("can't fullscreen the main Tk control widget");
                return(TCL_ERROR);
            } else {
                if(gomp_SetWindow(WinID)) return(TCL_ERROR);
                gomp_FullScreenWindow();
            }

            return(TCL_OK);
        }
        else {
            gomp_PrintERROR("'window' command not recognized");
            return(TCL_ERROR);
        }

    }

/*  E R R O R command not recognized         */
    return(TCL_OK);
}
#else /* ENABLE_GRAPHICS */
extern int i;
#endif
