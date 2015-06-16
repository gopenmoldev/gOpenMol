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

#include "gomwindow.h"
#include "gommain.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "hardcopy.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#endif /* ENABLE_GRAPHICS */

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS

/*********************************************************************/
int gomp_HardcopyCommand(ClientData clientdata, Tcl_Interp *interp,
                       int argc, const char **argv)
/*********************************************************************/
{
    static char   Text[BUFF_LEN];
    static char   Text1[BUFF_LEN];
    static int    ITemp;

/* #1   hardcopy  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "hard$copy")) {

        if(gomp_GetTermType() != GRAPHICS_AVAILABLE) {
            gomp_PrintERROR("no graphics display available");
            return(TCL_ERROR);
        }

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

        ITemp = atoi(Text);
        if(gomp_SetWindow(ITemp)) return(1);

        gomp_PopWindow();

#if defined(WIN32)
        (void)gomp_PeekMessageQueue();
#endif

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

/* bmp */
        if(gomp_StringMatch(Text , "bmp")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("hardcopy file name is missing");
                return(TCL_ERROR);
            }
            if(!gomp_HardcopyBMP(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "post$script")) {
            strncpy(Text, gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN-1);
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("hardcopy file name is missing");
                return(TCL_ERROR);
            }
            if(!gomp_HardcopyPostScript(Text , Text1))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "rgb")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("hardcopy file name is missing");
                return(TCL_ERROR);
            }
            if(!gomp_HardcopyRGB(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "xwd")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("hardcopy file name is missing");
                return(TCL_ERROR);
            }
            if(!gomp_HardcopyXWD(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "tga") ||
                gomp_StringMatch(Text , "targ$a")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("hardcopy file name is missing");
                return(TCL_ERROR);
            }
            if(!gomp_HardcopyTGA(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text , "jpg") ||
                gomp_StringMatch(Text , "jpeg")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("hardcopy file name is missing");
                return(TCL_ERROR);
            }
            if(!gomp_HardcopyJPEG(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else {
            gomp_PrintERROR("command 'hardcopy' not recognized");
            return(TCL_ERROR);
        }

    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'hardcopy' command not recognized");

    return(TCL_ERROR);

}
#else /* ENABLE_GRAPHICS */
extern int i;
#endif
