/*

Copyright (c) 1995 - 2004 by:
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

#if defined(WIN32)
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>

#if defined(GLUT)
#include <GL/glut.h>
#else
#include <GL/glaux.h>
#endif

#if !defined(WIN32)
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>
#include <X11/Intrinsic.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/glx.h>
#endif

#include <tcl.h>

#include "drawscene.h"
#include "gommain.h"
#include "gomproc.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#endif /* ENABLE_GRAPHICS */

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS

/*********************************************************************/
int gomp_DisplayCommand(ClientData clientdata, Tcl_Interp *interp,
                      int argc, const char **argv)
/*********************************************************************/
{
    static int   ITemp;
    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    static int   LoopDisplayState;
    static int   LoopDisplayValue;
    static float MsecsB,CsecsB;
    static float MsecsA,CsecsA;
    static int   Idiff;

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {

/* check to see if the event queue check is not to be done */
        strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);
        gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

/* an attemp to avoid some problems introduced in the complicated event handling */
        (void)gomp_SetMouseButtonState(OFF);

        LoopDisplayState = 0;

        ITemp = 0;
        if(Text[0] != (char)NULL) {
            if(gomp_StringMatch(Text , "noch$eck")) {
                gomp_DrawSceneCallback();
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text , "slee$p")) {
                
                LoopDisplayValue = atoi(Text1);
                if(LoopDisplayValue > 0)
                    LoopDisplayState = 1;
            }
            else if(gomp_StringMatch(Text , "chec$k")) {
            }
            else {
                ITemp = atoi(Text);
                if(ITemp > 0) {
                    gomp_DrawSceneCallback();
                    return(TCL_OK);
                }
            }
        }

        gomp_PeekMessageQueue();

        if(!LoopDisplayState) {
            gomp_DrawSceneCallback();
        } else {
            gomp_Get_cpu_secs(&MsecsB , &CsecsB);
            gomp_DrawSceneCallback();
            gomp_Get_cpu_secs(&MsecsA , &CsecsA);

            Idiff = (int)(1000.0 * (MsecsA + CsecsA - MsecsB - CsecsB));

            if(Idiff < LoopDisplayValue) 
                Tcl_Sleep(LoopDisplayValue - Idiff);
        }
    } else {
        gomp_PrintERROR("No display available");
        return(TCL_ERROR);
    }

    return(TCL_OK);
}
#else /* ENABLE_GRAPHICS */
extern int i;
#endif
