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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <tcl.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>

#if defined(GLUT)
#include <GL/glut.h>
#else
#include <GL/glaux.h>
#endif
#endif /* ENABLE_GRAPHICS */

#include "colouring.h"
#include "drawscene.h"
#include "g_status.h"
#include "gommain.h"
#include "model_file.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "selection.h"
#include "tclutils.h"

#include "stdafx.h"

#define OPENMOL_LABEL "*** gOpenMol ***\nVersion 0.2 r960401"

#define STRUCTURE_SELECTION_OFF   0
#define STRUCTURE_SELECTION_ON    1

#if 0
static int GetSelectedStructureList(void);
#endif

/*********************  OPENMOL label **********************************/
/*static const char *OPENMOL_label =
"\n"
"************************************************************************\n"
"                    OpenMol graphical interface\n"
"\n"
"                    Copyright (c) 1995 - 2000 by:\n"
"        Leif Laaksonen , Center for Scientific Computing, ESPOO, \n"
"                              FINLAND\n"
"            Confidential unpublished property of Leif Laaksonen\n"
"                        All rights reserved\n"
"\n"
"\n"
"                       in collaboration with\n"
"\n"
"                 OpenMol Molecular Astrophysics Group\n"
"                Max-Planck-Institut  fuer Astrophysik\n"
"                        Garching, GERMANY\n"
"************************************************************************\n"
"\n";*/
/*********************  OPENMOL label **********************************/

void MakeLOGOlabelWidget();
void MakeOPENMOLlabelWidget(void);

/***************************************************************************/
int gomp_UpdateMolecStructList()
/***************************************************************************/
{
    int         Code;
    
    (void)gomp_PushSelectedStructure(STRUCTURE_SELECTION_ON);

/* check if there is a graphics display available */

    if(gomp_GetTermType() != GRAPHICS_AVAILABLE) 
        return(1);

    Code = Tcl_GlobalEval(gomp_GetTclInterp(), "lulUpdateStructureList");

    if(Code != TCL_OK) {
        gomp_PrintERROR("can't update structure list");
        return(1);
    }

    return(0);
}

/****************************************************************************/
int gomp_ResetgOpenMol()
/****************************************************************************/
{
    int Temp;

    Temp = gomp_DeleteMolecStructs();
    /* set black as default background colour */
    (void)gomp_SetBGColor( 0.0 , 0.0 , 0.0);
#ifdef ENABLE_GRAPHICS
    gomp_DrawScene ();
#endif
    return(Temp);
}

/****************************************************************************/
int  gomp_PutText2Info1(const char *Text)
/****************************************************************************/
{
    int  Code;
    char OutText[BUFF_LEN];

/* check if graphics is available */
    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        sprintf(OutText,"lulPutText2Info1 {%s}",Text);
        Code = Tcl_GlobalEval(gomp_GetTclInterp(), OutText);
        if(Code != TCL_OK) {
            gomp_PrintERROR("can't put text to status box #1");
            return(1);
        }
    }

    return(0);
}
/****************************************************************************/
int  gomp_PutChars2Info2(int Value)
/****************************************************************************/
{
    static int  OldValue;
    static char OutText[BUFF_LEN];
    static int  Code;

/* check if graphics is available */
    if(gomp_GetTermType() != GRAPHICS_AVAILABLE) {
        return(0);
    }

    if(OldValue == Value) return(0);

    switch(Value) {

    case 0:
        sprintf(OutText,"lulPushText2TimeFlies {}");
        break;
    case 1:
        sprintf(OutText,"lulPushText2TimeFlies {[*]}");
        break;
    case 2:
        sprintf(OutText,"lulPushText2TimeFlies {[*][*]}");
        break;
    case 3:
        sprintf(OutText,"lulPushText2TimeFlies {[*][*][*]}");
        break;
    case 4:
        sprintf(OutText,"lulPushText2TimeFlies {[*][*][*][*]}");
        break;
    case 5:
        sprintf(OutText,"lulPushText2TimeFlies {[*][*][*][*][*]}");
        break;
    case 6:
        sprintf(OutText,"lulPushText2TimeFlies {[*][*][*][*][*][*]}");
        break;
    case 7:
        sprintf(OutText,"lulPushText2TimeFlies {[*][*][*][*][*][*][*]}");
        break;
    case 8:
        sprintf(OutText,"lulPushText2TimeFlies {[*][*][*][*][*][*][*][*]}");
        break;
    case 9:
        sprintf(OutText,"lulPushText2TimeFlies {[*][*][*][*][*][*][*][*][*]}");
        break;
    case 10:
        sprintf(OutText,"lulPushText2TimeFlies {[*][*][*][*][*][*][*][*][*][*]}");
        break;
    default:
        return(0);
    }

    Code = Tcl_GlobalEval(gomp_GetTclInterp(), OutText);
    if(Code != TCL_OK) {
        gomp_PrintERROR("can't put text to time flies field");
        return(1);
    }

    OldValue = Value;

    return(0);
}
/****************************************************************************/
int  gomp_PutText2Info3(const char *Text)
/****************************************************************************/
{
    int  Code;
    char OutText[BUFF_LEN];


/* check if graphics is available */
    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        sprintf(OutText,"lulPutText2Info3 {%s}",Text);
        Code = Tcl_GlobalEval(gomp_GetTclInterp(), OutText);
        if(Code != TCL_OK) {
            gomp_PrintERROR("can't put text to status box #3");
            return(1);
        }
    }

    return(0);
}
