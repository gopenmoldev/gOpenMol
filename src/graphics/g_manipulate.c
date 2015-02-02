/*

Copyright (c) 1993 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved


Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

#include <tcl.h>

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "selection.h"
#include "tclutils.h"

#include "stdafx.h"

static struct {
    int Rotation;
    int Translation;
} ActionState = { 1 , 0 };

#ifdef ENABLE_GRAPHICS
/***********************************************************************/
int gomp_Rotate( float Angle , float XAxis, float YAxis, float ZAxis)  
    /* own rotation driver */
/***********************************************************************/
{
    static float  RotMS[16];
    static const float *RotMP;
    static const float *Trans;
    static int    i;

    glMatrixMode(GL_MODELVIEW);

    for(i = 0 ; i < gomp_GetNumMolecStructs(); i++) {
        if(gomp_GetSelectedStructure(i) || !gomp_GetObjectCenterType()) {
            Trans = gomp_GetTranslateArrayMT(i);
            RotMP = gomp_GetSavedModelViewMatrixMT(i);

            glLoadIdentity();
            glTranslatef(Trans[0] , Trans[1] , Trans[2]);
            glRotatef(Angle , XAxis , YAxis , ZAxis);
            glTranslatef(-Trans[0] , -Trans[1] , -Trans[2]);
            glMultMatrixf(RotMP);
            glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);
            (void)gomp_SaveModelViewMatrixMT(i , RotMS);
        }
    }

    return(0);         /* end of rotation */

}
/***********************************************************************/
int gomp_Translate(float XAxis, float YAxis, float ZAxis)  
    /* own translation driver */
/***********************************************************************/
{
    static float  RotMS[16];
    static const float *RotMP;
    static const float *Trans;
    static float  NewTrans[3];
    static int    i;

    glMatrixMode(GL_MODELVIEW);

    for(i = 0 ; i < gomp_GetNumMolecStructs(); i++) {
        if(gomp_GetSelectedStructure(i) || !gomp_GetObjectCenterType()) {
            Trans       = gomp_GetTranslateArrayMT(i);
            NewTrans[0] = Trans[0] + XAxis;
            NewTrans[1] = Trans[1] + YAxis;
            NewTrans[2] = Trans[2] + ZAxis;
            (void)gomp_SaveTranslateArrayMT(i , NewTrans[0] , NewTrans[1] , NewTrans[2]);
            RotMP     = gomp_GetSavedModelViewMatrixMT(i);

            glLoadIdentity();
            glTranslatef(XAxis , YAxis , ZAxis);
            glMultMatrixf(RotMP);
            glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);
            (void)gomp_SaveModelViewMatrixMT(i , RotMS);
        }
    }

    return(0);         /* end of translation */


}
#endif /* ENABLE_GRAPHICS */

/***********************************************************************/
int gomp_GetRotationState()
/***********************************************************************/
{
    return(ActionState.Rotation);
}
/***********************************************************************/
int gomp_GetTranslationState()
/***********************************************************************/
{
    return(ActionState.Translation);
}

/***********************************************************************/
int gomp_SetRotationState(int Value)
/***********************************************************************/
{
    static const char *ITemp;

    ActionState.Rotation    = Value;

    if(Value) {
        ITemp = Tcl_SetVar(gomp_GetTclInterp(),"gomManipulationState","0",TCL_GLOBAL_ONLY);
        if(!ITemp) {
            gomp_PrintERROR("can't set tcl variable 'gomManipulationState'");
        }
        ActionState.Translation = 0;
    }
    else {
        ITemp = Tcl_SetVar(gomp_GetTclInterp(),"gomManipulationState","1",TCL_GLOBAL_ONLY);
        if(!ITemp) {
            gomp_PrintERROR("can't set tcl variable 'gomManipulationState'");
        }
        ActionState.Translation = 1;
    }

    return(0);
}
/***********************************************************************/
int gomp_SetTranslationState(int Value)
/***********************************************************************/
{
    static const char *ITemp;

    ActionState.Translation = Value;

    if(Value) {
        ITemp = Tcl_SetVar(gomp_GetTclInterp(),"gomManipulationState","1",TCL_GLOBAL_ONLY);
        if(!ITemp) {
            gomp_PrintERROR("can't set tcl variable 'gomManipulationState'");
        }
        ActionState.Rotation    = 0;
    }
    else {
        ITemp = Tcl_SetVar(gomp_GetTclInterp(),"gomManipulationState","0",TCL_GLOBAL_ONLY);
        if(!ITemp) {
            gomp_PrintERROR("can't set tcl variable 'gomManipulationState'");
        }
        ActionState.Rotation    = 1;
    }

    return(0);
}
