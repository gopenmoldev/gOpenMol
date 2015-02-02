/*

Copyright (c) 1993 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 by:
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

#include "drawscene.h"
#include "projview.h"
#include "stereo.h"

#include "stdafx.h"

#define YSTEREO        491    /* subfield hight in pixels */
#define YOFFSET_LEFT   532    /* YSTEREO + YBLANK */
#define XXMAXSCREEN     (glGetIntegerv(XXX_XPMAX, &gdtmp), gdtmp)
#define YYMAXSCREEN     (glGetIntegerv(XXX_YPMAX, &gdtmp), gdtmp)

#if 0
static struct {
    int  IsON;
    long MousePosition;
    long Monitor;
    short MatrixMode;
    long WinID;
    long SaveMainWindowID;
    float saveMV[4][4];
    float savePV[4][4];
    long  OLDXposition;
    long  OLDYposition;
    long  OLDXsize;
    long  OLDYsize;
    float dist;
    float eye;
} StereoSetting = { 0 , 0 , 0 , 0};
#endif

/*  structure to contain window information */

static struct {
    float Angle;       /* rotation angle */
    float Translate;   /* translation from the centre */
    int   Active;      /* == 0 off , != 0 on */
    int   Set;         /* == 0 reset , != do not reset */
} StereoPlot = {3.0 , 2.0 , 0 , 0}; 

/***********************************************************************/
int gomp_SetStereoPlotTranslate(float Translate)
/***********************************************************************/
{
    StereoPlot.Translate = Translate;

    return(0);
}
/***********************************************************************/
float gomp_GetStereoPlotTranslate()
/***********************************************************************/
{
    return(StereoPlot.Translate);

}

/***********************************************************************/
int gomp_SetStereoPlotAngle(float Angle)
/***********************************************************************/
{
    StereoPlot.Angle = Angle;

    return(0);
}
/***********************************************************************/
float gomp_GetStereoPlotAngle()
/***********************************************************************/
{
    return(StereoPlot.Angle);

}
/***********************************************************************/
int gomp_SetStereoPlotState(int State)
/***********************************************************************/
{
    StereoPlot.Active = State;

    return(0);
}
/***********************************************************************/
int gomp_GetStereoPlotState()
/***********************************************************************/
{
    return(StereoPlot.Active);

}

#if 0
/***********************************************************************/
int gomp_PlotStereoPair()
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS

    static float  RotMS[16];

    glMatrixMode(GL_MODELVIEW);

    glPushMatrix();
    glLoadMatrixf(gomp_GetSavedModelViewMatrix());
    glGetFloatv(GL_MODELVIEW_MATRIX, RotMS);
    glPopMatrix();


    gomp_Rotate(       StereoPlot.Angle     , 0.0 , 1.0 , 0.0);

    gomp_Translate(   -StereoPlot.Translate , 0.0 , 0.0);

    (void)gomp_UpdateScreen();

    glLoadMatrixf(RotMS);

    gomp_Rotate(       -StereoPlot.Angle     , 0.0 , 1.0 , 0.0 );
    gomp_Translate(     StereoPlot.Translate , 0.0 , 0.0);

    (void)gomp_UpdateScreen();

    glLoadMatrixf(RotMS);

    gomp_SaveModelViewMatrixMT(0 , RotMS);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
#endif
 
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/* Added by Sigismondo Boschi - CINECA - Italy:
   Stereo for High-end graphic machines, e.g. SGI Onyx2

   The high-end Stereo graphics is so called "quad buffer" or
   "stereo-in-a-window". It allows to have stereo and double-buffering
   without the need of splitting vertically the viewport - 
   and consequenty you do not need to have full-screen resolution.

   The basic steps are:

   1) to have a stereo-capable visual (GLUT_STEREO)
   2) draw the left eye in the BACK_LEFT drawing buffer:
   glDrawBuffer(GL_BACK_LEFT);
   3) draw the right eye in the BACK_RIGHT drawing buffer:
   glDrawBuffer(GL_BACK_RIGHT);
   3) put them in the FRONT_LEFT & FRONT_RIGHT by swapping once the buffers:
   glutSwapBuffers();
   If you need to plot "non-stereo" staff just plot it to the BACK buffer:
   glDrawBuffer(GL_BACK);
   It will be authomatically plotted in both the back buffers. In this way 
   you can have non-stereo images in a stereo-window.

*/
/***********************************************************************/
