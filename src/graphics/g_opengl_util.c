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
#include <math.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#ifndef WIN32
#include <GL/glx.h>
#endif
#endif /* ENABLE_GRAPHICS */

#if !defined(WIN32)
#include <X11/Intrinsic.h>
#endif

#include "gomfile.h"
#include "molecoord.h"
#include "molecule.h"
#include "printmsg.h"
#include "projview.h"

#include "stdafx.h"

/* If you change the value here do it also in 'maindefs.h'               */
#define HIT_RAD2  9.0   /* squared tolarence radius at hits (in pixels)  */

#define YES 1
#define NO  0
#define RABS(a)    ( ( a ) > 0.0 ? (a) : -(a))

#if 0
#if defined(WIN32)
static int QueryXExtension(void);
#else
static int QueryXExtension(Display *);
#endif

#ifdef ENABLE_GRAPHICS
static int OpenGLsupport(void);
#if defined(WIN32)
static struct {
    int Supported;
} OpenGLXavailable = { 1 };
#else
static struct {
    int Supported;
} OpenGLXavailable = { 0 };
#endif
   
static int MapObject2Window(int , int , int);
#endif
#endif

#ifdef ENABLE_GRAPHICS
static struct {
    float HitRadius2;                   /* hit radius (squared */
} Picking = { HIT_RAD2 };
#endif /* ENABLE_GRAPHICS */

/* ................................... */
static struct {
    float Near;
    float NearStep;
    float NearStepReset;
    float NearReset;
    float Far;
    float FarStep;
    float FarStepReset;
    float FarReset;
    float DistanceReset;
    float DistanceStep;
    float DistanceStepReset;
    float Angle;
    float Window;
    float WindowStep;
    float WindowStepReset;
} distPerspective = {10.0f        /* near                */,
                      0.1f        /* near step           */,
                      0.1f        /* near step reset     */,
                     10.0f        /* near reset          */,
                     11.0f        /* far                 */, 
                      0.1f        /* far step            */,
                      0.1f        /* far step reset      */,
                     11.0f        /* far reset           */,
                     10.0f        /* distance reset      */, 
                      0.1f        /* distance step       */, 
                      0.1f        /* distance step reset */,
                     45.0f        /* angle               */,
                      1.0f        /* window              */, 
                      0.1f        /* window step         */, 
                      0.1f        /* window step reset   */};

#if 0
static struct {
    float Xtrans;
    float Ytrans;
    float Ztrans;
} DisplayTranslation;

static int   SetDisplayTranslation(float , float , float);
static int   MoveDisplayTranslation(float , float , float);
static int   GetDisplayTranslation(float *, float *, float *);
static int   ResetDisplayTranslation(void);
static float GetDisplayTranslationX(void);
static float GetDisplayTranslationY(void);
static float GetDisplayTranslationZ(void);
#endif

#ifdef ENABLE_GRAPHICS
#if 0
/************************************************************************/
int MapObject2Window(int Wstr , int XC , int YC)
/************************************************************************/
{
    int    Atoms;
    int    i;
    int    Ret;
    const float *Xc;
    const float *Yc;
    const float *Zc;
    
    GLdouble       objx,objy,objz;
    GLdouble modelMatrix[16],projMatrix[16];
    GLint    viewport[4];
    GLdouble       winx,winy,winz; 

    glMatrixMode(GL_PROJECTION);
    glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
    glGetIntegerv(GL_VIEWPORT , viewport);
    
    Atoms   = gomp_GetNumAtomsInMolecStruct(Wstr);

    Xc      = gomp_GetAtomXCoordPointer(Wstr);
    Yc      = gomp_GetAtomYCoordPointer(Wstr);
    Zc      = gomp_GetAtomZCoordPointer(Wstr);

    for(i = 0 ; i < Atoms ; i++) {

        objx    = Xc[i];
        objy   = Yc[i];
        objz  = Zc[i];

        Ret = gluProject(objx , objy , objz ,
                         modelMatrix,
                         projMatrix,
                         viewport,
                         &winx,
                         &winy,
                         &winz);

        if(Ret == GL_FALSE) {
            gomp_PrintERROR("$Can't map the object coordinates to window ones");
            return(1);
        }

        {
            float Hit;
            Hit = sqrt((winx-(double)XC)*(winx-(double)XC) + (viewport[3] - winy-(double)YC)*(winy-(double)YC));
            printf("#%d %f\n",(i+1),Hit);
        }
/*
  printf(":%d: %f %f %f\n",(i+1),winx-(float)XC,winy-(float)YC,winz);
  printf("     %f %f %f %f %f %f\n",Xc[i],Yc[i],Zc[i],winx,winy,winz);
*/
    }

    return(0);
}
/************************************************************************/
#if defined(WIN32)
int QueryXExtension()
#else
int QueryXExtension(Display *CurrentDisplay)
#endif
/************************************************************************/
{
#if !defined(WIN32)
    Bool Answer;
    int  errorBase;
    int  eventBase;
#endif

#if defined(WIN32)
    OpenGLXavailable.Supported = 1;
    return(1);
#else
    Answer = glXQueryExtension(CurrentDisplay ,
                               &errorBase,
                               &eventBase);

    OpenGLXavailable.Supported = (int)Answer;
    return((int)Answer);
#endif
}

/************************************************************************/
int OpenGLsupport()
/************************************************************************/
{
    return(OpenGLXavailable.Supported);
}
#endif
/*

This function returns the atom index in the selected structure where
the hit occured.

If there is no hit the return value is negative (-1).

Leif Laaksonen 1994

*/

/************************************************************************/
int gomp_IdentifyAtomFromCoords(int Wstr , int XC , int YC)
/************************************************************************/
{
    static int    Atoms;
    static int    i;
    static int    Ret;
    static float  HitValue;
    static int    HitIndex;
    static const float *Xc;
    static const float *Yc;
    static const float *Zc;
    static float  Dx;
    static float  Dy;
    static float  Dist;
    static float  HitRad;
    static const float *RotMP;
    
    
    static GLdouble       objx,objy,objz;
    static GLdouble       modelMatrix[16],projMatrix[16];
    static GLint          viewport[4];
    static GLdouble       winx,winy,winz; 

    glPushMatrix();

    RotMP = gomp_GetSavedModelViewMatrixMT(Wstr);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {
    } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {
        gluLookAt(0.0 , 0.0 , 
                  (gomp_GetPerspectiveNear() + 0.5 * gomp_GetPerspectiveWindow()) ,
                  0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0);
    } else {
        gomp_PrintERROR("undefined Projection Transformation defined");
        glPopMatrix();
        return(1);
    }

    glMultMatrixf(RotMP);

    glMatrixMode(GL_PROJECTION);
    glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
    glGetIntegerv(GL_VIEWPORT , viewport);

    glPopMatrix();
    
    Atoms    = gomp_GetNumAtomsInMolecStruct(Wstr);

    Xc       = gomp_GetAtomXCoordPointer(Wstr);
    Yc      = gomp_GetAtomYCoordPointer(Wstr);
    Zc     = gomp_GetAtomZCoordPointer(Wstr);

    HitIndex = -1;
    HitValue = 1.e+30f;
    HitRad   = sqrt(Picking.HitRadius2);

    for(i = 0 ; i < Atoms ; i++) {

        objx    = (double)Xc[i];
        objy   = (double)Yc[i];
        objz  = (double)Zc[i];

        Ret = gluProject(objx , objy , objz ,
                         modelMatrix,
                         projMatrix,
                         viewport,
                         &winx,
                         &winy,
                         &winz);

        if(Ret == GL_FALSE) {
            gomp_PrintERROR("$Can't map the object coordinates to window ones");
            return(1);
        }

        Dx  = (float)(                        winx-(double)XC);
        Dy  = (float)((double)(viewport[3]) - winy-(double)YC);

        if((Dx > HitRad) || (Dy > HitRad)) continue;

        Dist = Dx * Dx + Dy * Dy;

        if(Dist < Picking.HitRadius2) {
            HitIndex = i;
            HitValue = sqrt(Dist);
            printf("You hit the atom: %5d (%f) '%s:%s(%d):%s(%d)'\n",
                   (i+1),HitValue,
                   gomp_GetAtomSegName(Wstr,i),
                   gomp_GetAtomResName(Wstr,i),gomp_GetAtomResNum1(Wstr,i),
                   gomp_GetAtomAtmName(Wstr,i),(i+1));
            return(i);
        }

    }

    return(-1);
}
#endif /* ENABLE_GRAPHICS */

/************************************************************************/
int gomp_SetPerspectiveNear(float Value)
/************************************************************************/
{
    distPerspective.Near = Value;

    return(0);
}
/************************************************************************/
int gomp_SetPerspectiveFar(float Value)
/************************************************************************/
{
    distPerspective.Far = Value;

    return(0);
}
/************************************************************************/
int gomp_SetPerspectiveStep(float Value)
/************************************************************************/
{
    distPerspective.DistanceStep = Value;

    return(0);
}

/************************************************************************/
int gomp_SetPerspectiveAngle(float Value)
/************************************************************************/
{
    distPerspective.Angle = Value;

    return(0);
}
/************************************************************************/
float gomp_GetPerspectiveNear(void)
/************************************************************************/
{
    return(distPerspective.Near);
}
/************************************************************************/
float gomp_GetPerspectiveFar(void)
/************************************************************************/
{
    return(distPerspective.Far);
}
/************************************************************************/
float gomp_GetPerspectiveAngle(void)
/************************************************************************/
{
    return(distPerspective.Angle);
}
/************************************************************************/
float gomp_GetPerspectiveStep(void)
/************************************************************************/
{
    return(distPerspective.DistanceStep);
}

/************************************************************************/
int gomp_SetPerspectiveWindow(float Value)
/************************************************************************/
{
    distPerspective.Window = Value;

    return(0);
}
/************************************************************************/
float gomp_GetPerspectiveWindow(void)
/************************************************************************/
{
    return(distPerspective.Window);

    return((float)0.0);
}
/************************************************************************/
int gomp_SetPerspectiveWindowStep(float Value)
/************************************************************************/
{
    distPerspective.WindowStep = Value;

    return(0);
}

/************************************************************************/
float gomp_GetPerspectiveWindowStep(void)
/************************************************************************/
{
    return(distPerspective.WindowStep);
}
/************************************************************************/
int gomp_SetPerspectiveDistanceStep(float Value)
/************************************************************************/
{
    distPerspective.DistanceStep = Value;

    return(0);
}

/************************************************************************/
float gomp_GetPerspectiveDistanceStep(void)
/************************************************************************/
{
    return(distPerspective.DistanceStep);
}

/************************************************************************/
int   gomp_ResetPerspectiveWindowAttributes()
/************************************************************************/
{
    distPerspective.DistanceStep = distPerspective.DistanceStepReset;
    distPerspective.WindowStep   = distPerspective.WindowStepReset;
    distPerspective.Near         = distPerspective.NearReset;
    distPerspective.Far          = distPerspective.FarReset;

    return(0);
}


/************************************************************************/
int   gomp_SetPerspectiveWindowAttributes(float Near , float Window)
/************************************************************************/
{
    (void)gomp_SetPerspectiveWindow(Window);
    (void)gomp_SetPerspectiveNear(Near);
/*      (void)gomp_SetPerspectiveFar(Near + Window + Window);*/
    (void)gomp_SetPerspectiveFar(Near + Window);

    return(0);
}

#ifdef ENABLE_GRAPHICS
/************************************************************************/
int   gomp_GetModelviewMatrix(float *Matrix)
/************************************************************************/
{
    glGetFloatv(GL_MODELVIEW_MATRIX, Matrix);

    return(0);
}
#endif /* ENABLE_GRAPHICS */

/************************************************************************/
int   gomp_PutModelviewMatrix(const float *Matrix)
/************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(Matrix);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
#ifdef ENABLE_GRAPHICS
/************************************************************************/
int   gomp_GetProjectionMatrix(float *Matrix)
/************************************************************************/
{
    glGetFloatv(GL_PROJECTION_MATRIX, Matrix);

    return(0);
}
#endif /* ENABLE_GRAPHICS */

/************************************************************************/
int   gomp_PutProjectionMatrix(const float *Matrix)
/************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    glMatrixMode(GL_PROJECTION);

    glLoadMatrixf(Matrix);

    glMatrixMode(GL_MODELVIEW);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
#if 0
/************************************************************************/
int   SetDisplayTranslation(float Xt, float Yt, float Zt)
/************************************************************************/
{
    DisplayTranslation.Xtrans = Xt;
    DisplayTranslation.Ytrans = Yt;
    DisplayTranslation.Ztrans = Zt;

    return(0);
}
/************************************************************************/
int   MoveDisplayTranslation(float Xt, float Yt, float Zt)
/************************************************************************/
{
    DisplayTranslation.Xtrans += Xt;
    DisplayTranslation.Ytrans += Yt;
    DisplayTranslation.Ztrans += Zt;

    return(0);
}

/************************************************************************/
int   GetDisplayTranslation(float *Xt, float *Yt, float *Zt)
/************************************************************************/
{
    *Xt = DisplayTranslation.Xtrans;
    *Yt = DisplayTranslation.Ytrans;
    *Zt = DisplayTranslation.Ztrans;

    return(0);
}

/************************************************************************/
float GetDisplayTranslationX()
/************************************************************************/
{
    return(DisplayTranslation.Xtrans);
}

/************************************************************************/
float GetDisplayTranslationY()
/************************************************************************/
{
    return(DisplayTranslation.Ytrans);
}

/************************************************************************/
float GetDisplayTranslationZ()
/************************************************************************/
{
    return(DisplayTranslation.Ztrans);
}

/************************************************************************/
int   ResetDisplayTranslation()
/************************************************************************/
{
    DisplayTranslation.Xtrans = 0.0;
    DisplayTranslation.Ytrans = 0.0;
    DisplayTranslation.Ztrans = 0.0;

    return(0);
}
#endif
/************************************************************************/
int   gomp_SetPerspectiveNearStep(float Value)
/************************************************************************/
{
    distPerspective.NearStep = Value;

    return(0);
}
/************************************************************************/
float gomp_GetPerspectiveNearStep()
/************************************************************************/
{

    return(distPerspective.NearStep);

}

/************************************************************************/
int   gomp_SetPerspectiveFarStep(float Value)
/************************************************************************/
{
    distPerspective.FarStep = Value;

    return(0);
}
/************************************************************************/
float gomp_GetPerspectiveFarStep()
/************************************************************************/
{

    return(distPerspective.FarStep);

}
/************************************************************************/
int   gomp_WriteDisplayAttributesGOM(FILE *Model_f)
/************************************************************************/
{

    fprintf(Model_f , "%f %f %f %f\n",
            distPerspective.Near,distPerspective.NearStep,
            distPerspective.NearStepReset,distPerspective.NearReset);

    fprintf(Model_f , "%f %f %f %f\n",
            distPerspective.Far,distPerspective.FarStep,
            distPerspective.FarStepReset,distPerspective.FarReset);

    fprintf(Model_f , "%f %f %f\n",
            distPerspective.DistanceReset,distPerspective.DistanceStep,
            distPerspective.DistanceStepReset);

    fprintf(Model_f , "%f %f %f\n",
            distPerspective.Window,distPerspective.WindowStep,
            distPerspective.WindowStepReset);
 
    return(0);
}
/************************************************************************/
int   gomp_ReadDisplayAttributesGOM(FILE *Model_f)
/************************************************************************/
{
    char InputText[BUFF_LEN];


    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%f %f %f %f",
           &distPerspective.Near,&distPerspective.NearStep,
           &distPerspective.NearStepReset,&distPerspective.NearReset);

    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%f %f %f %f",
           &distPerspective.Far,&distPerspective.FarStep,
           &distPerspective.FarStepReset,&distPerspective.FarReset);

    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%f %f %f",
           &distPerspective.DistanceReset,&distPerspective.DistanceStep,
           &distPerspective.DistanceStepReset);

    (void)gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%f %f %f",
           &distPerspective.Window,&distPerspective.WindowStep,
           &distPerspective.WindowStepReset);
 
    return(0);
}

