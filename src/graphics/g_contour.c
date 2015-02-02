/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <ctype.h>
#include <math.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#else
typedef double GLdouble;
#endif /* ENABLE_GRAPHICS */

#include "colouring.h"
#include "contour.h"
#include "gomfile.h"
#include "math_oper.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"

#include "stdafx.h"

static int ContourXClippingPlaneState     = 0;
static int ContourYClippingPlaneState     = 0;
static int ContourZClippingPlaneState     = 0;
static int ContourXYZClippingPlaneState1  = 0;
static int ContourXYZClippingPlaneState2  = 0;
static int ContourXYZClippingPlaneState3  = 0;

static GLdouble eqn1[4] = {1.0 , 0.0 , 0.0 , 0.0};
static GLdouble eqn2[4] = {0.0 , 1.0 , 0.0 , 0.0};
static GLdouble eqn3[4] = {0.0 , 0.0 , 1.0 , 0.0};
static GLdouble eqn4[4] = {0.0 , 0.0 , 0.0 , 0.0};
static GLdouble eqn5[4] = {0.0 , 0.0 , 0.0 , 0.0};
static GLdouble eqn6[4] = {0.0 , 0.0 , 0.0 , 0.0};

/* functions */

/************************************************************************/
/* NOTE: The main functionality of this function is to draw iso         */
/*       surfaces.                                                      */
/*       However, this function also stores surface triangles if        */
/*       surface method is set to saving state. Stored triangles are    */
/*       used by both the display function and the "show polygon"       */
/*       command.                                                       */
/************************************************************************/
/************************************************************************/
int gomp_PlotIsoSurf(void *userData,int Wstr,int drawFlags)
/************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int    i,k;
    static float  ColC[3];
    static float  ContVal;
    static float  SetRGBA[4];
    static const float *TransVec;
    static int    Which;
    static float  pA,pB,pC,pD;
    static float  ScaleMin,ScaleMax;

    if ( ! ( drawFlags & gom_PlotComplexElements) )
        return(-1);

    TransVec = gomp_GetTranslateArray();

    glEnable(GL_BLEND);
  
    glMatrixMode(GL_MODELVIEW); 

    glLineWidth((GLfloat)(gomp_GetContourLineWidth()+0.0));

    /* Delete stored polygon data if it can't be used anymore. */
    (void)gomp_DeleteInvalidatedPolygonData();

    for(k = 0 ; k < gomp_GetContoursDefined() ; k++) {

        if( gomp_GetContour2StructureMapping(k) != Wstr ) continue;

        if(gomp_GetCutPlanePlotStateXYZ(1)) {
            (void)gomp_GetCutPlaneXYZ(
                1 , &Which , &pA , &pB, &pC, &pD , &ScaleMin , &ScaleMax);
            if(k == Which)
                (void)gomp_CalculateJPRisosurface2(
                    k , gomp_GetContourProjection(k) , pA , pB , pC , pD , 
                    0.0 , 1.0   , 1 , ScaleMin , ScaleMax);
        }
        if(gomp_GetCutPlanePlotStateXYZ(2)) {
            (void)gomp_GetCutPlaneXYZ(
                2 , &Which , &pA , &pB, &pC, &pD, &ScaleMin , &ScaleMax);
            if(k == Which)
                (void)gomp_CalculateJPRisosurface2(
                    k , gomp_GetContourProjection(k) , pA , pB , pC , pD , 
                    0.0 , 1.0   , 1 , ScaleMin , ScaleMax);
        }
        if(gomp_GetCutPlanePlotStateXYZ(3)) {
            (void)gomp_GetCutPlaneXYZ(
                3 , &Which , &pA , &pB, &pC, &pD, &ScaleMin , &ScaleMax);
            if(k == Which)
                (void)gomp_CalculateJPRisosurface2(
                    k , gomp_GetContourProjection(k) , pA , pB , pC , pD , 
                    0.0 , 1.0   , 1 , ScaleMin , ScaleMax);
        }

        for(i = 0 ; i < gomp_GetContourLevels(k) ; i++) {

            if(!ContourInfo[k].levels[i].Display) continue;

            if(ContourInfo[k].levels[i].AlphaBlend < 1.0)
                glBlendFunc(GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA);
            else
                glBlendFunc(GL_ONE , GL_ZERO);

            if(gomp_GetContourCullFace(k , i)) {
                glEnable(GL_CULL_FACE);
                glCullFace(GL_FRONT);
            }

            ContVal = ContourInfo[k].levels[i].ColVal;

            ColC[0] = ContourInfo[k].levels[i].RedC;
            ColC[1] = ContourInfo[k].levels[i].GreenC;
            ColC[2] = ContourInfo[k].levels[i].BlueC;

            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&ColC[0] , &ColC[1] , &ColC[2]);

            SetRGBA[0] = ColC[0];
            SetRGBA[1] = ColC[1];
            SetRGBA[2] = ColC[2];
            SetRGBA[3] = ContourInfo[k].levels[i].AlphaBlend;
            glColor4fv(SetRGBA);

            if ( gomp_GetContourDisplayType(k,i) == CONTOUR_TYPE_LINE ) {
                gomp_CalculateIsoCurves(
                    k,
                    ContVal,
                    ContourInfo[k].levels[i].ClipAxis,
                    ContourInfo[k].levels[i].ClipPos);
            }
            else {
                if(gomp_GetContourXClippingPlaneState()) {
                    glClipPlane(GL_CLIP_PLANE0 , eqn1);
                    glEnable(GL_CLIP_PLANE0);
                }
                if(gomp_GetContourYClippingPlaneState()) {
                    glClipPlane(GL_CLIP_PLANE1 , eqn2);
                    glEnable(GL_CLIP_PLANE1);
                }
                if(gomp_GetContourZClippingPlaneState()) {
                    glClipPlane(GL_CLIP_PLANE2 , eqn3);
                    glEnable(GL_CLIP_PLANE2);
                } 
                if(gomp_GetContourXYZClippingPlaneState(1)) {
                    glClipPlane(GL_CLIP_PLANE3 , eqn4);
                    glEnable(GL_CLIP_PLANE3);
                }
                if(gomp_GetContourXYZClippingPlaneState(2)) {
                    glClipPlane(GL_CLIP_PLANE4 , eqn5);
                    glEnable(GL_CLIP_PLANE4);
                }
                if(gomp_GetContourXYZClippingPlaneState(3)) {
                    glClipPlane(GL_CLIP_PLANE5 , eqn6);
                    glEnable(GL_CLIP_PLANE5);
                }
                gomp_CalculateJPRisosurface1(
                    k , gomp_GetContourProjection(k) ,
                    ContVal ,
                    ContourInfo[k].levels[i].AlphaBlend ,
                    ContourInfo[k].levels[i].ContSmooth ,
                    gomp_GetContourDisplayType(k,i),
                    gomp_GetContourProjectionMin(k,i) ,
                    gomp_GetContourProjectionMax(k,i) ,
                    ! ( drawFlags & gom_PlotRealTimeRotation ) ,
                    gomp_GetSurfaceMethod() ?
                    &ContourInfo[k].levels[i].polygons : NULL );
            }
            
            glDisable(GL_CULL_FACE);
        }
    }

    glDisable(GL_BLEND);
    if(gomp_GetContourXClippingPlaneState())
        glDisable(GL_CLIP_PLANE0);
    if(gomp_GetContourYClippingPlaneState())
        glDisable(GL_CLIP_PLANE1);
    if(gomp_GetContourZClippingPlaneState())
        glDisable(GL_CLIP_PLANE2);
    if(gomp_GetContourXYZClippingPlaneState(1))
        glDisable(GL_CLIP_PLANE3);
    if(gomp_GetContourXYZClippingPlaneState(2))
        glDisable(GL_CLIP_PLANE4);
    if(gomp_GetContourXYZClippingPlaneState(3))
        glDisable(GL_CLIP_PLANE5);
#endif /* ENABLE_GRAPHICS */
    
    return(0);
}


/****************************************************************************/
int gomp_SetContourClippingPlaneState( int Axis, int What)
/****************************************************************************/
{
    if(Axis == 1) {
        switch(What) {
        case ON:
            ContourXClippingPlaneState = 1;
            break;
        case OFF:
            ContourXClippingPlaneState = 0;
            break;
        }
    } else if(Axis == 2) {
        switch(What) {
        case ON:
            ContourYClippingPlaneState = 1;
            break;
        case OFF:
            ContourYClippingPlaneState = 0;
            break;
        }
    } else if(Axis == 3) {
        switch(What) {
        case ON:
            ContourZClippingPlaneState = 1;
            break;
        case OFF:
            ContourZClippingPlaneState = 0;
            break;
        }
    } else {
        gomp_PrintERROR("wrong axis value! Has to be 1,2 or 3");
        return(1);
    }

    return(0);
}

/****************************************************************************/
int gomp_SetContourXYZClippingPlaneState(int Which , int What)
/****************************************************************************/
{

    switch(What) {
    case ON:
        switch(Which) {
        case 1:
            ContourXYZClippingPlaneState1 = 1;
            break;
        case 2:
            ContourXYZClippingPlaneState2 = 1;
            break;
        case 3:
            ContourXYZClippingPlaneState3 = 1;
            break;
        }
        break;
    case OFF:
        switch(Which) {
        case 1:
            ContourXYZClippingPlaneState1 = 0;
            break;
        case 2:
            ContourXYZClippingPlaneState2 = 0;
            break;
        case 3:
            ContourXYZClippingPlaneState3 = 0;
            break;
        }
        break;
    }

    return(0);
}

/****************************************************************************/
int gomp_GetContourXClippingPlaneState()
/****************************************************************************/
{
    return(ContourXClippingPlaneState);
}
/****************************************************************************/
int gomp_GetContourXYZClippingPlaneState(int Which)
/****************************************************************************/
{
    switch(Which) {
    case 1:
        return(ContourXYZClippingPlaneState1);
        break;
    case 2:
        return(ContourXYZClippingPlaneState2);
        break;
    case 3:
        return(ContourXYZClippingPlaneState3);
        break;
    }

    return(0);
}
/****************************************************************************/
int gomp_GetContourYClippingPlaneState()
/****************************************************************************/
{
    return(ContourYClippingPlaneState);
}
/****************************************************************************/
int gomp_GetContourZClippingPlaneState()
/****************************************************************************/
{
    return(ContourZClippingPlaneState);
}

/****************************************************************************/
int gomp_SetContourClippingPlaneParameters( int Axis, char Direction, float Value)
/****************************************************************************/
{
    const float *sumxyz;

    sumxyz         = gomp_GetTranslateArray();

    if(Axis == 1) {
        switch(Direction) {
        case '+':
            eqn1[0] = -1.0;
            eqn1[1] =  0.0;
            eqn1[2] =  0.0;
            eqn1[3] = (Value - sumxyz[0]);
            break;
        case '-':
            eqn1[0] = +1.0;
            eqn1[1] =  0.0;
            eqn1[2] =  0.0;
            eqn1[3] = -(Value - sumxyz[0]);
            break;
        }
    } else if(Axis == 2) {
        switch(Direction) {
        case '+':
            eqn2[0] =  0.0;
            eqn2[1] = -1.0;
            eqn2[2] =  0.0;
            eqn2[3] = (Value - sumxyz[1]);
            break;
        case '-':
            eqn2[0] =  0.0;
            eqn2[1] = +1.0;
            eqn2[2] =  0.0;
            eqn2[3] = -(Value - sumxyz[1]);
            break;
        }
    } else if(Axis == 3) {
        switch(Direction) {
        case '+':
            eqn3[0] =  0.0;
            eqn3[1] =  0.0;
            eqn3[2] = -1.0;
            eqn3[3] = (Value - sumxyz[2]);
            break;
        case '-':
            eqn3[0] =  0.0;
            eqn3[1] =  0.0;
            eqn3[2] = +1.0;
            eqn3[3] = -(Value - sumxyz[2]);
            break;
        }
    } else {
        gomp_PrintERROR("wrong axis parameter! Has to be 1, 2 or 3 (x,y,z)");
        return(1);
    }
    return(0);
}

/****************************************************************************/
int gomp_SetContourXYZClippingPlaneParameters(int Which, float x1 , float y1 , float z1 , 
                                            float x2 , float y2 , float z2 ,
                                            float x3 , float y3 , float z3)
/****************************************************************************/
{
    const float *sumxyz;
    float  ValueA;
    float  ValueB;
    float  ValueC;
    float  ValueD;
    float  p1[3];
    float  p2[3];
    float  p3[3];

    sumxyz         = gomp_GetTranslateArray();

    p1[0] = x1 - sumxyz[0];
    p1[1] = y1 - sumxyz[1];
    p1[2] = z1 - sumxyz[2];

    p2[0] = x2 - sumxyz[0];
    p2[1] = y2 - sumxyz[1];
    p2[2] = z2 - sumxyz[2];

    p3[0] = x3 - sumxyz[0];
    p3[1] = y3 - sumxyz[1];
    p3[2] = z3 - sumxyz[2];


    if(gomp_CalcPlaneFrom3Points( p1 , p2 , p3 , 
                                &ValueA ,&ValueB  ,&ValueC  ,&ValueD))
        return(1);

    switch(Which) {
    case 1:
        eqn4[0] =  ValueA;
        eqn4[1] =  ValueB;
        eqn4[2] =  ValueC;
        eqn4[3] =  ValueD;
        break;
    case 2:
        eqn5[0] =  ValueA;
        eqn5[1] =  ValueB;
        eqn5[2] =  ValueC;
        eqn5[3] =  ValueD;
        break;
    case 3:
        eqn6[0] =  ValueA;
        eqn6[1] =  ValueB;
        eqn6[2] =  ValueC;
        eqn6[3] =  ValueD;
        break;
    }

    return(0);
}
/************************************************************************/
int   gomp_WriteContourClipplaneInfo2ModelFile(FILE *Model_f)
/************************************************************************/
{
    if(!gomp_GetContoursDefined()) /* no contours defined */
        return(0);

/* CONTOUR CLIPPLANE- tag */
    fprintf(Model_f , "[Contour Clipplane]\n");
    fprintf(Model_f,"%d %f %f %f %f\n",
            ContourXClippingPlaneState,
            eqn1[0],eqn1[1],eqn1[2],eqn1[3]);
    fprintf(Model_f,"%d %f %f %f %f\n",
            ContourYClippingPlaneState,
            eqn2[0],eqn2[1],eqn2[2],eqn2[3]);
    fprintf(Model_f,"%d %f %f %f %f\n",
            ContourZClippingPlaneState,
            eqn3[0],eqn3[1],eqn3[2],eqn3[3]);

    fprintf(Model_f,"%d %f %f %f %f\n",
            ContourXYZClippingPlaneState1,
            eqn4[0],eqn4[1],eqn4[2],eqn4[3]);
    fprintf(Model_f,"%d %f %f %f %f\n",
            ContourXYZClippingPlaneState2,
            eqn5[0],eqn5[1],eqn5[2],eqn5[3]);
    fprintf(Model_f,"%d %f %f %f %f\n",
            ContourXYZClippingPlaneState3,
            eqn6[0],eqn6[1],eqn6[2],eqn6[3]);

    return(0);
}
/************************************************************************/
int   gomp_ReadContourClipplaneInfoFromModelFile(FILE *Model_f)
/************************************************************************/
{
    char  InputText[BUFF_LEN];
    double d1,d2,d3,d4;

/* CONTOUR CLIPPLANE- tag */
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%d %lf %lf %lf %lf",
           &ContourXClippingPlaneState,&d1,&d2,&d3,&d4);
    eqn1[0]=d1; eqn1[1]=d2; eqn1[2]=d3; eqn1[3]=d4;
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%d %lf %lf %lf %lf",
           &ContourYClippingPlaneState,&d1,&d2,&d3,&d4);
    eqn2[0]=d1; eqn2[1]=d2; eqn2[2]=d3; eqn2[3]=d4;
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%d %lf %lf %lf %lf",
           &ContourZClippingPlaneState,&d1,&d2,&d3,&d4);
    eqn3[0]=d1; eqn3[1]=d2; eqn3[2]=d3; eqn3[3]=d4;

    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%d %lf %lf %lf %lf",
           &ContourXYZClippingPlaneState1,&d1,&d2,&d3,&d4);
    eqn4[0]=d1; eqn4[1]=d2; eqn4[2]=d3; eqn4[3]=d4;
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%d %lf %lf %lf %lf",
           &ContourXYZClippingPlaneState2,&d1,&d2,&d3,&d4);
    eqn5[0]=d1; eqn5[1]=d2; eqn5[2]=d3; eqn5[3]=d4;
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText, "%d %lf %lf %lf %lf",
        &ContourXYZClippingPlaneState3,&d1,&d2,&d3,&d4);
    eqn6[0]=d1; eqn6[1]=d2; eqn6[2]=d3; eqn6[3]=d4;

    return(0);
}
