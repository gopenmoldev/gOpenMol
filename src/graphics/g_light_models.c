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

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif
#else
typedef float GLfloat;
#endif /* ENABLE_GRAPHICS */

#include "colouring.h"
#include "light_model.h"
#include "printmsg.h"

#include "stdafx.h"

#define  MAIN 1        /* main program */

#ifdef ENABLE_GRAPHICS
static GLfloat  mat_ambient[]    = { 0.2f, 0.2f, 0.2f, 1.0f};
static GLfloat  mat_diffuse[]    = { 0.8f, 0.8f, 0.8f, 1.0f};
static GLfloat  mat_emission[]   = { 0.0f, 0.0f, 0.0f, 1.0f};
#endif
static GLfloat  mat_specular[]   = { 0.8f, 0.8f, 0.8f, 1.0f};
static GLfloat  mat_shininess[]  = { 20.0f };
/*static GLfloat  no_mat[]         = { 0.0f , 0.0f , 0.0f , 1.0f};*/

static GLfloat  light_model_diffuse[]    = { 1.0f, 1.0f, 1.0f, 1.0f};
#ifdef ENABLE_GRAPHICS
static GLfloat  light_model_ambient[]    = { 0.0f, 0.0f, 0.0f, 1.0f};
static GLfloat  light_model_specular[]   = { 1.0f, 1.0f, 1.0f, 1.0f};
static GLfloat  light_model_position[]   = {  0.0f,  0.0f, 1.0f, 0.0};
static GLfloat  light_model_positionN[]  = {  0.0f,  1.0f, 1.0f, 0.0};
static GLfloat  light_model_positionNE[] = {  1.0f,  1.0f, 1.0f, 0.0};
static GLfloat  light_model_positionE[]  = {  1.0f,  0.0f, 1.0f, 0.0};
static GLfloat  light_model_positionSE[] = {  1.0f, -1.0f, 1.0f, 0.0};
static GLfloat  light_model_positionS[]  = {  0.0f, -1.0f, 1.0f, 0.0};
static GLfloat  light_model_positionSW[] = { -1.0f, -1.0f, 1.0f, 0.0};
static GLfloat  light_model_positionW[]  = { -1.0f,  0.0f, 1.0f, 0.0};
static GLfloat  light_model_positionNW[] = { -1.0f,  1.0f, 1.0f, 0.0};
#endif

/*static GLfloat  light_model_local_viewer = 1.0; */

static int      color_mapping_type = 0;   /* 0 == 1D texture mapping
                                             1 != Rainbow type of mapping by 
                                                  Eric Grosse */

/*
  Coding is:

  N (1)
  NW (8)           NE (2)
  W  (7)   C (0)   E  (3)
  SW (6)           SE (4)
  S (5)
*/
static struct {
    int Position;
    float DiffRed;
    float DiffGreen;
    float DiffBlue;
} LightProperties = { 0 , 1.0 , 1.0 , 1.0 };

#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int    gomp_SetUpLightMaterialModels()
/****************************************************************************/
{

/*------------------------------------
  Do the setting of lighting
  -------------------------------------*/

    glMaterialfv (GL_FRONT_AND_BACK,GL_AMBIENT,mat_ambient);
    glMaterialfv (GL_FRONT_AND_BACK,GL_DIFFUSE,mat_diffuse);
    glMaterialfv (GL_FRONT_AND_BACK,GL_EMISSION,mat_emission);
    glMaterialfv (GL_FRONT_AND_BACK,GL_SPECULAR,mat_specular);
    glMaterialfv (GL_FRONT_AND_BACK,GL_SHININESS,mat_shininess);

/*
  glMaterialfv (GL_FRONT_AND_BACK,GL_DIFFUSE,no_mat);
  glMaterialfv (GL_FRONT_AND_BACK,GL_SPECULAR,no_mat);
  glMaterialfv (GL_FRONT_AND_BACK,GL_EMISSION,no_mat);
*/

/*    glLightModelf (GL_LIGHT_MODEL_LOCAL_VIEWER,light_model_local_viewer); */

    glLightfv(GL_LIGHT0,GL_AMBIENT,light_model_ambient);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,light_model_diffuse);
    glLightfv(GL_LIGHT0,GL_SPECULAR,light_model_specular);
    glLightfv(GL_LIGHT0,GL_POSITION,light_model_position);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE , GL_TRUE);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);  

    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_DEPTH_TEST);
/*    glEnable(GL_DITHER);  */

    return(0);
}
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int    gomp_SetMaterialSpecularRed(float Value)
/****************************************************************************/
{
    mat_specular[0] = (GLfloat)Value;

#ifdef ENABLE_GRAPHICS
    glMaterialfv (GL_FRONT_AND_BACK,GL_SPECULAR,mat_specular);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/****************************************************************************/
int    gomp_SetMaterialSpecularGreen(float Value)
/****************************************************************************/
{
    mat_specular[1] = (GLfloat)Value;

#ifdef ENABLE_GRAPHICS
    glMaterialfv (GL_FRONT_AND_BACK,GL_SPECULAR,mat_specular);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/****************************************************************************/
int    gomp_SetMaterialSpecularBlue(float Value)
/****************************************************************************/
{
    mat_specular[2] = (GLfloat)Value;

#ifdef ENABLE_GRAPHICS
    glMaterialfv (GL_FRONT_AND_BACK,GL_SPECULAR,mat_specular);
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/****************************************************************************/
float gomp_GetMaterialSpecularRed()
/****************************************************************************/
{
    return(mat_specular[0]);
}
/****************************************************************************/
float gomp_GetMaterialSpecularGreen()
/****************************************************************************/
{
    return(mat_specular[1]);
}
/****************************************************************************/
float gomp_GetMaterialSpecularBlue()
/****************************************************************************/
{
    return(mat_specular[2]);
}

/****************************************************************************/
int gomp_SetMaterialShininess(float Value)
/****************************************************************************/
{
    mat_shininess[0] = Value;

#ifdef ENABLE_GRAPHICS
    glMaterialfv (GL_FRONT_AND_BACK,GL_SHININESS,mat_shininess);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/****************************************************************************/
float gomp_GetMaterialShininess()
/****************************************************************************/
{
    return(mat_shininess[0]);
}
/****************************************************************************/
int    gomp_SetColorMappingType(int Type)
/****************************************************************************/
{
    color_mapping_type = Type;

    return(0);

}
/****************************************************************************/
int    gomp_GetColorMappingType()
/****************************************************************************/
{
    return(color_mapping_type);
}

/****************************************************************************/
int    gomp_SetLightPosition(int Place)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    if(Place == 0) {
        glLightfv(GL_LIGHT0,GL_POSITION,light_model_position);
        LightProperties.Position = 0;
    } else if(Place == 1) {
        glLightfv(GL_LIGHT0,GL_POSITION,light_model_positionN);
        LightProperties.Position = 1;
    } else if(Place == 2) {
        glLightfv(GL_LIGHT0,GL_POSITION,light_model_positionNE);
        LightProperties.Position = 2;
    } else if(Place == 3) {
        glLightfv(GL_LIGHT0,GL_POSITION,light_model_positionE);
        LightProperties.Position = 3;
    } else if(Place == 4) {
        glLightfv(GL_LIGHT0,GL_POSITION,light_model_positionSE);
        LightProperties.Position = 4;
    } else if(Place == 5) {
        glLightfv(GL_LIGHT0,GL_POSITION,light_model_positionS);
        LightProperties.Position = 5;
    } else if(Place == 6) {
        glLightfv(GL_LIGHT0,GL_POSITION,light_model_positionSW);
        LightProperties.Position = 6;
    } else if(Place == 7) {
        glLightfv(GL_LIGHT0,GL_POSITION,light_model_positionW);
        LightProperties.Position = 7;
    } else if(Place == 8) {
        glLightfv(GL_LIGHT0,GL_POSITION,light_model_positionNW);
        LightProperties.Position = 8;
    }
#else
    if(Place >= 0 && Place <= 9) {
        LightProperties.Position = Place;
    }
#endif /* ! ENABLE_GRAPHICS */
    else {
        gomp_PrintERROR("wrong light position code (0 - 8)");
        return(1);
    }

    return(0);
}
/****************************************************************************/
int    gomp_GetLightPosition()
/****************************************************************************/
{
    return(LightProperties.Position);
}

/****************************************************************************/
int    gomp_SetLightDiffuseRed(float Value)
/****************************************************************************/
{
    LightProperties.DiffRed = Value;
    light_model_diffuse[0]  = Value;
#ifdef ENABLE_GRAPHICS
    glLightfv(GL_LIGHT0,GL_DIFFUSE,light_model_diffuse);
#endif /* ENABLE_GRAPHICS */
    return(0);
}
/****************************************************************************/
int    gomp_SetLightDiffuseGreen(float Value)
/****************************************************************************/
{
    LightProperties.DiffGreen = Value;
    light_model_diffuse[1]    = Value;
#ifdef ENABLE_GRAPHICS
    glLightfv(GL_LIGHT0,GL_DIFFUSE,light_model_diffuse);
#endif /* ENABLE_GRAPHICS */
    return(0);
}
/****************************************************************************/
int    gomp_SetLightDiffuseBlue(float Value)
/****************************************************************************/
{
    LightProperties.DiffBlue = Value;
    light_model_diffuse[2]   = Value;
#ifdef ENABLE_GRAPHICS
    glLightfv(GL_LIGHT0,GL_DIFFUSE,light_model_diffuse);
#endif /* ENABLE_GRAPHICS */
    return(0);
}
/****************************************************************************/
float  gomp_GetLightDiffuseRed()
/****************************************************************************/
{
    return(LightProperties.DiffRed);
}
/****************************************************************************/
float  gomp_GetLightDiffuseGreen()
/****************************************************************************/
{
    return(LightProperties.DiffGreen);
}
/****************************************************************************/
float  gomp_GetLightDiffuseBlue()
/****************************************************************************/
{
    return(LightProperties.DiffBlue);
}
