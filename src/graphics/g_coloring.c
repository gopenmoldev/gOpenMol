/*

Copyright (c) 1993 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved


Enhancements 2002 by:
Eero Häkkinen
*/

#include "maindefs.h"

#define YES 1
#define NO  0
#define RABS(a)    ( ( a ) > 0.0 ? (a) : -(a))

#include "gomstdio.h"
#include <math.h>

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

#include "stdafx.h"

typedef struct{
    float bgRED;         /* Back ground colour */
    float bgGREEN;
    float bgBLUE;
} Coloring;

static Coloring BGColoring  = { 0.0 , 0.0 , 0.0};
#if 0
static int      BGColoringActive = 0;
#endif

/***********************************************************************/
int gomp_SetBGColor(float CRed, float CGreen, float CBlue)
/***********************************************************************/
{
    BGColoring.bgRED    = CRed;
    BGColoring.bgGREEN  = CGreen;
    BGColoring.bgBLUE   = CBlue;

    return(0);
}

/***********************************************************************/
int gomp_GetBGColor(float *CRed, float *CGreen, float *CBlue)
/***********************************************************************/
{
    *CRed   = BGColoring.bgRED;
    *CGreen = BGColoring.bgGREEN;
    *CBlue  = BGColoring.bgBLUE;

    return(0);
}
#if 0
/***********************************************************************/
int gomp_SetBackgroundColour()
/***********************************************************************/
{
    BGColoringActive = 1;
/*
  (void)DisplayColourWidget(1);
*/
    return(0);
}
#endif
