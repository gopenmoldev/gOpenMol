/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved


Enhancements 2002 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
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

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include "cell.h"
#include "colouring.h"
#include "plot.h"

#include "stdafx.h"

static struct {
    float red;
    float green;
    float blue;
} CellColour = { 1.0 , 0.0 , 0.0};

static int PlotCell(float , float , float , float , float , float, 
                    float , float , float);

#ifdef ENABLE_GRAPHICS
/***********************************************************************/
int PlotCell(float SideX, float SideY , float SideZ ,
                float Xplace , float Yplace , float Zplace,
                float Red , float Green , float Blue)
/***********************************************************************/
{
    static float CR;
    static float CG;
    static float CB;


    glLineWidth((GLfloat)(gomp_GetCellLinewidth()+0.0));

    CR = CellColour.red;
    CG = CellColour.green;
    CB = CellColour.blue;

    if(!gomp_GetDisplayColourType())
        (void)gomp_RGB2Grayscale(&CR , &CG , &CB);

    glColor3f( CR , CG , CB);

    SideX = 0.5 * SideX;
    SideY = 0.5 * SideY;
    SideZ = 0.5 * SideZ;

/* set line width to at least 2 pixels */

/*     if(line_width < 2) linewidth(pixel2);*/

    glDisable(GL_LIGHTING);

    glPushMatrix();
    glTranslatef( Xplace , Yplace , Zplace);

    glBegin( GL_LINE_LOOP );
    glVertex3f(-SideX , -SideY , -SideZ);
    glVertex3f( SideX , -SideY , -SideZ);
    glVertex3f( SideX ,  SideY , -SideZ);
    glVertex3f(-SideX ,  SideY , -SideZ);
    glEnd();

    glBegin( GL_LINE_LOOP );
    glVertex3f(-SideX , -SideY ,  SideZ);
    glVertex3f( SideX , -SideY ,  SideZ);
    glVertex3f( SideX ,  SideY ,  SideZ);
    glVertex3f(-SideX ,  SideY ,  SideZ);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f( SideX ,  SideY , -SideZ);
    glVertex3f( SideX ,  SideY ,  SideZ);

    glVertex3f( SideX , -SideY , -SideZ);
    glVertex3f( SideX , -SideY ,  SideZ);


    glVertex3f(-SideX , -SideY , -SideZ);
    glVertex3f(-SideX , -SideY ,  SideZ);

    glVertex3f(-SideX ,  SideY , -SideZ);
    glVertex3f(-SideX ,  SideY ,  SideZ);
    glEnd();

/*     linewidth(line_width);*/

    glPopMatrix();

    return(0);
}
#endif /* ENABLE_GRAPHICS */
/***********************************************************************/
int gomp_PlotCellBox(void* userData,int Wstr,int drawFlags)
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS
    if ( ! ( drawFlags & gom_PlotSimpleElements ) || Wstr != 0 )
        return(-1);

    (void)PlotCell(gomp_GetCellA() , gomp_GetCellB() , gomp_GetCellC() ,
                   gomp_GetCellXtrans() , gomp_GetCellYtrans() , gomp_GetCellZtrans() ,
                      1.0 , 0.0 , 0.0);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/***********************************************************************/
int gomp_SetCellColour(float red , float green , float blue )
/***********************************************************************/
{
    CellColour.red   = red;
    CellColour.green = green;
    CellColour.blue  = blue;

    return(0);
}
/***********************************************************************/
int gomp_GetCellColour(float *red, float *green, float *blue)
/***********************************************************************/
{

    *red   = CellColour.red;
    *green = CellColour.green;
    *blue  = CellColour.blue;

    return(0);
}
