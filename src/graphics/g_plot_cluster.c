/*

Copyright (c) 1990 - 2004 by:
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
#include <string.h>
#include <stdlib.h>

#ifdef ENABLE_GRAPHICS
#if !defined(WIN32)
#include <X11/Intrinsic.h>
#endif

#if defined(IRIX)
#include <X11/StringDefs.h>
#endif
#include <X11/keysym.h>
#include <X11/cursorfont.h>

#if defined(WIN32)
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include <X11/X.h>
#include <X11/Xlib.h>
#endif /* ENABLE_GRAPHICS */

#include "cluster.h"
#include "colors.h"
#include "gomtext.h"
#include "printmsg.h"

#include "stdafx.h"

#define CLUSTER_ON  1
#define CLUSTER_OFF 0

#define MOUSEXMAP(x)  ((ClusterData.NumSets*((x)-Wxorg))/(Wxsize))
#define MOUSEYMAP(y)  ((ClusterData.NumSets*((y)-Wyorg))/(Wysize))

#ifdef ENABLE_GRAPHICS
/*  drawit handles the 2-d plotting of the cluster matrix  */
/***********************************************************************/
int gomp_PlotClusterMatrix(int num1, int num2,
                         float min1, float min2,float min3,float min4,
                         float max1, float max2,float max3,float max4,
                         float rest)
/***********************************************************************/
{
    static float  vec[2],dist,st;
    static char   text[BUFF_LEN];
    static int    i,j,k;
    static int    mm;
    static const float *ClusterArray;
    static char   PropertyFont[] = "*";


    if(!gomp_GetClusterStatus()) {
        (void)gomp_PrintERROR("No cluster data available for plotting");
        (void)gomp_SetDisplayCLUSTERmatrix(CLUSTER_OFF);
        return(1);
    }

    ClusterArray = gomp_GetClusterData();

/* Plot the distance matrix */ 

    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);

    glDisable(GL_LIGHTING);

    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(1, num1, 1, num2 );
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(mm);


    st=0.5;

    for(i =0   ; i < num1 - 1 ; i++) {
        for(j=i+1 ; j < num2     ; j++)   {

            dist = ClusterArray[i + j * (j - 1) / 2];

            glColor3fv(gomp_WHITEv);

            if(dist > min1 && dist < max1 ) glColor3fv(gomp_GREENv);

            if(dist >= min2 && dist < max2) glColor3fv(gomp_BLUEv);

            if(dist >= min3 && dist < max3) glColor3fv(gomp_REDv);

            if(dist >= min4 && dist < max4) glColor3fv(gomp_YELLOWv);

            if(dist >= rest ) glColor3fv(gomp_CYANv);

            glBegin(GL_QUADS);

            vec[0]=(float)(i);
            vec[1]=(float)(j); 
            glVertex2fv(vec);
            vec[0]=(float)(i+1);
            vec[1]=(float)(j  );
            glVertex2fv(vec);
            vec[0]=(float)(i+1);
            vec[1]=(float)(j+1);
            glVertex2fv(vec);
            vec[0]=(float)(i  );
            vec[1]=(float)(j+1);
            glVertex2fv(vec);

            glEnd();

        }
    }

    glColor3fv(gomp_BLACKv);

    k=num1/50;

    if( k > 0) {

        for(i=1 ; i <= k ; i++) {

            j=50*i;

            glBegin(GL_LINES); 
            glVertex2i(j, 1);
            glVertex2i(j, num1);
            glEnd();
        }
    }

    k=num2/50;

    if( k > 0) {

        for(i=1 ; i <= k ; i++) {

            j=50*i;

            glBegin(GL_LINES); 
            glVertex2i(1, j);
            glVertex2i(num2, j);
            glEnd();

        }
    }

    glColor3fv(gomp_WHITEv);
    glRasterPos2f(num1/2.1, 4.0*num2/10.);
    sprintf(text,"         RMSD <  %4.2f  WHITE",min1);
    gomp_PrintString(text , PropertyFont);
    /* OGLXXX charstr: check list numbering */
    glColor3fv(gomp_GREENv);
    glRasterPos2f(num1/2.1, 3.5*num2/10.);
    sprintf(text," %4.2f  < RMSD <  %4.2f  GREEN",min1,max1);
    gomp_PrintString(text , PropertyFont);
    /* OGLXXX charstr: check list numbering */
    glColor3fv(gomp_BLUEv);
    glRasterPos2f(num1/2.1, 3.0*num2/10.);
    sprintf(text," %4.2f  < RMSD <  %4.2f  BLUE ",min2,max2);
    gomp_PrintString(text , PropertyFont);
    /* OGLXXX charstr: check list numbering */
    glColor3fv(gomp_REDv);
    glRasterPos2f(num1/2.1, 2.5*num2/10.);
    sprintf(text," %4.2f  < RMSD <  %4.2f  RED  ",min3,max3);
    gomp_PrintString(text , PropertyFont);
    /* OGLXXX charstr: check list numbering */
    glColor3fv(gomp_YELLOWv);
    glRasterPos2f(num1/2.1, 2.0*num2/10.);
    sprintf(text," %4.2f  < RMSD <  %4.2f  YELLOW",min4,max4);
    gomp_PrintString(text , PropertyFont);
    /* OGLXXX charstr: check list numbering */
    glColor3fv(gomp_CYANv);
    glRasterPos2f(num1/2.1, 1.5*num2/10.);
    sprintf(text," %4.2f  < RMSD           CYAN",rest);
    gomp_PrintString(text , PropertyFont);

    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(mm);

    glEnable(GL_LIGHTING);

    return(0);
}
#endif /* ENABLE_GRAPHICS */
