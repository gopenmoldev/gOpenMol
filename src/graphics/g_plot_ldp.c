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
#include "gomstdio.h"
#include "gommath.h"
#include <string.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

#if defined(GLUT)
#include <GL/glut.h>
#endif

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include "cell.h"
#include "colors.h"
#include "drawscene.h"
#include "gomtext.h"
#include "gomwindow.h"
#include "ldp.h"
#include "molecoord.h"
#include "objseg.h"
#include "plot_prop.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"

#include "stdafx.h"

#define MOUSEXMAP(x)  ((gomp_GetNumLdpAtoms1() * ((x)))/(WSize[2]))
#define MOUSEYMAP(y)  ((gomp_GetNumLdpAtoms1() * (WSize[3] - (y)))/(WSize[3]))

#if defined(GLUT)
static void LDPMouseFunc(int , int , int , int );
#if 0
static void LDPMouseMotion(int, int);
static void LDPDisplayFunc(void);
#endif
#endif

#ifdef ENABLE_GRAPHICS
/*  drawit handles the 2-d plotting of the distance matrix  */
/***********************************************************************/
int gomp_CalculateLDParraydraw_ldp(int num1, int num2, float min1, float min2, 
                                 float min3, float min4,
                                 float max1, float max2,
                                 float max3, float max4,
                                 float rest)
/***********************************************************************/
{
    float vec[3],dist,st;
    char text[BUFF_LEN];
    static float tmp1,tmp2,tmp3;
    static int i,ii,j,jj,k;
    static int SaveWindowID;
    static const float *x;
    static const float *y;
    static const float *z;
    static const int *AtomList1;
    static const int *AtomList2;
    static int   Wstr;
    static int   mm;
    static float xboxl,yboxl,zboxl;
    static float boxl1,boxl2,boxl3;

#if defined(GLUT)
    SaveWindowID = glutGetWindow();
    if(gomp_GetWindowingStyle() == MULTI_WINDOWING) {
        j = 0;
        for(i = 0 ; i < gomp_GetNumDefinedWindows() ; i++) {
            if(gomp_GetWindowTypeFromStack(i) == LDP_WINDOW) {
                j = 1;
                break;
            }
        }
        if(j) {
            glutSetWindow(gomp_GetWindowIDFromStack(i));
        } else {
            glutInitWindowPosition( 100 , 100);
            glutInitWindowSize( 400 , 400 );
            glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
            sprintf(text,"(%d) : %s",gomp_GetNumDefinedWindows() + 1,
                    "LDP Window");
            i = glutCreateWindow(text);
            (void)gomp_PushToWindowStack(i , LDP_WINDOW , "LDP Window");
            glutDisplayFunc(gomp_DrawSceneCallback);
            glutReshapeFunc(gomp_Reshape);
            glutMouseFunc(LDPMouseFunc);
            gomp_WindowInit();
        }
    }
#endif

    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);

    glDisable(GL_LIGHTING);

    glPushMatrix();
    glLoadIdentity();
    glOrtho((GLdouble)1   , (GLdouble)num1, 
            (GLdouble)1   , (GLdouble)num2,
            (GLdouble)0.0 , (GLdouble)num1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(mm);


    st=0.5;

    Wstr      = 0;
    AtomList1 = gomp_GetLdpAtomList1();
    AtomList2 = gomp_GetLdpAtomList2();
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);

/* take the periodic boundaries into account */

    boxl1 = gomp_GetCellA();
    boxl2 = gomp_GetCellB();
    boxl3 = gomp_GetCellC();
    xboxl = 1. / boxl1;
    yboxl = 1. / boxl2;
    zboxl = 1. / boxl3;

    for( i=0   ; i < num1 ; i++) {
        for(j=i+1 ; j < num2 ; j++)   {

            jj = AtomList2[j];
            ii = AtomList1[i];

            tmp1 = (x[ii]-x[jj]);
            tmp2 = (y[ii]-y[jj]);
            tmp3 = (z[ii]-z[jj]);

            tmp1 = tmp1 - boxl1 * nearbyint(xboxl * tmp1);
            tmp2 = tmp2 - boxl2 * nearbyint(yboxl * tmp2);
            tmp3 = tmp3 - boxl3 * nearbyint(zboxl * tmp3);

            dist=sqrt( tmp1 * tmp1 + 
                       tmp2 * tmp2 + 
                       tmp3 * tmp3);


            glColor3fv(gomp_WHITEv);

            if(dist > min1 && dist < max1 ) glColor3fv(gomp_GREENv);

            if(dist >= min2 && dist < max2) glColor3fv(gomp_BLUEv);

            if(dist >= min3 && dist < max3) glColor3fv(gomp_REDv);

            if(dist >= min4 && dist < max4) glColor3fv(gomp_YELLOWv);

            if(dist >= rest ) glColor3fv(gomp_CYANv);


            glBegin(GL_QUADS);

            glNormal3f(0.0 , 0.0 ,1.0);
            vec[0]=(float)(i);
            vec[1]=(float)(j); 
            vec[2]= 0.0;
            glVertex3fv(vec);
            vec[0]=(float)(i+1);
            vec[1]=(float)(j  );
            glVertex3fv(vec);
            vec[0]=(float)(i+1);
            vec[1]=(float)(j+1);
            glVertex3fv(vec);
            vec[0]=(float)(i  );
            vec[1]=(float)(j+1);
            glVertex3fv(vec);

            glEnd();

        }
    }

    glColor3fv(gomp_BLACKv);

    k=num1/50;

    if( k > 0) {

        for(i=1 ; i <= k ; i++) {

            j=50*i;

            glBegin(GL_LINES); 
            glVertex3i(j, 1 , 1);
            glVertex3i(j, num1 , 1);
            glEnd();
        }
    }

    k=num2/50;

    if( k > 0) {

        for(i=1 ; i <= k ; i++) {

            j=50*i;

            /* OGLXXX glBegin: Use GL_LINES if only one line segment is desired. */
            glBegin(GL_LINES); 
            glVertex3i(1, j , 1);
            glVertex3i(num2, j , 1);
            glEnd();
        }
    }


    glColor3fv(gomp_WHITEv);
    glRasterPos3f(num1/2.1, 4.0*num2/10. , 0.0);
    sprintf(text,"         dist <  %4.2f  WHITE",min1);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_GREENv);
    glRasterPos3f(num1/2.1, 3.5*num2/10. , 0.0);
    sprintf(text," %4.2f  < dist <  %4.2f  GREEN",min1,max1);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_BLUEv);
    glRasterPos3f(num1/2.1, 3.0*num2/10. , 0.0);
    sprintf(text," %4.2f  < dist <  %4.2f  BLUE ",min2,max2);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_REDv);
    glRasterPos3f(num1/2.1, 2.5*num2/10. , 0.0);
    sprintf(text," %4.2f  < dist <  %4.2f  RED  ",min3,max3);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_YELLOWv);
    glRasterPos3f(num1/2.1, 2.0*num2/10. , 0.0);
    sprintf(text," %4.2f  < dist <  %4.2f  YELLOW",min4,max4);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_CYANv);
    glRasterPos3f(num1/2.1, 1.5*num2/10. , 0.0);
    sprintf(text," %4.2f  < dist           CYAN",rest);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);


    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(mm);

    glEnable(GL_LIGHTING);

#if defined(GLUT)
    glutSetWindow(SaveWindowID);
#endif

    return(0);
}
#if 0
/***********************************************************************/
int CalculateLDParraydraw_ldp3(int num1, int num2, float min1, float min2, 
                               float min3, float min4,
                               float max1, float max2,
                               float max3, float max4,
                               float rest)
/***********************************************************************/
{
    float vec[3],dist,st;
    char text[BUFF_LEN];
    static float tmp1,tmp2,tmp3;
    static int i,ii,j,jj,k;

    static const float *x;
    static const float *y;
    static const float *z;
    static const int *AtomList1;
    static const int *AtomList2;
    static int   Wstr;
    static int   mm;
    static float xboxl,yboxl,zboxl;
    static float boxl1,boxl2,boxl3;

    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);

    glDisable(GL_LIGHTING);

    glPushMatrix();
    glLoadIdentity();
    glOrtho( (GLdouble)1   , (GLdouble)num1, 
             (GLdouble)1   , (GLdouble)num2,
             -(GLdouble)num1, (GLdouble)num1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(mm);


    st=0.5;

    Wstr      = 0;
    AtomList1 = gomp_GetLdpAtomList1();
    AtomList2 = gomp_GetLdpAtomList2();
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);

/* take the periodic boundaries into account */

    boxl1 = gomp_GetCellA();
    boxl2 = gomp_GetCellB();
    boxl3 = gomp_GetCellC();
    xboxl = 1. / boxl1;
    yboxl = 1. / boxl2;
    zboxl = 1. / boxl3;

    for( i=0   ; i < num1 ; i++) {
        for(j=i+1 ; j < num2 ; j++)   {

            jj = AtomList2[j];
            ii = AtomList1[i];

            tmp1 = (x[ii]-x[jj]);
            tmp2 = (y[ii]-y[jj]);
            tmp3 = (z[ii]-z[jj]);

            dist=sqrt( tmp1 * tmp1 + 
                       tmp2 * tmp2 + 
                       tmp3 * tmp3);


            glColor3fv(gomp_WHITEv);

            if(dist > min1 && dist < max1 ) glColor3fv(gomp_GREENv);

            if(dist >= min2 && dist < max2) glColor3fv(gomp_BLUEv);

            if(dist >= min3 && dist < max3) glColor3fv(gomp_REDv);

            if(dist >= min4 && dist < max4) glColor3fv(gomp_YELLOWv);

            if(dist >= rest ) glColor3fv(gomp_CYANv);

/*
  {
  float col[4];
  glGetFloatv(GL_CURRENT_COLOR,col);
  printf("%d %d %f %f %f %f\n",ii,jj,dist,col[0],col[1],col[2]);
  }
*/
            /* OGLXXX
             * special cases for polygons:
             *  independant quads: use GL_QUADS
             *  independent triangles: use GL_TRIANGLES
             */

            glBegin(GL_QUADS);

            glNormal3f(0.0 , 0.0 ,1.0);
            vec[0]=(float)(i);
            vec[1]=(float)(j);
            vec[2]=(float)0.0;
            glVertex3fv(vec);
            vec[0]=(float)(i+1);
            vec[1]=(float)(j  );
            vec[2]=(float)0.0;
            glVertex3fv(vec);
            vec[0]=(float)(i+1);
            vec[1]=(float)(j+1);
            vec[2]=(float)0.0;
            glVertex3fv(vec);
            vec[0]=(float)(i  );
            vec[1]=(float)(j+1);
            vec[2]=(float)0.0;
            glVertex3fv(vec);

            glEnd();

        }
    }

    glColor3fv(gomp_BLACKv);

    k=num1/50;

    if( k > 0) {

        for(i=1 ; i <= k ; i++) {

            j=50*i;

            glBegin(GL_LINES); 
            glVertex3i(j, 1   , 0);
            glVertex3i(j, num1, 0);
            glEnd();
        }
    }

    k=num2/50;

    if( k > 0) {

        for(i=1 ; i <= k ; i++) {

            j=50*i;

            /* OGLXXX glBegin: Use GL_LINES if only one line segment is desired. */
            glBegin(GL_LINES); 
            glVertex3i(1,    j, 0);
            glVertex3i(num2, j, 0);
            glEnd();
        }
    }


    glColor3fv(gomp_WHITEv);
    glRasterPos3f(num1/2.1, 4.0*num2/10., 0.1);
    sprintf(text,"         dist <  %4.2f  WHITE",min1);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_GREENv);
    glRasterPos3f(num1/2.1, 3.5*num2/10., 0.1);
    sprintf(text," %4.2f  < dist <  %4.2f  GREEN",min1,max1);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_BLUEv);
    glRasterPos3f(num1/2.1, 3.0*num2/10., 0.1);
    sprintf(text," %4.2f  < dist <  %4.2f  BLUE ",min2,max2);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_REDv);
    glRasterPos3f(num1/2.1, 2.5*num2/10., 0.1);
    sprintf(text," %4.2f  < dist <  %4.2f  RED  ",min3,max3);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_YELLOWv);
    glRasterPos3f(num1/2.1, 2.0*num2/10., 0.1);
    sprintf(text," %4.2f  < dist <  %4.2f  YELLOW",min4,max4);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);

    glColor3fv(gomp_CYANv);
    glRasterPos3f(num1/2.1, 1.5*num2/10., 0.1);
    sprintf(text," %4.2f  < dist           CYAN",rest);
    /* OGLXXX charstr: check list numbering */
    /*
      glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
    */
    gomp_PrintString(text , gomp_PropertyFont);


    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(mm);

    glEnable(GL_LIGHTING);

    return(0);
}
#endif
#if defined(GLUT)
#if 0
/***********************************************************************/
void LDPDisplayFunc()
/***********************************************************************/
{
    if(gomp_DrawScene()) {
        gomp_PrintERROR("problems in making the LDP matrix");
    }
}
#endif
/****************************************************************************/
void LDPMouseFunc(int button , int state , int x, int y)
/****************************************************************************/
{
    int  WSize[4];
    const int *Aindex;
    int  Xm;
    int  Ym;
    const float *Xc;
    const float *Yc;
    const float *Zc;
    const float *sumxyz;
    char   Temp1[BUFF_LEN];
    char   Temp2[BUFF_LEN];
    char   Temp3[BUFF_LEN];
    char   Temp4[BUFF_LEN];
    char   Temp5[BUFF_LEN];
    char   Temp6[BUFF_LEN];

    if(button == GLUT_LEFT_BUTTON) {
        if(state == GLUT_DOWN) {
            glGetIntegerv(GL_VIEWPORT , WSize);
            Xm = MOUSEXMAP(x);
            Ym = MOUSEYMAP(y);
            if(Xm >= Ym) {
                (void)gomp_DelArrowSeg();
                gomp_DrawScene();
                return;
            }
            Xc = gomp_GetAtomXCoordPointer(0);
            Yc = gomp_GetAtomYCoordPointer(0);
            Zc = gomp_GetAtomZCoordPointer(0);
            sumxyz = gomp_GetTranslateArray();
            Aindex = gomp_GetLdpAtomList1();
            printf("%d %d (%d %d)\n",MOUSEXMAP(x),MOUSEYMAP(y),Aindex[MOUSEXMAP(x)],Aindex[MOUSEYMAP(y)]);

            sprintf(Temp1,"%f",Xc[Aindex[MOUSEXMAP(x)]] + sumxyz[0]);
            sprintf(Temp2,"%f",Yc[Aindex[MOUSEXMAP(x)]] + sumxyz[1]);
            sprintf(Temp3,"%f",Zc[Aindex[MOUSEXMAP(x)]] + sumxyz[2]);

            sprintf(Temp4,"%f",Xc[Aindex[MOUSEYMAP(y)]] + sumxyz[0]);
            sprintf(Temp5,"%f",Yc[Aindex[MOUSEYMAP(y)]] + sumxyz[1]);
            sprintf(Temp6,"%f",Zc[Aindex[MOUSEYMAP(y)]] + sumxyz[2]);

            (void)gomp_PushArrowStack(Temp1,Temp2,Temp3,
                                 Temp4,Temp5,Temp6,
                                 "0.5" , "red" , "",
                                 "");
            gomp_DrawScene();
        }
    }
}
#if 0
/****************************************************************************/
void LDPMouseMotion(int x, int y)
/****************************************************************************/
{
    printf("%d %d\n",x,y);
}
#endif

#endif
#endif /* ENABLE_GRAPHICS */
