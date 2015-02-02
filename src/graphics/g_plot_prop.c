/*

Copyright (c) 1995 - 2004 by:
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
#include "gommath.h"
#include <tcl.h>

#if defined(WIN32)
#include <windows.h>
#undef FILE_EXECUTE
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include "gommonitor.h"
#include "gomstring.h"
#include "gomtext.h"
#include "measure.h"
#include "molecoord.h"
#include "plot_prop.h"
#include "plot.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

static int DefineDrawLineType(int);

char gomp_PropertyFont[BUFF_LEN] = "*";

/**************************************************************************/
int gomp_PlotMonitorDistance(void* userData,int Wstr,int drawFlags)
/**************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static const float *x,*y,*z;
    static int    i,j,k;
    static float  col_vec[3];
    static float  vec[3];
    static float  dist_lengt,xpart,ypart,zpart;
    static int    Loop;
    static int    Loop1;
    static char   text[BUFF_LEN];
    static const char *TValue;
    static char   font[BUFF_LEN];
    static char   prec[BUFF_LEN];
    static int    dist_len;
    static const int *dist_vec;

    if ( ! ( drawFlags & gom_PlotSimpleElements ) || Wstr != 0 )
        return(-1);

    dist_len = gomp_GetDistMonitorSamples();
    dist_vec = gomp_GetDistMonitorSamplesList();

    TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomAtomLabelFont",TCL_GLOBAL_ONLY);

    if(TValue) {
        gomp_CopyString(font,TValue,BUFF_LEN);
    } else {
        gomp_CopyString(font,"BITMAP_HELVETICA_12",BUFF_LEN);
    }

    TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomMonitorDistancePrecision",TCL_GLOBAL_ONLY);

    if(TValue) {
        sprintf(prec,"%s",TValue);
    } else {
        strncpy(prec,"%.2f",strlen("%.2f"));
    }

/* put dotted line on  */
    glEnable(GL_LINE_STIPPLE);
    glDisable(GL_LIGHTING);
    glLineWidth((GLfloat)1.0);
/* ..................  */

    x = gomp_GetAtomXCoordPointer(Wstr);
    y = gomp_GetAtomYCoordPointer(Wstr);
    z = gomp_GetAtomZCoordPointer(Wstr);

    Loop  = 0;
    Loop1 = 0;

    for(i = 0 ; i < dist_len ; i += 2) {

        if(DefineDrawLineType(gomp_GetDistMonitorType(Loop1))) {
            return(1);
        }

        if(gomp_GetDistMonitorColor(Loop1 , &col_vec[0],
                                  &col_vec[1],
                                  &col_vec[2])) {
            return(1);
        }
        glColor3fv(col_vec);

/*
  if(line_width < 2) linewidth(2);
  else
  linewidth(line_width);
*/
        glBegin(GL_LINES);

        j = dist_vec[Loop];
        vec[0] = x[j];
        vec[1] = y[j];
        vec[2] = z[j];
        glVertex3fv(vec);

        k = dist_vec[Loop + 1];
        vec[0] = x[k];
        vec[1] = y[k];
        vec[2] = z[k];
        glVertex3fv(vec);

        glEnd();

        Loop += 2;
        Loop1++;

/*    linewidth(line_width); */
    
        xpart = (x[j]-x[k]);
        ypart = (y[j]-y[k]);
        zpart = (z[j]-z[k]);
        dist_lengt = sqrt(xpart*xpart + ypart*ypart +zpart*zpart);

        sprintf(text,prec,dist_lengt);

        glRasterPos3f(xpart * 0.5 + x[k] + 0.1,
                      ypart * 0.5 + y[k] , 
                      zpart * 0.5 + z[k]);
        gomp_PrintString(text , font);

    }

/* put default line style back */
    glDisable(GL_LINE_STIPPLE);
/* ........................... */
#endif /* ENABLE_GRAPHICS */

    return(0);
}


/**************************************************************************/
int gomp_PlotMonitorAngle(void* userData,int Wstr,int drawFlags)
/**************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static const float *x,*y,*z;
    static int   i,Loop,Loop1;
    static int   angi,angii,angiii;
    static float angle;
    static float ang_col[3];
    static float ang_c[9];
    static char  text[BUFF_LEN];
    static float xpart,ypart,zpart;
    static const char *TValue;
    static char  font[BUFF_LEN];
    static char  prec[BUFF_LEN];
    static int   ang_len;
    static const int *ang_vec;

    if ( ! ( drawFlags & gom_PlotSimpleElements ) || Wstr != 0 )
        return(-1);

    ang_len = gomp_GetAngMonitorSamples();
    ang_vec = gomp_GetAngMonitorSamplesList();

    TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomAtomLabelFont",TCL_GLOBAL_ONLY);

    if(TValue) {
        gomp_CopyString(font,TValue,BUFF_LEN);
    } else {
        gomp_CopyString(font,"BITMAP_HELVETICA_12",BUFF_LEN);
    }

    TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomMonitorAnglePrecision",TCL_GLOBAL_ONLY);

    if(TValue) {
        sprintf(prec,"%s",TValue);
    } else {
        strncpy(prec,"%.2f",strlen("%.2f"));
    }

    x = gomp_GetAtomXCoordPointer(Wstr);
    y = gomp_GetAtomYCoordPointer(Wstr);
    z = gomp_GetAtomZCoordPointer(Wstr);

/* put dashed line on  */
    glEnable(GL_LINE_STIPPLE);
    glDisable(GL_LIGHTING);
    glLineWidth((GLfloat)1.0);
/* ..................  */

    Loop  = 0;
    Loop1 = 0;

    for(i = 0 ; i < ang_len ; i += 3) {

        if(DefineDrawLineType(gomp_GetAngMonitorType(Loop1))) {
            return(1);
        }

        if(gomp_GetAngMonitorColor(Loop1 , &ang_col[0],
                                 &ang_col[1],
                                 &ang_col[2])) {
            return(1);
        }
 
        glColor3fv(ang_col);

        angi     = ang_vec[Loop  ];
        angii   = ang_vec[Loop+1];
        angiii = ang_vec[Loop+2];

        glBegin(GL_LINE_STRIP);

        ang_c[0] = x[angi];
        ang_c[1] = y[angi];
        ang_c[2] = z[angi];

        glVertex3fv(&ang_c[0]);

        ang_c[3] = x[angii];
        ang_c[4] = y[angii];
        ang_c[5] = z[angii];

        glVertex3fv(&ang_c[3]);

        ang_c[6] = x[angiii];
        ang_c[7] = y[angiii];
        ang_c[8] = z[angiii];

        glVertex3fv(&ang_c[6]);

        glEnd();

        Loop += 3;
        Loop1++;

/*     linewidth(line_width);*/

        gomp_BondAngle(ang_c[0] , ang_c[1] , ang_c[2] ,
                     ang_c[3] , ang_c[4] , ang_c[5] ,
                     ang_c[6] , ang_c[7] , ang_c[8] , &angle);

        angle = angle * 180.0/M_PI;


        xpart = (x[angi]-x[angiii]);
        ypart = (y[angi]-y[angiii]);
        zpart = (z[angi]-z[angiii]);

        sprintf(text,prec,angle); 

        glRasterPos3f(xpart * 0.5 + x[angiii] , 
                      ypart * 0.5 + y[angiii] ,
                      zpart * 0.5 + z[angiii]);
        gomp_PrintString(text , font);

    }

    glDisable(GL_LINE_STIPPLE);
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/**************************************************************************/
int gomp_PlotMonitorTorsion(void *userData,int Wstr,int drawFlags)
/**************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static const float *x,*y,*z;
    static int   i,Loop,Loop1;
    static int   torsi,torsii,torsiii,torsiv;
    static float ang_col[3];
    static float angle;
    static float tors_c[12];
    static float xpart,ypart,zpart;
    static char  text[BUFF_LEN];
    static const char *TValue;
    static char  font[BUFF_LEN];
    static char  prec[BUFF_LEN];
    static int   tors_len;
    static const int *tors_vec;

    if ( ! ( drawFlags & gom_PlotSimpleElements ) || Wstr != 0 )
        return(-1);

    tors_len = gomp_GetTorsMonitorSamples();
    tors_vec = gomp_GetTorsMonitorSamplesList();

    TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomAtomLabelFont",TCL_GLOBAL_ONLY);

    if(TValue) {
        gomp_CopyString(font,TValue,BUFF_LEN);
    } else {
        gomp_CopyString(font,"BITMAP_HELVETICA_12",BUFF_LEN);
    }

    TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomMonitorTorsionPrecision",TCL_GLOBAL_ONLY);

    if(TValue) {
        sprintf(prec,"%s",TValue);
    } else {
        strncpy(prec,"%.2f",strlen("%.2f"));
    }

    x = gomp_GetAtomXCoordPointer(Wstr);
    y = gomp_GetAtomYCoordPointer(Wstr);
    z = gomp_GetAtomZCoordPointer(Wstr);

/* put dash/dot/dash line on  */
    glEnable(GL_LINE_STIPPLE);
    glDisable(GL_LIGHTING);
    glLineWidth((GLfloat)1.0);
/* ..................  */

    Loop  = 0;
    Loop1 = 0;

    for(i = 0 ; i < tors_len ; i += 4) {

        if(DefineDrawLineType(gomp_GetTorsMonitorType(Loop1))) {
            return(1);
        }

        if(gomp_GetTorsMonitorColor(Loop1 , &ang_col[0],
                                  &ang_col[1],
                                  &ang_col[2])) {
            return(1);
        }

        glColor3fv(ang_col); 


        torsi     = tors_vec[Loop  ];
        torsii   = tors_vec[Loop+1];
        torsiii = tors_vec[Loop+2];
        torsiv = tors_vec[Loop+3];

        glBegin(GL_LINE_STRIP);

        tors_c[0] = x[torsi];
        tors_c[1] = y[torsi];
        tors_c[2] = z[torsi];

        glVertex3fv(&tors_c[0]);

        tors_c[3] = x[torsii];
        tors_c[4] = y[torsii];
        tors_c[5] = z[torsii];

        glVertex3fv(&tors_c[3]);

        tors_c[6] = x[torsiii];
        tors_c[7] = y[torsiii];
        tors_c[8] = z[torsiii];

        glVertex3fv(&tors_c[6]);

        tors_c[9] = x[torsiv];
        tors_c[10] = y[torsiv];
        tors_c[11] = z[torsiv];

        glVertex3fv(&tors_c[9]);

        glEnd();

        Loop += 4;
        Loop1++;

/*    linewidth(line_width);*/

        gomp_floDihedAngle(tors_c[0] , tors_c[1 ] , tors_c[2 ] , 
                      tors_c[3] , tors_c[4 ] , tors_c[5 ] ,
                      tors_c[6] , tors_c[7 ] , tors_c[8 ] ,
                      tors_c[9] , tors_c[10] , tors_c[11] , &angle);

        angle = angle * 180./M_PI;


        xpart = (x[torsii]-x[torsiii]);
        ypart = (y[torsii]-y[torsiii]);
        zpart = (z[torsii]-z[torsiii]);

        sprintf(text,prec,angle); 

        glRasterPos3f(xpart * 0.5 + x[torsiii], 
                      ypart * 0.5 + y[torsiii], 
                      zpart * 0.5 + z[torsiii]);
        gomp_PrintString(text , font);
/*
  cmov(xpart/2.+x[torsiii] , ypart/2.+y[torsiii] , zpart/2.+z[torsiii]);
  charstr(text);
*/
/*
  retv = show_text_world(text,txt_param.txfont,txt_param.txcolor,
  txt_param.txtsize,
  xpart/2.+x[torsiii] , 
  ypart/2.+y[torsiii] , 
  zpart/2.+z[torsiii] ,
  BACKBUFFER);
  }
*/
    }

    glDisable(GL_LINE_STIPPLE);
#endif /* ENABLE_GRAPHICS */

    return(0);
     
}

#ifdef ENABLE_GRAPHICS
/**************************************************************************/
int DefineDrawLineType(int type)
/**************************************************************************/
{

    switch(type) {

    case 1: /* distance type  '* * * *' */
        glLineStipple(1 , 0x0101);
        break;
    case 2: /* angle type     '- - - -' */
        glLineStipple(1 , 0x00FF);
        break;
    case 3: /* torsion type   '* - - * - -'*/
        glLineStipple(1 , 0x1C47);
        break;
    case 4: /* disable line stiple */
        glDisable(GL_LINE_STIPPLE);
        break;
    default:
        gomp_PrintERROR("undefinied line type defined");
        return(1);
        break;
    }
    return(0);
}
#endif /* ENABLE_GRAPHICS */
