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

#define RABS(a)    ( ( a ) > 0.0 ? (a) : -(a))

#include "gomstdio.h"
#include "gommath.h"

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include "colouring.h"
#include "molecoord.h"
#include "molecule.h"
#include "plot_molec.h"

#include "stdafx.h"

#if 0
static void Cylinder( float , float , float ,
                      float , float , float , float , float);
#endif

static void DrawCylinder(GLdouble , GLdouble , GLdouble , 
                         int      , int,
                         float    , float    , float,
                         float    , float    , float,
                         float);
static float ArrowAD2CD = 1.3f;
static float ArrowL2H   = 0.7f;

/***********************************************************************/
int gomp_PlotMoleculeLicorice(int Wstr) /* display the molecule as liquorice */
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS
    register int  i,j,k,l;
    static float  xh,yh,zh;
    static float  diffx,diffy,diffz;
    static float  xi,yi,zi;
    static float  xk,yk,zk;
    static float  ax,ay,az;
    
    static float  rad_c,rad_s;
    static float  CylinderLong;
    register const float *ColRED,*ColGREEN,*ColBLUE;
    register const float *XC,*YC,*ZC;
    static int from,to;
    static const char *DispList;
    static const char *LicoList;
    static int b_disp_a;
    static const int *Connect;
    static int   CylinderQ;
    static int   CylinderQ2;
    static int   SphereQ;
    static int   SphereQ2;
 
    static GLUquadricObj    *plainQuad;
    static GLUquadricObj    *SphereQuad;
    static int               QuadLoop = 0;

    static float  CRi;
    static float  CGi;
    static float  CBi;

    static float  CRk;
    static float  CGk;
    static float  CBk;

    if(!QuadLoop) {
        plainQuad = gluNewQuadric();
        gluQuadricDrawStyle(plainQuad, GLU_FILL);
        SphereQuad = gluNewQuadric();
        gluQuadricDrawStyle(SphereQuad, GLU_FILL);
        QuadLoop = 1;
    }

    CylinderQ  = gomp_GetCylinderQuality();
    CylinderQ2 = 2 * CylinderQ;
    SphereQ    = gomp_GetSphereQuality();
    SphereQ2   = 2 * SphereQ;

    b_disp_a = gomp_GetBondDisplayStyle();

    glEnable(GL_LIGHTING);

    from  = 0;
    to    = gomp_GetNumAtomsInMolecStruct(Wstr);

    ColRED   = gomp_GetAtomColourRedPointer(Wstr);
    ColGREEN = gomp_GetAtomColourGreenPointer(Wstr);
    ColBLUE  = gomp_GetAtomColourBluePointer(Wstr);
    XC       = gomp_GetAtomXCoordPointer(Wstr);
    YC       = gomp_GetAtomYCoordPointer(Wstr);
    ZC       = gomp_GetAtomZCoordPointer(Wstr);
    DispList = gomp_GetAtomDisplayStatePointer(Wstr);
    LicoList = gomp_GetAtomLicoDisplayStatePointer(Wstr);

    /* throw out the cylinders */

    rad_c = gomp_GetAtomLicoRadC(Wstr); /* radius of the cylinder */
    
    for(i = from ; i < to ; i++ ) {    /* main loop */
        
        if(LicoList[i] != 1) continue;  /* check if there should be a surface    */

        if(DispList[i] != 1) continue;  /* check if the atom should be displayed */

        xi   = XC[i];
        yi  = YC[i];
        zi = ZC[i];

        Connect = gomp_GetAtomConnection(Wstr, i);

        l = Connect[0];

        for( j = 1 ; j <= l ; j ++) {
  
/* the points are x[i],y[i],z[i] and x[k],y[k],z[k] */
            k = Connect[j];

            if(LicoList[k] != 1) continue; /* check if there should be a surface */

            if(DispList[k] != 1) continue; /* check if the atom should be displayed */

            if(k > i) continue;       /* avoid counting the bond twice */

            xk   = XC[k];
            yk  = YC[k];
            zk = ZC[k];

            if(b_disp_a) {

                diffx   = xk-xi;
                diffy  = yk-yi;
                diffz = zk-zi;

                CylinderLong = 0.5 * sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

                xh   = (xi + xk) * 0.5;
                yh  = (yi + yk) * 0.5;
                zh = (zi + zk) * 0.5;

                CRi   = ColRED[i];
                CGi  = ColGREEN[i];
                CBi = ColBLUE[i];

                if(!gomp_GetDisplayColourType())
                    (void)gomp_RGB2Grayscale(&CRi , &CGi , &CBi);

                glColor4f(CRi , 
                          CGi , 
                          CBi , 1.0); 

                ax   = xk;
                ay  = yk;
                az = zk;

                gomp_Vector2Angles(xi,yi,zi,&ax,&ay,&az);

                glColor4f(CRi , 
                          CGi , 
                          CBi , 1.0); 

                glPushMatrix();
                glTranslatef(xi,yi,zi);
                glRotatef(ax   ,  0.0 ,  0.0 , 1.0);
                glRotatef(ay   ,  0.0 ,  1.0 , 0.0);
                gluCylinder(plainQuad, rad_c , rad_c , CylinderLong , 
                            CylinderQ2 , CylinderQ);
                glPopMatrix();
     
                CRk   = ColRED[k];
                CGk  = ColGREEN[k];
                CBk = ColBLUE[k];

                if(!gomp_GetDisplayColourType())
                    (void)gomp_RGB2Grayscale(&CRk , &CGk , &CBk);

                glColor4f(CRk , 
                          CGk , 
                          CBk , 1.0); 

                glPushMatrix();
                glTranslatef(xh,yh,zh);
                glRotatef(ax   ,  0.0 ,  0.0 , 1.0);
                glRotatef(ay   ,  0.0 ,  1.0 , 0.0);
                gluCylinder(plainQuad, rad_c , rad_c , CylinderLong , 
                            CylinderQ2 , CylinderQ);
                glPopMatrix();
              
            } else {

                diffx   = xk-xi;
                diffy  = yk-yi;
                diffz = zk-zi;

                CylinderLong = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

                CRi   = ColRED[i];
                CGi  = ColGREEN[i];
                CBi = ColBLUE[i];

                CRk   = ColRED[k];
                CGk  = ColGREEN[k];
                CBk = ColBLUE[k];

                if(!gomp_GetDisplayColourType()) {
                    (void)gomp_RGB2Grayscale(&CRi , &CGi , &CBi);
                    (void)gomp_RGB2Grayscale(&CRk , &CGk , &CBk);
                }

                ax   = xk;
                ay  = yk;
                az = zk;

                gomp_Vector2Angles(xi,yi,zi,&ax,&ay,&az);

                glPushMatrix();
                glTranslatef(xi,yi,zi);
                glRotatef(ax   ,  0.0 ,  0.0 , 1.0);
                glRotatef(ay   ,  0.0 ,  1.0 , 0.0);
                DrawCylinder(rad_c , rad_c , CylinderLong , 
                             CylinderQ2 , CylinderQ,
                             CRi,CGi,CBi,
                             CRk,CGk,CBk, 1.0);
                glPopMatrix();
            }
        }
    }

/* throw out the balls */

    rad_s = gomp_GetAtomLicoRadS(Wstr);   /* radius of the sphere */

    for(k = from ; k < to ; k++) {
        
        if(LicoList[k] != 1) continue;  /* check if there should be a surface    */

        if(DispList[k] != 1) continue;  /* check if the atom should be displayed */

        CRk   = ColRED[k];
        CGk  = ColGREEN[k];
        CBk = ColBLUE[k];

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&CRk , &CGk , &CBk);

        glColor4f(CRk , 
                  CGk , 
                  CBk , 1.0); 

        glPushMatrix();
        glTranslatef(XC[k],YC[k],ZC[k]);
        gluSphere(SphereQuad, (double)rad_s , SphereQ2 , SphereQ);
        glPopMatrix();
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/***********************************************************************/
void gomp_Vector2Angles(float vx,float vy,float vz,float *px,float *py,float *pz)
/***********************************************************************/
{
    static float r,x,y,z;
    static float alpha, beta;

    x = *px - vx;
    y = *py - vy;
    z = *pz - vz;
    r = sqrt(x*x+y*y+z*z);
    x /= r;
    y /= r;
    z /= r;

    beta = atan2(y,x);
    r = y*sin(beta) + x*cos(beta);
    alpha  = atan2(r,z);

    *px = 180.0 * beta/M_PI;
    *py = 180.0 * alpha/M_PI;
    *pz = 0.0;
}

#ifdef ENABLE_GRAPHICS
#if 0
/***********************************************************************/
void Cylinder(float x0,float y0,float z0,
              float x1,float y1,float z1,
              float Rad1 , float Rad2)
/***********************************************************************/
{
    static float ax,ay,az,r;
    static GLUquadricObj    *plainQuad;
    static int               NoMore=0;
    static int   CylinderQ;
    static int   CylinderQ2;

    if(!NoMore) {
        plainQuad = gluNewQuadric();
        gluQuadricDrawStyle(plainQuad, GLU_FILL);
        NoMore = 1;
    }

    CylinderQ  = gomp_GetCylinderQuality();
    CylinderQ2 = 2 * CylinderQ;

    ax = x1;
    ay = y1;
    az = z1;

    x1 = ax - x0;
    y1 = ay - y0;
    z1 = az - z0;
    r = sqrt(x1*x1+y1*y1+z1*z1);

    gomp_Vector2Angles(x0,y0,z0,&ax,&ay,&az);

    glPushMatrix();

    glTranslatef(x0,y0,z0);

    glRotatef(ax   ,  0.0 ,  0.0 , 1.0);
    glRotatef(ay   ,  0.0 ,  1.0 , 0.0);
    gluCylinder(plainQuad, Rad1 , Rad2 , r , CylinderQ2 , CylinderQ);
    glPopMatrix();
}
#endif

#define ARROW_SCALE 0.7

/***********************************************************************/
void gomp_Arrow(float x0,float y0,float z0,float x1,float y1,float z1,float Rad)
/***********************************************************************/
{
    static float ax,ay,az,r;
    static GLUquadricObj    *plainQuad;
    static int               NoMore=0;
    static float x2,y2,z2;
    static int   CylinderQ;
    static int   CylinderQ2;

    if(!NoMore) {
        plainQuad = gluNewQuadric();
        gluQuadricDrawStyle(plainQuad, GLU_FILL);
        NoMore = 1;
    }

    CylinderQ  = gomp_GetCylinderQuality();
    CylinderQ2 = 2 * CylinderQ;

    ax = x1;
    ay = y1;
    az = z1;

    x1 = ax - x0;
    y1 = ay - y0;
    z1 = az - z0;
    r = sqrt(x1*x1+y1*y1+z1*z1);

    gomp_Vector2Angles(x0,y0,z0,&ax,&ay,&az);

    glPushMatrix();

    glTranslatef(x0,y0,z0);

    glRotatef(ax   ,  0.0 ,  0.0 , 1.0);
    glRotatef(ay   ,  0.0 ,  1.0 , 0.0);
    gluCylinder(plainQuad, Rad , Rad , ArrowL2H * r , 
                CylinderQ2 , CylinderQ);
    glPopMatrix();

    x2 = x0 + ArrowL2H * x1;
    y2 = y0 + ArrowL2H * y1;
    z2 = z0 + ArrowL2H * z1;

    glPushMatrix();

    glTranslatef(x2,y2,z2);

    glRotatef(ax   ,  0.0 ,  0.0 , 1.0);
    glRotatef(ay   ,  0.0 ,  1.0 , 0.0);
    gluCylinder(plainQuad, ArrowAD2CD * Rad , 0.001 , (1.0 - ArrowL2H) * r , 
                CylinderQ2 , CylinderQ);
    glPopMatrix();

}

/***********************************************************************/
void DrawCylinder(GLdouble rad_c1 , GLdouble rad_c2 , GLdouble CylinderLong , 
                  int CylinderQ2 , int CylinderQ,
                  float ColREDi  , float ColGREENi, float ColBLUEi,
                  float ColREDk  , float ColGREENk, float ColBLUEk,
                  float Alpha)
/***********************************************************************/
{
    int i;
    float Delta;
    float x1,y1,x2,y2;

    glEnable(GL_LIGHTING);

    Delta = 2.0 * M_PI/(float)CylinderQ2;

    for(i = 0 ; i < CylinderQ2 ; i++) {

        x1 = cos(i * Delta);
        y1 = sin(i * Delta);

        x2 = cos(i * Delta + Delta);
        y2 = sin(i * Delta + Delta);

        glBegin( GL_QUADS );

        glColor4f(ColREDk , ColGREENk , ColBLUEk , Alpha);
        glNormal3f(x1 , y1 , 0.0);
        glVertex3f( rad_c2 * x1 , rad_c2 * y1 , CylinderLong );

        glColor4f(ColREDi , ColGREENi , ColBLUEi , Alpha);
        glNormal3f(x1 , y1 , 0.0);
        glVertex3f( rad_c1 * x1 , rad_c1 * y1 , 0.0 );

        glColor4f(ColREDi , ColGREENi , ColBLUEi , Alpha);
        glNormal3f(x2 , y2 , 0.0);
        glVertex3f( rad_c1 * x2 , rad_c1 * y2 , 0.0 );

        glColor4f(ColREDk , ColGREENk , ColBLUEk , Alpha);
        glNormal3f(x2 , y2 , 0.0);
        glVertex3f( rad_c2 * x2 , rad_c2 * y2 , CylinderLong );

        glEnd();

    }

}

/***********************************************************************/
void gomp_PlotCylinder(float xi     , float yi     , float zi    ,
                     float xk     , float yk     , float zk    ,
                     float redi   , float greeni , float bluei ,
                     float redk   , float greenk , float bluek ,
                     float alpha   ,
                     float rad_ci , float rad_ck)
/***********************************************************************/
{
    static int   CylinderQ;
    static int   CylinderQ2;
    static float CylinderLong;
    static float ax,ay,az;
    static float diffx;
    static float diffy;
    static float diffz;

    diffx = xk-xi;
    diffy = yk-yi;
    diffz = zk-zi;

    CylinderLong = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    CylinderQ    = gomp_GetCylinderQuality();
    CylinderQ2   = 2 * CylinderQ;


    ax = xk;
    ay = yk;
    az = zk;

    gomp_Vector2Angles(xi,yi,zi,&ax,&ay,&az);

    glPushMatrix();
    glTranslatef(xi,yi,zi);
    glRotatef(ax   ,  0.0 ,  0.0 , 1.0);
    glRotatef(ay   ,  0.0 ,  1.0 , 0.0);
    DrawCylinder((GLdouble)rad_ci , (GLdouble)rad_ck , (GLdouble)CylinderLong ,
                 CylinderQ2 , CylinderQ,
                 redi,greeni,bluei,
                 redk,greenk,bluek,
                 alpha);
    glPopMatrix();

}
#endif /* ENABLE_GRAPHICS */
/***********************************************************************/
int   gomp_SetArrowAD2CD(float value)
/***********************************************************************/
{
    ArrowAD2CD = value;

    return(0);
}
/***********************************************************************/
int   gomp_SetArrowL2H(float value)
/***********************************************************************/
{
    ArrowL2H = value;

    return(0);
}
