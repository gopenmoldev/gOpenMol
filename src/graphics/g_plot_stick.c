/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2001, 2002 by:
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

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include "bond.h"
#include "colouring.h"
#include "molecoord.h"
#include "molecule.h"
#include "plot_molec.h"

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS
static int plot_cross(float,float,float);
#endif
static float CrossLen = 0.2f;

/***********************************************************************/
int gomp_PlotMoleculeStick(int Wstr)             /* Plot molecule as sticks */
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS
    register int i,j,k,l,m;
    register const float *ColRED,*ColGREEN,*ColBLUE;
    register const float *XC,*YC,*ZC;
    static const int *cnct_run;
    static const char *disp_run;
    static const char *DispList;
    static int first,last;
    static int plotCross;
    static float CRi,CGi,CBi;
    static float CRk,CGk,CBk;
    static float CRh,CGh,CBh;
    static float xh,yh,zh;
    static float xi,yi,zi;
    static float xk,yk,zk;

    static int b_disp_a;

    (void)gomp_GetHBondColour(&CRh , &CGh, &CBh);
    if(!gomp_GetDisplayColourType())
        (void)gomp_RGB2Grayscale(&CRh , &CGh , &CBh);

    glDisable(GL_LIGHTING);

    glLineWidth((GLfloat)(gomp_GetMoleculeLineWidth()+0.0));

    first = 0;
    last  = gomp_GetNumAtomsInMolecStruct(Wstr);

    ColRED   = gomp_GetAtomColourRedPointer(Wstr);
    ColGREEN = gomp_GetAtomColourGreenPointer(Wstr);
    ColBLUE  = gomp_GetAtomColourBluePointer(Wstr);
    XC       = gomp_GetAtomXCoordPointer(Wstr);
    YC       = gomp_GetAtomYCoordPointer(Wstr);
    ZC       = gomp_GetAtomZCoordPointer(Wstr);
    disp_run = DispList = gomp_GetAtomDisplayStatePointer(Wstr);

    b_disp_a = gomp_GetBondDisplayStyle();

/* 1 */
    if(!b_disp_a) {
            
        for(i = first ; i < last ; i++ ) { /* main loop over atoms */

            if(*disp_run++ != 1) continue;  /* look into display list */

            CRi = ColRED[i];
            CGi = ColGREEN[i];
            CBi = ColBLUE[i];

            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&CRi , &CGi , &CBi);

            xi   = XC[i];
            yi  = YC[i];
            zi = ZC[i];

            /* Plot bonds */
            plotCross = 0;

            cnct_run = gomp_GetAtomConnection(Wstr, i);

            l = *cnct_run++;

            if( l >= 1 ) {
                for( j = 1 ; j <= l ; j ++) {

                    k = *cnct_run++;
                    if(k > i) continue;     /* avoid counting the bond twice */

                    if(DispList[k]) {
                        glBegin(GL_LINES);
                        glColor3f( CRi , CGi , CBi); /* default colour */
                        glVertex3f(xi , yi , zi);

                        CRk = ColRED[k];
                        CGk = ColGREEN[k];
                        CBk = ColBLUE[k];

                        if(!gomp_GetDisplayColourType())
                            (void)gomp_RGB2Grayscale(&CRk , &CGk , &CBk);

                        glColor3f( CRk , CGk , CBk); /* default colour */
                        glVertex3f( XC[k] , YC[k] , ZC[k]);
                        glEnd();
                    }
                    else
                        plotCross = 1;
                }
            }
            else
                plotCross = 1;

            if( plotCross ) {
                glColor3f( CRi , CGi , CBi); /* default colour */
                (void)plot_cross(xi,yi,zi);
            }
            
            /* Plot hydrogen bonds */
            cnct_run = gomp_GetAtomHydrogenBond(Wstr, i);

            l = *cnct_run++;

            for( j = 1; j <= l; ++j) {
                k = *cnct_run++; 
                
                if(k > i) continue;     /* avoid counting the bond twice */

                if(!DispList[k]) continue;

                xk   = XC[k];
                yk  = YC[k];
                zk = ZC[k];

                glBegin(GL_LINES);
                glColor3f( CRh , CGh , CBh ); /* hbond colour */

                for( m = 1; m < 9; m++)
                    glVertex3f( ( (9-m)*xi + m*xk )/9,
                                ( (9-m)*yi + m*yk )/9,
                                ( (9-m)*zi + m*zk )/9 );

                glEnd();
            }
        }
    }
    else {
        
        for(i = first ; i < last ; i++ ) { /* main loop over atoms */

            if(*disp_run++ != 1) continue;  /* look into display list */

            CRi = ColRED[i];
            CGi = ColGREEN[i];
            CBi = ColBLUE[i];
            
            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&CRi , &CGi , &CBi);

            xi  = XC[i];
            yi  = YC[i];
            zi  = ZC[i];

            /* Plot bonds */
            plotCross = 0;

            cnct_run = gomp_GetAtomConnection(Wstr, i);
            
            l = *cnct_run++;

            if( l >= 1 ) {
                for( j = 1 ; j <= l ; j ++) {

                    k = *cnct_run++; 
                    
                    if(DispList[k]) {

                        if(k > i) continue;     /* avoid counting the bond twice */

                        xh     = (XC[k] + xi) / 2.; /* second point */
                        yh    = (YC[k] + yi) / 2.;
                        zh   = (ZC[k] + zi) / 2.;
                          
                        glBegin(GL_LINES);
                        /* OGLXXX color values need to be scaled */
                        glColor3f( CRi , CGi , CBi); /* default colour */
                        glVertex3f(xi , yi , zi);
                        glVertex3f(xh , yh , zh);
                        glEnd();

                        CRk = ColRED[k];
                        CGk = ColGREEN[k];
                        CBk = ColBLUE[k];

                        if(!gomp_GetDisplayColourType())
                            (void)gomp_RGB2Grayscale(&CRk , &CGk , &CBk);
                        
                        /* OGLXXX for multiple, independent line segments: use GL_LINES */
                        glBegin(GL_LINES);
                        /* OGLXXX color values need to be scaled */
                        glColor3f( CRk , CGk , CBk); /* default colour */
                        glVertex3f(xh , yh , zh);
                        glVertex3f(XC[k] , YC[k] , ZC[k]);
                        glEnd();
                    }
                    else
                        plotCross = 1;
                }
            }
            else
                plotCross = 1;

            if( plotCross ) {
                glColor3f( CRi , CGi , CBi); /* default colour */
                (void)plot_cross(xi,yi,zi);
            }

            /* Plot hydrogen bonds */
            cnct_run = gomp_GetAtomHydrogenBond(Wstr, i);

            l = *cnct_run++;
            
            for( j = 1; j <= l; ++j) {
                k = *cnct_run++; 
                
                if(k > i) continue;     /* avoid counting the bond twice */

                if(!DispList[k]) continue;

                xk   = XC[k];
                yk  = YC[k];
                zk = ZC[k];

                glBegin(GL_LINES);
                glColor3f( CRh , CGh , CBh ); /* hbond colour */

                for( m = 1; m < 9; m++)
                    glVertex3f( ( (9-m)*xi + m*xk )/9,
                                ( (9-m)*yi + m*yk )/9,
                                ( (9-m)*zi + m*zk )/9 );

                glEnd();
            }
        }
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}

#ifdef ENABLE_GRAPHICS
/***********************************************************************/
int plot_cross(float xpt,float ypt,float zpt)  
    /* plot a cross at xpt,ypt,zpt */
/***********************************************************************/
{
    static float cross_len; /* size of the cross */

    cross_len = CrossLen;

    glBegin(GL_LINES); 
    glVertex3f(xpt-cross_len , ypt , zpt);
    glVertex3f(xpt+cross_len , ypt , zpt);
    glEnd();

    glBegin(GL_LINES); /* y - direction */
    glVertex3f(xpt , ypt-cross_len , zpt);
    glVertex3f(xpt , ypt+cross_len , zpt);
    glEnd();

    glBegin(GL_LINES); /* z - direction */
    glVertex3f(xpt , ypt , zpt-cross_len);
    glVertex3f(xpt , ypt , zpt+cross_len);
    glEnd();

    return(0);
}
#endif /* ENABLE_GRAPHICS */

/***********************************************************************/
int   gomp_SetCrossLen(float value)
/***********************************************************************/
{
    CrossLen = value;

    return(0);
}
/***********************************************************************/
float gomp_GetCrossLen(void)
/***********************************************************************/
{
    return(CrossLen);
}

