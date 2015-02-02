/*

Copyright (c) 2004 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by:
Eero Häkkinen
*/

#include "maindefs.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <tcl.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

#include "colouring.h"
#include "contour.h"
#include "gomstdio.h"
#include "plot_molec.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS

#define VVAdd(u,v) do { \
    u[0] += v[0]; \
    u[1] += v[1]; \
    u[2] += v[2]; \
} while(0)

/****************************************************************************/
int gomp_CalculateIsoCurves(int Which1, float IsoValue,
                            char Axis,  float Location)
/****************************************************************************/
{
    int PointsX, PointsY, PointsZ, PointsXY, base, i, j, k, l, m;
    int limits[3][2];
    unsigned short int tests;
    GLfloat C[3], D[3], M[3], V[3], dd[4];
    const float *Data;
    struct {
        int ind;
        GLfloat D[3];
    } corners[5];

    const float *sumxyz;

    if( Which1 <  0 ||
        Which1 >= gomp_GetContoursDefined()) {
        gomp_FormatERROR(
            "contour data index 1 (%d) is out of allowed range (1 , %d)",
            Which1+1,gomp_GetContoursDefined());
        return(1);
    }

    if( IsoValue < ContourInfo[Which1].min ||
        IsoValue > ContourInfo[Which1].max) {
        gomp_FormatERROR(
            "isovalue (%f) is outside defined range (%f , %f)",
            IsoValue,ContourInfo[Which1].min , ContourInfo[Which1].max); 
        return(1);
    }

    D[0]     = (ContourInfo[Which1].Xmax - ContourInfo[Which1].Xmin) / 
        (ContourInfo[Which1].xdim - 1);

    D[1]     = (ContourInfo[Which1].Ymax - ContourInfo[Which1].Ymin) / 
        (ContourInfo[Which1].ydim - 1);

    D[2]     = (ContourInfo[Which1].Zmax - ContourInfo[Which1].Zmin) / 
        (ContourInfo[Which1].zdim - 1);

    sumxyz   = gomp_GetTranslateArray();

    PointsX  = ContourInfo[Which1].xdim;
    PointsY  = ContourInfo[Which1].ydim;
    PointsZ  = ContourInfo[Which1].zdim;
    PointsXY = PointsX * PointsY;

    Data     = ContourInfo[Which1].data;
    
    limits[0][0] = limits[1][0] = limits[2][0] = 0;
    limits[0][1] = PointsX - 2;
    limits[1][1] = PointsY - 2;
    limits[2][1] = PointsZ - 2;

    memset(corners, 0, sizeof(corners));
    switch ( Axis ) {
    case 'x':
        corners[1].ind  = PointsX;
        corners[2].ind  = PointsX + PointsXY;
        corners[3].ind  = PointsXY;
        corners[1].D[1] = D[1];
        corners[2].D[1] = D[1];
        corners[2].D[2] = D[2];
        corners[3].D[2] = D[2];
        limits[0][0] = limits[0][1] = (int)(
            ( Location - ContourInfo[Which1].Xmin + sumxyz[0] ) / D[0] + 0.5);
        if ( limits[0][0] < 0 || limits[0][0] >= PointsX )
            return(0);
        break;
    case 'y':
        corners[1].ind  = 1;
        corners[2].ind  = 1 + PointsXY;
        corners[3].ind  = PointsXY;
        corners[1].D[0] = D[0];
        corners[2].D[0] = D[0];
        corners[2].D[2] = D[2];
        corners[3].D[2] = D[2];
        limits[1][0] = limits[1][1] = (int)(
            ( Location - ContourInfo[Which1].Ymin + sumxyz[1] ) / D[1] + 0.5);
        if ( limits[1][0] < 0 || limits[1][0] >= PointsY )
            return(0);
        break;
    case 'z':
        corners[1].ind  = 1;
        corners[2].ind  = 1 + PointsX;
        corners[3].ind  = PointsX;
        corners[1].D[0] = D[0];
        corners[2].D[0] = D[0];
        corners[2].D[1] = D[1];
        corners[3].D[1] = D[1];
        limits[2][0] = limits[2][1] = (int)(
            ( Location - ContourInfo[Which1].Zmin + sumxyz[2] ) / D[2] + 0.5);
        if ( limits[2][0] < 0 || limits[2][0] >= PointsZ )
            return(0);
        break;
    default:
        return(0);
    }
    memcpy(&corners[4], &corners[0], sizeof(corners[0]));
    

    glDisable(GL_LIGHTING);

    glLineWidth((GLfloat)(gomp_GetContourLineWidth()+0.0));

    glBegin(GL_LINES);
    for ( i = limits[0][0] ; i <= limits[0][1] ; i++ ) {
        for ( j = limits[1][0] ; j <= limits[1][1] ; j++ ) {
            for ( k = limits[2][0] ; k <= limits[2][1] ; k++ ) {
                base = i + j * PointsX + k * PointsXY;
                switch (
                    ( ( Data[base + corners[0].ind] < IsoValue ) ? 0x1 : 0 ) |
                    ( ( Data[base + corners[1].ind] < IsoValue ) ? 0x2 : 0 ) |
                    ( ( Data[base + corners[2].ind] < IsoValue ) ? 0x4 : 0 ) |
                    ( ( Data[base + corners[3].ind] < IsoValue ) ? 0x8 : 0 ) ){
                case 0:
                case 0x1 | 0x2 | 0x4 | 0x8:
                    /* No intersections. */
                    continue;
                }
                
                C[0] = i * D[0] + ContourInfo[Which1].Xmin - sumxyz[0];
                C[1] = j * D[1] + ContourInfo[Which1].Ymin - sumxyz[1];
                C[2] = k * D[2] + ContourInfo[Which1].Zmin - sumxyz[2];

                tests = 0;
                for ( l = 0 ; l < 4 ; l++ ) {
                    float diff =
                        Data[base + corners[l  ].ind] -
                        Data[base + corners[l+1].ind];
                    if ( diff >= 0 && diff <= 0 )
                        continue;
                    dd[l] =
                        ( Data[base + corners[l  ].ind] - IsoValue ) / diff;
                    if ( 0 <= dd[l] && dd[l] <= 1 )
                        tests |= 1 << l;
                }

                switch ( tests ) {
                case 0x1 | 0x2:
                case 0x1 | 0x4:
                case 0x1 | 0x8:
                case 0x2 | 0x4:
                case 0x2 | 0x8:
                case 0x4 | 0x8:
                    /**
                     * Two intersections. Draw a posect line.
                     * Two of the following if statements will be executed.
                     */
                    for ( l = 0 ; l < 4 ; l++ ) {
                        if ( ! ( tests & ( 1 << l ) ) )
                            continue;
                        for ( m = 0 ; m < 3 ; m++ )
                            V[m] =
                                C[m] +  dd[l] * corners[l+1].D[m] +
                                ( 1 - dd[l] ) * corners[l  ].D[m];
                        glVertex3fv(V);
                    }
                    break;
                default:
                    M[0] = C[0] + corners[2].D[0] / 2;
                    M[1] = C[1] + corners[2].D[1] / 2;
                    M[2] = C[2] + corners[2].D[2] / 2;
                    /**
                     * One, three or all of the following if
                     * statements will be executed.
                     * Connect lines in the middle of the square.
                     */
                    for ( l = 0 ; l < 4 ; l++ ) {
                        if ( ! ( tests & ( 1 << l ) ) )
                            continue;
                        for ( m = 0 ; m < 3 ; m++ )
                            V[m] =
                                C[m] +  dd[l] * corners[l+1].D[m] +
                                ( 1 - dd[l] ) * corners[l  ].D[m];
                        glVertex3fv(V);
                        glVertex3fv(M);
                    }
                }
            }
        }
    }
    glEnd();

    glEnable(GL_LIGHTING);

    return(0);
}
#endif /* ENABLE_GRAPHICS */
