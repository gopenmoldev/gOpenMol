/*

Copyright (c) 1996 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#if defined(GLUT)
#include <GL/glut.h>
#else
#include <GL/glaux.h>
#endif

#include <GL/gl.h>
#endif /* ENABLE_GRAPHICS */

#include <tcl.h>

#include "colors.h"
#include "colouring.h"
#include "contour.h"
#include "gomfile.h"
#include "gomtext.h"
#include "math_oper.h"
#include "memalloc.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#include "stdafx.h"

#define NUMBER_BINS 41

static struct {
    float *DataX;
    float *DataY;
    float *DataZ;
    float  MinDX;
    float  MaxDX;
    float  MinDY;
    float  MaxDY;
    float  MinDZ;
    float  MaxDZ;
    int    pxX;
    int    pyX;
    int    pzX;
    int    pxY;
    int    pyY;
    int    pzY;
    int    pxZ;
    int    pyZ;
    int    pzZ;

    float  MinXX;
    float  MaxXX;
    float  MinYX;
    float  MaxYX;
    float  MinZX;
    float  MaxZX; 

    float  MinXY;
    float  MaxXY;
    float  MinYY;
    float  MaxYY;
    float  MinZY;
    float  MaxZY; 

    float  MinXZ;
    float  MaxXZ;
    float  MinYZ;
    float  MaxYZ;
    float  MinZZ;
    float  MaxZZ; 

    float  DXX;
    float  DYX;
    float  DZX;

    float  DXY;
    float  DYY;
    float  DZY;

    float  DXZ;
    float  DYZ;
    float  DZZ;

    float *xyzCoordsX;
    float *xyzCoordsY;
    float *xyzCoordsZ;
    int    InUseX;
    int    InUseY;
    int    InUseZ;
    int    PlotStateX;
    int    PlotStateY;
    int    PlotStateZ;
    int    PlotStateXYZ1;
    int    PlotStateXYZ2;
    int    PlotStateXYZ3;
    int    PlotWhichXYZ1;
    int    PlotWhichXYZ2;
    int    PlotWhichXYZ3;
    float  PlotScaleMinXYZ1;
    float  PlotScaleMinXYZ2;
    float  PlotScaleMinXYZ3;
    float  PlotScaleMaxXYZ1;
    float  PlotScaleMaxXYZ2;
    float  PlotScaleMaxXYZ3;
    float  PlotPlaneXYZ1[4];
    float  PlotPlaneXYZ2[4];
    float  PlotPlaneXYZ3[4];
    int    PlotManX;
    int    PlotManY;
    int    PlotManZ;
} CutPlane;

#ifdef ENABLE_GRAPHICS
static struct {
    float *Data;
    float  Min;
    float  Max;
    int    px;
    int    py;
    int    pz;
    int    NumBins;
    float  BinMin;
    float  BinMax;
    float  BinCut;

    float  MinX;
    float  MaxX;
    float  MinY;
    float  MaxY;
    float  MinZ;
    float  MaxZ; 

    int    InUse;
    int    PlotState;
} CutPlaneSpectrum;
#endif

static int    CutPlaneType_3D_X = 0;
static int    CutPlaneType_3D_Y = 0;
static int    CutPlaneType_3D_Z = 0;
static float  CutPlaneDampingConstant = 1.0;

static gom_Plotter *CutPlaneCallbackHandle;

static int Check = 0;
static int    PlotCutPlaneX(void);
static int    PlotCutPlane_3D_X(void);
static int    GetCutPlanePlotStateX(void);
static int    EnableCutPlanePlotStateX(void);
static int    DisableCutPlanePlotStateX(void);

static int    PlotCutPlaneY(void);
static int    PlotCutPlane_3D_Y(void);
static int    GetCutPlanePlotStateY(void);
static int    EnableCutPlanePlotStateY(void);
static int    DisableCutPlanePlotStateY(void);

static int    PlotCutPlaneZ(void);
static int    PlotCutPlane_3D_Z(void);
static int    GetCutPlanePlotStateZ(void);
static int    EnableCutPlanePlotStateZ(void);
static int    DisableCutPlanePlotStateZ(void);

#if 0
static float  GetCutPlaneXcoord(void);
static float  GetCutPlaneYcoord(void);
static float  GetCutPlaneZcoord(void);
static float  GetCutPlaneYZmin(int);
static float  GetCutPlaneYZmax(int);
static float  GetCutPlaneXZmin(int);
static float  GetCutPlaneXZmax(int);
static float  GetCutPlaneXYmin(int);
static float  GetCutPlaneXYmax(int);
static int    EnableCutPlanePlotStateXYZ(int);
#endif

static int    PlotCutPlaneSpectrumX(void);
static int    PlotCutPlaneSpectrumY(void);
static int    PlotCutPlaneSpectrumZ(void);

/****************************************************************************/
int gomp_PrepareCutPlaneZ(int Which , float Zcoord , 
                        const char *Text1 , const char *Text2, const char *Text3)
/****************************************************************************/
{
    int    i;
    int    j;
    int    k1;
    int    k2;
    int    Loop;
    int    Jump1;
    int    Jump2;
    float  Temp1;
    float  Temp2;
    float  Temp3;
    const float *sumxyz;
    char   Text[BUFF_LEN];

    if(Which <  0 ||
       Which >= gomp_GetContoursDefined()) {
        sprintf(Text,"contour data index (%d) out of allowed range (1 <=> %d)",
                Which+1,gomp_GetContoursDefined());
        gomp_PrintERROR(Text);
        return(1);
    }

    if(Zcoord < ContourInfo[Which].Zmin ||
       Zcoord > ContourInfo[Which].Zmax) {
        sprintf(Text,"z-coordinate is outside range (%f , %f)",
                ContourInfo[Which].Zmin , ContourInfo[Which].Zmax); 
        gomp_PrintERROR(Text);
        return(1);
    }

#ifdef ENABLE_GRAPHICS
    (void)gomp_Prepare1DTexture();
#endif /* ENABLE_GRAPHICS */

    if(CutPlane.InUseZ) {
        gomp_DeleteCutPlaneDataZ();
    }

    CutPlane.InUseZ = 0;

    sumxyz         = gomp_GetTranslateArray();

    if(strlen(Text1) != 0)
        CutPlane.MinDZ = atof(Text1);
    else
        CutPlane.MinDZ = ContourInfo[Which].min;

    if(strlen(Text2) != 0)
        CutPlane.MaxDZ = atof(Text2);
    else
        CutPlane.MaxDZ = ContourInfo[Which].max;

    CutPlane.PlotManZ = CUTPLANE_NO_ACTION;
    if(strlen(Text3) != 0) {
        if(      gomp_StringMatch(Text3 , "sqrt")) {
            CutPlane.PlotManZ = CUTPLANE_SQRT_ACTION;
        } 
        else if( gomp_StringMatch(Text3 , "log10"))   {
            CutPlane.PlotManZ = CUTPLANE_LOG10_ACTION;
        } 
    }

    sprintf(Text,"Min value: %e (blue)",CutPlane.MinDZ);
    gomp_PrintMessage(Text);
    sprintf(Text,"Max value: %e (red)", CutPlane.MaxDZ);
    gomp_PrintMessage(Text);

    CutPlane.pxZ   = ContourInfo[Which].xdim;
    CutPlane.pyZ   = ContourInfo[Which].ydim;
    CutPlane.pzZ   = ContourInfo[Which].zdim;

    CutPlane.MinXZ = ContourInfo[Which].Xmin;
    CutPlane.MaxXZ = ContourInfo[Which].Xmax;

    CutPlane.MinYZ = ContourInfo[Which].Ymin;
    CutPlane.MaxYZ = ContourInfo[Which].Ymax;

    CutPlane.MinZZ = ContourInfo[Which].Zmin;
    CutPlane.MaxZZ = ContourInfo[Which].Zmax;

    CutPlane.DXZ   = (ContourInfo[Which].Xmax - ContourInfo[Which].Xmin) / 
        (float)(ContourInfo[Which].xdim - 1);

    CutPlane.DYZ   = (ContourInfo[Which].Ymax - ContourInfo[Which].Ymin) / 
        (float)(ContourInfo[Which].ydim - 1);

    CutPlane.DZZ   = (ContourInfo[Which].Zmax - ContourInfo[Which].Zmin) / 
        (float)(ContourInfo[Which].zdim - 1);

    k1 = (int)((Zcoord - CutPlane.MinZZ)/CutPlane.DZZ);
    k2 = k1 + 1;
    Temp3 = (Zcoord - (k1 * CutPlane.DZZ + CutPlane.MinZZ)) / CutPlane.DZZ;

    CutPlane.DataZ      = gomp_AllocateFloatVector(    CutPlane.pxZ * CutPlane.pyZ);
    CutPlane.xyzCoordsZ = gomp_AllocateFloatVector(3 * CutPlane.pxZ * CutPlane.pyZ); 

    Loop   = 0;
    Temp1  = CutPlane.MaxDZ - CutPlane.MinDZ;
    Jump1  = k1 * CutPlane.pxZ * CutPlane.pyZ;
    Jump2  = k2 * CutPlane.pxZ * CutPlane.pyZ;

    for(j = 0 ; j < CutPlane.pyZ ; j++)    {

        for( i = 0 ; i < CutPlane.pxZ ; i++) {

            Temp2 = ContourInfo[Which].data[Jump1 + j * CutPlane.pxZ + i] +
                Temp3 * 
                (ContourInfo[Which].data[Jump2 + j * CutPlane.pxZ + i] -
                 ContourInfo[Which].data[Jump1 + j * CutPlane.pxZ + i]);


            switch(CutPlane.PlotManZ) {

            case CUTPLANE_NO_ACTION:
                CutPlane.DataZ[Loop] = (Temp2 - CutPlane.MinDZ) / Temp1;
                break;
            case CUTPLANE_SQRT_ACTION:
                CutPlane.DataZ[Loop] = sqrt(fabs(Temp2 - CutPlane.MinDZ) 
                                            / Temp1);
                break;
            case CUTPLANE_LOG10_ACTION:
                CutPlane.DataZ[Loop] = 1 - log(fabs(Temp2 - CutPlane.MinDZ)) 
                    / log(Temp1);
                break;
            }

            CutPlane.xyzCoordsZ[3 * Loop    ] = 
                i * CutPlane.DXZ + CutPlane.MinXZ - sumxyz[0];
            CutPlane.xyzCoordsZ[3 * Loop + 1] = 
                j * CutPlane.DYZ + CutPlane.MinYZ - sumxyz[1];
            CutPlane.xyzCoordsZ[3 * Loop + 2] = Zcoord - sumxyz[2];

            /*
              printf("%f %f %f %f %f\n",Temp2,Temp1,CutPlane.xyzCoordsZ[3 * Loop    ],
              CutPlane.xyzCoordsZ[3 * Loop + 1],
              CutPlane.xyzCoordsZ[3 * Loop + 2]);
            */
            Loop++;
        }
    }

    CutPlane.InUseZ = Which + 1;

    (void)EnableCutPlanePlotStateZ();

    return(0);
}
#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int PlotCutPlaneZ()
/****************************************************************************/
{
    static int   i;
    static int   j;
    static int   Loop;
    static int   Jump1;
    static int   Jump2;
    static int   Jump3;
    static int   Jump4;
    static int   Jump5;
    static int   Jump6;
    static int   ColorMapping;
    static float Temp1;
    static float Temp2;
    static float Temp3;
    static float Temp4;
    static float Xc;
    static float Yc;
    static float Zc;
    static float Red;
    static float Green;
    static float Blue;

    if(!GetCutPlanePlotStateZ())
        return(0);

    ColorMapping = gomp_GetColorMappingType();
/*    glDisable(GL_LIGHTING);*/
    glEnable(GL_LIGHTING);
    glColor3f(  1.0 , 1.0 , 1.0);
    glNormal3f( 0.0 , 0.0 , 1.0);

    if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
        (void)gomp_Prepare1DTexture();
        glEnable(GL_TEXTURE_1D);
    }
    else {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
    }


    for(j = 0 ; j < CutPlane.pyZ - 1 ; j++)    {

        Jump1 =  j     * CutPlane.pxZ;
        Jump2 =  Jump1 + CutPlane.pxZ;

        for( i = 0 ; i < CutPlane.pxZ - 1; i++) {

            Jump3 = 3 * (Jump1 + i);
            Jump4 =      Jump3 + 3;
            Jump6 = 3 * (Jump2 + i);
            Jump5 =      Jump6 + 3;

            Temp1 = CutPlane.DataZ[Jump1 + i    ];
            Temp2 = CutPlane.DataZ[Jump1 + i + 1];
            Temp3 = CutPlane.DataZ[Jump2 + i + 1];
            Temp4 = CutPlane.DataZ[Jump2 + i    ];

            glPushMatrix();

        
            if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {

/* strange but this does not work on my HP OmniBook XE3
   glBegin( GL_QUADS );
*/
                glBegin( GL_POLYGON );

                Xc =  CutPlane.xyzCoordsZ[ Jump3    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsZ[ Jump3 + 2];

                glTexCoord1f( Temp1 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump4    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsZ[ Jump4 + 2];

                glTexCoord1f( Temp2 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump5    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsZ[ Jump5 + 2];

                glTexCoord1f( Temp3 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump6    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsZ[ Jump6 + 2];

                glTexCoord1f( Temp4 );
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            } else {

                glBegin( GL_QUADS );

                Xc =  CutPlane.xyzCoordsZ[ Jump3    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsZ[ Jump3 + 2];

                gomp_PreRainbow((double)Temp1 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump4    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsZ[ Jump4 + 2];

                gomp_PreRainbow((double)Temp2 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump5    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsZ[ Jump5 + 2];

                gomp_PreRainbow((double)Temp3 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump6    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsZ[ Jump6 + 2];

                gomp_PreRainbow((double)Temp4 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            }

            glPopMatrix();

            Loop++;
        }
    }

    glDisable(GL_TEXTURE_1D);

    return(0);
}

/****************************************************************************/
int PlotCutPlane_3D_Z()
/****************************************************************************/
{
    static int   i;
    static int   j;
    static int   Loop;
    static int   Jump1;
    static int   Jump2;
    static int   Jump3;
    static int   Jump4;
    static int   Jump5;
    static int   Jump6;
    static int   ColorMapping;
    static float Temp1;
    static float Temp2;
    static float Temp3;
    static float Temp4;
    static float Temp5;
    static float Xc;
    static float Yc;
    static float Zc;
    static float Red;
    static float Green;
    static float Blue;
    static float AmpDamping = 1.0;

    if(!GetCutPlanePlotStateZ())
        return(0);

    AmpDamping = gomp_GetCutPlaneDamping();

    ColorMapping = gomp_GetColorMappingType();
/*    glDisable(GL_LIGHTING);*/
    glEnable(GL_LIGHTING);
    glColor3f(  1.0 , 1.0 , 1.0);
    glNormal3f( 0.0 , 0.0 , 1.0);

    if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
        (void)gomp_Prepare1DTexture();
        glEnable(GL_TEXTURE_1D);
    }
    else {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
    }

    Temp5  = CutPlane.MinDZ / (CutPlane.MaxDZ - CutPlane.MinDZ);

    for(j = 0 ; j < CutPlane.pyZ - 1   ; j++)    {

        Jump1 =  j     * CutPlane.pxZ;
        Jump2 =  Jump1 + CutPlane.pxZ;

        for( i = 0 ; i < CutPlane.pxZ - 1; i++) {

            Jump3 = 3 * (Jump1 + i);
            Jump4 =      Jump3 + 3;
            Jump6 = 3 * (Jump2 + i);
            Jump5 =      Jump6 + 3;

            Temp1 = CutPlane.DataZ[Jump1 + i    ];
            Temp2 = CutPlane.DataZ[Jump1 + i + 1];
            Temp3 = CutPlane.DataZ[Jump2 + i + 1];
            Temp4 = CutPlane.DataZ[Jump2 + i    ];

            glPushMatrix();

            if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {

/* strange but this does not work on my HP OmniBook XE3
   glBegin( GL_QUADS );
*/
                glBegin( GL_POLYGON );

                Xc =  CutPlane.xyzCoordsZ[ Jump3    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsZ[ Jump3 + 2] + AmpDamping * (Temp1 + Temp5);

                glTexCoord1f( Temp1 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump4    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsZ[ Jump4 + 2] + AmpDamping * (Temp2 + Temp5);

                glTexCoord1f( Temp2 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump5    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsZ[ Jump5 + 2] + AmpDamping * (Temp3 + Temp5);

                glTexCoord1f( Temp3 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump6    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsZ[ Jump6 + 2] + AmpDamping * (Temp4 + Temp5);

                glTexCoord1f( Temp4 );
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            } else {

                glBegin( GL_QUADS );

                Xc =  CutPlane.xyzCoordsZ[ Jump3    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsZ[ Jump3 + 2] + AmpDamping * (Temp1 + Temp5);

                gomp_PreRainbow((double)Temp1 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump4    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsZ[ Jump4 + 2] + AmpDamping * (Temp2 + Temp5);

                gomp_PreRainbow((double)Temp2 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump5    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsZ[ Jump5 + 2] + AmpDamping * (Temp3 + Temp5);

                gomp_PreRainbow((double)Temp3 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsZ[ Jump6    ];
                Yc =  CutPlane.xyzCoordsZ[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsZ[ Jump6 + 2] + AmpDamping * (Temp4 + Temp5);

                gomp_PreRainbow((double)Temp4 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            }

            glPopMatrix();

            Loop++;
        }
    }

    glDisable( GL_TEXTURE_1D );

    return(0);
}
#endif /* ENABLE_GRAPHICS */


/****************************************************************************/
int gomp_PrepareCutPlaneY(int Which , float Ycoord , 
                        const char *Text1 , const char *Text2, const char *Text3)
/****************************************************************************/
{
    int    i;
    int    j;
    int    k1;
    int    k2;
    int    Loop;
    float  Temp1;
    float  Temp2;
    float  Temp3;
    const float *sumxyz;
    char   Text[BUFF_LEN];

    if(Which <  0  ||
       Which >= gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour data index out of range");
        return(1);
    }

    if(Ycoord < ContourInfo[Which].Ymin ||
       Ycoord > ContourInfo[Which].Ymax) {
        sprintf(Text,"y-coordinate is outside range (%f , %f)",
                ContourInfo[Which].Ymin , ContourInfo[Which].Ymax); 
        gomp_PrintERROR(Text);
        return(1);
    }

#ifdef ENABLE_GRAPHICS
    (void)gomp_Prepare1DTexture();
#endif /* ENABLE_GRAPHICS */

    if(CutPlane.InUseY) {
        (void)gomp_DeleteCutPlaneDataY();
    }

    CutPlane.InUseY = 0;


    sumxyz         = gomp_GetTranslateArray();

    if(strlen(Text1) != 0)
        CutPlane.MinDY = atof(Text1);
    else
        CutPlane.MinDY = ContourInfo[Which].min;

    if(strlen(Text2) != 0)
        CutPlane.MaxDY = atof(Text2);
    else
        CutPlane.MaxDY = ContourInfo[Which].max;

    CutPlane.PlotManY = CUTPLANE_NO_ACTION;
    if(strlen(Text3) != 0) {
        if(      gomp_StringMatch(Text3 , "sqrt")) {
            CutPlane.PlotManY = CUTPLANE_SQRT_ACTION;
        } 
        else if( gomp_StringMatch(Text3 , "log10"))   {
            CutPlane.PlotManY = CUTPLANE_LOG10_ACTION;
        } 
    }

    sprintf(Text,"Min value: %e (blue)",CutPlane.MinDY);
    gomp_PrintMessage(Text);
    sprintf(Text,"Max value: %e (red)", CutPlane.MaxDY);
    gomp_PrintMessage(Text);

    CutPlane.pxY   = ContourInfo[Which].xdim;
    CutPlane.pyY   = ContourInfo[Which].ydim;
    CutPlane.pzY   = ContourInfo[Which].zdim;

    CutPlane.MinXY = ContourInfo[Which].Xmin;
    CutPlane.MaxXY = ContourInfo[Which].Xmax;

    CutPlane.MinYY = ContourInfo[Which].Ymin;
    CutPlane.MaxYY = ContourInfo[Which].Ymax;

    CutPlane.MinZY = ContourInfo[Which].Zmin;
    CutPlane.MaxZY = ContourInfo[Which].Zmax;

    CutPlane.DXY   = (ContourInfo[Which].Xmax - ContourInfo[Which].Xmin) / 
        (float)(ContourInfo[Which].xdim - 1);

    CutPlane.DYY   = (ContourInfo[Which].Ymax - ContourInfo[Which].Ymin) / 
        (float)(ContourInfo[Which].ydim - 1);

    CutPlane.DZY   = (ContourInfo[Which].Zmax - ContourInfo[Which].Zmin) / 
        (float)(ContourInfo[Which].zdim - 1);

    k1 = (int)((Ycoord - CutPlane.MinYY)/CutPlane.DYY);
    k2 = k1 + 1;
    Temp3 = (Ycoord - (k1 * CutPlane.DYY + CutPlane.MinYY)) / CutPlane.DYY;

    CutPlane.DataY      = gomp_AllocateFloatVector(    CutPlane.pxY * CutPlane.pzY);
    CutPlane.xyzCoordsY = gomp_AllocateFloatVector(3 * CutPlane.pxY * CutPlane.pzY); 

    Loop   = 0;
    Temp1  = CutPlane.MaxDY - CutPlane.MinDY;

    for(j = 0 ; j < CutPlane.pzY ; j++)    {

        for( i = 0 ; i < CutPlane.pxY ; i++) {

            Temp2 = ContourInfo[Which].data[
                j * CutPlane.pxY * CutPlane.pyY + k1 * CutPlane.pxY + i] +
                Temp3 * 
                (ContourInfo[Which].data[
                    j * CutPlane.pxY * CutPlane.pyY + k2 * CutPlane.pxY + i] -
                 ContourInfo[Which].data[
                     j * CutPlane.pxY * CutPlane.pyY + k1 * CutPlane.pxY + i]);

            switch(CutPlane.PlotManY) {

            case CUTPLANE_NO_ACTION:

                CutPlane.DataY[Loop] = (Temp2 - CutPlane.MinDY) / Temp1;
                break;
            case CUTPLANE_SQRT_ACTION:
                CutPlane.DataY[Loop] = sqrt(fabs(Temp2 - CutPlane.MinDY) 
                                            / Temp1);
                break;
            case CUTPLANE_LOG10_ACTION:
                CutPlane.DataY[Loop] = 1. - log10(fabs(Temp2 - CutPlane.MinDY))
                    / log10(Temp1);
                break;
            }

            CutPlane.xyzCoordsY[3 * Loop    ] = 
                i * CutPlane.DXY + CutPlane.MinXY - sumxyz[0];
            CutPlane.xyzCoordsY[3 * Loop + 1] = Ycoord - sumxyz[1];
            CutPlane.xyzCoordsY[3 * Loop + 2] = 
                j * CutPlane.DZY + CutPlane.MinZY - sumxyz[2];

            Loop++;
        }
    }

    CutPlane.InUseY = Which + 1;

    (void)EnableCutPlanePlotStateY();

    return(0);
}
#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int PlotCutPlaneY()
/****************************************************************************/
{
    static int   i;
    static int   j;
    static int   Loop;
    static int   Jump1;
    static int   Jump2;
    static int   Jump3;
    static int   Jump4;
    static int   Jump5;
    static int   Jump6;
    static int   ColorMapping;
    static float Temp1;
    static float Temp2;
    static float Temp3;
    static float Temp4;
    static float Xc;
    static float Yc;
    static float Zc;
    static float Red;
    static float Green;
    static float Blue;

    if(!GetCutPlanePlotStateY())
        return(0);

    ColorMapping = gomp_GetColorMappingType();
    /*    glDisable(GL_LIGHTING);*/
    glEnable(GL_LIGHTING);
    glColor3f(  1.0 ,  1.0 , 1.0);
    glNormal3f( 0.0 , -1.0 , 0.0);

    if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
        (void)gomp_Prepare1DTexture();
        glEnable(GL_TEXTURE_1D);
    }
    else {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
    }

    for(j = 0 ; j < CutPlane.pzY - 1   ; j++)    {

        Jump1 =  j      * CutPlane.pxY;
        Jump2 =  Jump1  + CutPlane.pxY;

        for( i = 0 ; i < CutPlane.pxY - 1; i++) {

            Jump3 = 3 * (Jump1 + i);
            Jump4 =      Jump3 + 3;
            Jump6 = 3 * (Jump2 + i);
            Jump5 =      Jump6 + 3;

            Temp1 = CutPlane.DataY[Jump1 + i    ];
            Temp2 = CutPlane.DataY[Jump1 + i + 1];
            Temp3 = CutPlane.DataY[Jump2 + i + 1];
            Temp4 = CutPlane.DataY[Jump2 + i    ];

            glPushMatrix();
 
            if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {

/* strange but this does not work on my HP OmniBook XE3
   glBegin( GL_QUADS );
*/
                glBegin( GL_POLYGON );

                Xc =  CutPlane.xyzCoordsY[ Jump3    ];
                Yc =  CutPlane.xyzCoordsY[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsY[ Jump3 + 2];

                glTexCoord1f( Temp1 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump4    ];
                Yc =  CutPlane.xyzCoordsY[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsY[ Jump4 + 2];

                glTexCoord1f( Temp2 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump5    ];
                Yc =  CutPlane.xyzCoordsY[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsY[ Jump5 + 2];

                glTexCoord1f( Temp3 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump6    ];
                Yc =  CutPlane.xyzCoordsY[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsY[ Jump6 + 2];

                glTexCoord1f( Temp4 );
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            } else {

                glBegin( GL_QUADS );

                Xc =  CutPlane.xyzCoordsY[ Jump3    ];
                Yc =  CutPlane.xyzCoordsY[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsY[ Jump3 + 2];

                gomp_PreRainbow((double)Temp1 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump4    ];
                Yc =  CutPlane.xyzCoordsY[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsY[ Jump4 + 2];

                gomp_PreRainbow((double)Temp2 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump5    ];
                Yc =  CutPlane.xyzCoordsY[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsY[ Jump5 + 2];

                gomp_PreRainbow((double)Temp3 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump6    ];
                Yc =  CutPlane.xyzCoordsY[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsY[ Jump6 + 2];

                gomp_PreRainbow((double)Temp4 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            }

            glPopMatrix();

            Loop++;
        }
    }

    glDisable( GL_TEXTURE_1D );

    return(0);
}

/****************************************************************************/
int PlotCutPlane_3D_Y()
/****************************************************************************/
{
    static int   i;
    static int   j;
    static int   Loop;
    static int   Jump1;
    static int   Jump2;
    static int   Jump3;
    static int   Jump4;
    static int   Jump5;
    static int   Jump6;
    static int   ColorMapping;
    static float Temp1;
    static float Temp2;
    static float Temp3;
    static float Temp4;
    static float Temp5;
    static float Xc;
    static float Yc;
    static float Zc;
    static float Red;
    static float Green;
    static float Blue;
    static float AmpDamping = 1.0;

    if(!GetCutPlanePlotStateY())
        return(0);

    AmpDamping = gomp_GetCutPlaneDamping();

    ColorMapping = gomp_GetColorMappingType();
    /*    glDisable(GL_LIGHTING);*/
    glEnable(GL_LIGHTING);
    glColor3f(  1.0 ,  1.0 , 1.0);
    glNormal3f( 0.0 , -1.0 , 0.0);

    if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
        (void)gomp_Prepare1DTexture();
        glEnable(GL_TEXTURE_1D);
    }
    else {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
    }

    Temp5  = CutPlane.MinDY / (CutPlane.MaxDY - CutPlane.MinDY);

    for(j = 0 ; j < CutPlane.pzY - 1   ; j++)    {

        Jump1 =  j      * CutPlane.pxY;
        Jump2 =  Jump1  + CutPlane.pxY;

        for( i = 0 ; i < CutPlane.pxY - 1; i++) {

            Jump3 = 3 * (Jump1 + i);
            Jump4 =      Jump3 + 3;
            Jump6 = 3 * (Jump2 + i);
            Jump5 =      Jump6 + 3;

            Temp1 = CutPlane.DataY[Jump1 + i    ];
            Temp2 = CutPlane.DataY[Jump1 + i + 1];
            Temp3 = CutPlane.DataY[Jump2 + i + 1];
            Temp4 = CutPlane.DataY[Jump2 + i    ];

            glPushMatrix();
 
            if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {

/* strange but this does not work on my HP OmniBook XE3
   glBegin( GL_QUADS );
*/
                glBegin( GL_POLYGON );

                Xc =  CutPlane.xyzCoordsY[ Jump3    ];
                Yc =  CutPlane.xyzCoordsY[ Jump3 + 1] + AmpDamping * (Temp1 + Temp5);
                Zc =  CutPlane.xyzCoordsY[ Jump3 + 2];

                glTexCoord1f( Temp1 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump4    ];
                Yc =  CutPlane.xyzCoordsY[ Jump4 + 1] + AmpDamping * (Temp2 + Temp5);
                Zc =  CutPlane.xyzCoordsY[ Jump4 + 2];

                glTexCoord1f( Temp2 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump5    ];
                Yc =  CutPlane.xyzCoordsY[ Jump5 + 1] + AmpDamping * (Temp3 + Temp5); 
                Zc =  CutPlane.xyzCoordsY[ Jump5 + 2];

                glTexCoord1f( Temp3 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump6    ];
                Yc =  CutPlane.xyzCoordsY[ Jump6 + 1] + AmpDamping * (Temp4 + Temp5); 
                Zc =  CutPlane.xyzCoordsY[ Jump6 + 2];

                glTexCoord1f( Temp4 );
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            } else {

                glBegin( GL_QUADS );

                Xc =  CutPlane.xyzCoordsY[ Jump3    ];
                Yc =  CutPlane.xyzCoordsY[ Jump3 + 1] + AmpDamping * (Temp1 + Temp5);
                Zc =  CutPlane.xyzCoordsY[ Jump3 + 2];

                gomp_PreRainbow((double)Temp1 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump4    ];
                Yc =  CutPlane.xyzCoordsY[ Jump4 + 1] + AmpDamping * (Temp2 + Temp5);
                Zc =  CutPlane.xyzCoordsY[ Jump4 + 2];

                gomp_PreRainbow((double)Temp2 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump5    ];
                Yc =  CutPlane.xyzCoordsY[ Jump5 + 1] + AmpDamping * (Temp3 + Temp5); 
                Zc =  CutPlane.xyzCoordsY[ Jump5 + 2];

                gomp_PreRainbow((double)Temp3 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsY[ Jump6    ];
                Yc =  CutPlane.xyzCoordsY[ Jump6 + 1] + AmpDamping * (Temp4 + Temp5); 
                Zc =  CutPlane.xyzCoordsY[ Jump6 + 2];

                gomp_PreRainbow((double)Temp4 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            }

            glPopMatrix();

            Loop++;
        }
    }

    glDisable( GL_TEXTURE_1D );

    return(0);
}
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int gomp_PrepareCutPlaneX(int Which , float Xcoord , 
                        const char *Text1, const char *Text2 , const char *Text3)
/****************************************************************************/
{
    int    i;
    int    j;
    int    k1;
    int    k2;
    int    Loop;
    float  Temp1;
    float  Temp2;
    float  Temp3;
    const float *sumxyz;
    char   Text[BUFF_LEN];

    if(Which <  0 ||
       Which >= gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour data index out of range");
        return(1);
    }

    if(Xcoord < ContourInfo[Which].Xmin ||
       Xcoord > ContourInfo[Which].Xmax) {
        sprintf(Text,"x-coordinate is outside range (%f , %f)",
                ContourInfo[Which].Xmin , ContourInfo[Which].Xmax); 
        gomp_PrintERROR(Text);
        return(1);
    }

#ifdef ENABLE_GRAPHICS
    (void)gomp_Prepare1DTexture();
#endif /* ENABLE_GRAPHICS */

    if(CutPlane.InUseX) {
        (void)gomp_DeleteCutPlaneDataX();
    }

    CutPlane.InUseX = 0;

    sumxyz         = gomp_GetTranslateArray();

    if(strlen(Text1) != 0)
        CutPlane.MinDX = atof(Text1);
    else
        CutPlane.MinDX = ContourInfo[Which].min;

    if(strlen(Text2) != 0)
        CutPlane.MaxDX = atof(Text2);
    else
        CutPlane.MaxDX = ContourInfo[Which].max;

    CutPlane.PlotManX = CUTPLANE_NO_ACTION;
    if(strlen(Text3) != 0) {
        if(      gomp_StringMatch(Text3 , "sqrt")) {
            CutPlane.PlotManX = CUTPLANE_SQRT_ACTION;
        } 
        else if( gomp_StringMatch(Text3 , "log10"))   {
            CutPlane.PlotManX = CUTPLANE_LOG10_ACTION;
        } 
    }

    sprintf(Text,"Min value: %e (blue)",CutPlane.MinDX);
    gomp_PrintMessage(Text);
    sprintf(Text,"Max value: %e (red)", CutPlane.MaxDX);
    gomp_PrintMessage(Text);

    CutPlane.pxX   = ContourInfo[Which].xdim;
    CutPlane.pyX   = ContourInfo[Which].ydim;
    CutPlane.pzX   = ContourInfo[Which].zdim;

    CutPlane.MinXX = ContourInfo[Which].Xmin;
    CutPlane.MaxXX = ContourInfo[Which].Xmax;

    CutPlane.MinYX = ContourInfo[Which].Ymin;
    CutPlane.MaxYX = ContourInfo[Which].Ymax;

    CutPlane.MinZX = ContourInfo[Which].Zmin;
    CutPlane.MaxZX = ContourInfo[Which].Zmax;

    CutPlane.DXX   = (ContourInfo[Which].Xmax - ContourInfo[Which].Xmin) / 
        (float)(ContourInfo[Which].xdim - 1);

    CutPlane.DYX   = (ContourInfo[Which].Ymax - ContourInfo[Which].Ymin) / 
        (float)(ContourInfo[Which].ydim - 1);

    CutPlane.DZX   = (ContourInfo[Which].Zmax - ContourInfo[Which].Zmin) / 
        (float)(ContourInfo[Which].zdim - 1);


    k1 = (int)((Xcoord - CutPlane.MinXX)/CutPlane.DXX);
    k2 = k1 + 1;
    Temp3 = (Xcoord - (k1 * CutPlane.DXX + CutPlane.MinXX)) / CutPlane.DXX;

    CutPlane.DataX      = gomp_AllocateFloatVector(    CutPlane.pyX * CutPlane.pzX);
    CutPlane.xyzCoordsX = gomp_AllocateFloatVector(3 * CutPlane.pyX * CutPlane.pzX); 

    Loop   = 0;
    Temp1  = CutPlane.MaxDX - CutPlane.MinDX;

    for(j = 0 ; j < CutPlane.pzX ; j++)    {

        for( i = 0 ; i < CutPlane.pyX ; i++) {

            Temp2 = ContourInfo[Which].data[
                j * CutPlane.pxX * CutPlane.pyX + i * CutPlane.pxX + k1] +
                Temp3 * 
                (ContourInfo[Which].data[
                    j * CutPlane.pxX * CutPlane.pyX + i * CutPlane.pxX + k2] -
                 ContourInfo[Which].data[
                     j * CutPlane.pxX * CutPlane.pyX + i * CutPlane.pxX + k1]);

            switch(CutPlane.PlotManX) {
            case CUTPLANE_NO_ACTION:
                CutPlane.DataX[Loop] = (Temp2 - CutPlane.MinDX) / Temp1;
                break;
            case CUTPLANE_SQRT_ACTION:
                CutPlane.DataX[Loop] = sqrt(fabs(Temp2 - CutPlane.MinDX)
                                            / Temp1);
                break;
            case CUTPLANE_LOG10_ACTION:
                CutPlane.DataX[Loop] = 1. - log10(fabs(Temp2 - CutPlane.MinDX))
                    / log10(Temp1);
                break;
            }

            CutPlane.xyzCoordsX[3 * Loop    ] = Xcoord - sumxyz[0];
            CutPlane.xyzCoordsX[3 * Loop + 1] = 
                i * CutPlane.DYX + CutPlane.MinYX - sumxyz[1];
            CutPlane.xyzCoordsX[3 * Loop + 2] = 
                j * CutPlane.DZX + CutPlane.MinZX - sumxyz[2];
/*
  printf("%f %f %f %d %d\n",CutPlane.xyzCoordsX[3 * Loop    ],
  CutPlane.xyzCoordsX[3 * Loop  +1],
  CutPlane.xyzCoordsX[3 * Loop  +2],i,j);
*/
            Loop++;
        }
    }

    CutPlane.InUseX = Which + 1;

    (void)EnableCutPlanePlotStateX();

    return(0);
}
#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int PlotCutPlaneX()
/****************************************************************************/
{
    static int   i;
    static int   j;
    static int   Loop;
    static int   Jump1;
    static int   Jump2;
    static int   Jump3;
    static int   Jump4;
    static int   Jump5;
    static int   Jump6;
    static int   ColorMapping;
    static float Temp1;
    static float Temp2;
    static float Temp3;
    static float Temp4;
    static float Xc;
    static float Yc;
    static float Zc;
    static float Red;
    static float Green;
    static float Blue;

    if(!GetCutPlanePlotStateX())
        return(0);

    ColorMapping = gomp_GetColorMappingType();
    /*    glDisable(GL_LIGHTING);*/
    glEnable(GL_LIGHTING);
    glColor3f(  1.0 , 1.0 , 1.0);
    glNormal3f( 1.0 , 0.0 , 0.0);

    if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
        (void)gomp_Prepare1DTexture();
        glEnable(GL_TEXTURE_1D);
    }
    else {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
    }

    for(j = 0 ; j < CutPlane.pzX - 1   ; j++)    {

        Jump1 =  j      * CutPlane.pyX;
        Jump2 =  Jump1  + CutPlane.pyX;

        for( i = 0 ; i < CutPlane.pyX - 1; i++) {

            Jump3 = 3 * (Jump1 + i);
            Jump4 =      Jump3 + 3;
            Jump6 = 3 * (Jump2 + i);
            Jump5 =      Jump6 + 3;

            Temp1 = CutPlane.DataX[Jump1 + i    ];
            Temp2 = CutPlane.DataX[Jump1 + i + 1];
            Temp3 = CutPlane.DataX[Jump2 + i + 1];
            Temp4 = CutPlane.DataX[Jump2 + i    ];

            glPushMatrix();
 
            if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {

/* strange but this does not work on my HP OmniBook XE3
   glBegin( GL_QUADS );
*/
                glBegin( GL_POLYGON );

                Xc =  CutPlane.xyzCoordsX[ Jump3    ];
                Yc =  CutPlane.xyzCoordsX[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsX[ Jump3 + 2];

                glTexCoord1f( Temp1 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump4    ];
                Yc =  CutPlane.xyzCoordsX[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsX[ Jump4 + 2];

                glTexCoord1f( Temp2 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump5    ];
                Yc =  CutPlane.xyzCoordsX[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsX[ Jump5 + 2];

                glTexCoord1f( Temp3 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump6    ];
                Yc =  CutPlane.xyzCoordsX[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsX[ Jump6 + 2];

                glTexCoord1f( Temp4 );
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            } else {

                glBegin( GL_QUADS );

                Xc =  CutPlane.xyzCoordsX[ Jump3    ];
                Yc =  CutPlane.xyzCoordsX[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsX[ Jump3 + 2];

                gomp_PreRainbow((double)Temp1 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump4    ];
                Yc =  CutPlane.xyzCoordsX[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsX[ Jump4 + 2];

                gomp_PreRainbow((double)Temp2 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump5    ];
                Yc =  CutPlane.xyzCoordsX[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsX[ Jump5 + 2];

                gomp_PreRainbow((double)Temp3 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump6    ];
                Yc =  CutPlane.xyzCoordsX[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsX[ Jump6 + 2];

                gomp_PreRainbow((double)Temp4 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            }

            glPopMatrix();

            Loop++;
        }
    }

    glDisable( GL_TEXTURE_1D );

    return(0);
}

/****************************************************************************/
int PlotCutPlane_3D_X()
/****************************************************************************/
{
    static int   i;
    static int   j;
    static int   Loop;
    static int   Jump1;
    static int   Jump2;
    static int   Jump3;
    static int   Jump4;
    static int   Jump5;
    static int   Jump6;
    static int   ColorMapping;
    static float Temp1;
    static float Temp2;
    static float Temp3;
    static float Temp4;
    static float Temp5;
    static float Xc;
    static float Yc;
    static float Zc;
    static float Red;
    static float Green;
    static float Blue;
    static float AmpDamping = 1.0;

    if(!GetCutPlanePlotStateX())
        return(0);

    AmpDamping = gomp_GetCutPlaneDamping();

    ColorMapping = gomp_GetColorMappingType();
    /*    glDisable(GL_LIGHTING);*/
    glEnable(GL_LIGHTING);
    glColor3f(  1.0 , 1.0 , 1.0);
    glNormal3f( 1.0 , 0.0 , 0.0);

    if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
        (void)gomp_Prepare1DTexture();
        glEnable(GL_TEXTURE_1D);
    }
    else {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
    }

    Temp5  = CutPlane.MinDX / (CutPlane.MaxDX - CutPlane.MinDX);

    for(j = 0 ; j < CutPlane.pzX - 1   ; j++)    {

        Jump1 =  j      * CutPlane.pyX;
        Jump2 =  Jump1  + CutPlane.pyX;

        for( i = 0 ; i < CutPlane.pyX - 1; i++) {

            Jump3 = 3 * (Jump1 + i);
            Jump4 =      Jump3 + 3;
            Jump6 = 3 * (Jump2 + i);
            Jump5 =      Jump6 + 3;

            Temp1 = CutPlane.DataX[Jump1 + i    ];
            Temp2 = CutPlane.DataX[Jump1 + i + 1];
            Temp3 = CutPlane.DataX[Jump2 + i + 1];
            Temp4 = CutPlane.DataX[Jump2 + i    ];

            glPushMatrix();
 
            if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {

/* strange but this does not work on my HP OmniBook XE3
   glBegin( GL_QUADS );
*/
                glBegin( GL_POLYGON );

                Xc =  CutPlane.xyzCoordsX[ Jump3    ] + AmpDamping * (Temp1 + Temp5);
                Yc =  CutPlane.xyzCoordsX[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsX[ Jump3 + 2];

                glTexCoord1f( Temp1 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump4    ] + AmpDamping * (Temp2 + Temp5);
                Yc =  CutPlane.xyzCoordsX[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsX[ Jump4 + 2];

                glTexCoord1f( Temp2 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump5    ] + AmpDamping * (Temp3 + Temp5);
                Yc =  CutPlane.xyzCoordsX[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsX[ Jump5 + 2];

                glTexCoord1f( Temp3 );
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump6    ] + AmpDamping * (Temp4 + Temp5);
                Yc =  CutPlane.xyzCoordsX[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsX[ Jump6 + 2];

                glTexCoord1f( Temp4 );
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            } else {

                glBegin( GL_QUADS );

                Xc =  CutPlane.xyzCoordsX[ Jump3    ] + AmpDamping * (Temp1 + Temp5);
                Yc =  CutPlane.xyzCoordsX[ Jump3 + 1];
                Zc =  CutPlane.xyzCoordsX[ Jump3 + 2];

                gomp_PreRainbow((double)Temp1 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump4    ] + AmpDamping * (Temp2 + Temp5);
                Yc =  CutPlane.xyzCoordsX[ Jump4 + 1];
                Zc =  CutPlane.xyzCoordsX[ Jump4 + 2];

                gomp_PreRainbow((double)Temp2 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump5    ] + AmpDamping * (Temp3 + Temp5);
                Yc =  CutPlane.xyzCoordsX[ Jump5 + 1]; 
                Zc =  CutPlane.xyzCoordsX[ Jump5 + 2];

                gomp_PreRainbow((double)Temp3 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                Xc =  CutPlane.xyzCoordsX[ Jump6    ] + AmpDamping * (Temp4 + Temp5);
                Yc =  CutPlane.xyzCoordsX[ Jump6 + 1]; 
                Zc =  CutPlane.xyzCoordsX[ Jump6 + 2];

                gomp_PreRainbow((double)Temp4 , &Red , &Green , &Blue);
                glColor3f(Red , Green , Blue);
                glVertex3f( Xc , Yc , Zc );

                glEnd();
            }

            glPopMatrix();

            Loop++;
        }
    }

    glDisable( GL_TEXTURE_1D );

    return(0);
}

/* 64 256 */
#define LEVELS    64
#define LEVELS4  256

/****************************************************************************/
int gomp_Prepare1DTexture()
/****************************************************************************/ 
{
    static int i;
    static int Levels;
    static float gs;

    static unsigned char ColorMap[LEVELS4][3];

    if(Check) {
        glTexImage1D( GL_TEXTURE_1D,0,3,256,0,GL_RGB,GL_UNSIGNED_BYTE,ColorMap );
        return(0);
    }

    Check = 1;

    Levels = LEVELS;

    if(gomp_GetDisplayColourType()) { /* colour type */
#if defined(JUNK)
/* from blue to red through white ... */
        for( i = 0 ; i < Levels; i++ )
        {
            ColorMap[i][2]      = 255;
            ColorMap[i][1]      = 4*i;
            ColorMap[i][0]      = 0;

            ColorMap[Levels + i][2] = 255;
            ColorMap[Levels + i][1] = 255;
            ColorMap[Levels + i][0] = 4*i;

            ColorMap[2 * Levels + i][2] = 255 - 4*i;
            ColorMap[2 * Levels + i][1] = 255;
            ColorMap[2 * Levels + i][0] = 255;

            ColorMap[3 * Levels + i][2] = 0;
            ColorMap[3 * Levels + i][1] = 255 - 4*i;
            ColorMap[3 * Levels + i][0] = 255;
        }

        ColorMap[255][2] = 0;
        ColorMap[255][1] = 0;
        ColorMap[255][0] = 255;
#endif
/* from blue to red through green ... */
        for( i = 0 ; i < Levels; i++ )
        {
            ColorMap[i][2]      = 255;
            ColorMap[i][1]      = 4*i;
            ColorMap[i][0]      = 0;

            ColorMap[Levels + i][2] = 255 - 4*i;
            ColorMap[Levels + i][1] = 255;
            ColorMap[Levels + i][0] = 0;

            ColorMap[2 * Levels + i][2] = 0;
            ColorMap[2 * Levels + i][1] = 255;
            ColorMap[2 * Levels + i][0] = 4*i;

            ColorMap[3 * Levels + i][2] = 0;
            ColorMap[3 * Levels + i][1] = 255 - 4*i;
            ColorMap[3 * Levels + i][0] = 255;
        }

        ColorMap[255][2] = 0;
        ColorMap[255][1] = 0;
        ColorMap[255][0] = 255;
    } else { /* grayscale */
        for( i = 0 ; i < Levels; i++ )
        {
            gs  = (0.299 * (float)0       + 
                   0.587 * (float)(4 * i) + 
                   0.114 * (float)255);
            ColorMap[i][2]      = (int)gs;
            ColorMap[i][1]      = (int)gs;
            ColorMap[i][0]      = (int)gs;

            gs  = (0.299 * (float)0       + 
                   0.587 * (float)(255)   + 
                   0.114 * (float)(255 - 4 * i));
            ColorMap[Levels + i][2] = (int)gs;
            ColorMap[Levels + i][1] = (int)gs;
            ColorMap[Levels + i][0] = (int)gs;

            gs  = (0.299 * (float)(4 * i) + 
                   0.587 * (float)(255)   + 
                   0.114 * (float)(0));
            ColorMap[2 * Levels + i][2] = (int)gs;
            ColorMap[2 * Levels + i][1] = (int)gs;
            ColorMap[2 * Levels + i][0] = (int)gs;

            gs  = (0.299 * (float)(255)          + 
                   0.587 * (float)(255 - 4 * i)  + 
                   0.114 * (float)(0));
            ColorMap[3 * Levels + i][2] = (int)gs;
            ColorMap[3 * Levels + i][1] = (int)gs;
            ColorMap[3 * Levels + i][0] = (int)gs;
        }

        ColorMap[255][2] = (int)(0.299 * 255.);
        ColorMap[255][1] = (int)(0.299 * 255.);
        ColorMap[255][0] = (int)(0.299 * 255.);
    }

    glTexParameterf( GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
    glTexParameterf( GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameterf( GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP );

    glTexEnvf( GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE );
    glTexImage1D( GL_TEXTURE_1D,0,3,256,0,GL_RGB,GL_UNSIGNED_BYTE,ColorMap );

    /*    glEnable( GL_TEXTURE_1D );*/

    return(1);
}
#endif /* ENABLE_GRAPHICS */
/****************************************************************************/
static int PlotCutPlane(void* userData,int Wstr,int drawFlags)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    if ( ! ( drawFlags & gom_PlotComplexElements ) || Wstr != 0 )
        return(-1);

    if(GetCutPlanePlotStateX()) {
        if(gomp_GetCutPlaneType_3D_X())
            PlotCutPlane_3D_X();
        else
            PlotCutPlaneX();
    }
    if(GetCutPlanePlotStateY()) {
        if(gomp_GetCutPlaneType_3D_Y())
            PlotCutPlane_3D_Y();
        else
            PlotCutPlaneY();
    }
    if(GetCutPlanePlotStateZ()) {
        if(gomp_GetCutPlaneType_3D_Z())
            PlotCutPlane_3D_Z();
        else
            PlotCutPlaneZ();
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/****************************************************************************/
int    GetCutPlanePlotStateX()
/****************************************************************************/
{
    return(CutPlane.PlotStateX);
}

/****************************************************************************/
int    EnableCutPlanePlotStateX()
/****************************************************************************/
{
    CutPlane.PlotStateX = 1;
    
    if( !CutPlaneCallbackHandle )
        CutPlaneCallbackHandle = gomp_RegisterPlotter(
            PlotCutPlane,NULL,PLOTTER_NAME_CUTPLANE,PLOTTER_ORDER_CUTPLANE);

    return(0);
}

/****************************************************************************/
int    DisableCutPlanePlotStateX()
/****************************************************************************/
{
    CutPlane.PlotStateX = 0;
    
    if( CutPlaneCallbackHandle &&
        !GetCutPlanePlotStateY() &&
        !GetCutPlanePlotStateZ() ) {
        gomp_UnregisterPlotter(CutPlaneCallbackHandle);
        CutPlaneCallbackHandle = NULL;
    }
    
    return(0);
}
/****************************************************************************/
int    gomp_DeleteCutPlaneDataX()
/****************************************************************************/
{

    if(CutPlane.InUseX) {
        free(CutPlane.DataX);
        free(CutPlane.xyzCoordsX);
    }
    /*
      else {
      gomp_PrintERROR("no (X) cut plane data available");
      }
    */

    CutPlane.PlotStateX  = 0;
    CutPlane.InUseX      = 0;
    CutPlane.PlotManX    = 0;

    DisableCutPlanePlotStateX();

    return(0);
}

/****************************************************************************/
int    GetCutPlanePlotStateY()
/****************************************************************************/
{
    return(CutPlane.PlotStateY);
}

/****************************************************************************/
int    EnableCutPlanePlotStateY()
/****************************************************************************/
{
    CutPlane.PlotStateY = 1;

    if( !CutPlaneCallbackHandle )
        CutPlaneCallbackHandle = gomp_RegisterPlotter(
            PlotCutPlane,NULL,PLOTTER_NAME_CUTPLANE,PLOTTER_ORDER_CUTPLANE);

    return(0);
}

/****************************************************************************/
int    DisableCutPlanePlotStateY()
/****************************************************************************/
{
    CutPlane.PlotStateY = 0;

    if( CutPlaneCallbackHandle &&
        !GetCutPlanePlotStateX() &&
        !GetCutPlanePlotStateZ() ) {
        gomp_UnregisterPlotter(CutPlaneCallbackHandle);
        CutPlaneCallbackHandle = NULL;
    }

    return(0);
}

/****************************************************************************/
int    gomp_DeleteCutPlaneDataY()
/****************************************************************************/
{

    if(CutPlane.InUseY) {
        free(CutPlane.DataY);
        free(CutPlane.xyzCoordsY);
    }
    /*
      else {
      gomp_PrintERROR("no (Y) cut plane data available");
      }
    */

    CutPlane.PlotStateY  = 0;
    CutPlane.InUseY      = 0;
    CutPlane.PlotManY    = 0;

    DisableCutPlanePlotStateY();

    return(0);
}

/****************************************************************************/
int    GetCutPlanePlotStateZ()
/****************************************************************************/
{
    return(CutPlane.PlotStateZ);
}

/****************************************************************************/
int    EnableCutPlanePlotStateZ()
/****************************************************************************/
{
    CutPlane.PlotStateZ = 1;

    if( !CutPlaneCallbackHandle )
        CutPlaneCallbackHandle = gomp_RegisterPlotter(
            PlotCutPlane,NULL,PLOTTER_NAME_CUTPLANE,PLOTTER_ORDER_CUTPLANE);

    return(0);
}

/****************************************************************************/
int    DisableCutPlanePlotStateZ()
/****************************************************************************/
{
    CutPlane.PlotStateZ = 0;

    if( CutPlaneCallbackHandle &&
        !GetCutPlanePlotStateX() &&
        !GetCutPlanePlotStateY() ) {
        gomp_UnregisterPlotter(CutPlaneCallbackHandle);
        CutPlaneCallbackHandle = NULL;
    }

    return(0);
}

/****************************************************************************/
int    gomp_DeleteCutPlaneDataZ()
/****************************************************************************/
{

    if(CutPlane.InUseZ) {
        free(CutPlane.DataZ);
        free(CutPlane.xyzCoordsZ);
    }
    /*
      else {
      gomp_PrintERROR("no (Z) cut plane data available");
      }
    */

    CutPlane.PlotStateZ  = 0;
    CutPlane.InUseZ      = 0;
    CutPlane.PlotManZ    = 0;

    DisableCutPlanePlotStateZ();

    return(0);
}

/****************************************************************************/
int    gomp_GetCutPlanePlotStateXYZ(int Which)
/****************************************************************************/
{
    switch(Which) {
    case 1:
        return(CutPlane.PlotStateXYZ1);
    case 2:
        return(CutPlane.PlotStateXYZ2);
    case 3:
        return(CutPlane.PlotStateXYZ3);
    }

    return(0);
}
#if 0
/****************************************************************************/
int    EnableCutPlanePlotStateXYZ(int Which)
/****************************************************************************/
{
    switch(Which) {
    case 1:
        CutPlane.PlotStateXYZ1 = 1;
        break;
    case 2:
        CutPlane.PlotStateXYZ2 = 1;
        break;
    case 3:
        CutPlane.PlotStateXYZ3 = 1;
        break;
    }
    return(0);
}
#endif
/****************************************************************************/
int    gomp_DisableCutPlanePlotStateXYZ(int Which)
/****************************************************************************/
{
    switch(Which) {
    case 1:
        CutPlane.PlotStateXYZ1 = 0;
        break;
    case 2:
        CutPlane.PlotStateXYZ2 = 0;
        break;
    case 3:
        CutPlane.PlotStateXYZ3 = 0;
        break;
    }
    return(0);
}

/****************************************************************************/
int    gomp_WriteCutPlane2ModelFile(FILE *Model_f)
/****************************************************************************/
{
    const float *sumxyz;

    sumxyz = gomp_GetTranslateArray();

    if(CutPlane.InUseX) {
/* X-cut plane - tag */
        fprintf(Model_f , "[Cut plane X]\n");
        fprintf(Model_f , "%d %d %d %f %f %f\n",CutPlane.InUseX,
                CutPlane.PlotStateX,
                CutPlane.PlotManX,
                CutPlane.MinDX,
                CutPlane.MaxDX,
                CutPlane.xyzCoordsX[0] +
                sumxyz[0]);
    } /* end of x-plane */

    if(CutPlane.InUseY) {
/* Y-cut plane - tag */
        fprintf(Model_f , "[Cut plane Y]\n");
        fprintf(Model_f , "%d %d %d %f %f %f\n",CutPlane.InUseY,
                CutPlane.PlotStateY,
                CutPlane.PlotManY,
                CutPlane.MinDY,
                CutPlane.MaxDY,
                CutPlane.xyzCoordsY[1] +
                sumxyz[1]);
    } /* end of y-cut plane */

    if(CutPlane.InUseZ) {
/* Z-cut plane - tag */
        fprintf(Model_f , "[Cut plane Z]\n");
        fprintf(Model_f , "%d %d %d %f %f %f\n",CutPlane.InUseZ,
                CutPlane.PlotStateZ,
                CutPlane.PlotManZ,
                CutPlane.MinDZ,
                CutPlane.MaxDZ,
                CutPlane.xyzCoordsZ[2] +
                sumxyz[2]);
    } /* end of z-cut plane */
/* XYZ -plane */
    if(CutPlane.PlotStateXYZ1 || CutPlane.PlotStateXYZ2 || CutPlane.PlotStateXYZ3) {
        fprintf(Model_f , "[Cut plane XYZ]\n");
        fprintf(Model_f , "%d %d %f %f %f %f %f %f\n",CutPlane.PlotStateXYZ1,
                CutPlane.PlotWhichXYZ1,
                CutPlane.PlotPlaneXYZ1[0],
                CutPlane.PlotPlaneXYZ1[1],
                CutPlane.PlotPlaneXYZ1[2],
                CutPlane.PlotPlaneXYZ1[3],
                CutPlane.PlotScaleMinXYZ1,
                CutPlane.PlotScaleMaxXYZ1);
        fprintf(Model_f , "%d %d %f %f %f %f %f %f\n",CutPlane.PlotStateXYZ2,
                CutPlane.PlotWhichXYZ2,
                CutPlane.PlotPlaneXYZ2[0],
                CutPlane.PlotPlaneXYZ2[1],
                CutPlane.PlotPlaneXYZ2[2],
                CutPlane.PlotPlaneXYZ2[3],
                CutPlane.PlotScaleMinXYZ2,
                CutPlane.PlotScaleMaxXYZ2);
        fprintf(Model_f , "%d %d %f %f %f %f %f %f\n",CutPlane.PlotStateXYZ3,
                CutPlane.PlotWhichXYZ3,
                CutPlane.PlotPlaneXYZ3[0],
                CutPlane.PlotPlaneXYZ3[1],
                CutPlane.PlotPlaneXYZ3[2],
                CutPlane.PlotPlaneXYZ3[3],
                CutPlane.PlotScaleMinXYZ3,
                CutPlane.PlotScaleMaxXYZ3);
    }

    return(0);
}
/****************************************************************************/
int    gomp_ReadCutPlaneXFromModelFile(FILE *Model_f)
/****************************************************************************/
{
    char  InputText[BUFF_LEN];
    char  Text1[BUFF_LEN];
    char  Text2[BUFF_LEN];
    char  Text3[BUFF_LEN];
    float Temp1;
    int   InUseX,PlotStateX,ManX;
    float MinDX,MaxDX;

/* start the process ... */

    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText ,"%d %d %d %f %f %f", &InUseX,
           &PlotStateX,
           &ManX,
           &MinDX,
           &MaxDX,
           &Temp1);

    sprintf(Text1,"%f",MinDX);
    sprintf(Text2,"%f",MaxDX);

    switch(ManX) {
    case CUTPLANE_NO_ACTION:
        sprintf(Text3,"%s"," ");
        break;
    case CUTPLANE_SQRT_ACTION:
        sprintf(Text3,"%s","sqrt");
        break;
    case CUTPLANE_LOG10_ACTION:
        sprintf(Text3,"%s","log10");
        break;
    }

    return(gomp_PrepareCutPlaneX(InUseX - 1 , Temp1 , Text1 , Text2 , Text3));
}

/****************************************************************************/
int    gomp_ReadCutPlaneYFromModelFile(FILE *Model_f)
/****************************************************************************/
{
    char  InputText[BUFF_LEN];
    char  Text1[BUFF_LEN];
    char  Text2[BUFF_LEN];
    char  Text3[BUFF_LEN];
    float Temp1;
    int   InUseY,PlotStateY,ManY;
    float MinDY,MaxDY;

/* start the process ... */

    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText ,"%d %d %d %f %f %f", &InUseY,
           &PlotStateY,
           &ManY,
           &MinDY,
           &MaxDY,
           &Temp1);

    sprintf(Text1,"%f",MinDY);
    sprintf(Text2,"%f",MaxDY);

    switch(ManY) {
    case CUTPLANE_NO_ACTION:
        sprintf(Text3,"%s"," ");
        break;
    case CUTPLANE_SQRT_ACTION:
        sprintf(Text3,"%s","sqrt");
        break;
    case CUTPLANE_LOG10_ACTION:
        sprintf(Text3,"%s","log10");
        break;
    }

    return(gomp_PrepareCutPlaneY(InUseY - 1 , Temp1 , Text1 , Text2 , Text3));
}

/****************************************************************************/
int    gomp_ReadCutPlaneZFromModelFile(FILE *Model_f)
/****************************************************************************/
{
    char  InputText[BUFF_LEN];
    char  Text1[BUFF_LEN];
    char  Text2[BUFF_LEN];
    char  Text3[BUFF_LEN];
    float Temp1;
    int   InUseZ,PlotStateZ,ManZ;
    float MinDZ,MaxDZ;

/* start the process ... */

    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText ,"%d %d %d %f %f %f", &InUseZ,
           &PlotStateZ,
           &ManZ,
           &MinDZ,
           &MaxDZ,
           &Temp1);

    sprintf(Text1,"%f",MinDZ);
    sprintf(Text2,"%f",MaxDZ);

    switch(ManZ) {
    case CUTPLANE_NO_ACTION:
        sprintf(Text3,"%s"," ");
        break;
    case CUTPLANE_SQRT_ACTION:
        sprintf(Text3,"%s","sqrt");
        break;
    case CUTPLANE_LOG10_ACTION:
        sprintf(Text3,"%s","log10");
        break;
    }

    return(gomp_PrepareCutPlaneZ(InUseZ - 1 , Temp1 , Text1 , Text2 , Text3));
}
#if 0
/****************************************************************************/
float  GetCutPlaneYZmin(int Which)
/****************************************************************************/
{
    if(CutPlane.InUseX) {
        return(CutPlane.MinDX);
    }
    else {
        return(ContourInfo[Which].min);
    }
}
/****************************************************************************/
float  GetCutPlaneYZmax(int Which)
/****************************************************************************/
{
    if(CutPlane.InUseX) {
        return(CutPlane.MaxDX);
    }
    else {
        return(ContourInfo[Which].max);
    }
}

/****************************************************************************/
float  GetCutPlaneXZmin(int Which)
/****************************************************************************/
{
    if(CutPlane.InUseY) {
        return(CutPlane.MinDY);
    }
    else {
        return(ContourInfo[Which].min);
    }
}
/****************************************************************************/
float  GetCutPlaneXZmax(int Which)
/****************************************************************************/
{
    if(CutPlane.InUseY) {
        return(CutPlane.MaxDY);
    }
    else {
        return(ContourInfo[Which].max);
    }
}

/****************************************************************************/
float  GetCutPlaneXYmin(int Which)
/****************************************************************************/
{
    if(CutPlane.InUseZ) {
        return(CutPlane.MinDZ);
    }
    else {
        return(ContourInfo[Which].min);
    }
}
/****************************************************************************/
float  GetCutPlaneXYmax(int Which)
/****************************************************************************/
{
    if(CutPlane.InUseZ) {
        return(CutPlane.MaxDZ);
    }
    else {
        return(ContourInfo[Which].max);
    }
}

/****************************************************************************/
float  GetCutPlaneXcoord()
/****************************************************************************/
{
    if(CutPlane.InUseX) {
        return(CutPlane.xyzCoordsX[0]);
    }
    else {
        return(0.0);
    }
}

/****************************************************************************/
float  GetCutPlaneYcoord()
/****************************************************************************/
{
    if(CutPlane.InUseY) {
        return(CutPlane.xyzCoordsX[1]);
    }
    else {
        return(0.0);
    }
}

/****************************************************************************/
float  GetCutPlaneZcoord()
/****************************************************************************/
{
    if(CutPlane.InUseZ) {
        return(CutPlane.xyzCoordsX[2]);
    }
    else {
        return(0.0);
    }
}
#endif
/****************************************************************************/
int    gomp_GetCutplaneSmoothTypeX()
/****************************************************************************/
{

    return(CutPlane.PlotManX);
}
/****************************************************************************/
int    gomp_GetCutplaneSmoothTypeY()
/****************************************************************************/
{

    return(CutPlane.PlotManY);
}
/****************************************************************************/
int    gomp_GetCutplaneSmoothTypeZ()
/****************************************************************************/
{

    return(CutPlane.PlotManZ);
}

/****************************************************************************/
int gomp_PrepareCutPlaneSpectrumZ(int Which   , const char *Text0,
                                const char *Text1 , const char *Text2,
                                const char *Text3)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    int    i;
    int    j;
    int    k;
    int    id;
    int    Jump1;
    float  Temp1;
    float  Temp2;
    float  Temp3;
    float  Delta;
    char   Text[BUFF_LEN];


    if(Which <  0 ||
       Which >= gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour data index out of range");
        return(1);
    }

    (void)gomp_Prepare1DTexture();

    if(CutPlaneSpectrum.InUse) {
        free(CutPlaneSpectrum.Data);
    }

    CutPlaneSpectrum.InUse = 0;

    if(strlen(Text0) != 0)
        CutPlaneSpectrum.NumBins = atoi(Text0);
    else
        CutPlaneSpectrum.NumBins = NUMBER_BINS;

    if(strlen(Text1) != 0)
        CutPlaneSpectrum.Min = atof(Text1);
    else
        CutPlaneSpectrum.Min = ContourInfo[Which].min;

    if(strlen(Text2) != 0)
        CutPlaneSpectrum.Max = atof(Text2);
    else
        CutPlaneSpectrum.Max = ContourInfo[Which].max;

    if(strlen(Text3) != 0)
        CutPlaneSpectrum.BinCut = atof(Text3);
    else
        CutPlaneSpectrum.BinCut = 1.0;

    sprintf(Text,"Number of bins: %d",CutPlaneSpectrum.NumBins);
    gomp_PrintMessage(Text);
    sprintf(Text,"Cutting at bin value: %f",CutPlaneSpectrum.BinCut);
    gomp_PrintMessage(Text);
    sprintf(Text,"Min value: %e (blue)",CutPlaneSpectrum.Min);
    gomp_PrintMessage(Text);
    sprintf(Text,"Max value: %e (red)", CutPlaneSpectrum.Max);
    gomp_PrintMessage(Text);

    CutPlaneSpectrum.px   = ContourInfo[Which].xdim;
    CutPlaneSpectrum.py   = ContourInfo[Which].ydim;
    CutPlaneSpectrum.pz   = ContourInfo[Which].zdim;

    CutPlaneSpectrum.MinX = ContourInfo[Which].Xmin;
    CutPlaneSpectrum.MaxX = ContourInfo[Which].Xmax;

    CutPlaneSpectrum.MinY = ContourInfo[Which].Ymin;
    CutPlaneSpectrum.MaxY = ContourInfo[Which].Ymax;

    CutPlaneSpectrum.MinZ = ContourInfo[Which].Zmin;
    CutPlaneSpectrum.MaxZ = ContourInfo[Which].Zmax;


    CutPlaneSpectrum.Data      = gomp_AllocateFloatVector(CutPlaneSpectrum.pz * 
                                           CutPlaneSpectrum.NumBins);
    for(i = 0 ; i < CutPlaneSpectrum.pz * CutPlaneSpectrum.NumBins ; i++)
        CutPlaneSpectrum.Data[i] = 0.0;

    Temp1  = CutPlaneSpectrum.Max - CutPlaneSpectrum.Min;
    Delta  = 1.0 / (float)(CutPlaneSpectrum.NumBins - 1);

    CutPlaneSpectrum.BinMin  = 0.0;
    CutPlaneSpectrum.BinMax  = 0.0;

    for(k = 0 ; k < CutPlaneSpectrum.pz ; k++)     {

        Jump1  = k * CutPlaneSpectrum.px * CutPlaneSpectrum.py;
 
        for(j = 0 ; j < CutPlaneSpectrum.py ; j++)     {

            for( i = 0 ; i < CutPlaneSpectrum.px ; i++)    {

                Temp2 = ContourInfo[Which].data[Jump1 + j * CutPlaneSpectrum.px + i];

                Temp3 = (Temp2 - CutPlaneSpectrum.Min) / Temp1;

                id = (int)(Temp3/Delta);

                if(id < (CutPlaneSpectrum.NumBins - 1)) {
                    CutPlaneSpectrum.Data[id + k * CutPlaneSpectrum.NumBins] += 1.0;
                    Temp2 = CutPlaneSpectrum.Data[id + k * CutPlaneSpectrum.NumBins];
                    if(Temp2 > CutPlaneSpectrum.BinMax)
                        CutPlaneSpectrum.BinMax = Temp2;
                }
            }
        }
    }

    CutPlaneSpectrum.InUse = Which + 1;

    if(CutPlaneSpectrum.BinCut > CutPlaneSpectrum.BinMax) {
        gomp_PrintERROR("Bin cut value > Bin max value, will put cut value = max value");
        CutPlaneSpectrum.BinCut = CutPlaneSpectrum.BinMax;
    }

    (void)PlotCutPlaneSpectrumZ();
#endif /* ENABLE_GRAPHICS */

    return(0);
}
#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int PlotCutPlaneSpectrumZ()
/****************************************************************************/
{
    static int   i;
    static int   j;
    static int   Jump1;
    static int   Jump2;
    static float Temp1;
    static float Temp2;
    static float Temp3;
    static float Temp4;
    static float Temp5;
    static float Xc;
    static float Yc;
    static int   mm;
    static char  Text[BUFF_LEN];

    /*    glDisable(GL_LIGHTING);*/
 
    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);

    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D((GLdouble)0.0, (GLdouble)CutPlaneSpectrum.NumBins - 1.0, 
               (GLdouble)0.0, (GLdouble)CutPlaneSpectrum.pz      - 1.0);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(mm);

    glColor3f(  1.0 , 1.0 , 1.0);
    glNormal3f( 0.0 , 0.0 , 1.0);

/* clear background */
    glClearColor( 1.0 , 1.0 , 1.0 , 0.0);

    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
/* ....... done ............. */

/* add some text */
    glColor3fv(gomp_BLACKv);

    glPushMatrix();
    glRasterPos3f(0.5, (GLfloat)CutPlaneSpectrum.pz - 2. , 0.1);
    sprintf(Text,"<Range: %.4e ... %.4e>",
            CutPlaneSpectrum.MinZ,CutPlaneSpectrum.MaxZ);
    gomp_PrintString(Text , "*");
    glPopMatrix();

    glRasterPos3f((GLfloat)CutPlaneSpectrum.NumBins/3., 0.5 , 0.1);
    sprintf(Text,"<Range: %.4e ... %.4e>",
            CutPlaneSpectrum.Min,CutPlaneSpectrum.Max);
    gomp_PrintString(Text , "*");

    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_1D);

    Temp5 =  CutPlaneSpectrum.BinCut - CutPlaneSpectrum.BinMin;
    
    for(j = 0 ; j < CutPlaneSpectrum.pz - 1   ; j++)    {

        Jump1 =  j     * CutPlaneSpectrum.NumBins;
        Jump2 =  Jump1 + CutPlaneSpectrum.NumBins;

        for( i = 0 ; i < CutPlaneSpectrum.NumBins - 1; i++) {

            Temp1 = 
                (CutPlaneSpectrum.Data[Jump1 + i    ] - CutPlaneSpectrum.BinMin) / Temp5;
            Temp2 = 
                (CutPlaneSpectrum.Data[Jump1 + i + 1] - CutPlaneSpectrum.BinMin) / Temp5;
            Temp3 = 
                (CutPlaneSpectrum.Data[Jump2 + i + 1] - CutPlaneSpectrum.BinMin) / Temp5;
            Temp4 = 
                (CutPlaneSpectrum.Data[Jump2 + i    ] - CutPlaneSpectrum.BinMin) / Temp5;

            glPushMatrix();
 
            glBegin( GL_QUADS );

            Xc =  (float)i;
            Yc =  (float)j;

            glTexCoord1f( Temp1 );
            glVertex2f( Xc , Yc );

            Xc =  (float)(i+1);
            Yc =  (float)j;

            glTexCoord1f( Temp2 );
            glVertex2f( Xc , Yc );

            Xc =  (float)(i+1);
            Yc =  (float)(j+1); 

            glTexCoord1f( Temp3 );
            glVertex2f( Xc , Yc );

            Xc =  (float)i;
            Yc =  (float)(j+1); 

            glTexCoord1f( Temp4 );
            glVertex2f( Xc , Yc );

            glEnd();

            glPopMatrix();

        }
    }

    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_1D);

    /*    glEnable(GL_LIGHTING);*/
    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(mm);

#if defined(GLUT)
    glutSwapBuffers();
#else
    auxSwapBuffers();
#endif

    return(0);
}
#endif /* ENABLE_GRAPHICS */
/*                  */
/****************************************************************************/
int gomp_PrepareCutPlaneSpectrumY(int Which   , const char *Text0,
                                const char *Text1 , const char *Text2,
                                const char *Text3)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    int    i;
    int    j;
    int    k;
    int    id;
    float  Temp1;
    float  Temp2;
    float  Temp3;
    float  Delta;
    char   Text[BUFF_LEN];


    if(Which <  0 ||
       Which >= gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour data index out of range");
        return(1);
    }

    (void)gomp_Prepare1DTexture();

    if(CutPlaneSpectrum.InUse) {
        free(CutPlaneSpectrum.Data);
    }

    CutPlaneSpectrum.InUse = 0;

    if(strlen(Text0) != 0)
        CutPlaneSpectrum.NumBins = atoi(Text0);
    else
        CutPlaneSpectrum.NumBins = NUMBER_BINS;

    if(strlen(Text1) != 0)
        CutPlaneSpectrum.Min = atof(Text1);
    else
        CutPlaneSpectrum.Min = ContourInfo[Which].min;

    if(strlen(Text2) != 0)
        CutPlaneSpectrum.Max = atof(Text2);
    else
        CutPlaneSpectrum.Max = ContourInfo[Which].max;

    if(strlen(Text3) != 0)
        CutPlaneSpectrum.BinCut = atof(Text3);
    else
        CutPlaneSpectrum.BinCut = 1.0;

    sprintf(Text,"Number of bins: %d",CutPlaneSpectrum.NumBins);
    gomp_PrintMessage(Text);
    sprintf(Text,"Cutting at bin value: %f",CutPlaneSpectrum.BinCut);
    gomp_PrintMessage(Text);
    sprintf(Text,"Min value: %e (blue)",CutPlaneSpectrum.Min);
    gomp_PrintMessage(Text);
    sprintf(Text,"Max value: %e (red)", CutPlaneSpectrum.Max);
    gomp_PrintMessage(Text);

    CutPlaneSpectrum.px   = ContourInfo[Which].xdim;
    CutPlaneSpectrum.py   = ContourInfo[Which].ydim;
    CutPlaneSpectrum.pz   = ContourInfo[Which].zdim;

    CutPlaneSpectrum.MinX = ContourInfo[Which].Xmin;
    CutPlaneSpectrum.MaxX = ContourInfo[Which].Xmax;

    CutPlaneSpectrum.MinY = ContourInfo[Which].Ymin;
    CutPlaneSpectrum.MaxY = ContourInfo[Which].Ymax;

    CutPlaneSpectrum.MinZ = ContourInfo[Which].Zmin;
    CutPlaneSpectrum.MaxZ = ContourInfo[Which].Zmax;


    CutPlaneSpectrum.Data      = gomp_AllocateFloatVector(CutPlaneSpectrum.py * 
                                           CutPlaneSpectrum.NumBins);
    for(i = 0 ; i < CutPlaneSpectrum.py * CutPlaneSpectrum.NumBins ; i++)
        CutPlaneSpectrum.Data[i] = 0.0;

    Temp1  = CutPlaneSpectrum.Max - CutPlaneSpectrum.Min;
    Delta  = 1.0 / (float)(CutPlaneSpectrum.NumBins - 1);

    CutPlaneSpectrum.BinMin  = 0.0;
    CutPlaneSpectrum.BinMax  = 0.0;

    for(k = 0 ; k < CutPlaneSpectrum.py ; k++)     {
 
        for(j = 0 ; j < CutPlaneSpectrum.pz ; j++)     {

            for( i = 0 ; i < CutPlaneSpectrum.px ; i++)    {

                Temp2 = ContourInfo[Which].data[
                    j * CutPlaneSpectrum.py * CutPlaneSpectrum.px + 
                    k * CutPlaneSpectrum.px + i];

                Temp3 = (Temp2 - CutPlaneSpectrum.Min) / Temp1;

                id = (int)(Temp3/Delta);

                if(id < (CutPlaneSpectrum.NumBins - 1)) {
                    CutPlaneSpectrum.Data[id + k * CutPlaneSpectrum.NumBins] += 1.0;
                    Temp2 = CutPlaneSpectrum.Data[id + k * CutPlaneSpectrum.NumBins];
                    if(Temp2 > CutPlaneSpectrum.BinMax)
                        CutPlaneSpectrum.BinMax = Temp2;
                }
            }
        }
    }

    CutPlaneSpectrum.InUse = Which + 1;

    if(CutPlaneSpectrum.BinCut > CutPlaneSpectrum.BinMax) {
        gomp_PrintERROR("Bin cut value > Bin max value, will put cut value = max value");
        CutPlaneSpectrum.BinCut = CutPlaneSpectrum.BinMax;
    }

    (void)PlotCutPlaneSpectrumY();
#endif /* ENABLE_GRAPHICS */

    return(0);
}
#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int PlotCutPlaneSpectrumY()
/****************************************************************************/
{
    static int   i;
    static int   j;
    static int   Jump1;
    static int   Jump2;
    static float Temp1;
    static float Temp2;
    static float Temp3;
    static float Temp4;
    static float Temp5;
    static float Xc;
    static float Yc;
    static int   mm;
    static char Text[BUFF_LEN];

    /*    glDisable(GL_LIGHTING);*/
 
    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);

    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D((GLdouble)0.0, (GLdouble)CutPlaneSpectrum.NumBins - 1.0, 
               (GLdouble)0.0, (GLdouble)CutPlaneSpectrum.py      - 1.0);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(mm);

    glColor3f(  1.0 , 1.0 , 1.0);
    glNormal3f( 0.0 , 0.0 , 1.0);

/* clear background */
    glClearColor( 1.0 , 1.0 , 1.0 , 0.0);

    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
/* ....... done ............. */

/* add some text */
    glColor3fv(gomp_BLACKv);

    glPushMatrix();
    glRasterPos3f(0.5, (GLfloat)CutPlaneSpectrum.py - 2. , 0.1);
    sprintf(Text,"<Range: %.4e ... %.4e>",
            CutPlaneSpectrum.MinY,CutPlaneSpectrum.MaxY);
    gomp_PrintString(Text , "*");
    glPopMatrix();

    glRasterPos3f((GLfloat)CutPlaneSpectrum.NumBins/3., 0.5 , 0.1);
    sprintf(Text,"<Range: %.4e ... %.4e>",
            CutPlaneSpectrum.Min,CutPlaneSpectrum.Max);
    gomp_PrintString(Text , "*");

    glEnable(GL_LIGHTING);
    glEnable( GL_TEXTURE_1D );

    Temp5 =  CutPlaneSpectrum.BinCut - CutPlaneSpectrum.BinMin;

    for(j = 0 ; j < CutPlaneSpectrum.py - 1   ; j++)    {

        Jump1 =  j     * CutPlaneSpectrum.NumBins;
        Jump2 =  Jump1 + CutPlaneSpectrum.NumBins;

        for( i = 0 ; i < CutPlaneSpectrum.NumBins - 1; i++) {

            Temp1 = 
                (CutPlaneSpectrum.Data[Jump1 + i    ] - CutPlaneSpectrum.BinMin) / Temp5;
            Temp2 = 
                (CutPlaneSpectrum.Data[Jump1 + i + 1] - CutPlaneSpectrum.BinMin) / Temp5;
            Temp3 = 
                (CutPlaneSpectrum.Data[Jump2 + i + 1] - CutPlaneSpectrum.BinMin) / Temp5;
            Temp4 = 
                (CutPlaneSpectrum.Data[Jump2 + i    ] - CutPlaneSpectrum.BinMin) / Temp5;

            glPushMatrix();
 
            glBegin( GL_QUADS );

            Xc =  (float)i;
            Yc =  (float)j;

            glTexCoord1f( Temp1 );
            glVertex2f( Xc , Yc );

            Xc =  (float)(i+1);
            Yc =  (float)j;

            glTexCoord1f( Temp2 );
            glVertex2f( Xc , Yc );

            Xc =  (float)(i+1);
            Yc =  (float)(j+1); 

            glTexCoord1f( Temp3 );
            glVertex2f( Xc , Yc );

            Xc =  (float)i;
            Yc =  (float)(j+1); 

            glTexCoord1f( Temp4 );
            glVertex2f( Xc , Yc );

            glEnd();

            glPopMatrix();

        }
    }

    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_1D);
    /*    glEnable(GL_LIGHTING);*/
    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(mm);

#if defined(GLUT)
    glutSwapBuffers();
#else
    auxSwapBuffers();
#endif
    return(0);
}
#endif /* ENABLE_GRAPHICS */
/*        */
/****************************************************************************/
int gomp_PrepareCutPlaneSpectrumX(int Which   , const char *Text0,
                                const char *Text1 , const char *Text2,
                                const char *Text3)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    int    i;
    int    j;
    int    k;
    int    id;
    float  Temp1;
    float  Temp2;
    float  Temp3;
    float  Delta;
    char   Text[BUFF_LEN];


    if(Which <  0 ||
       Which >= gomp_GetContoursDefined()) {
        gomp_PrintERROR("contour data index out of range");
        return(1);
    }

    (void)gomp_Prepare1DTexture();

    if(CutPlaneSpectrum.InUse) {
        free(CutPlaneSpectrum.Data);
    }

    CutPlaneSpectrum.InUse = 0;

    if(strlen(Text0) != 0)
        CutPlaneSpectrum.NumBins = atoi(Text0);
    else
        CutPlaneSpectrum.NumBins = NUMBER_BINS;

    if(strlen(Text1) != 0)
        CutPlaneSpectrum.Min = atof(Text1);
    else
        CutPlaneSpectrum.Min = ContourInfo[Which].min;

    if(strlen(Text2) != 0)
        CutPlaneSpectrum.Max = atof(Text2);
    else
        CutPlaneSpectrum.Max = ContourInfo[Which].max;

    if(strlen(Text3) != 0)
        CutPlaneSpectrum.BinCut = atof(Text3);
    else
        CutPlaneSpectrum.BinCut = 1.0;

    sprintf(Text,"Number of bins: %d",CutPlaneSpectrum.NumBins);
    gomp_PrintMessage(Text);
    sprintf(Text,"Cutting at bin value: %f",CutPlaneSpectrum.BinCut);
    gomp_PrintMessage(Text);
    sprintf(Text,"Min value: %e (blue)",CutPlaneSpectrum.Min);
    gomp_PrintMessage(Text);
    sprintf(Text,"Max value: %e (red)", CutPlaneSpectrum.Max);
    gomp_PrintMessage(Text);

    CutPlaneSpectrum.px   = ContourInfo[Which].xdim;
    CutPlaneSpectrum.py   = ContourInfo[Which].ydim;
    CutPlaneSpectrum.pz   = ContourInfo[Which].zdim;

    CutPlaneSpectrum.MinX = ContourInfo[Which].Xmin;
    CutPlaneSpectrum.MaxX = ContourInfo[Which].Xmax;

    CutPlaneSpectrum.MinY = ContourInfo[Which].Ymin;
    CutPlaneSpectrum.MaxY = ContourInfo[Which].Ymax;

    CutPlaneSpectrum.MinZ = ContourInfo[Which].Zmin;
    CutPlaneSpectrum.MaxZ = ContourInfo[Which].Zmax;


    CutPlaneSpectrum.Data      = gomp_AllocateFloatVector(CutPlaneSpectrum.px * 
                                           CutPlaneSpectrum.NumBins);
    for(i = 0 ; i < CutPlaneSpectrum.px * CutPlaneSpectrum.NumBins ; i++)
        CutPlaneSpectrum.Data[i] = 0.0;

    Temp1  = CutPlaneSpectrum.Max - CutPlaneSpectrum.Min;
    Delta  = 1.0 / (float)(CutPlaneSpectrum.NumBins - 1);

    CutPlaneSpectrum.BinMin  = 0.0;
    CutPlaneSpectrum.BinMax  = 0.0;

    for(k = 0 ; k < CutPlaneSpectrum.px ; k++)     {
 
        for(j = 0 ; j < CutPlaneSpectrum.pz ; j++)     {

            for( i = 0 ; i < CutPlaneSpectrum.py ; i++)    {
 
                Temp2 = ContourInfo[Which].data[
                    j * CutPlaneSpectrum.py * CutPlaneSpectrum.px + 
                    i * CutPlaneSpectrum.px + k];

                Temp3 = (Temp2 - CutPlaneSpectrum.Min) / Temp1;

                id = (int)(Temp3/Delta);

                if(id < (CutPlaneSpectrum.NumBins - 1)) {
                    CutPlaneSpectrum.Data[id + k * CutPlaneSpectrum.NumBins] += 1.0;
                    Temp2 = CutPlaneSpectrum.Data[id + k * CutPlaneSpectrum.NumBins];
                    if(Temp2 > CutPlaneSpectrum.BinMax)
                        CutPlaneSpectrum.BinMax = Temp2;
                }

            }
        }
    }

    CutPlaneSpectrum.InUse = Which + 1;

    if(CutPlaneSpectrum.BinCut > CutPlaneSpectrum.BinMax) {
        gomp_PrintERROR("Bin cut value > Bin max value, will put cut value = max value");
        CutPlaneSpectrum.BinCut = CutPlaneSpectrum.BinMax;
    }

    (void)PlotCutPlaneSpectrumX();
#endif /* ENABLE_GRAPHICS */

    return(0);
}
#ifdef ENABLE_GRAPHICS
/****************************************************************************/
int PlotCutPlaneSpectrumX()
/****************************************************************************/
{
    static int   i;
    static int   j;
    static int   Jump1;
    static int   Jump2;
    static float Temp1;
    static float Temp2;
    static float Temp3;
    static float Temp4;
    static float Temp5;
    static float Xc;
    static float Yc;
    static int   mm;
    static char  Text[BUFF_LEN];

    /*    glDisable(GL_LIGHTING);*/
 
    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);

    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D((GLdouble)0.0, (GLdouble)CutPlaneSpectrum.NumBins - 1.0, 
               (GLdouble)0.0, (GLdouble)CutPlaneSpectrum.px      - 1.0);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(mm);

    glColor3f(  1.0 , 1.0 , 1.0);
    glNormal3f( 0.0 , 0.0 , 1.0);

/* clear background */
    glClearColor( 1.0 , 1.0 , 1.0 , 0.0);

    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
/* ....... done ............. */
    glColor3fv(gomp_BLACKv);

    glPushMatrix();
    glRasterPos3f(0.5, (GLfloat)CutPlaneSpectrum.px - 2. , 0.1);
    sprintf(Text,"<Range: %.4e ... %.4e>",
            CutPlaneSpectrum.MinX,CutPlaneSpectrum.MaxX);
    gomp_PrintString(Text , "*");
    glPopMatrix();

    glRasterPos3f((GLfloat)CutPlaneSpectrum.NumBins/3., 0.5 , 0.1);
    sprintf(Text,"<Range: %.4e ... %.4e>",
            CutPlaneSpectrum.Min,CutPlaneSpectrum.Max);
    gomp_PrintString(Text , "*");

    glEnable(GL_LIGHTING);
    glEnable( GL_TEXTURE_1D );

    Temp5 =  CutPlaneSpectrum.BinCut - CutPlaneSpectrum.BinMin;

    for(j = 0 ; j < CutPlaneSpectrum.px - 1   ; j++)    {

        Jump1 =  j     * CutPlaneSpectrum.NumBins;
        Jump2 =  Jump1 + CutPlaneSpectrum.NumBins;

        for( i = 0 ; i < CutPlaneSpectrum.NumBins - 1; i++) {

            Temp1 = 
                (CutPlaneSpectrum.Data[Jump1 + i    ] - CutPlaneSpectrum.BinMin) / Temp5;
            Temp2 = 
                (CutPlaneSpectrum.Data[Jump1 + i + 1] - CutPlaneSpectrum.BinMin) / Temp5;
            Temp3 = 
                (CutPlaneSpectrum.Data[Jump2 + i + 1] - CutPlaneSpectrum.BinMin) / Temp5;
            Temp4 = 
                (CutPlaneSpectrum.Data[Jump2 + i    ] - CutPlaneSpectrum.BinMin) / Temp5;

            glPushMatrix();
 
            glBegin( GL_QUADS );

            Xc =  (float)i;
            Yc =  (float)j;

            glTexCoord1f( Temp1 );
            glVertex2f( Xc , Yc );

            Xc =  (float)(i+1);
            Yc =  (float)j;

            glTexCoord1f( Temp2 );
            glVertex2f( Xc , Yc );

            Xc =  (float)(i+1);
            Yc =  (float)(j+1); 

            glTexCoord1f( Temp3 );
            glVertex2f( Xc , Yc );

            Xc =  (float)i;
            Yc =  (float)(j+1); 

            glTexCoord1f( Temp4 );
            glVertex2f( Xc , Yc );

            glEnd();

            glPopMatrix();

        }
    }

    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_1D);
    /*    glEnable(GL_LIGHTING);*/
    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(mm);

    glEnable(GL_LIGHTING);

#if defined(GLUT)
    glutSwapBuffers();
#else
    auxSwapBuffers();
#endif

    return(0);
}
#endif /* ENABLE_GRAPHICS */

/****************************************************************************/
int    gomp_GetCutPlaneType_3D_X()
/****************************************************************************/
{
    return(CutPlaneType_3D_X);
}
/****************************************************************************/
int    gomp_GetCutPlaneType_3D_Y()
/****************************************************************************/
{
    return(CutPlaneType_3D_Y);
}
/****************************************************************************/
int    gomp_GetCutPlaneType_3D_Z()
/****************************************************************************/
{
    return(CutPlaneType_3D_Z);
}

/****************************************************************************/
int    gomp_SetCutPlaneType_3D_X(int Value)
/****************************************************************************/
{
    CutPlaneType_3D_X = Value;

    return(0);
}
/****************************************************************************/
int    gomp_SetCutPlaneType_3D_Y(int Value)
/****************************************************************************/
{
    CutPlaneType_3D_Y = Value;

    return(0);
}
/****************************************************************************/
int    gomp_SetCutPlaneType_3D_Z(int Value)
/****************************************************************************/
{
    CutPlaneType_3D_Z = Value;

    return(0);
}

/****************************************************************************/
int     gomp_SetCutPlaneDamping(float Value)
/****************************************************************************/
{

    CutPlaneDampingConstant = Value;

    return(0);
}

/****************************************************************************/
float   gomp_GetCutPlaneDamping()
/****************************************************************************/
{

    return(CutPlaneDampingConstant);
}
/****************************************************************************/
int gomp_ResetPrepare1DTexture() 
/****************************************************************************/
{
    Check = 0;

    return(0);
}
/****************************************************************************/
int gomp_PrepareCutPlaneXYZ(int Alt , int Which , 
                          const char *Text1 , const char *Text2, const char *Text3,
                          const char *Text4 , const char *Text5, const char *Text6,
                          const char *Text7 , const char *Text8, const char *Text9,
                          const char *Text10, const char *Text11)
/****************************************************************************/
{
    float p1[3];
    float p2[3];
    float p3[3];
    float A = 0.0,B = 0.0,C = 0.0,D = 0.0;
    float ScaleMin;
    float ScaleMax;

    p1[0] = atof(Text1);
    p1[1] = atof(Text2);
    p1[2] = atof(Text3);

    p2[0] = atof(Text4);
    p2[1] = atof(Text5);
    p2[2] = atof(Text6);

    p3[0] = atof(Text7);
    p3[1] = atof(Text8);
    p3[2] = atof(Text9);

/* scale min */
    if(strlen(Text10) != 0) {
        ScaleMin = atof(Text10);
    } else {
        ScaleMin = ContourInfo[Which].min;
    }
/* scale max */
    if(strlen(Text11) != 0) {
        ScaleMax = atof(Text11);
    } else {
        ScaleMax = ContourInfo[Which].max;
    }

    if(gomp_CalcPlaneFrom3Points( p1 , p2 , p3 , 
                                &A ,&B  ,&C  ,&D))
        return(1);

    switch(Alt) {
    case 1:
        CutPlane.PlotStateXYZ1    = 1;
        CutPlane.PlotWhichXYZ1    = Which;
        CutPlane.PlotPlaneXYZ1[0] = A;
        CutPlane.PlotPlaneXYZ1[1] = B;
        CutPlane.PlotPlaneXYZ1[2] = C;
        CutPlane.PlotPlaneXYZ1[3] = D;
        CutPlane.PlotScaleMinXYZ1 = ScaleMin;
        CutPlane.PlotScaleMaxXYZ1 = ScaleMax;
        break;
    case 2:
        CutPlane.PlotStateXYZ2    = 1;
        CutPlane.PlotWhichXYZ2    = Which;
        CutPlane.PlotPlaneXYZ2[0] = A;
        CutPlane.PlotPlaneXYZ2[1] = B;
        CutPlane.PlotPlaneXYZ2[2] = C;
        CutPlane.PlotPlaneXYZ2[3] = D;
        CutPlane.PlotScaleMinXYZ2 = ScaleMin;
        CutPlane.PlotScaleMaxXYZ2 = ScaleMax;
        break;
    case 3:
        CutPlane.PlotStateXYZ3    = 1;
        CutPlane.PlotWhichXYZ3    = Which;
        CutPlane.PlotPlaneXYZ3[0] = A;
        CutPlane.PlotPlaneXYZ3[1] = B;
        CutPlane.PlotPlaneXYZ3[2] = C;
        CutPlane.PlotPlaneXYZ3[3] = D;
        CutPlane.PlotScaleMinXYZ3 = ScaleMin;
        CutPlane.PlotScaleMaxXYZ3 = ScaleMax;
        break;
    }
    return(0);
}

/****************************************************************************/
int gomp_GetCutPlaneXYZ(int Alt , int *Which , 
                      float *A , float *B, float *C, float *D , 
                      float *ScaleMin , float *ScaleMax)
/****************************************************************************/
{

    switch(Alt) {
    case 1:
        *Which = CutPlane.PlotWhichXYZ1;
        *A     = CutPlane.PlotPlaneXYZ1[0];
        *B     = CutPlane.PlotPlaneXYZ1[1];
        *C     = CutPlane.PlotPlaneXYZ1[2];
        *D     = CutPlane.PlotPlaneXYZ1[3];
        *ScaleMin = CutPlane.PlotScaleMinXYZ1;
        *ScaleMax = CutPlane.PlotScaleMaxXYZ1;
        break;
    case 2:
        *Which = CutPlane.PlotWhichXYZ2;
        *A     = CutPlane.PlotPlaneXYZ2[0];
        *B     = CutPlane.PlotPlaneXYZ2[1];
        *C     = CutPlane.PlotPlaneXYZ2[2];
        *D     = CutPlane.PlotPlaneXYZ2[3];
        *ScaleMin = CutPlane.PlotScaleMinXYZ2;
        *ScaleMax = CutPlane.PlotScaleMaxXYZ2;
        break;
    case 3:
        *Which = CutPlane.PlotWhichXYZ3;
        *A     = CutPlane.PlotPlaneXYZ3[0];
        *B     = CutPlane.PlotPlaneXYZ3[1];
        *C     = CutPlane.PlotPlaneXYZ3[2];
        *D     = CutPlane.PlotPlaneXYZ3[3];
        *ScaleMin = CutPlane.PlotScaleMinXYZ3;
        *ScaleMax = CutPlane.PlotScaleMaxXYZ3;
        break;
    }

    return(0);
}
/****************************************************************************/
int    gomp_ReadCutPlaneXYZFromModelFile(FILE *Model_f)
/****************************************************************************/
{
    char  InputText[BUFF_LEN];
    int   State,Which;
    float XYZ[4];
    float Min,Max;

/* start the process ... */
/* plane #1*/
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText ,"%d %d %f %f %f %f %f %f", &State,
           &Which,
           &XYZ[0],
           &XYZ[1],
           &XYZ[2],
           &XYZ[3],
           &Min,
           &Max);

    CutPlane.PlotStateXYZ1        = State;
    CutPlane.PlotWhichXYZ1        = Which;
    CutPlane.PlotPlaneXYZ1[0]     = XYZ[0];
    CutPlane.PlotPlaneXYZ1[1]     = XYZ[1];
    CutPlane.PlotPlaneXYZ1[2]     = XYZ[2];
    CutPlane.PlotPlaneXYZ1[3]     = XYZ[3];
    CutPlane.PlotScaleMinXYZ1     = Min;
    CutPlane.PlotScaleMaxXYZ1     = Max;
/* plane #2*/
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText ,"%d %d %f %f %f %f %f %f", &State,
           &Which,
           &XYZ[0],
           &XYZ[1],
           &XYZ[2],
           &XYZ[3],
           &Min,
           &Max);

    CutPlane.PlotStateXYZ2        = State;
    CutPlane.PlotWhichXYZ2        = Which;
    CutPlane.PlotPlaneXYZ2[0]     = XYZ[0];
    CutPlane.PlotPlaneXYZ2[1]     = XYZ[1];
    CutPlane.PlotPlaneXYZ2[2]     = XYZ[2];
    CutPlane.PlotPlaneXYZ2[3]     = XYZ[3];
    CutPlane.PlotScaleMinXYZ2     = Min;
    CutPlane.PlotScaleMaxXYZ2     = Max;
/* plane #3*/
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText ,"%d %d %f %f %f %f %f %f", &State,
           &Which,
           &XYZ[0],
           &XYZ[1],
           &XYZ[2],
           &XYZ[3],
           &Min,
           &Max);

    CutPlane.PlotStateXYZ3        = State;
    CutPlane.PlotWhichXYZ3        = Which;
    CutPlane.PlotPlaneXYZ3[0]     = XYZ[0];
    CutPlane.PlotPlaneXYZ3[1]     = XYZ[1];
    CutPlane.PlotPlaneXYZ3[2]     = XYZ[2];
    CutPlane.PlotPlaneXYZ3[3]     = XYZ[3];
    CutPlane.PlotScaleMinXYZ3     = Min;
    CutPlane.PlotScaleMaxXYZ3     = Max;

    return(0);
}

