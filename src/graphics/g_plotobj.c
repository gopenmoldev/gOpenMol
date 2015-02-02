/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include "gommath.h"
#include <ctype.h>
#include <tcl.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

#include "axis.h"
#include "colouring.h"
#include "gomcast.h"
#include "gomfile.h"
#include "gomstring.h"
#include "gomtext.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "objseg.h"
#include "plot_cpk.h"
#include "plot_molec.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "rforce.h"
#include "tclutils.h"

#include "stdafx.h"

#define SURFACE_NORM_TRIGGER 1.0e-10

typedef struct {
    int (*DelSeg)(void);
} DelSegFunc;
static int ResetMolecStructureSpecificData(void *, int, int);
#define ADD_RESET_LISTENER(Type) \
    if ( ! Plot##Type.MolecStructDeleteListener ) { \
        static DelSegFunc func = {gomp_Del##Type##Seg}; \
        Plot##Type.MolecStructDeleteListener = \
            gomp_AddMolecStructDeleteListener( \
                ResetMolecStructureSpecificData,&func); \
    }
#define CANCEL_RESET_LISTENER(Type) \
    if ( Plot##Type.MolecStructDeleteListener ) { \
        gomp_CancelMolecStructDeleteListener( \
            Plot##Type.MolecStructDeleteListener); \
        Plot##Type.MolecStructDeleteListener = NULL; \
    }
    

/* line .............................*/
static struct LineSeg {
    gom_PlotterData Plotter;
    gom_MolecStructDeleteListener *MolecStructDeleteListener;
    int    Lines;          /* number of lines            */
    float *LineCoord;      /* pointer to the coordinates */
    float *Red;            /* pointer to contain the colour */
    float *Green;
    float *Blue;
} PlotLine = { { NULL , NULL } , NULL ,
               0 , NULL , NULL , NULL , NULL};
#define PlotLineSegIsChanging() \
    gomp_InvalidatePlotterDelayed(&PlotLine.Plotter)
/* line .............................*/

/* sphere .............................*/
static struct SphereSeg {
    gom_PlotterData Plotter;
    gom_MolecStructDeleteListener *MolecStructDeleteListener;
    int    Spheres;          /* number of spheres            */
    float *SphereCoord;      /* pointer to the coordinates   */
    float *Radius;           /* pointer to sphere radius     */
    float *XYZScale;
    float *Red;              /* pointer to contain the colour */
    float *Green;
    float *Blue;
    float *Alpha;
} PlotSphere = { { NULL , NULL } , NULL ,
                 0 , NULL , NULL , NULL , NULL , NULL , NULL};
#define PlotSphereSegIsChanging() \
    gomp_InvalidatePlotterDelayed(&PlotSphere.Plotter)
/* sphere .............................*/

/* plane .............................*/
static struct PlaneSeg {
    gom_PlotterData Plotter;
    gom_MolecStructDeleteListener *MolecStructDeleteListener;
    int    Planes;           /* number of planes              */
    int *PlaneAxis;          /* plane on axis (x,y or z)      */
    float *PlaneCoord;       /* pointer to the coordinates    */
    float *Red;              /* pointer to contain the colour */
    float *Green;
    float *Blue;
    float *Alpha;
} PlotPlane = { { NULL , NULL } , NULL ,
                0 , NULL , NULL , NULL , NULL , NULL};
#define PlotPlaneSegIsChanging() \
    gomp_InvalidatePlotterDelayed(&PlotPlane.Plotter)
/* plane .............................*/

#if 0
static struct {
    int NumCurves;
    int Step;
    int *NumObs;
    float *Xvalues;
    float *Yvalues;
    float *Zvalues;
} PlotCurves = { 0 , 100 , NULL , NULL , NULL ,NULL};
#endif

/* cyliner .............................*/
static struct CylinderSeg {
    gom_PlotterData Plotter;
    gom_MolecStructDeleteListener *MolecStructDeleteListener;
    int    Cylinders;          /* number of cylinders           */
    float *CylinderCoord;      /* pointer to the coordinates    */
    float *CylinderRad;        /* cylinder radius               */
    float *Red;            /* pointer to contain the colour */
    float *Green;
    float *Blue;
    float *Alpha;
} PlotCylinder = { { NULL , NULL } , NULL ,
                   0 , NULL , NULL , NULL , NULL , NULL , NULL};
#define PlotCylinderSegIsChanging() \
    gomp_InvalidatePlotterDelayed(&PlotCylinder.Plotter)
/* cyliner .............................*/

/* arrow .............................. */
static struct ArrowSeg {
    gom_PlotterData Plotter;
    gom_MolecStructDeleteListener *MolecStructDeleteListener;
    int    Arrows;          /* number of arrows              */
    float *ArrowCoord;      /* pointer to the coordinates    */
    float *ArrowRad;        /* arrow radius                  */
    float *Red;             /* pointer to contain the colour */
    float *Green;
    float *Blue;
    float *Alpha;
} PlotArrow = { { NULL , NULL } , NULL ,
                0 , NULL , NULL , NULL , NULL , NULL , NULL};
#define PlotArrowSegIsChanging() \
    gomp_InvalidatePlotterDelayed(&PlotArrow.Plotter)

static struct TriangleMesh {
    gom_PlotterData Plotter;
    gom_MolecStructDeleteListener *MolecStructDeleteListener;
    int   Triangles;
    float *x;
    float *y;
    float *z;
    float *u;
    float *v;
    float *w;
    float *c;
    float *Alpha;
} PlotTriangle = { { NULL , NULL } , NULL ,
                        0 , NULL , NULL , NULL , NULL, NULL, NULL, NULL};
#define PlotTriangeMeshIsChanging() \
    gomp_InvalidatePlotterDelayed(&PlotTriangle.Plotter)

static int PlotArrowSeg(void*,int,int);
static int PlotSphereSeg(void*,int,int);
static int PlotCylinderSeg(void*,int,int);
static int PlotPlaneSeg(void*,int,int);
static int PlotLineSeg(void*,int,int);
#if 0
static int SplitColourLine(const char *, float *, float *, float *);
#endif

static int PlotTriangleSeg(void*,int,int);

/* arrow .............................. */
static int PlotOneLineArrow(float,float,float,float,float,float);  
static struct LineArrowSeg {
    int    Arrows;          /* number of arrows              */
    int    Style;           /* display style                 */
    float *ArrowCoord;      /* pointer to the coordinates    */
    float *Red;             /* pointer to contain the colour */
    float *Green;
    float *Blue;
} PlotLineArrow = { 0 , 0 , NULL , NULL , NULL , NULL };

#ifdef ENABLE_GRAPHICS
static int PlotLineArrowSeg(int,int);
#endif
#if 0
static int GetLineArrowEntry(int , float *, float *, float *,
                             float *, float *, float *,
                             float *, float *, float *);
static int SetLineArrowEntryColor(int , float, float, float);
#endif
#ifdef ENABLE_GRAPHICS
static int GetGradientDisplayStyle(void);
#endif

/****************************************************************************/
int ResetMolecStructureSpecificData(void *Func, int Wstr, int Dstr)
/****************************************************************************/
{
    DelSegFunc *dsf = Func;
    if ( Wstr < 0 || /* full reset */
         ( Wstr == 0 && Dstr != 1 ) )
        /* Structure Wstr is about to be deleted or merged to
         * structure which won't be the first structure
         */
        dsf->DelSeg();
    return(0);
}
/****************************************************************************/
int gomp_PushLineStack(const char *InputP1,const char *InputP2,const char *InputP3,
                  const char *InputP4,const char *InputP5,const char *InputP6,
                  const char *InputP7,const char *InputP8)
/****************************************************************************/
{
    float  Xc1,Yc1,Zc1;
    float  Xc2,Yc2,Zc2;
    int    iappend = 0;
    int    Location;
    float  LineRed;
    float  LineGreen;
    float  LineBlue;
    const float *sumxyz;

    Xc1 = atof(InputP1);
    Yc1 = atof(InputP2);
    Zc1 = atof(InputP3);

    Xc2 = atof(InputP4);
    Yc2 = atof(InputP5);
    Zc2 = atof(InputP6);

    if(gomp_ColourName2RGB(InputP7 , &LineRed , &LineGreen , &LineBlue)) {
        gomp_PrintMessage("?ERROR - can't resolve the colour");
        return(1);
    }

    if(gomp_StringMatch(InputP8,"appe$nd")) iappend = 1;

    sumxyz = gomp_GetTranslateArray();

    PlotLineSegIsChanging();

    if(!iappend)
        gomp_DelLineSeg();

    if(!PlotLine.Lines) { /* no old lines */

        PlotLine.LineCoord = gomp_AllocateFloatVector(6);
        PlotLine.LineCoord[0] = Xc1 - sumxyz[0];
        PlotLine.LineCoord[1] = Yc1 - sumxyz[1];
        PlotLine.LineCoord[2] = Zc1 - sumxyz[2];
        PlotLine.LineCoord[3] = Xc2 - sumxyz[0];
        PlotLine.LineCoord[4] = Yc2 - sumxyz[1];
        PlotLine.LineCoord[5] = Zc2 - sumxyz[2];

        PlotLine.Red    = gomp_AllocateFloatVector(1);
        PlotLine.Green  = gomp_AllocateFloatVector(1);
        PlotLine.Blue   = gomp_AllocateFloatVector(1);

        PlotLine.Red[0]    = (float)LineRed;
        PlotLine.Green[0]  = (float)LineGreen;
        PlotLine.Blue[0]   = (float)LineBlue;
      
        PlotLine.Lines     =         1;
    } else {
        PlotLine.LineCoord = gomp_ReallocateFloatVector(PlotLine.LineCoord,
                                           (PlotLine.Lines + 1) * 6);
        Location = 6 * PlotLine.Lines;
        PlotLine.LineCoord[Location]     = Xc1 - sumxyz[0];
        PlotLine.LineCoord[Location + 1] = Yc1 - sumxyz[1];
        PlotLine.LineCoord[Location + 2] = Zc1 - sumxyz[2];
        PlotLine.LineCoord[Location + 3] = Xc2 - sumxyz[0];
        PlotLine.LineCoord[Location + 4] = Yc2 - sumxyz[1];
        PlotLine.LineCoord[Location + 5] = Zc2 - sumxyz[2];

        PlotLine.Red    =
            gomp_ReallocateFloatVector(PlotLine.Red   , (PlotLine.Lines + 1));
        PlotLine.Green  =
            gomp_ReallocateFloatVector(PlotLine.Green , (PlotLine.Lines + 1));
        PlotLine.Blue   =
            gomp_ReallocateFloatVector(PlotLine.Blue  , (PlotLine.Lines + 1));

        PlotLine.Red[PlotLine.Lines]    = (float)LineRed;
        PlotLine.Green[PlotLine.Lines]  = (float)LineGreen;
        PlotLine.Blue[PlotLine.Lines]   = (float)LineBlue;

        PlotLine.Lines    +=         1;
    }

    gomp_SetPlotterRegistrationState(
        1, &PlotLine.Plotter, PlotLineSeg, NULL,
        PLOTTER_NAME_PLOT_LINE, PLOTTER_ORDER_PLOT_LINE);

    ADD_RESET_LISTENER(Line);

    return(0);
}
/****************************************************************************/
int PlotLineSeg(void *userData,int Wstr,int drawFlags)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int i;
    static float LineCol[3];

    if ( ! ( drawFlags & gom_PlotSimpleElements ) || Wstr != 0 )
        return(-1);
    
    glLineWidth((GLfloat)(gomp_GetMoleculeLineWidth()+0.0));

    glDisable(GL_LIGHTING);

    for(i = 0 ; i < PlotLine.Lines ; i++) {
        LineCol[0] = PlotLine.Red[i];
        LineCol[1] = PlotLine.Green[i];
        LineCol[2] = PlotLine.Blue[i];

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&LineCol[0] , &LineCol[1] , &LineCol[2]);

        glColor3f(LineCol[0] , LineCol[1] , LineCol[2]);

        glBegin(GL_LINES);
        glVertex3fv(&PlotLine.LineCoord[6*i]);       
        glVertex3fv(&PlotLine.LineCoord[6*i + 3]);
        glEnd();
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/****************************************************************************/
int gomp_DelLineSeg()
/****************************************************************************/
{
    PlotLineSegIsChanging();
    
    if(PlotLine.Lines) {
        free(PlotLine.LineCoord);
        free(PlotLine.Red);
        free(PlotLine.Green);
        free(PlotLine.Blue);
    }
    PlotLine.Lines = 0;

    gomp_UnregisterPlotter(PlotLine.Plotter.plotter);
    PlotLine.Plotter.plotter = NULL;

    CANCEL_RESET_LISTENER(Line);

    return(0);
}

/****************************************************************************/
int gomp_PushSphereStack(const char *InputP1,const char *InputP2,const char *InputP3,
                    const char *InputP4,const char *InputP5,const char *InputP6,
                    const char *InputP7,const char *InputP8,const char *InputP9,
                    const char *InputP10)
/****************************************************************************/
{
    float  Xc1,Yc1,Zc1;
    float  radius;
    int    iappend = 0;
    int    Location;
    float  SphereRed;
    float  SphereGreen;
    float  SphereBlue;
    float  XScale;
    float  YScale;
    float  ZScale;
    float  Alpha;
    const float *sumxyz;

    Xc1 = atof(InputP1);
    Yc1 = atof(InputP2);
    Zc1 = atof(InputP3);

    radius = atof(InputP4);

    if(gomp_ColourName2RGB(InputP5, &SphereRed ,&SphereGreen ,&SphereBlue)) {
        gomp_PrintMessage("?ERROR - can't resolve the colour");
        return(1);
    }
    if(gomp_StringMatch(InputP6,"appe$nd") ||
       gomp_StringMatch(InputP7,"appe$nd") ||
       gomp_StringMatch(InputP8,"appe$nd") ||
       gomp_StringMatch(InputP9,"appe$nd") ||
       gomp_StringMatch(InputP10,"appe$nd"))
        iappend = 1;

    if(gomp_IsStringAFloat(InputP6)) {
        XScale  = atof(InputP6);
        YScale  = atof(InputP7);
        ZScale  = atof(InputP8);
    }
    else {
        XScale  = 1.0;
        YScale  = 1.0;
        ZScale  = 1.0;
    }
       
    if(XScale < 0.001) XScale = 1.0;
    if(YScale < 0.001) YScale = 1.0;
    if(ZScale < 0.001) ZScale = 1.0;

    if(gomp_IsStringAFloat(InputP9)) {
        Alpha = atof(InputP9);
    } else {
        Alpha = 1.0;
    }

    sumxyz = gomp_GetTranslateArray();

    PlotSphereSegIsChanging();

    if(!iappend)
        gomp_DelSphereSeg();

    if(!PlotSphere.Spheres) { /* no old lines */

        PlotSphere.SphereCoord = gomp_AllocateFloatVector(3);
        PlotSphere.XYZScale    = gomp_AllocateFloatVector(3);
        PlotSphere.Radius      = gomp_AllocateFloatVector(1);
        PlotSphere.Alpha       = gomp_AllocateFloatVector(1);
        PlotSphere.SphereCoord[0] = Xc1 - sumxyz[0];
        PlotSphere.SphereCoord[1] = Yc1 - sumxyz[1];
        PlotSphere.SphereCoord[2] = Zc1 - sumxyz[2];

        PlotSphere.XYZScale[0]    = XScale * radius;
        PlotSphere.XYZScale[1]    = YScale * radius;
        PlotSphere.XYZScale[2]    = ZScale * radius;

        PlotSphere.Red    = gomp_AllocateFloatVector(1);
        PlotSphere.Green  = gomp_AllocateFloatVector(1);
        PlotSphere.Blue   = gomp_AllocateFloatVector(1);

        PlotSphere.Red[0]    = (float)SphereRed;
        PlotSphere.Green[0]  = (float)SphereGreen;
        PlotSphere.Blue[0]   = (float)SphereBlue;
      
        PlotSphere.Radius[0] =    radius; 
        PlotSphere.Alpha[0]  =     Alpha;
        PlotSphere.Spheres   =         1;
    } else {
        PlotSphere.SphereCoord = gomp_ReallocateFloatVector(PlotSphere.SphereCoord,
                                               (PlotSphere.Spheres + 1) * 3);
        PlotSphere.XYZScale    = gomp_ReallocateFloatVector(PlotSphere.XYZScale,
                                               (PlotSphere.Spheres + 1) * 3);
        PlotSphere.Radius      = gomp_ReallocateFloatVector(PlotSphere.Radius,
                                               (PlotSphere.Spheres + 1));
        PlotSphere.Alpha       = gomp_ReallocateFloatVector(PlotSphere.Alpha,
                                               (PlotSphere.Spheres + 1));
        Location = 3 * PlotSphere.Spheres;
        PlotSphere.SphereCoord[Location]     = Xc1 - sumxyz[0];
        PlotSphere.SphereCoord[Location + 1] = Yc1 - sumxyz[1];
        PlotSphere.SphereCoord[Location + 2] = Zc1 - sumxyz[2];

        PlotSphere.XYZScale[Location    ]   = XScale * radius;
        PlotSphere.XYZScale[Location + 1]   = YScale * radius;
        PlotSphere.XYZScale[Location + 2]   = ZScale * radius;

        PlotSphere.Radius[PlotSphere.Spheres]  = radius; 
        PlotSphere.Alpha[PlotSphere.Spheres]   = Alpha; 

        PlotSphere.Red    = gomp_ReallocateFloatVector(PlotSphere.Red   , 
                                          (PlotSphere.Spheres + 1));
        PlotSphere.Green  = gomp_ReallocateFloatVector(PlotSphere.Green , 
                                          (PlotSphere.Spheres + 1));
        PlotSphere.Blue   = gomp_ReallocateFloatVector(PlotSphere.Blue  , 
                                          (PlotSphere.Spheres + 1));

        PlotSphere.Red[PlotSphere.Spheres]    = (float)SphereRed;
        PlotSphere.Green[PlotSphere.Spheres]  = (float)SphereGreen;
        PlotSphere.Blue[PlotSphere.Spheres]   = (float)SphereBlue;

        PlotSphere.Spheres  +=         1;
    }

    gomp_SetPlotterRegistrationState(
        1, &PlotSphere.Plotter, PlotSphereSeg, NULL,
        PLOTTER_NAME_PLOT_SPHERE, PLOTTER_ORDER_PLOT_SPHERE);

    ADD_RESET_LISTENER(Sphere);

    return(0);
}
/****************************************************************************/
int PlotSphereSeg(void *userData,int Wstr,int drawFlags)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int i;
    static float param[4];
    static int SphereQ;
    static int SphereQ2;

    if ( ! ( drawFlags & gom_PlotComplexElements ) || Wstr != 0 )
        return(-1);

    glEnable(GL_LIGHTING);

    for(i = 0 ; i < PlotSphere.Spheres ; i++) {

        param[0] = PlotSphere.Red[i];
        param[1] = PlotSphere.Green[i];
        param[2] = PlotSphere.Blue[i];
        param[3] = PlotSphere.Alpha[i];
 
        SphereQ  = gomp_GetSphereQuality();
        SphereQ2 = 2 * SphereQ;

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&param[0] , &param[1] , &param[2]);

        if(param[3] < 1.0) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA);
        } else {
            glBlendFunc(GL_ONE , GL_ZERO);
            glDisable(GL_BLEND);
        }

        glColor4f(param[0] , param[1] , param[2] , param[3]);

        glPushMatrix();
        param[0]   = PlotSphere.SphereCoord[3 * i    ];
        param[1]   = PlotSphere.SphereCoord[3 * i + 1];
        param[2]   = PlotSphere.SphereCoord[3 * i + 2];
        glTranslatef(param[0],param[1],param[2]);
        glScalef(PlotSphere.XYZScale[3 * i    ] ,
                 PlotSphere.XYZScale[3 * i + 1] ,
                 PlotSphere.XYZScale[3 * i + 2]);
/*   gluSphere(SphereQuad,(double)PlotSphere.Radius[i], SphereQ2 , SphereQ);*/
        gluSphere(gomp_SphereQuad,1.0, SphereQ2 , SphereQ);
        glPopMatrix();

        glBlendFunc(GL_ONE , GL_ZERO);
        glDisable(GL_BLEND);
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/****************************************************************************/
int gomp_DelSphereSeg()
/****************************************************************************/
{
    PlotSphereSegIsChanging();
    
    if(PlotSphere.Spheres) { 
        free(PlotSphere.Red);
        free(PlotSphere.Green);
        free(PlotSphere.Blue);
        free(PlotSphere.SphereCoord);
        free(PlotSphere.XYZScale);
        free(PlotSphere.Radius);
        free(PlotSphere.Alpha);
    }
    PlotSphere.Spheres = 0;

    gomp_UnregisterPlotter(PlotSphere.Plotter.plotter);
    PlotSphere.Plotter.plotter = NULL;

    CANCEL_RESET_LISTENER(Sphere);

    return(0);
}

/****************************************************************************/
int gomp_PushPlaneStack(
    const char *InputP1,
    const char *InputP2,
    const char *InputP3,
    const char *InputP4,
    const char *InputP5,
    const char *InputP6,
    const char *InputP7)
/****************************************************************************/
{

    float Xc1,Yc1,Zc1;
    float Xc2,Yc2,Zc2;
    float Xc3,Yc3,Zc3;
    float Xc4,Yc4,Zc4;
    float Da;
    float Db;
    float Dc;
    int   iappend = 0;
    int   Location;
    int   Axis;
    float PlaneRed;
    float PlaneGreen;
    float PlaneBlue;
    const float *sumxyz;
    float Alpha;

    Axis = -1;

    if(gomp_StringMatch(InputP1,"x") ||
       gomp_StringMatch(InputP1,"X")) {
        Axis = 0;
    }
    else if(gomp_StringMatch(InputP1,"y") ||
            gomp_StringMatch(InputP1,"Y")) {
        Axis = 1;
    }
    else if(gomp_StringMatch(InputP1,"z") ||
            gomp_StringMatch(InputP1,"Z")) {
        Axis = 2;
    }

    Da = atof(InputP2);
    Db = atof(InputP3);
    Dc = atof(InputP4);

    switch(Axis) {
    case 0: /* x */
        Xc1 = Da;
        Yc1 = Db;
        Zc1 = Dc;

        Xc2 = Da;
        Yc2 = Db;
        Zc2 = -Dc;

        Xc3 = Da;
        Yc3 = -Db;
        Zc3 = -Dc;

        Xc4 = Da;
        Yc4 = -Db;
        Zc4 = Dc;
        break;
    case 1: /* y */
        Xc1 = Db;
        Yc1 = Da;
        Zc1 = Dc;

        Xc2 = Db;
        Yc2 = Da;
        Zc2 = -Dc;

        Xc3 = -Db;
        Yc3 = Da;
        Zc3 = -Dc;

        Xc4 = -Db;
        Yc4 = Da;
        Zc4 = Dc;
        break;
    case 2: /* z */
        Xc1 = Db;
        Yc1 = Dc;
        Zc1 = Da;

        Xc2 = Db;
        Yc2 = -Dc;
        Zc2 = Da;

        Xc3 = -Db;
        Yc3 = -Dc;
        Zc3 = Da;

        Xc4 = -Db;
        Yc4 = Dc;
        Zc4 = Da;
        break;
    default:
        return(1);
    }

    if(Axis < 0) {
        gomp_PrintERROR("undefined axis, must be one of x,y or z");
        return(1);
    }
    if(gomp_ColourName2RGB(InputP5, &PlaneRed ,&PlaneGreen ,&PlaneBlue)) {
        gomp_PrintMessage("?ERROR - can't resolve the colour");
        return(1);
    }

    if(gomp_IsStringAFloat(InputP6)) {
        Alpha = atof(InputP6);
    } else {
        Alpha = 1.0;
    }

    if(gomp_StringMatch(InputP6,"appe$nd") ||
       gomp_StringMatch(InputP6,"appe$nd")) iappend = 1;

    sumxyz = gomp_GetTranslateArray();

    PlotPlaneSegIsChanging();
    
    if(!iappend)
        gomp_DelPlaneSeg();

    if(!PlotPlane.Planes) { /* no old Planes */

        PlotPlane.PlaneCoord     = gomp_AllocateFloatVector(12);
        PlotPlane.PlaneCoord[0]  = Xc1 - sumxyz[0];
        PlotPlane.PlaneCoord[1]  = Yc1 - sumxyz[1];
        PlotPlane.PlaneCoord[2]  = Zc1 - sumxyz[2];
        PlotPlane.PlaneCoord[3]  = Xc2 - sumxyz[0];
        PlotPlane.PlaneCoord[4]  = Yc2 - sumxyz[1];
        PlotPlane.PlaneCoord[5]  = Zc2 - sumxyz[2];
        PlotPlane.PlaneCoord[6]  = Xc3 - sumxyz[0];
        PlotPlane.PlaneCoord[7]  = Yc3 - sumxyz[1];
        PlotPlane.PlaneCoord[8]  = Zc3 - sumxyz[2];
        PlotPlane.PlaneCoord[9]  = Xc4 - sumxyz[0];
        PlotPlane.PlaneCoord[10] = Yc4 - sumxyz[1];
        PlotPlane.PlaneCoord[11] = Zc4 - sumxyz[2];

        PlotPlane.PlaneAxis    = gomp_AllocateIntVector(1);
        PlotPlane.PlaneAxis[0] = Axis;

        PlotPlane.Red        = gomp_AllocateFloatVector(1);
        PlotPlane.Green      = gomp_AllocateFloatVector(1);
        PlotPlane.Blue       = gomp_AllocateFloatVector(1);

        PlotPlane.Alpha      = gomp_AllocateFloatVector(1);

        PlotPlane.Red[0]     = PlaneRed;
        PlotPlane.Green[0]   = PlaneGreen;
        PlotPlane.Blue[0]    = PlaneBlue;

        PlotPlane.Alpha[0]   = Alpha;
        
        PlotPlane.Planes     =         1;
    } else {
        PlotPlane.PlaneCoord = gomp_ReallocateFloatVector(PlotPlane.PlaneCoord,
                                             (PlotPlane.Planes + 1) * 12);
        Location = 12 * PlotPlane.Planes;
        PlotPlane.PlaneCoord[Location]      = Xc1 - sumxyz[0];
        PlotPlane.PlaneCoord[Location + 1]  = Yc1 - sumxyz[1];
        PlotPlane.PlaneCoord[Location + 2]  = Zc1 - sumxyz[2];
        PlotPlane.PlaneCoord[Location + 3]  = Xc2 - sumxyz[0];
        PlotPlane.PlaneCoord[Location + 4]  = Yc2 - sumxyz[1];
        PlotPlane.PlaneCoord[Location + 5]  = Zc2 - sumxyz[2];
        PlotPlane.PlaneCoord[Location + 6]  = Xc3 - sumxyz[0];
        PlotPlane.PlaneCoord[Location + 7]  = Yc3 - sumxyz[1];
        PlotPlane.PlaneCoord[Location + 8]  = Zc3 - sumxyz[2];
        PlotPlane.PlaneCoord[Location + 9]  = Xc4 - sumxyz[0];
        PlotPlane.PlaneCoord[Location + 10] = Yc4 - sumxyz[1];
        PlotPlane.PlaneCoord[Location + 11] = Zc4 - sumxyz[2];

        PlotPlane.PlaneAxis                   = 
            gomp_ReallocateIntVector(PlotPlane.PlaneAxis   , (PlotPlane.Planes + 1));
        PlotPlane.PlaneAxis[PlotPlane.Planes] = Axis;

        PlotPlane.Red    =
            gomp_ReallocateFloatVector(PlotPlane.Red   , (PlotPlane.Planes + 1));
        PlotPlane.Green  =
            gomp_ReallocateFloatVector(PlotPlane.Green , (PlotPlane.Planes + 1));
        PlotPlane.Blue   =
            gomp_ReallocateFloatVector(PlotPlane.Blue  , (PlotPlane.Planes + 1));
        PlotPlane.Alpha  =
            gomp_ReallocateFloatVector(PlotPlane.Alpha , (PlotPlane.Planes + 1));

        PlotPlane.Red[PlotPlane.Planes]    = PlaneRed;
        PlotPlane.Green[PlotPlane.Planes]  = PlaneGreen;
        PlotPlane.Blue[PlotPlane.Planes]   = PlaneBlue;

        PlotPlane.Alpha[PlotPlane.Planes]  = Alpha;

        PlotPlane.Planes    +=         1;
    }

    gomp_SetPlotterRegistrationState(
        1, &PlotPlane.Plotter, PlotPlaneSeg, NULL,
        PLOTTER_NAME_PLOT_PLANE, PLOTTER_ORDER_PLOT_PLANE);

    ADD_RESET_LISTENER(Plane);

    return(0);
}
/****************************************************************************/
int PlotPlaneSeg(void *userData,int Wstr,int drawFlags)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    int i;
    float  CRi;
    float  CGi;
    float  CBi;

    if ( ! ( drawFlags & gom_PlotComplexElements ) || Wstr != 0 )
        return(-1);

    glEnable(GL_LIGHTING);

    for(i = 0 ; i < PlotPlane.Planes ; i++) {

        if(PlotPlane.Alpha[i] < 1.0) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA);
        } else {
            glBlendFunc(GL_ONE , GL_ZERO);
            glDisable(GL_BLEND);
        }

        CRi = PlotPlane.Red[i];
        CGi = PlotPlane.Green[i];
        CBi = PlotPlane.Blue[i];

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&CRi , &CGi , &CBi);

        glColor4f(CRi ,
                  CGi ,
                  CBi ,
                  PlotPlane.Alpha[i]);

        glBegin(GL_QUADS);

        switch(PlotPlane.PlaneAxis[i]) {
        case 0:
            glNormal3f( -1.0 ,  0.0 ,  0.0);
            break;
        case 1:
            glNormal3f( 0.0 ,  1.0 ,  0.0);
            break;
        case 2:
            glNormal3f( 0.0 ,  0.0 , -1.0);
            break;
        }

        glVertex3fv(&PlotPlane.PlaneCoord[12*i]);       
        glVertex3fv(&PlotPlane.PlaneCoord[12*i + 3]);
        glVertex3fv(&PlotPlane.PlaneCoord[12*i + 6]);
        glVertex3fv(&PlotPlane.PlaneCoord[12*i + 9]);
        glEnd();
        
        glBlendFunc(GL_ONE , GL_ZERO);
        glDisable(GL_BLEND);

    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/****************************************************************************/
int gomp_DelPlaneSeg()
/****************************************************************************/
{
    PlotPlaneSegIsChanging();
    
    if(PlotPlane.Planes) {
        free(PlotPlane.PlaneAxis);
        free(PlotPlane.PlaneCoord);
        free(PlotPlane.Alpha);
        free(PlotPlane.Red);
        free(PlotPlane.Green);
        free(PlotPlane.Blue);
    }
    PlotPlane.Planes = 0;

    gomp_UnregisterPlotter(PlotPlane.Plotter.plotter);
    PlotPlane.Plotter.plotter = NULL;

    CANCEL_RESET_LISTENER(Plane);

    return(0);
}

/****************************************************************************/
int gomp_PushCylinderStack(
    const char *InputP1,const char *InputP2,const char *InputP3,
    const char *InputP4,const char *InputP5,const char *InputP6,
    const char *InputP7,const char *InputP8,const char *InputP9,
    const char *InputP10)
/****************************************************************************/
{

    float  Xc1,Yc1,Zc1;
    float  Xc2,Yc2,Zc2;
    float  Radius1;
    float  Radius2;
    int    iappend = 0;
    int    Location;
    float  CylinderRed1;
    float  CylinderGreen1;
    float  CylinderBlue1;
    float  CylinderRed2;
    float  CylinderGreen2;
    float  CylinderBlue2;
    const float *sumxyz;
    int    argc, code;
    const char **argv;
    float  Alpha;

    Xc1 = atof(InputP1);
    Yc1 = atof(InputP2);
    Zc1 = atof(InputP3);

    Xc2 = atof(InputP4);
    Yc2 = atof(InputP5);
    Zc2 = atof(InputP6);

    code = Tcl_SplitList(gomp_GetTclInterp(), InputP7, &argc, &argv);

    if(argc == 1) {  
        Radius1 = atof(InputP7);
        Radius2 = Radius1;
        Tcl_Free((char *)CONST_CAST(char **, argv));
    } else if(argc == 2) {
        Radius1 = atof(argv[0]);
        Radius2 = atof(argv[1]);
        Tcl_Free((char *)CONST_CAST(char **, argv));
    } else {
        Tcl_Free((char *)CONST_CAST(char **, argv));
        gomp_PrintERROR("wrong number of entries in the radius field");
        return(1);
    }

    sumxyz = gomp_GetTranslateArray();

    code = Tcl_SplitList(gomp_GetTclInterp(), InputP8, &argc, &argv);

    if(argc == 1) {  
        if(gomp_ColourName2RGB(InputP8 ,  &CylinderRed1 , 
                          &CylinderGreen1 , 
                          &CylinderBlue1)) {
            gomp_PrintMessage("?ERROR - can't resolve the colour");
            Tcl_Free((char *)CONST_CAST(char **, argv));
            return(1);
        }
        CylinderRed2   = CylinderRed1;
        CylinderGreen2 = CylinderGreen1;
        CylinderBlue2  = CylinderBlue1;
        Tcl_Free((char *)CONST_CAST(char **, argv));
    } else if(argc == 2) {
/* start point */
        if(gomp_ColourName2RGB(argv[0] ,  &CylinderRed1 , 
                          &CylinderGreen1 , 
                          &CylinderBlue1)) {
            gomp_PrintMessage("?ERROR - can't resolve the colour #1");
            Tcl_Free((char *)CONST_CAST(char **, argv));
            return(1);
        }
/* end point */
        if(gomp_ColourName2RGB(argv[1] ,  &CylinderRed2 , 
                          &CylinderGreen2 , 
                          &CylinderBlue2)) {
            gomp_PrintMessage("?ERROR - can't resolve the colour #2");
            Tcl_Free((char *)CONST_CAST(char **, argv));
            return(1);
        }
        Tcl_Free((char *)CONST_CAST(char **, argv));
    } else {
        Tcl_Free((char *)CONST_CAST(char **, argv));
        gomp_PrintERROR("wrong number of entries in the colour field");
        return(1);
    }

    if(gomp_IsStringAFloat(InputP9)) {
        Alpha = atof(InputP9);
    } else {
        Alpha = 1.0;
    }

    if(gomp_StringMatch(InputP9,"appe$nd") ||
       gomp_StringMatch(InputP10,"appe$nd") ) iappend = 1;

    PlotCylinderSegIsChanging();
    
    if(!iappend)
        gomp_DelCylinderSeg();

    if(!PlotCylinder.Cylinders) { /* no old Cylinders */

        PlotCylinder.CylinderCoord = gomp_AllocateFloatVector(6);
        PlotCylinder.CylinderCoord[0]   = Xc1 - sumxyz[0];
        PlotCylinder.CylinderCoord[1]  = Yc1 - sumxyz[1];
        PlotCylinder.CylinderCoord[2] = Zc1 - sumxyz[2];
        PlotCylinder.CylinderCoord[3]   = Xc2 - sumxyz[0];
        PlotCylinder.CylinderCoord[4]  = Yc2 - sumxyz[1];
        PlotCylinder.CylinderCoord[5] = Zc2 - sumxyz[2];
        PlotCylinder.CylinderRad        = gomp_AllocateFloatVector(2);
        PlotCylinder.CylinderRad[0]     = Radius1;
        PlotCylinder.CylinderRad[1]     = Radius2;

        PlotCylinder.Alpha     =  gomp_AllocateFloatVector(1);
        PlotCylinder.Alpha[0]  =  Alpha;

        PlotCylinder.Red    = gomp_AllocateFloatVector(2);
        PlotCylinder.Green  = gomp_AllocateFloatVector(2);
        PlotCylinder.Blue   = gomp_AllocateFloatVector(2);

        PlotCylinder.Red[0]    = (float)CylinderRed1;
        PlotCylinder.Green[0]  = (float)CylinderGreen1;
        PlotCylinder.Blue[0]   = (float)CylinderBlue1;
        PlotCylinder.Red[1]    = (float)CylinderRed2;
        PlotCylinder.Green[1]  = (float)CylinderGreen2;
        PlotCylinder.Blue[1]   = (float)CylinderBlue2;
      
        PlotCylinder.Cylinders     =         1;
    } else {
        PlotCylinder.CylinderCoord =
            gomp_ReallocateFloatVector(PlotCylinder.CylinderCoord,
                             (PlotCylinder.Cylinders + 1) * 6);
        Location = 6 * PlotCylinder.Cylinders;
        PlotCylinder.CylinderCoord[Location]     = Xc1 - sumxyz[0];
        PlotCylinder.CylinderCoord[Location + 1] = Yc1 - sumxyz[1];
        PlotCylinder.CylinderCoord[Location + 2] = Zc1 - sumxyz[2];
        PlotCylinder.CylinderCoord[Location + 3] = Xc2 - sumxyz[0];
        PlotCylinder.CylinderCoord[Location + 4] = Yc2 - sumxyz[1];
        PlotCylinder.CylinderCoord[Location + 5] = Zc2 - sumxyz[2];
        PlotCylinder.CylinderRad = gomp_ReallocateFloatVector(PlotCylinder.CylinderRad,
                                                 2 * (PlotCylinder.Cylinders + 1));
        PlotCylinder.CylinderRad[2 * PlotCylinder.Cylinders]    = Radius1;
        PlotCylinder.CylinderRad[2 * PlotCylinder.Cylinders + 1]= Radius2;

        PlotCylinder.Alpha = gomp_ReallocateFloatVector(PlotCylinder.Alpha,
                                           (PlotCylinder.Cylinders + 1));
        PlotCylinder.Alpha[PlotCylinder.Cylinders]= Alpha;

        PlotCylinder.Red    =
            gomp_ReallocateFloatVector(PlotCylinder.Red   , 
                             2 * (PlotCylinder.Cylinders + 1));
        PlotCylinder.Green  =
            gomp_ReallocateFloatVector(PlotCylinder.Green , 
                             2 * (PlotCylinder.Cylinders + 1));
        PlotCylinder.Blue   =
            gomp_ReallocateFloatVector(PlotCylinder.Blue  , 
                             2 * (PlotCylinder.Cylinders + 1));

        PlotCylinder.Red[  2 * PlotCylinder.Cylinders]     =
            (float)CylinderRed1;
        PlotCylinder.Green[2 * PlotCylinder.Cylinders]     =
            (float)CylinderGreen1;
        PlotCylinder.Blue[ 2 * PlotCylinder.Cylinders]     =
            (float)CylinderBlue1;
        PlotCylinder.Red[  2 * PlotCylinder.Cylinders + 1] =
            (float)CylinderRed2;
        PlotCylinder.Green[2 * PlotCylinder.Cylinders + 1] =
            (float)CylinderGreen2;
        PlotCylinder.Blue[ 2 * PlotCylinder.Cylinders + 1] =
            (float)CylinderBlue2;

        PlotCylinder.Cylinders    +=         1;
    }

    gomp_SetPlotterRegistrationState(
        1, &PlotCylinder.Plotter, PlotCylinderSeg, NULL,
        PLOTTER_NAME_PLOT_CYLINDER, PLOTTER_ORDER_PLOT_CYLINDER);

    ADD_RESET_LISTENER(Cylinder);
    
    return(0);
}
/****************************************************************************/
int gomp_DelCylinderSeg()
/****************************************************************************/
{
    PlotCylinderSegIsChanging();
    
    if(PlotCylinder.Cylinders) {
        free(PlotCylinder.Red);
        free(PlotCylinder.Green);
        free(PlotCylinder.Blue);
        free(PlotCylinder.CylinderCoord);
        free(PlotCylinder.Alpha);
        free(PlotCylinder.CylinderRad);
    }
    PlotCylinder.Cylinders = 0;

    gomp_UnregisterPlotter(PlotCylinder.Plotter.plotter);
    PlotCylinder.Plotter.plotter = NULL;

    CANCEL_RESET_LISTENER(Cylinder);

    return(0);
}

/***********************************************************************/
int PlotCylinderSeg(void *userData,int Wstr,int drawFlags)
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int i;
    static int i2;
    static int i6;
    static float radd1;
    static float radd2;
    static float xi,yi,zi;
    static float xk,yk,zk;
    static int from,to;
    static float  CRi;
    static float  CGi;
    static float  CBi;
    static float  CRj;
    static float  CGj;
    static float  CBj;

    if ( ! ( drawFlags & gom_PlotComplexElements ) || Wstr != 0 )
        return(-1);

    glEnable(GL_LIGHTING);

    from = 0;
    to   = PlotCylinder.Cylinders;


/* throw out the cylinders */

    for(i = from ; i < to ; i++ ) {

        i2 = i + i;
        i6 = i2 + i2 + i2;
    
        radd1 = PlotCylinder.CylinderRad[i2];      /* radius of the cylinder */
        radd2 = PlotCylinder.CylinderRad[i2 + 1];  /* radius of the cylinder */

        xi = PlotCylinder.CylinderCoord[i6];
        yi = PlotCylinder.CylinderCoord[i6 + 1];
        zi = PlotCylinder.CylinderCoord[i6 + 2];

        xk = PlotCylinder.CylinderCoord[i6 + 3];
        yk = PlotCylinder.CylinderCoord[i6 + 4];
        zk = PlotCylinder.CylinderCoord[i6 + 5];

        CRi = PlotCylinder.Red[i2];
        CGi = PlotCylinder.Green[i2];
        CBi = PlotCylinder.Blue[i2];
        CRj = PlotCylinder.Red[i2 + 1];
        CGj = PlotCylinder.Green[i2 + 1];
        CBj = PlotCylinder.Blue[i2 + 1];

        if(PlotCylinder.Alpha[i] < 1.0) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA);
        } else {
            glBlendFunc(GL_ONE , GL_ZERO);
            glDisable(GL_BLEND);
        }

        if(!gomp_GetDisplayColourType()) {
            (void)gomp_RGB2Grayscale(&CRi , &CGi , &CBi);
            (void)gomp_RGB2Grayscale(&CRj , &CGj , &CBj);
        }

        gomp_PlotCylinder(xi , yi , zi ,
                        xk , yk , zk ,
                        CRi     , CGi     , CBi     ,
                        CRj     , CGj     , CBj ,
                        PlotCylinder.Alpha[i]    ,
                        radd1 , radd2);

        glBlendFunc(GL_ONE , GL_ZERO);
        glDisable(GL_BLEND);

    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}

#if 0
/***********************************************************************/
int  SplitColourLine(const char *InputP , float *LineRed, float *LineGreen , 
                     float *LineBlue) 
/***********************************************************************/
{
    sscanf(InputP,"%f %f %f",LineRed,LineGreen,LineBlue);

    return(0);
}
#endif



/****************************************************************************/
int gomp_PushArrowStack(
    const char *InputP1,const char *InputP2,const char *InputP3,
    const char *InputP4,const char *InputP5,const char *InputP6,
    const char *InputP7,const char *InputP8,const char *InputP9,
    const char *InputP10)
/****************************************************************************/
{
    float  Xc1,Yc1,Zc1;
    float  Xc2,Yc2,Zc2;
    float  Radius;
    int    iappend = 0;
    int    Location;
    float  ArrowRed;
    float  ArrowGreen;
    float  ArrowBlue;
    float  Alpha;
    const float *sumxyz;

    Xc1 = atof(InputP1);
    Yc1 = atof(InputP2);
    Zc1 = atof(InputP3);

    Xc2 = atof(InputP4);
    Yc2 = atof(InputP5);
    Zc2 = atof(InputP6);

    Radius = atof(InputP7);

    sumxyz = gomp_GetTranslateArray();

    if(gomp_ColourName2RGB(InputP8 , &ArrowRed , &ArrowGreen , &ArrowBlue)) {
        gomp_PrintMessage("?ERROR - can't resolve the colour");
        return(1);
    }

    if(gomp_IsStringAFloat(InputP9)) {
        Alpha = atof(InputP9);
    } else {
        Alpha = 1.0;
    }

    if(gomp_StringMatch(InputP9,"appe$nd") ||
       gomp_StringMatch(InputP10,"appe$nd") ) iappend = 1;

    PlotArrowSegIsChanging();

    if(!iappend)
        gomp_DelArrowSeg();

    if(!PlotArrow.Arrows) { /* no old Arrows */

        PlotArrow.ArrowCoord    = gomp_AllocateFloatVector(6);
        PlotArrow.ArrowCoord[0] = Xc1 - sumxyz[0];
        PlotArrow.ArrowCoord[1] = Yc1 - sumxyz[1];
        PlotArrow.ArrowCoord[2] = Zc1 - sumxyz[2];
        PlotArrow.ArrowCoord[3] = Xc2 - sumxyz[0];
        PlotArrow.ArrowCoord[4] = Yc2 - sumxyz[1];
        PlotArrow.ArrowCoord[5] = Zc2 - sumxyz[2];
        PlotArrow.ArrowRad      = gomp_AllocateFloatVector(1);
        PlotArrow.ArrowRad[0]   = Radius;

        PlotArrow.Alpha         =  gomp_AllocateFloatVector(1);
        PlotArrow.Alpha[0]      =  Alpha;

        PlotArrow.Red           = gomp_AllocateFloatVector(1);
        PlotArrow.Green         = gomp_AllocateFloatVector(1);
        PlotArrow.Blue          = gomp_AllocateFloatVector(1);

        PlotArrow.Red[0]        = (float)ArrowRed;
        PlotArrow.Green[0]      = (float)ArrowGreen;
        PlotArrow.Blue[0]       = (float)ArrowBlue;
      
        PlotArrow.Arrows        = 1;
    } else {
        PlotArrow.ArrowCoord = gomp_ReallocateFloatVector(PlotArrow.ArrowCoord,
                                             (PlotArrow.Arrows + 1) * 6);
        Location = 6 * PlotArrow.Arrows;
        PlotArrow.ArrowCoord[Location]     = Xc1 - sumxyz[0];
        PlotArrow.ArrowCoord[Location + 1] = Yc1 - sumxyz[1];
        PlotArrow.ArrowCoord[Location + 2] = Zc1 - sumxyz[2];
        PlotArrow.ArrowCoord[Location + 3] = Xc2 - sumxyz[0];
        PlotArrow.ArrowCoord[Location + 4] = Yc2 - sumxyz[1];
        PlotArrow.ArrowCoord[Location + 5] = Zc2 - sumxyz[2];
        PlotArrow.ArrowRad = gomp_ReallocateFloatVector(PlotArrow.ArrowRad,
                                           (PlotArrow.Arrows + 1));
        PlotArrow.ArrowRad[PlotArrow.Arrows]= Radius;

        PlotArrow.Alpha = gomp_ReallocateFloatVector(PlotArrow.Alpha,
                                        (PlotArrow.Arrows + 1));
        PlotArrow.Alpha[PlotArrow.Arrows]= Alpha;

        PlotArrow.Red    = gomp_ReallocateFloatVector(PlotArrow.Red   , 
                                         (PlotArrow.Arrows + 1));
        PlotArrow.Green  = gomp_ReallocateFloatVector(PlotArrow.Green , 
                                         (PlotArrow.Arrows + 1));
        PlotArrow.Blue   = gomp_ReallocateFloatVector(PlotArrow.Blue  , 
                                         (PlotArrow.Arrows + 1));

        PlotArrow.Red[PlotArrow.Arrows]    = (float)ArrowRed;
        PlotArrow.Green[PlotArrow.Arrows]  = (float)ArrowGreen;
        PlotArrow.Blue[PlotArrow.Arrows]   = (float)ArrowBlue;

        PlotArrow.Arrows    +=         1;
    }

    gomp_SetPlotterRegistrationState(
        1, &PlotArrow.Plotter, PlotArrowSeg, NULL,
        PLOTTER_NAME_PLOT_ARROW, PLOTTER_ORDER_PLOT_ARROW);

    ADD_RESET_LISTENER(Arrow);

    return(0);
}
/****************************************************************************/
int gomp_DelArrowSeg()
/****************************************************************************/
{
    PlotArrowSegIsChanging();
    
    if(PlotArrow.Arrows) {
        free(PlotArrow.Red);
        free(PlotArrow.Green);
        free(PlotArrow.Blue);
        free(PlotArrow.ArrowCoord);
        free(PlotArrow.Alpha);
        free(PlotArrow.ArrowRad);
    }
    PlotArrow.Arrows = 0;

    gomp_UnregisterPlotter(PlotArrow.Plotter.plotter);
    PlotArrow.Plotter.plotter = NULL;

    CANCEL_RESET_LISTENER(Arrow);

    return(0);
}
/***********************************************************************/
int PlotArrowSeg(void *userData,int Wstr,int drawFlags)
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int i;
    static float rad;
    static float xi,yi,zi;
    static float xk,yk,zk;
    static int from,to;
    static float param[4];
    static const char *Value;
    static float value1,value2;

    if ( ! ( drawFlags & gom_PlotComplexElements ) || Wstr != 0 )
        return(-1);

    Value = Tcl_GetVar(gomp_GetTclInterp(),"gomCylinderArrowControl",TCL_GLOBAL_ONLY);

    if(Value) {
        sscanf(Value,"%f %f",&value1,&value2);
        (void)gomp_SetArrowAD2CD(value1);
        (void)gomp_SetArrowL2H(value2);
    } 

    glEnable(GL_LIGHTING);

    glLineWidth((GLfloat)(gomp_GetMoleculeLineWidth()+0.0));

    from = 0;
    to   = PlotArrow.Arrows;


/* throw out the Arrows */

    for(i = from ; i < to ; i++ ) {

        param[0]    = PlotArrow.Red[i];
        param[1]   = PlotArrow.Green[i];
        param[2]  = PlotArrow.Blue[i];
        param[3] = PlotArrow.Alpha[i];

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&param[0] , &param[1] , &param[2]);

        if(param[3] < 1.0) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA);
        } else {
            glBlendFunc(GL_ONE , GL_ZERO);
            glDisable(GL_BLEND);
        }

        glColor4f(param[0] , param[1] , param[2] , param[3]);
    
        rad = PlotArrow.ArrowRad[i];  /* radius of the Arrow */

        xi = PlotArrow.ArrowCoord[6 * i];
        yi = PlotArrow.ArrowCoord[6 * i + 1];
        zi = PlotArrow.ArrowCoord[6 * i + 2];

        xk = PlotArrow.ArrowCoord[6 * i + 3];
        yk = PlotArrow.ArrowCoord[6 * i + 4];
        zk = PlotArrow.ArrowCoord[6 * i + 5];

        gomp_Arrow(xi,yi,zi,xk,yk,zk,rad);

        glBlendFunc(GL_ONE , GL_ZERO);
        glDisable(GL_BLEND);

    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/*************************************************************************/
int gomp_PlotVectorData(void *userData,int Wstr,int drawFlags)
/*************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int    LineArrowRvalue;
    static int    i,j;
    static float  scale = 1.0;
    static float  rad   = 1.0;
    static const float *x;
    static const float *y;
    static const float *z;
    static const float *rc;
    static const float *gc;
    static const float *bc;
    static float  xt;
    static float  yt;
    static float  zt;
    static const int *Aindex;
    static const float *fx;
    static const float *fy;
    static const float *fz;
    static float  CRi;
    static float  CGi;
    static float  CBi;

    LineArrowRvalue = PlotLineArrowSeg( Wstr , drawFlags );

    if( ! ( drawFlags & gom_PlotComplexElements ) &&
        ! ( drawFlags & gom_PlotSimplifiedElements ) )
        return(LineArrowRvalue);

    if( Wstr != gomp_GetVectorStructureIndex() )
        return(LineArrowRvalue);

/*
  if(gomp_GetVectorListLength() != 
  gomp_GetNumAtomsInMolecStruct(k)) {
  gomp_PrintERROR("number of atoms in vector array != atom list");
  return(1);
  }
*/
    scale  = gomp_GetVectorScale();
    rad    = gomp_GetVectorRadius();
    fx     = gomp_GetVectorForceXp();
    fy     = gomp_GetVectorForceYp();
    fz     = gomp_GetVectorForceZp();
    x      = gomp_GetAtomXCoordPointer(Wstr);
    y      = gomp_GetAtomYCoordPointer(Wstr);
    z      = gomp_GetAtomZCoordPointer(Wstr);
    rc     = gomp_GetAtomColourRedPointer(Wstr);
    gc     = gomp_GetAtomColourGreenPointer(Wstr);
    bc     = gomp_GetAtomColourBluePointer(Wstr);
    Aindex = gomp_GetVectorListArray();

/* real tube based arrows */
    if ( drawFlags & gom_PlotComplexElements ) {
        glEnable(GL_LIGHTING);

        for(i = 0 ; i < gomp_GetVectorListLength() ; i++) {

            j  = Aindex[i];
            xt = x[j];
            yt = y[j];
            zt = z[j];

            CRi = rc[j];
            CGi = gc[j];
            CBi = bc[j];

            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&CRi , &CGi , &CBi);

            glColor4f(CRi , CGi , CBi , 1.0);

            (void)gomp_Arrow(xt , yt , zt , 
                        xt + scale *  fx[j],
                        yt + scale *  fy[j] ,
                        zt + scale *  fz[j] ,
                        rad);
        }
    }
/* just lines to get an idea about the directions and magnitudes */
    else {
        glDisable(GL_LIGHTING);

        for(i = 0 ; i < gomp_GetVectorListLength() ; i++) {

            j  = Aindex[i];
            xt = x[j];
            yt = y[j];
            zt = z[j];

            glColor4f(rc[j] , gc[j] , bc[j] , 1.0);

            (void)PlotOneLineArrow(xt,yt,zt,
                                   scale *  fx[j],
                                   scale *  fy[j],
                                   scale *  fz[j]);
        }
    }

    if( LineArrowRvalue > 0 )
        return(LineArrowRvalue);
#endif /* ENABLE_GRAPHICS */

    return(0);
}

#define PLOT_X_SIZE       30
#define PLOT_Y_TOP        40
#define PLOT_Y_BOT        40
#define PLOT_SIZE_LIMIT   100

static struct {
    int   plot;    /* on/off switch */
    float mmin;
    float mmax;
    float Range;
    int   num_level;
} ColourScalePlot = { 0 , 0.0 , 1.0f , 0.1f , 11};
#define ColourScalePlotIsChanging() /* do nothing */

int   gomp_DisplayColourScale(void);

int   gomp_SetColourScalePlotMin(float);
float gomp_GetColourScalePlotMin(void);
int   gomp_SetColourScalePlotMax(float);
float gomp_GetColourScalePlotMax(void);
int   gomp_SetColourScalePlotStep(float);
float gomp_GetColourScalePlotStep(void);
int   gomp_SetColourScalePlotLevels(int);
int   gomp_GetColourScalePlotLevels(void);
int   gomp_SetColourScalePlotStatus(int);
int   gomp_GetColourScalePlotStatus(void);

/*************************************************************************/
int   gomp_SetColourScalePlotMin(float Min)
/*************************************************************************/
{
    ColourScalePlotIsChanging();
    ColourScalePlot.mmin = Min;
    return(0);
} 
/*************************************************************************/
float gomp_GetColourScalePlotMin()
/*************************************************************************/
{
    return(ColourScalePlot.mmin);
} 
/*************************************************************************/
int   gomp_SetColourScalePlotMax(float Max)
/*************************************************************************/
{
    ColourScalePlotIsChanging();
    ColourScalePlot.mmax = Max;
    return(0);
} 
/*************************************************************************/
float gomp_GetColourScalePlotMax()
/*************************************************************************/
{
    return(ColourScalePlot.mmax);
} 

/*************************************************************************/
int   gomp_SetColourScalePlotStep(float Step)
/*************************************************************************/
{
    ColourScalePlotIsChanging();
    ColourScalePlot.Range = Step;
    return(0);
} 

/*************************************************************************/
float gomp_GetColourScalePlotStep()
/*************************************************************************/
{
    return(ColourScalePlot.Range);
} 

/*************************************************************************/
int   gomp_SetColourScalePlotLevels(int Levels)
/*************************************************************************/
{
    ColourScalePlotIsChanging();
    ColourScalePlot.num_level = Levels;
    return(0);
} 

/*************************************************************************/
int   gomp_GetColourScalePlotLevels()
/*************************************************************************/
{
    return(ColourScalePlot.num_level);
} 
/*************************************************************************/
int   gomp_SetColourScalePlotStatus(int StatusValue)
/*************************************************************************/
{
    ColourScalePlotIsChanging();
    ColourScalePlot.plot = StatusValue;
    return(0);
} 

/*************************************************************************/
int   gomp_GetColourScalePlotStatus()
/*************************************************************************/
{
    return(ColourScalePlot.plot);
} 

/*************************************************************************/
int gomp_DisplayColourScale()
/*************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int    x_step;
    static long   ysz;
    static char   text[BUFF_LEN];
    static float  co[2];

    static float  delta,delta1;
    static float  value;
    static float  rrggbb[3];
    static double step;
    static int    i;
    static int    WSize[4];
    static int    mm;
    static int    ColorMapping;
    static const char *TValue;
    static char   font[BUFF_LEN];

    TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomAtomLabelFont",TCL_GLOBAL_ONLY);

    if(TValue) {
        gomp_CopyString(font,TValue,BUFF_LEN);
    } else {
        gomp_CopyString(font,"BITMAP_HELVETICA_12",BUFF_LEN);
    }

    x_step = PLOT_X_SIZE;

    glGetIntegerv(GL_VIEWPORT , WSize);


/* check that there is space for the scale ... */
    if((WSize[3] - PLOT_Y_BOT - PLOT_Y_TOP) < PLOT_SIZE_LIMIT) {
        gomp_PrintERROR("the window is too small for the scale plot ");
        ColourScalePlot.plot = 0;
        return(1);
    }

    glDisable(GL_LIGHTING);
/*     glEnable(GL_LIGHTING);*/
    glColor3f(  1.0 , 1.0 , 1.0);
    glNormal3f( 0.0 , 0.0 , 1.0);

    glGetIntegerv(GL_MATRIX_MODE, &mm);
    glMatrixMode(GL_PROJECTION);
        
    glPushMatrix();
    glLoadIdentity();

    gluOrtho2D( 1 , WSize[2] , 1 , WSize[3] );

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
/*     glMatrixMode(mm); */

    glColor3f(1.0 , 1.0 , 1.0);

    ColorMapping = gomp_GetColorMappingType(); 
    if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
        (void)gomp_Prepare1DTexture();
        glEnable(GL_TEXTURE_1D);
    } else {
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
    }

    ysz = (WSize[3] - PLOT_Y_BOT - PLOT_Y_TOP) / 
        ColourScalePlot.num_level;

    if(ColourScalePlot.num_level - 1 <= 0) {
        gomp_PrintERROR("number of bins has to be > 1");
        ColourScalePlot.plot = 0;
        return(1);
    }

    delta  = 1.0 / (float)(ColourScalePlot.num_level - 1);

    delta1 = (ColourScalePlot.mmax - ColourScalePlot.mmin) /
        ((float)ColourScalePlot.num_level - 1.0);

    for(i = 0 ; i < ColourScalePlot.num_level ; i++) {

        step = (double)(i) * delta;

        if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
            glTexCoord1f( (float)step );
        } else {
            gomp_PreRainbow(step,&rrggbb[0],&rrggbb[1],&rrggbb[2]);
            glColor3fv(rrggbb);
        }

        glBegin(GL_QUADS);

        co[0] = (float)(x_step);
        co[1] = (float)(i * ysz) + PLOT_Y_BOT;

        glVertex2fv(co);

        co[0] = (float)(x_step);
        co[1] = (float)(i * ysz + ysz) + PLOT_Y_BOT;

        glVertex2fv(co);

        co[0] = 0.0;
        co[1] = (float)(i * ysz + ysz) + PLOT_Y_BOT;

        glVertex2fv(co);

        co[0] = 0.0;
        co[1] = (float)(i * ysz)     + PLOT_Y_BOT;

        glVertex2fv(co);

        glEnd();
    }

#if defined(WIN32)
    for(i = 0 ; i < ColourScalePlot.num_level ; i++) {
        co[0] = (float)(x_step);
        co[1] = (float)((i + i + 1) * ysz)/2. /*+ 1*/ + PLOT_Y_BOT;

        step = (double)(i) * delta ;

        if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {
            glTexCoord1f( (float)step );
        } else {
            gomp_PreRainbow(step,&rrggbb[0],&rrggbb[1],&rrggbb[2]);
            glColor3fv(rrggbb);
        }

        value = i * delta1 + ColourScalePlot.mmin;
        sprintf(text," %6.2f ",value);
        glRasterPos2f(co[0],co[1]);
        (void)gomp_PrintString(text , font);
    }
#else
/*
  Strange thing but I can't get the text to be colored correctly on Linux
  so I'm switching from 1D texture to normal.
*/
    glDisable(GL_TEXTURE_1D);
    glDisable(GL_LIGHTING);
    for(i = 0 ; i < ColourScalePlot.num_level ; i++) {
        co[0] = (float)(x_step);
        co[1] = (float)((i + i + 1) * ysz)/2. /*+ 1*/ + PLOT_Y_BOT;

        step = (double)(i) * delta ;

        gomp_PreRainbow(step,&rrggbb[0],&rrggbb[1],&rrggbb[2]);
        glColor3fv(rrggbb);

        value = i * delta1 + ColourScalePlot.mmin;
        sprintf(text," %6.2f ",value);
        glRasterPos2f(co[0],co[1]);
        (void)gomp_PrintString(text , font);
    }
#endif

/*      glGetIntegerv(GL_MATRIX_MODE, &mm); */
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(mm);

    glDisable(GL_TEXTURE_1D);
    glEnable(GL_LIGHTING);
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/*************************************************************************/
int gomp_WriteSphereSeg2ModelFile(FILE *Model_f)
/*************************************************************************/
{
    int i,i3,Show;

/* here we go ... */
    if(!PlotSphere.Spheres) return(0);

/* SPHERE - tag */
    fprintf(Model_f , "[Sphere]\n");
    Show = PlotSphere.Plotter.plotter != NULL;
    fprintf(Model_f , "%d %d\n",PlotSphere.Spheres , Show);

    i3    = 0;
    for(i = 0 ; i < PlotSphere.Spheres ; i++) {
        fprintf(Model_f , "%f %f %f %f %f %f\n",
                PlotSphere.SphereCoord[i3]     , 
                PlotSphere.SphereCoord[i3 + 1] ,
                PlotSphere.SphereCoord[i3 + 2] ,
                PlotSphere.XYZScale[i3]        ,
                PlotSphere.XYZScale[i3 + 1]    ,
                PlotSphere.XYZScale[i3 + 2]);
        fprintf(Model_f , "%f %f %f %f %f\n",
                PlotSphere.Radius[i]           ,
                PlotSphere.Alpha[i]            ,
                PlotSphere.Red[i]              ,
                PlotSphere.Green[i]            ,
                PlotSphere.Blue[i]);

        i3 += 3;
    }

    return(0);
}
/*************************************************************************/
int gomp_ReadSphereSegFromModelFile(FILE *Model_f)
/*************************************************************************/
{
    int   i,i3,Show;
    char  InputText[BUFF_LEN];

/* here we go ... */
    PlotSphereSegIsChanging();
 
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText , "%d %d",&PlotSphere.Spheres , &Show);
/* memory handling ....   */
    gomp_DelSphereSeg();
    PlotSphere.SphereCoord = gomp_AllocateFloatVector(3 * PlotSphere.Spheres);
    PlotSphere.XYZScale    = gomp_AllocateFloatVector(3 * PlotSphere.Spheres);
    PlotSphere.Radius      = gomp_AllocateFloatVector(PlotSphere.Spheres);
    PlotSphere.Alpha       = gomp_AllocateFloatVector(PlotSphere.Spheres);

    PlotSphere.Red         = gomp_AllocateFloatVector(PlotSphere.Spheres);
    PlotSphere.Green       = gomp_AllocateFloatVector(PlotSphere.Spheres);
    PlotSphere.Blue        = gomp_AllocateFloatVector(PlotSphere.Spheres);
/* end of memory handling */

    gomp_SetPlotterRegistrationState(
        Show, &PlotSphere.Plotter, PlotSphereSeg, NULL,
        PLOTTER_NAME_PLOT_SPHERE, PLOTTER_ORDER_PLOT_SPHERE);

    ADD_RESET_LISTENER(Sphere);

/* read the rest ... */
    i3    = 0;
    for(i = 0 ; i < PlotSphere.Spheres ; i++) {
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText , "%f %f %f %f %f %f\n",
               &PlotSphere.SphereCoord[i3    ] , 
               &PlotSphere.SphereCoord[i3 + 1] ,
               &PlotSphere.SphereCoord[i3 + 2] ,
               &PlotSphere.XYZScale[i3    ]    ,
               &PlotSphere.XYZScale[i3 + 1]    ,
               &PlotSphere.XYZScale[i3 + 2]);
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText , "%f %f %f %f %f\n",
               &PlotSphere.Radius[i]           ,
               &PlotSphere.Alpha[i]            ,
               &PlotSphere.Red[i]              ,
               &PlotSphere.Green[i]            ,
               &PlotSphere.Blue[i]);

        i3 += 3;
    }

    return(0);
}

/* ........ */

/*************************************************************************/
int gomp_WriteCylinderSeg2ModelFile(FILE *Model_f)
/*************************************************************************/
{
    int i,i6,Show;

/* here we go ... */
    if(!PlotCylinder.Cylinders) return(0);

/* Cylinder - tag */
    fprintf(Model_f , "[Cylinder]\n");
    Show = PlotCylinder.Plotter.plotter != 0;
    fprintf(Model_f , "%d %d\n",PlotCylinder.Cylinders , Show);

    i6    = 0;
    for(i = 0 ; i < PlotCylinder.Cylinders ; i++) {
        fprintf(Model_f , "%f %f %f %f %f %f\n",
                PlotCylinder.CylinderCoord[i6]     , 
                PlotCylinder.CylinderCoord[i6 + 1] ,
                PlotCylinder.CylinderCoord[i6 + 2] ,
                PlotCylinder.CylinderCoord[i6 + 3] ,
                PlotCylinder.CylinderCoord[i6 + 4] ,
                PlotCylinder.CylinderCoord[i6 + 5]);
        fprintf(Model_f , "%f %f %f %f %f\n",
                PlotCylinder.CylinderRad[i]      ,
                PlotCylinder.Alpha[i]            ,
                PlotCylinder.Red[i]              ,
                PlotCylinder.Green[i]            ,
                PlotCylinder.Blue[i]);

        i6 += 6;
    }

    return(0);
}
/*************************************************************************/
int gomp_ReadCylinderSegFromModelFile(FILE *Model_f)
/*************************************************************************/
{
    int   i,i6,Show;
    char  InputText[BUFF_LEN];

/* here we go ... */
    PlotCylinderSegIsChanging();
    
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText , "%d %d",&PlotCylinder.Cylinders , &Show);
/* memory handling ....   */
    gomp_DelCylinderSeg();
    PlotCylinder.CylinderCoord = gomp_AllocateFloatVector(6 * PlotCylinder.Cylinders);
    PlotCylinder.CylinderRad   = gomp_AllocateFloatVector(PlotCylinder.Cylinders);
    PlotCylinder.Alpha         = gomp_AllocateFloatVector(PlotCylinder.Cylinders);

    PlotCylinder.Red           = gomp_AllocateFloatVector(PlotCylinder.Cylinders);
    PlotCylinder.Green         = gomp_AllocateFloatVector(PlotCylinder.Cylinders);
    PlotCylinder.Blue          = gomp_AllocateFloatVector(PlotCylinder.Cylinders);
/* end of memory handling */

    gomp_SetPlotterRegistrationState(
        Show, &PlotCylinder.Plotter, PlotCylinderSeg, NULL,
        PLOTTER_NAME_PLOT_CYLINDER, PLOTTER_ORDER_PLOT_CYLINDER);

    ADD_RESET_LISTENER(Cylinder);

/* read the rest ... */
    i6    = 0;
    for(i = 0 ; i < PlotCylinder.Cylinders ; i++) {
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText , "%f %f %f %f %f %f\n",
               &PlotCylinder.CylinderCoord[i6    ] , 
               &PlotCylinder.CylinderCoord[i6 + 1] ,
               &PlotCylinder.CylinderCoord[i6 + 2] ,
               &PlotCylinder.CylinderCoord[i6 + 3] ,
               &PlotCylinder.CylinderCoord[i6 + 4] ,
               &PlotCylinder.CylinderCoord[i6 + 5]);
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText , "%f %f %f %f %f\n",
               &PlotCylinder.CylinderRad[i]      ,
               &PlotCylinder.Alpha[i]            ,
               &PlotCylinder.Red[i]              ,
               &PlotCylinder.Green[i]            ,
               &PlotCylinder.Blue[i]);

        i6 += 6;
    }

    return(0);
}


/* ........ */

/*************************************************************************/
int gomp_WriteArrowSeg2ModelFile(FILE *Model_f)
/*************************************************************************/
{
    int i,i6,Show;

/* here we go ... */
    if(!PlotArrow.Arrows) return(0);

/* Arrow - tag */
    fprintf(Model_f , "[Arrow]\n");
    Show = PlotArrow.Plotter.plotter != 0;
    fprintf(Model_f , "%d %d\n",PlotArrow.Arrows , Show);

    i6    = 0;
    for(i = 0 ; i < PlotArrow.Arrows ; i++) {
        fprintf(Model_f , "%f %f %f %f %f %f\n",
                PlotArrow.ArrowCoord[i6]     , 
                PlotArrow.ArrowCoord[i6 + 1] ,
                PlotArrow.ArrowCoord[i6 + 2] ,
                PlotArrow.ArrowCoord[i6 + 3] ,
                PlotArrow.ArrowCoord[i6 + 4] ,
                PlotArrow.ArrowCoord[i6 + 5]);
        fprintf(Model_f , "%f %f %f %f %f\n",
                PlotArrow.ArrowRad[i]           ,
                PlotArrow.Alpha[i]            ,
                PlotArrow.Red[i]              ,
                PlotArrow.Green[i]            ,
                PlotArrow.Blue[i]);

        i6 += 6;
    }

    return(0);
}
/*************************************************************************/
int gomp_ReadArrowSegFromModelFile(FILE *Model_f)
/*************************************************************************/
{
    int   i,i6,Show;
    char  InputText[BUFF_LEN];

/* here we go ... */
    PlotArrowSegIsChanging();
    
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText , "%d %d",&PlotArrow.Arrows , &Show);
/* memory handling ....   */
    gomp_DelArrowSeg();
    PlotArrow.ArrowCoord    = gomp_AllocateFloatVector(6 * PlotArrow.Arrows);
    PlotArrow.ArrowRad      = gomp_AllocateFloatVector(PlotArrow.Arrows);
    PlotArrow.Alpha         = gomp_AllocateFloatVector(PlotArrow.Arrows);

    PlotArrow.Red           = gomp_AllocateFloatVector(PlotArrow.Arrows);
    PlotArrow.Green         = gomp_AllocateFloatVector(PlotArrow.Arrows);
    PlotArrow.Blue          = gomp_AllocateFloatVector(PlotArrow.Arrows);
/* end of memory handling */

    gomp_SetPlotterRegistrationState(
        Show, &PlotArrow.Plotter, PlotArrowSeg, NULL,
        PLOTTER_NAME_PLOT_ARROW, PLOTTER_ORDER_PLOT_ARROW);

    ADD_RESET_LISTENER(Arrow);

/* read the rest ... */
    i6    = 0;
    for(i = 0 ; i < PlotArrow.Arrows ; i++) {
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText , "%f %f %f %f %f %f\n",
               &PlotArrow.ArrowCoord[i6    ] , 
               &PlotArrow.ArrowCoord[i6 + 1] ,
               &PlotArrow.ArrowCoord[i6 + 2] ,
               &PlotArrow.ArrowCoord[i6 + 3] ,
               &PlotArrow.ArrowCoord[i6 + 4] ,
               &PlotArrow.ArrowCoord[i6 + 5]);
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText , "%f %f %f %f %f\n",
               &PlotArrow.ArrowRad[i]           ,
               &PlotArrow.Alpha[i]            ,
               &PlotArrow.Red[i]              ,
               &PlotArrow.Green[i]            ,
               &PlotArrow.Blue[i]);

        i6 += 6;
    }

    return(0);
}
/* ........ */

/*************************************************************************/
int gomp_WritePlaneSeg2ModelFile(FILE *Model_f)
/*************************************************************************/
{
    int i,i12,Show;

/* here we go ... */
    if(!PlotPlane.Planes) return(0);

/* Plane - tag */
    fprintf(Model_f , "[Plane]\n");
    Show = PlotPlane.Plotter.plotter != 0;
    fprintf(Model_f , "%d %d\n",PlotPlane.Planes , Show);

    i12   = 0;
    for(i = 0 ; i < PlotPlane.Planes ; i++) {
        fprintf(Model_f , "%f %f %f %f %f %f %f %f %f %f %f %f\n",
                PlotPlane.PlaneCoord[i12]     , 
                PlotPlane.PlaneCoord[i12 + 1] ,
                PlotPlane.PlaneCoord[i12 + 2] ,
                PlotPlane.PlaneCoord[i12 + 3] ,
                PlotPlane.PlaneCoord[i12 + 4] ,
                PlotPlane.PlaneCoord[i12 + 5] ,
                PlotPlane.PlaneCoord[i12 + 6] ,
                PlotPlane.PlaneCoord[i12 + 7] ,
                PlotPlane.PlaneCoord[i12 + 8] ,
                PlotPlane.PlaneCoord[i12 + 9] ,
                PlotPlane.PlaneCoord[i12 + 10] ,
                PlotPlane.PlaneCoord[i12 + 11]);
        fprintf(Model_f , "%f %f %f %f\n",
                PlotPlane.Alpha[i]            ,
                PlotPlane.Red[i]              ,
                PlotPlane.Green[i]            ,
                PlotPlane.Blue[i]);

        i12 += 12;
    }

    return(0);
}
/*************************************************************************/
int gomp_ReadPlaneSegFromModelFile(FILE *Model_f)
/*************************************************************************/
{
    int   i,i12,Show;
    char  InputText[BUFF_LEN];

/* here we go ... */
    PlotPlaneSegIsChanging();
 
    gomp_Fgets(InputText,BUFF_LEN,Model_f);
    sscanf(InputText , "%d %d",&PlotPlane.Planes , &Show);
/* memory handling ....   */
    gomp_DelPlaneSeg();
    PlotPlane.PlaneCoord    = gomp_AllocateFloatVector(12 * PlotPlane.Planes);
    PlotPlane.Alpha         = gomp_AllocateFloatVector(PlotPlane.Planes);

    PlotPlane.Red           = gomp_AllocateFloatVector(PlotPlane.Planes);
    PlotPlane.Green         = gomp_AllocateFloatVector(PlotPlane.Planes);
    PlotPlane.Blue          = gomp_AllocateFloatVector(PlotPlane.Planes);
/* end of memory handling */

    gomp_SetPlotterRegistrationState(
        Show, &PlotPlane.Plotter, PlotPlaneSeg, NULL,
        PLOTTER_NAME_PLOT_PLANE, PLOTTER_ORDER_PLOT_PLANE);

    ADD_RESET_LISTENER(Plane);

/* read the rest ... */
    i12    = 0;
    for(i = 0 ; i < PlotPlane.Planes ; i++) {
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText , "%f %f %f %f %f %f %f %f %f %f %f %f\n",
               &PlotPlane.PlaneCoord[i12    ] , 
               &PlotPlane.PlaneCoord[i12 + 1] ,
               &PlotPlane.PlaneCoord[i12 + 2] ,
               &PlotPlane.PlaneCoord[i12 + 3] ,
               &PlotPlane.PlaneCoord[i12 + 4] ,
               &PlotPlane.PlaneCoord[i12 + 5] ,
               &PlotPlane.PlaneCoord[i12 + 6] ,
               &PlotPlane.PlaneCoord[i12 + 7] ,
               &PlotPlane.PlaneCoord[i12 + 8] ,
               &PlotPlane.PlaneCoord[i12 + 9] ,
               &PlotPlane.PlaneCoord[i12 + 10] ,
               &PlotPlane.PlaneCoord[i12 + 11]);
        gomp_Fgets(InputText,BUFF_LEN,Model_f);
        sscanf(InputText , "%f %f %f %f\n",
               &PlotPlane.Alpha[i]            ,
               &PlotPlane.Red[i]              ,
               &PlotPlane.Green[i]            ,
               &PlotPlane.Blue[i]);

        i12 += 12;
    }

    return(0);
}

/*************************************************************************/
int gomp_PlotCoordAxis(void *userData,int Wstr,int drawFlags)
/*************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int    i,j;
    static int    Switch;
    static float  Xaxis;
    static float  Yaxis;
    static float  Zaxis;
    static const int *ListP;
    static const int *LengthP;
    static const int *StructureP;
    static int    Structure;
    static const float *XcP;
    static const float *YcP;
    static const float *ZcP;
    static float  CR;
    static float  CG;
    static float  CB;
    static const char *TValue;
    static char   font[BUFF_LEN];

    if ( ! ( drawFlags & gom_PlotSimpleElements ) )
        return(-1);

    TValue = Tcl_GetVar(gomp_GetTclInterp(),"gomAtomLabelFont",TCL_GLOBAL_ONLY);

    if(TValue) {
        gomp_CopyString(font,TValue,BUFF_LEN);
    } else {
        gomp_CopyString(font,"BITMAP_HELVETICA_12",BUFF_LEN);
    }

    glLineWidth((GLfloat)(gomp_GetMoleculeLineWidth()+0.0));

    glDisable(GL_LIGHTING);

    if ( Wstr==0 && gomp_GetPlotAxisCoordPoints() ) {

        (void)gomp_GetPlotAxisLength(&Xaxis,&Yaxis,&Zaxis);

        XcP = gomp_GetPlotAxisXCoord();
        YcP = gomp_GetPlotAxisYCoord();
        ZcP = gomp_GetPlotAxisZCoord();

        for ( j = 0 ; j < gomp_GetPlotAxisCoordPoints() ; j++ ) {

/* X axis */
            CR = 1.0;
            CG = 0.0;
            CB = 0.0;
            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&CR , &CG , &CB);

            glColor3f(CR , CG , CB);
            glRasterPos3f(XcP[j] + Xaxis, YcP[j] , ZcP[j]);
            gomp_PrintString("X" , font);

            glBegin(GL_LINES);
            glVertex3f(XcP[j], YcP[j] , ZcP[j]);       
            glVertex3f(XcP[j] + Xaxis, YcP[j] , ZcP[j]);
            glEnd();
/* Y axis */
            CR = 0.0;
            CG = 1.0;
            CB = 0.0;
            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&CR , &CG , &CB);

            glColor3f(CR , CG , CB);
            glRasterPos3f(XcP[j], YcP[j]  + Yaxis, ZcP[j]);
            gomp_PrintString("Y" , font);

            glBegin(GL_LINES);
            glVertex3f(XcP[j] , YcP[j] , ZcP[j]);       
            glVertex3f(XcP[j] , YcP[j] + Yaxis, ZcP[j]);
            glEnd();
/* Z axis */
            CR = 0.0;
            CG = 0.0;
            CB = 1.0;
            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&CR , &CG , &CB);

            glColor3f(CR , CG , CB);
            glRasterPos3f(XcP[j], YcP[j] , ZcP[j]  + Yaxis);
            gomp_PrintString("Z" , font);

            glBegin(GL_LINES);
            glVertex3f(XcP[j] , YcP[j] , ZcP[j]);       
            glVertex3f(XcP[j] , YcP[j] , ZcP[j] + Zaxis);
            glEnd();
 
        }
    }

    if(!gomp_GetPlotAxisListLength())
        return(0);

    Switch     = 0;
    ListP      = gomp_GetPlotAxisListP();
    LengthP    = gomp_GetPlotAxisLengthP();
    StructureP = gomp_GetPlotAxisStructureP();
    (void)gomp_GetPlotAxisLength(&Xaxis,&Yaxis,&Zaxis);

    for ( i = 0 ; i < gomp_GetPlotAxisListLength() ; i++ ) {

        Structure = StructureP[i];

        if( Structure != Wstr ) continue;

        XcP       = gomp_GetAtomXCoordPointer(Structure);
        YcP       = gomp_GetAtomYCoordPointer(Structure);
        ZcP       = gomp_GetAtomZCoordPointer(Structure);

        for(j = 0 ; j < LengthP[i] ; j++) {

/* X axis */
            CR = 1.0;
            CG = 0.0;
            CB = 0.0;
            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&CR , &CG , &CB);

            glColor3f(CR , CG , CB);
            glRasterPos3f(XcP[ListP[Switch]] + Xaxis, YcP[ListP[Switch]] , 
                          ZcP[ListP[Switch]]);
            gomp_PrintString("X" , font);

            glBegin(GL_LINES);
            glVertex3f(XcP[ListP[Switch]], YcP[ListP[Switch]] , 
                       ZcP[ListP[Switch]]);       
            glVertex3f(XcP[ListP[Switch]] + Xaxis, YcP[ListP[Switch]] , 
                       ZcP[ListP[Switch]]);
            glEnd();
/* Y axis */
            CR = 0.0;
            CG = 1.0;
            CB = 0.0;
            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&CR , &CG , &CB);

            glColor3f(CR , CG , CB);
            glRasterPos3f(XcP[ListP[Switch]], YcP[ListP[Switch]]  + Yaxis, 
                          ZcP[ListP[Switch]]);
            gomp_PrintString("Y" , font);

            glBegin(GL_LINES);
            glVertex3f(XcP[ListP[Switch]] , YcP[ListP[Switch]] , 
                       ZcP[ListP[Switch]]);       
            glVertex3f(XcP[ListP[Switch]] , YcP[ListP[Switch]] + Yaxis, 
                       ZcP[ListP[Switch]]);
            glEnd();
/* Z axis */
            CR = 0.0;
            CG = 0.0;
            CB = 1.0;
            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&CR , &CG , &CB);

            glColor3f(CR , CG , CB);
            glRasterPos3f(XcP[ListP[Switch]], YcP[ListP[Switch]] , 
                          ZcP[ListP[Switch]]  + Yaxis);
            gomp_PrintString("Z" , font);

            glBegin(GL_LINES);
            glVertex3f(XcP[ListP[Switch]] , YcP[ListP[Switch]] , 
                       ZcP[ListP[Switch]]);       
            glVertex3f(XcP[ListP[Switch]] , YcP[ListP[Switch]] , 
                       ZcP[ListP[Switch]] + Zaxis);
            glEnd();
 
            Switch++;
        }
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}


/****************************************************************************/
int gomp_PushTriangleStack(
    const char *InputP1, const char *InputP2, const char *InputP3,
    const char *InputP4, const char *InputP5, const char *InputP6,
    const char *InputP7, const char *InputP8, const char *InputP9,
    const char *InputP10,const char *InputP11,const char *InputP12,
    const char *InputP13,const char *InputP14,const char *InputP15,
    const char *InputP16,const char *InputP17,const char *InputP18,
    const char *InputP19,const char *InputP20,const char *InputP21,
    const char *InputP22,const char *InputP23)
/****************************************************************************/
{
    float  Xc1,Yc1,Zc1;
    float  Xc2,Yc2,Zc2;
    float  Xc3,Yc3,Zc3;
    float  Un1,Un2,Un3;
    float  Vn1,Vn2,Vn3;
    float  Wn1,Wn2,Wn3;
    float  Cc1,Cc2,Cc3;
    float  Sum;
    int    iappend = 0;
    int    Location;
    float  Alpha;
    const float *sumxyz;

    Xc1 = atof(InputP1);
    Yc1 = atof(InputP2);
    Zc1 = atof(InputP3);

    Xc2 = atof(InputP7);
    Yc2 = atof(InputP8);
    Zc2 = atof(InputP9);

    Xc3 = atof(InputP13);
    Yc3 = atof(InputP14);
    Zc3 = atof(InputP15);

    Un1 = atof(InputP4);
    Vn1 = atof(InputP5);
    Wn1 = atof(InputP6);

    Un2 = atof(InputP10);
    Vn2 = atof(InputP11);
    Wn2 = atof(InputP12);

    Un3 = atof(InputP16);
    Vn3 = atof(InputP17);
    Wn3 = atof(InputP18);

    Cc1 = atof(InputP19);
    Cc2 = atof(InputP20);
    Cc3 = atof(InputP21);

    Alpha = 1.0;
    if(!gomp_StringMatch(InputP22,"appe$nd")) {
        Alpha = atof(InputP22);
    }

    iappend = 0;
    if(gomp_StringMatch(InputP23,"appe$nd") ||
       gomp_StringMatch(InputP22,"appe$nd") ) iappend = 1;
       
    sumxyz = gomp_GetTranslateArray();

/* check if the normals are set */

    if(sqrt(Un1 * Un1 + Vn1 * Vn1 + Wn1 * Wn1) +
       sqrt(Un2 * Un2 + Vn2 * Vn2 + Wn2 * Wn2) +
       sqrt(Un3 * Un3 + Vn3 * Vn3 + Wn3 * Wn3) < SURFACE_NORM_TRIGGER) {
        float x1,y1,z1;
        float x2,y2,z2;
        float tx,ty,tz;

        x1 = Xc2 - Xc1 ;
        y1 = Yc2 - Yc1;
        z1 = Zc2 - Zc1;

        x2 = Xc3 - Xc1;
        y2 = Yc3 - Yc1;
        z2 = Zc3 - Zc1;

        tx = y1 * z2 - z1 * y2;
        ty = z1 * x2 - x1 * z2;
        tz = x1 * y2 - y1 * x2;

        Sum = sqrt(tx*tx + ty*ty + tz*tz);
        tx  = tx/Sum;
        ty  = ty/Sum;
        tz  = tz/Sum;

        Un1  = tx;
        Vn1  = ty;
        Wn1  = tz;
      
        Un2  = tx;
        Vn2  = ty;
        Wn2  = tz;

        Un3  = tx;
        Vn3  = ty;
        Wn3  = tz;
    }

    PlotTriangeMeshIsChanging();

    if(!iappend)
        gomp_DelTriangleSeg();

    if(!PlotTriangle.Triangles) { /* no old triangles */

        PlotTriangle.x = gomp_AllocateFloatVector(3);
        PlotTriangle.y = gomp_AllocateFloatVector(3);
        PlotTriangle.z = gomp_AllocateFloatVector(3);

        PlotTriangle.u = gomp_AllocateFloatVector(3);
        PlotTriangle.v = gomp_AllocateFloatVector(3);
        PlotTriangle.w = gomp_AllocateFloatVector(3);

        PlotTriangle.c = gomp_AllocateFloatVector(3);

        PlotTriangle.Alpha = gomp_AllocateFloatVector(1);

        PlotTriangle.x[0]  = Xc1 - sumxyz[0];
        PlotTriangle.y[0]  = Yc1 - sumxyz[1];
        PlotTriangle.z[0]  = Zc1 - sumxyz[2];

        PlotTriangle.x[1]  = Xc2 - sumxyz[0];
        PlotTriangle.y[1]  = Yc2 - sumxyz[1];
        PlotTriangle.z[1]  = Zc2 - sumxyz[2];

        PlotTriangle.x[2]  = Xc3 - sumxyz[0];
        PlotTriangle.y[2]  = Yc3 - sumxyz[1];
        PlotTriangle.z[2]  = Zc3 - sumxyz[2];

        PlotTriangle.u[0]  = Un1;
        PlotTriangle.v[0]  = Vn1;
        PlotTriangle.w[0]  = Wn1;
      
        PlotTriangle.u[1]  = Un2;
        PlotTriangle.v[1]  = Vn2;
        PlotTriangle.w[1]  = Wn2;

        PlotTriangle.u[2]  = Un3;
        PlotTriangle.v[2]  = Vn3;
        PlotTriangle.w[2]  = Wn3;

        PlotTriangle.c[0]  = Cc1;
        PlotTriangle.c[1]  = Cc2;
        PlotTriangle.c[2]  = Cc3;

        PlotTriangle.Alpha[0]    =     Alpha;
        PlotTriangle.Triangles   =         1;
    } else {
        PlotTriangle.x = gomp_ReallocateFloatVector(PlotTriangle.x,
                                            (PlotTriangle.Triangles + 1) * 3);
        PlotTriangle.y = gomp_ReallocateFloatVector(PlotTriangle.y,
                                            (PlotTriangle.Triangles + 1) * 3);
        PlotTriangle.z = gomp_ReallocateFloatVector(PlotTriangle.z,
                                            (PlotTriangle.Triangles + 1) * 3);

        PlotTriangle.u = gomp_ReallocateFloatVector(PlotTriangle.u,
                                            (PlotTriangle.Triangles + 1) * 3);
        PlotTriangle.v = gomp_ReallocateFloatVector(PlotTriangle.v,
                                            (PlotTriangle.Triangles + 1) * 3);
        PlotTriangle.w = gomp_ReallocateFloatVector(PlotTriangle.w,
                                            (PlotTriangle.Triangles + 1) * 3);

        PlotTriangle.c = gomp_ReallocateFloatVector(PlotTriangle.c,
                                            (PlotTriangle.Triangles + 1) * 3);

        PlotTriangle.Alpha = gomp_ReallocateFloatVector(PlotTriangle.Alpha,
                                                (PlotTriangle.Triangles + 1));

        Location = 3 * PlotTriangle.Triangles;

        PlotTriangle.x[Location]  = Xc1 - sumxyz[0];
        PlotTriangle.y[Location]  = Yc1 - sumxyz[1];
        PlotTriangle.z[Location]  = Zc1 - sumxyz[2];

        PlotTriangle.x[Location + 1]  = Xc2 - sumxyz[0];
        PlotTriangle.y[Location + 1]  = Yc2 - sumxyz[1];
        PlotTriangle.z[Location + 1]  = Zc2 - sumxyz[2];

        PlotTriangle.x[Location + 2]  = Xc3 - sumxyz[0];
        PlotTriangle.y[Location + 2]  = Yc3 - sumxyz[1];
        PlotTriangle.z[Location + 2]  = Zc3 - sumxyz[2];

        PlotTriangle.u[Location]      = Un1;
        PlotTriangle.v[Location]      = Vn1;
        PlotTriangle.w[Location]      = Wn1;

        PlotTriangle.u[Location + 1]  = Un2;
        PlotTriangle.v[Location + 1]  = Vn2;
        PlotTriangle.w[Location + 1]  = Wn2;

        PlotTriangle.u[Location + 2]  = Un3;
        PlotTriangle.v[Location + 2]  = Vn3;
        PlotTriangle.w[Location + 2]  = Wn3;

        PlotTriangle.c[Location    ]  = Cc1;
        PlotTriangle.c[Location + 1]  = Cc2;
        PlotTriangle.c[Location + 2]  = Cc3;

        PlotTriangle.Alpha[PlotTriangle.Triangles]   = Alpha; 
 
        PlotTriangle.Triangles  +=         1;
    }

    gomp_SetPlotterRegistrationState(
        1, &PlotTriangle.Plotter, PlotTriangleSeg, NULL,
        PLOTTER_NAME_PLOT_TRIANGLE, PLOTTER_ORDER_PLOT_TRIANGLE);

    ADD_RESET_LISTENER(Triangle);
        
    return(0);
}
/****************************************************************************/
int PlotTriangleSeg(void *userData,int Wstr,int drawFlags)
/****************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int i,tri;
    static int ColorMapping;
    static float Red;
    static float Green;
    static float Blue;

    if ( ! ( drawFlags & gom_PlotComplexElements ) || Wstr != 0 )
        return(-1);

    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    (void)gomp_Prepare1DTexture();
    glDisable(GL_TEXTURE_1D);

    ColorMapping = gomp_GetColorMappingType();

    if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE)
        glEnable(GL_TEXTURE_1D);
    else
        glDisable(GL_TEXTURE_1D);

    for (tri = 0; tri < PlotTriangle.Triangles; tri++) {

        i = 3 * tri;

        if(ColorMapping == COLOR_MAPPING_TYPE_TEXTURE) {

            if(PlotTriangle.Alpha[tri] < 1.0) {
                glBlendFunc(GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA);
                glColor4f(1.0 , 1.0 , 1.0 , PlotTriangle.Alpha[tri]);
            }
            else {
                glBlendFunc(GL_ONE , GL_ZERO);
                glColor4f(1.0 , 1.0 , 1.0 , 1.0);
            }

            glBegin(GL_POLYGON);
            glTexCoord1f(PlotTriangle.c[i]);
            glNormal3f(PlotTriangle.u[i] , 
                       PlotTriangle.v[i] , 
                       PlotTriangle.w[i]);
            glVertex3f(PlotTriangle.x[i] , 
                       PlotTriangle.y[i] , 
                       PlotTriangle.z[i]);

            glTexCoord1f(PlotTriangle.c[i + 1]);
            glNormal3f(PlotTriangle.u[i + 1] , 
                       PlotTriangle.v[i + 1] , 
                       PlotTriangle.w[i + 1]);
            glVertex3f(PlotTriangle.x[i + 1] , 
                       PlotTriangle.y[i + 1] , 
                       PlotTriangle.z[i + 1]);

            glTexCoord1f(PlotTriangle.c[i + 2]);
            glNormal3f(PlotTriangle.u[i + 2] , 
                       PlotTriangle.v[i + 2] , 
                       PlotTriangle.w[i + 2]);
            glVertex3f(PlotTriangle.x[i + 2] , 
                       PlotTriangle.y[i + 2] , 
                       PlotTriangle.z[i + 2]);
            glEnd();
        } else {

            if(PlotTriangle.Alpha[tri] < 1.0) {
                glBlendFunc(GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA);
                glColor4f(1.0 , 1.0 , 1.0 , PlotTriangle.Alpha[tri]);
            }
            else {
                glBlendFunc(GL_ONE , GL_ZERO);
                glColor4f(1.0 , 1.0 , 1.0 , 1.0);
            }

            glBegin(GL_POLYGON);
            gomp_PreRainbow((double)PlotTriangle.c[i] , &Red , &Green , &Blue);
            glColor4f(Red , Green , Blue , PlotTriangle.Alpha[tri]);
            glNormal3f(PlotTriangle.u[i] , 
                       PlotTriangle.v[i] , 
                       PlotTriangle.w[i]);
            glVertex3f(PlotTriangle.x[i] , 
                       PlotTriangle.y[i] , 
                       PlotTriangle.z[i]);

            gomp_PreRainbow((double)PlotTriangle.c[i + 1] , &Red , &Green , &Blue);
            glColor4f(Red , Green , Blue , PlotTriangle.Alpha[tri]);
            glNormal3f(PlotTriangle.u[i + 1] , 
                       PlotTriangle.v[i + 1] , 
                       PlotTriangle.w[i + 1]);
            glVertex3f(PlotTriangle.x[i + 1] , 
                       PlotTriangle.y[i + 1] , 
                       PlotTriangle.z[i + 1]);

            gomp_PreRainbow((double)PlotTriangle.c[i + 2] , &Red , &Green , &Blue);
            glColor4f(Red , Green , Blue , PlotTriangle.Alpha[tri]);
            glNormal3f(PlotTriangle.u[i + 2] , 
                       PlotTriangle.v[i + 2] , 
                       PlotTriangle.w[i + 2]);
            glVertex3f(PlotTriangle.x[i + 2] , 
                       PlotTriangle.y[i + 2] , 
                       PlotTriangle.z[i + 2]);
            glEnd();
        }
    }
#if defined(DEBUG)
    printf("1: %f %f %f | %f %f %f\n",PlotTriangle.x[i] , 
           PlotTriangle.y[i] , 
           PlotTriangle.z[i], PlotTriangle.u[i] , 
           PlotTriangle.v[i] , 
           PlotTriangle.w[i]);
    printf("2: %f %f %f | %f %f %f\n",PlotTriangle.x[i+1] , 
           PlotTriangle.y[i+1] , 
           PlotTriangle.z[i+1], PlotTriangle.u[i+1] , 
           PlotTriangle.v[i+1] , 
           PlotTriangle.w[i+1]);
    printf("3: %f %f %f | %f %f %f\n",PlotTriangle.x[i+2] , 
           PlotTriangle.y[i+2] , 
           PlotTriangle.z[i+2], PlotTriangle.u[i+2] , 
           PlotTriangle.v[i+2] , 
           PlotTriangle.w[i+2]);
#endif

    glDisable( GL_TEXTURE_1D );
    glDisable( GL_NORMALIZE  );
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/****************************************************************************/
int gomp_DelTriangleSeg()
/****************************************************************************/
{
    PlotTriangeMeshIsChanging();
    
    if(PlotTriangle.Triangles) {
        free(PlotTriangle.x);
        free(PlotTriangle.y);
        free(PlotTriangle.z);
        free(PlotTriangle.u);
        free(PlotTriangle.v);
        free(PlotTriangle.w);
        free(PlotTriangle.c);
    }
    PlotTriangle.Triangles = 0;

    gomp_UnregisterPlotter(PlotTriangle.Plotter.plotter);
    PlotTriangle.Plotter.plotter = NULL;

    CANCEL_RESET_LISTENER(Triangle);

    return(0);
}
/***********************************************************************/
int PlotOneLineArrow(
    float xpt1,float ypt1,float zpt1,float xpt2,float ypt2,float zpt2)  
    /* plot a line arrow */
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static float length;
    static float xp;
    static float yp;
    static float zp;
    static float cangle = 0.7071f;
    static float xpt3, xp3, x3;
    static float ypt3, yp3, y3;
    static float zpt3, zp3, z3;
    static float xpt4, xp4, x4;
    static float ypt4, yp4, y4;
    static float zpt4, zp4, z4;
    static float xpt5, xp5, x5;
    static float ypt5, yp5, y5;
    static float zpt5, zp5, z5;
    static float xpt6, xp6, x6;
    static float ypt6, yp6, y6;
    static float zpt6, zp6, z6;
    static float cosa,cosb,sina,sinb,temp;

    length = sqrt(xpt2 * xpt2 + ypt2 * ypt2 + zpt2 * zpt2);
    if(length < 1.0e-06) return(0);
    temp   = (length / 5.0) * cangle;
    
    xp = xpt2;
    yp = ypt2;
    zp = zpt2;

    gomp_Vector2Angles(0.0,0.0,0.0,&xp,&yp,&zp);

    cosa = cos(M_PI * yp /180.); 
    sina = sin(M_PI * yp /180.); 

    cosb = cos(M_PI * xp /180.); 
    sinb = sin(M_PI * xp /180.); 

    xpt3 = 0.0;
    ypt3 =  temp;
    zpt3 = length - temp;
    xpt4 = 0.0;
    ypt4 = -temp;
    zpt4 = length - temp;
    xpt5 = temp;
    ypt5 = 0.0;
    zpt5 = length - temp;
    xpt6 = -temp;
    ypt6 = 0.0;
    zpt6 = length - temp;

    xp3  =  cosa * xpt3 + sina * zpt3;
    yp3  =  ypt3;
    zp3  = -sina * xpt3 + cosa * zpt3;

    xp4  =  cosa * xpt4 + sina * zpt4;
    yp4  =  ypt4;
    zp4  = -sina * xpt4 + cosa * zpt4;

    xp5  =  cosa * xpt5 + sina * zpt5;
    yp5  =  ypt5;
    zp5  = -sina * xpt5 + cosa * zpt5;

    xp6  =  cosa * xpt6 + sina * zpt6;
    yp6  =  ypt6;
    zp6  = -sina * xpt6 + cosa * zpt6;

    x3   =  cosb * xp3 - sinb * yp3;
    y3   =  sinb * xp3 + cosb * yp3;
    z3   =  zp3;

    x4   =  cosb * xp4 - sinb * yp4;
    y4   =  sinb * xp4 + cosb * yp4;
    z4   =  zp4;

    x5   =  cosb * xp5 - sinb * yp5;
    y5   =  sinb * xp5 + cosb * yp5;
    z5   =  zp5;

    x6   =  cosb * xp6 - sinb * yp6;
    y6   =  sinb * xp6 + cosb * yp6;
    z6   =  zp6;

    xp = xpt2 + xpt1;
    yp = ypt2 + ypt1;
    zp = zpt2 + zpt1;

    glBegin(GL_LINES); 
    glVertex3f(xpt1       , ypt1      , zpt1);
    glVertex3f(xp, yp, zp);
    glEnd();

    glBegin(GL_LINE_STRIP);
    glVertex3f(xp , yp , zp);
    glVertex3f(x3 + xpt1 , y3 +ypt1 , z3 + zpt1);
    glVertex3f(x4 + xpt1 , y4 +ypt1 , z4 + zpt1);
    glVertex3f(xp , yp , zp);
    glEnd();

    glBegin(GL_LINE_STRIP);
    glVertex3f(xp , yp , zp);
    glVertex3f(x5 + xpt1 , y5 + ypt1 , z5 + zpt1);
    glVertex3f(x6 + xpt1 , y6 + ypt1 , z6 + zpt1);
    glVertex3f(xp , yp , zp);
    glEnd();
#endif /* ENABLE_GRAPHICS */

    return(0);
}

/****************************************************************************/
int gomp_PushLineArrowStack(float Xc1,float Yc1,float Zc1,
                            float Xc2,float Yc2,float Zc2,
                            const char *InputP7,const char *InputP8)
/****************************************************************************/
{
    int    iappend = 0;
    int    Location;
    float  ArrowRed;
    float  ArrowGreen;
    float  ArrowBlue;
    const float *sumxyz;

    ArrowRed = 1.0;
    ArrowGreen = ArrowBlue = 0.0;

    sumxyz = gomp_GetTranslateArray();

    if(gomp_StringMatch(InputP7,"appe$nd")) {
        iappend = 1;
    } else { 

        if(gomp_ColourName2RGB(InputP7 , &ArrowRed , 
                          &ArrowGreen , 
                          &ArrowBlue)) {
            gomp_PrintMessage("?ERROR - can't resolve the colour");
            return(1);
        }
    }

    if(gomp_StringMatch(InputP8,"appe$nd")) iappend = 1;

    if(!iappend)
        gomp_DelLineArrowSeg();


    if(!PlotLineArrow.Arrows) { /* no old Arrows */

        PlotLineArrow.ArrowCoord = gomp_AllocateFloatVector(6);
        PlotLineArrow.ArrowCoord[0]   = Xc1 - sumxyz[0];
        PlotLineArrow.ArrowCoord[1]  = Yc1 - sumxyz[1];
        PlotLineArrow.ArrowCoord[2] = Zc1 - sumxyz[2];
        PlotLineArrow.ArrowCoord[3]   = Xc2;
        PlotLineArrow.ArrowCoord[4]  = Yc2;
        PlotLineArrow.ArrowCoord[5] = Zc2;

        PlotLineArrow.Red    = gomp_AllocateFloatVector(1);
        PlotLineArrow.Green  = gomp_AllocateFloatVector(1);
        PlotLineArrow.Blue   = gomp_AllocateFloatVector(1);

        PlotLineArrow.Red[0]    = (float)ArrowRed;
        PlotLineArrow.Green[0]  = (float)ArrowGreen;
        PlotLineArrow.Blue[0]   = (float)ArrowBlue;
      
        PlotLineArrow.Arrows     =         1;
      
        return(0);
    } else {
        PlotLineArrow.ArrowCoord = gomp_ReallocateFloatVector(PlotLineArrow.ArrowCoord,
                                                 (PlotLineArrow.Arrows + 1) * 6);
        Location = 6 * PlotLineArrow.Arrows;
        PlotLineArrow.ArrowCoord[Location]       = Xc1 - sumxyz[0];
        PlotLineArrow.ArrowCoord[Location + 1]  = Yc1 - sumxyz[1];
        PlotLineArrow.ArrowCoord[Location + 2] = Zc1 - sumxyz[2];
        PlotLineArrow.ArrowCoord[Location + 3]   = Xc2;
        PlotLineArrow.ArrowCoord[Location + 4]  = Yc2;
        PlotLineArrow.ArrowCoord[Location + 5] = Zc2;

        PlotLineArrow.Red    = gomp_ReallocateFloatVector(PlotLineArrow.Red   , 
                                             (PlotLineArrow.Arrows + 1));
        PlotLineArrow.Green  = gomp_ReallocateFloatVector(PlotLineArrow.Green , 
                                             (PlotLineArrow.Arrows + 1));
        PlotLineArrow.Blue   = gomp_ReallocateFloatVector(PlotLineArrow.Blue  , 
                                             (PlotLineArrow.Arrows + 1));

        PlotLineArrow.Red[PlotLineArrow.Arrows]    = (float)ArrowRed;
        PlotLineArrow.Green[PlotLineArrow.Arrows]  = (float)ArrowGreen;
        PlotLineArrow.Blue[PlotLineArrow.Arrows]   = (float)ArrowBlue;

        PlotLineArrow.Arrows    +=         1;

        return(0);
    }
}
/****************************************************************************/
int gomp_DelLineArrowSeg()
/****************************************************************************/
{
    if(PlotLineArrow.Arrows) {
        free(PlotLineArrow.Red);
        free(PlotLineArrow.Green);
        free(PlotLineArrow.Blue);
        free(PlotLineArrow.ArrowCoord);
    }
    PlotLineArrow.Arrows = 0;

    return(0);
}

#ifdef ENABLE_GRAPHICS
/***********************************************************************/
int PlotLineArrowSeg(int Wstr,int drawFlags)   /* plot line arrow(s) */
/***********************************************************************/
{
    static int i;
    static float xi,yi,zi;
    static float xk,yk,zk;
    static int from,to;
    static float param[3];
    static float vlen;
    static float scale;
    static float range_min;
    static float range_max;
    static float delta;
    static float step;
    static const char *value;
    static int   ArrowType = 0;
    static float ArrowRad  = 0.1f;

    if( Wstr != 0 ) return(-1);

    ArrowType = 0;
    ArrowRad  = 0.1f;

    value  = Tcl_GetVar(gomp_GetTclInterp() , "gomVectorDefaultType", TCL_GLOBAL_ONLY);
    if(value) {
        sscanf(value,"%d %f",&ArrowType,&ArrowRad);
        if(ArrowRad < 0.0001) ArrowRad = 0.1f;
    } 

    if(ArrowType) {
        if( !( drawFlags & gom_PlotComplexElements ) &&
            !( drawFlags & gom_PlotSimplifiedElements ) )
            return(0);

        if ( drawFlags & gom_PlotComplexElements )
            glEnable(GL_LIGHTING);
        else
            glDisable(GL_LIGHTING);
    } else {
        if ( ! ( drawFlags & gom_PlotSimpleElements ) )
            return(0);

        glDisable(GL_LIGHTING);
    }

    glLineWidth((GLfloat)(gomp_GetMoleculeLineWidth()+0.0));

    from = 0;
    to   = PlotLineArrow.Arrows;
    scale  = gomp_GetVectorScale();
    (void)gomp_GetVectorDisplayRange(&range_min , &range_max);
    delta = range_max - range_min;
/* throw out the Arrows */

    for(i = from ; i < to ; i++ ) {

        xi = PlotLineArrow.ArrowCoord[6 * i];
        yi = PlotLineArrow.ArrowCoord[6 * i + 1];
        zi = PlotLineArrow.ArrowCoord[6 * i + 2];

        xk = PlotLineArrow.ArrowCoord[6 * i + 3];
        yk = PlotLineArrow.ArrowCoord[6 * i + 4];
        zk = PlotLineArrow.ArrowCoord[6 * i + 5];

        xk *= scale;
        yk *= scale;
        zk *= scale;

        vlen = sqrt(xk * xk + yk * yk + zk * zk);

        if((vlen < range_min) || (vlen > range_max)) continue;

        if(GetGradientDisplayStyle() == 0) {
            param[0]    = PlotLineArrow.Red[i];
            param[1]   = PlotLineArrow.Green[i];
            param[2]  = PlotLineArrow.Blue[i];
        } else if (GetGradientDisplayStyle() == 1) {
            step = (double)((vlen - range_min) / delta);
 
            gomp_PreRainbow(step,&param[0],&param[1],&param[2]);

        }

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&param[0] , &param[1] , &param[2]);

        glColor3f(param[0] , param[1] , param[2] );

        if(ArrowType) {
            if ( drawFlags & gom_PlotComplexElements )
                gomp_Arrow(xi,yi,zi,(xi + xk),(yi + yk),(zi + zk),ArrowRad);
            else
                (void)PlotOneLineArrow(xi,yi,zi,xk,yk,zk);
        } else
            (void)PlotOneLineArrow(xi,yi,zi,xk,yk,zk);
    }

    return(0);
}
#endif /* ENABLE_GRAPHICS */
/***********************************************************************/
int gomp_GetTotalLineArrowEntries() 
/***********************************************************************/
{
    return(PlotLineArrow.Arrows);
}
#if 0
/***********************************************************************/
int GetLineArrowEntry(int Which,
                      float *xc1 , float *yc1, float *zc1,
                      float *xc2 , float *yc2, float *zc2,
                      float *redc, float *greenc, float *bluec) 
/***********************************************************************/
{
    if(!PlotLineArrow.Arrows) {
        gomp_PrintERROR("no line arrow entries available");
        return(1);
    }

    if(Which < 1 || Which > PlotLineArrow.Arrows) {
        gomp_PrintERROR("line arrow index is outside allowed range");
        return(1);
    }

    *xc1 = PlotLineArrow.ArrowCoord[6 * Which];
    *yc1 = PlotLineArrow.ArrowCoord[6 * Which + 1];
    *zc1 = PlotLineArrow.ArrowCoord[6 * Which + 2];

    *xc2 = PlotLineArrow.ArrowCoord[6 * Which + 3];
    *yc2 = PlotLineArrow.ArrowCoord[6 * Which + 4];
    *zc2 = PlotLineArrow.ArrowCoord[6 * Which + 5];

    *redc   = PlotLineArrow.Red[Which];
    *greenc = PlotLineArrow.Green[Which];
    *bluec  = PlotLineArrow.Blue[Which];

    return(0);
}

/***********************************************************************/
int SetLineArrowEntryColor(int Which, float redc, float greenc, float bluec) 
/***********************************************************************/
{
    if(!PlotLineArrow.Arrows) {
        gomp_PrintERROR("no line arrow entries available");
        return(1);
    }

    if(Which < 1 || Which > PlotLineArrow.Arrows) {
        gomp_PrintERROR("line arrow index is outside allowed range");
        return(1);
    }

    PlotLineArrow.Red[Which] = redc;
    PlotLineArrow.Green[Which] = greenc;
    PlotLineArrow.Blue[Which]  = bluec;

    return(0);
}
#endif
/***********************************************************************/
int gomp_SetGradientDisplayStyle(int Type) 
/***********************************************************************/
{
    PlotLineArrow.Style = Type;

    return(0);
}
/***********************************************************************/
int GetGradientDisplayStyle() 
/***********************************************************************/
{
    return(PlotLineArrow.Style);
}
