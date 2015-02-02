/*

Copyright (c) 1994 - 2005 by:
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
#include "plot_cpk.h"
#include "plot_molec.h"

#include "stdafx.h"

#ifdef ENABLE_GRAPHICS
GLUquadricObj    *gomp_SphereQuad = NULL;
#endif /* ENABLE_GRAPHICS */


/***********************************************************************/
int gomp_PlotMoleculeCPK(int Wstr)          /* plot molecule as cpk spheres */
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS
    register int k;
    static float rad;
    register const float *ColRED,*ColGREEN,*ColBLUE;
    register const float *XC,*YC,*ZC;
    static int from,to;
    static const char *DispList;
    static const char *CPKList;
    static const float *Scale_CPK;
    static int   SphereQ;
    static int   SphereQ2;
    static float  CRk;
    static float  CGk;
    static float  CBk;

    glEnable(GL_LIGHTING);

    SphereQ   = gomp_GetSphereQuality();
    SphereQ2  = 2 * SphereQ;

    from      = 0;
    to        = gomp_GetNumAtomsInMolecStruct(Wstr);

    ColRED    = gomp_GetAtomColourRedPointer(Wstr);
    ColGREEN = gomp_GetAtomColourGreenPointer(Wstr);
    ColBLUE = gomp_GetAtomColourBluePointer(Wstr);
    XC        = gomp_GetAtomXCoordPointer(Wstr);
    YC       = gomp_GetAtomYCoordPointer(Wstr);
    ZC      = gomp_GetAtomZCoordPointer(Wstr);
    DispList  = gomp_GetAtomDisplayStatePointer(Wstr);
    CPKList   = gomp_GetAtomCPKDisplayStatePointer(Wstr);
    Scale_CPK = gomp_GetAtomCPKScalePointer(Wstr);
    
    for(k = from ; k < to ; k++) {      /* main loop */

        if(DispList[k] != 1) continue;  /* check if the atom should be displayed */

        if(CPKList[k]  != 1) continue; /* check if there should be a surface */

        CRk   = ColRED[k];
        CGk  = ColGREEN[k];
        CBk = ColBLUE[k];

        if(!gomp_GetDisplayColourType())
            (void)gomp_RGB2Grayscale(&CRk , &CGk , &CBk);

        glColor4f(CRk , CGk , CBk , 1.0);

        rad = Scale_CPK[k] * gomp_GetAtomVdwRad(Wstr , k);

        glPushMatrix();
        glTranslatef(XC[k],YC[k],ZC[k]); 
        gluSphere(gomp_SphereQuad, (double)rad , SphereQ2 , SphereQ);
        glPopMatrix();
    }
#endif /* ENABLE_GRAPHICS */

    return(0);    
}

/***********************************************************************/
int gomp_SetSphereQuad()
/***********************************************************************/
{
#ifdef ENABLE_GRAPHICS
    gomp_SphereQuad = gluNewQuadric();
    gluQuadricDrawStyle(gomp_SphereQuad, GLU_FILL);
    gluQuadricOrientation(gomp_SphereQuad , GLU_OUTSIDE);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
