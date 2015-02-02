/*

Copyright (c) 1990 - 2005 by:
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
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <ctype.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>

#if defined(MOTIF)
#include <GL/GLwMDrawA.h>
#endif

#include <X11/X.h>
#include <X11/Xlib.h>

#if !defined(WIN32)
#include <X11/Intrinsic.h>
#endif

#if defined(IRIX)
#include <X11/StringDefs.h>
#endif
#include <X11/keysym.h>
#include <X11/cursorfont.h>
#endif /* ENABLE_GRAPHICS */

#include "molecule.h"
#include "plot_molec.h"
#include "plot.h"
#include "printmsg.h"
#include "trace.h"
#include "trajectory.h"

#include "stdafx.h"

/************************************************************************/
int gomp_PlotAtomTrace(void* userData,int Wstr,int drawFlags)
/************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int    i,j,k,jj,loop,loopI,in_set,in_set1;
    static float  vec[3];
    static int    first_frame,last_frame,delta_frame;
    static const int *trace_atoms;
    static const int *trace_list;
    static const float *trcx;
    static const float *trcy;
    static const float *trcz;
    static const float *red;
    static const float *green;
    static const float *blue;
    static int    temp;

    if ( ! ( drawFlags & gom_PlotSimpleElements ) || Wstr != 0 )
        return(-1);

    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintERROR("Number of frames is not defined ");
        return(1);
    }

    if(!gomp_GetTraceState()) {
        gomp_PrintERROR("No trace defined to be shown");
        return(1);
    }

    (void)gomp_GetTrajectoryDisplayParams(&first_frame ,
                                        &last_frame  ,
                                        &delta_frame);

    in_set  = 0;
    in_set1 = 0;

    trace_atoms = gomp_GetTraceAtomsInSet();
    trace_list  = gomp_GetTraceAtomList();
    trcx        = gomp_GetTraceAtomXCoord();
    trcy        = gomp_GetTraceAtomYCoord();
    trcz        = gomp_GetTraceAtomZCoord();
    red         = gomp_GetAtomColourRedPointer(Wstr);
    green       = gomp_GetAtomColourGreenPointer(Wstr);
    blue        = gomp_GetAtomColourBluePointer(Wstr);

    glDisable(GL_LIGHTING);

    glLineWidth((GLfloat)(gomp_GetMoleculeLineWidth()+0.0));

    for(k = 0 ; k < gomp_GetTraceSets() ; k++) {

        for(j = 0 ; j < trace_atoms[k] ; j++) {

            glBegin(GL_LINE_STRIP);
 
            loopI = 0;
            for(i = (first_frame - 1) ; 
                i <  last_frame       ; 
                i += delta_frame) {

                jj = in_set + loopI * trace_atoms[k];

                loopI++;

                vec[0] = trcx[j + jj]; 
                vec[1] = trcy[j + jj]; 
                vec[2] = trcz[j + jj]; 

                temp = trace_list[j + in_set1];
                glColor3f(red[temp] , green[temp] , blue[temp]);
                glVertex3fv(vec);
                loop++;

                if(loop > 250) {
                    glEnd();
                    glBegin(GL_LINE_STRIP);
                    loop = 0;
                }
            }

            glEnd();

        }
        in_set  += (last_frame - first_frame + 1) * trace_atoms[k];
        in_set1 += trace_atoms[k];
    }
#endif /* ENABLE_GRAPHICS */

    return(0);
}

