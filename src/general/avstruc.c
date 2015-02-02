/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2003 - 2005 by:
Eero Häkkinen

*/

#include "maindefs.h"

#include "gomstdio.h"
#include "gommath.h"
#include <stdlib.h>

#include "cell.h"
#include "memalloc.h"
#include "molecule.h"
#include "molecoord.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

/* ................... */

/* The original structure stored in x,y and z will be replaced by the
   average structure. The average is for the FIRST atom set */
/**************************************************************************/
int gomp_TrajAvStructure()  
    /* average structure from the trajectory */
/**************************************************************************/
{
    static int    i,j,atom_max;
    static float *avx,*avy,*avz,fstep;
    static char   OutText[BUFF_LEN];
    static float  cellaa,cellbb,cellcc;
    static float  cella,cellb,cellc;
    static float  work1,work2,work3;
    static float  temp1,temp2,temp3;
    static FILE  *File_p;
    static float *x;
    static float *y;
    static float *z;
    static const float *sumxyz;
    static int    Wstr;
    static int    FirstF;
    static int    LastF;
    static int    StepF;

/* check that the info in the trajectory file matches the current ones */
    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintERROR("Number of frames is not defined ");
        return(1);
    }

    (void)gomp_GetTrajectoryDisplayParams(&FirstF ,
                                        &LastF  ,
                                        &StepF);

    if(*gomp_GetTrajectoryFileName() == '\0') {
        gomp_PrintERROR("?ERROR - trajectory file name is unknown");
        return(1);
    }

    gomp_PrintMessage("Calculating average structure");
    sprintf(OutText,"First frame: %d, Last frame: %d, Step: %d",
            FirstF,LastF,StepF);
    gomp_PrintMessage(OutText);

    /* open trajectory file */

    if(gomp_GetTrajectoryFileType() == XMOL_TRAJ      ||
       gomp_GetTrajectoryFileType() == GROMOS96A_TRAJ ||
       gomp_GetTrajectoryFileType() == TINKER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FAMBER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FDL_POLY_TRAJ) {
        File_p = fopen(gomp_GetTrajectoryFileName(),"r");
    } else {
        File_p = fopen(gomp_GetTrajectoryFileName(),"rb");
    }
    if(File_p == NULL) {
        sprintf(OutText,">>> Can't open input file: %s",
                gomp_GetTrajectoryFileName());
        gomp_PrintERROR(OutText);
        return(1); }

    Wstr     = 0;
    x        = gomp_GetModifiableAtomXCoordPointer(Wstr);
    y       = gomp_GetModifiableAtomYCoordPointer(Wstr);
    z      = gomp_GetModifiableAtomZCoordPointer(Wstr);
    sumxyz   = gomp_GetTranslateArray();

/* get some gomp_ratch space */
    atom_max = gomp_GetTotalNumberOfAtoms();

    avx = gomp_AllocateFloatVector(atom_max);
    avy = gomp_AllocateFloatVector(atom_max);
    avz = gomp_AllocateFloatVector(atom_max);

/* cell dimensions */
    cella = gomp_GetCellA();
    cellb = gomp_GetCellB();
    cellc = gomp_GetCellC();

    cellaa = (float)(1./gomp_GetCellA());
    cellbb = (float)(1./gomp_GetCellB());
    cellcc = (float)(1./gomp_GetCellC());

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
        avx[i] = avy[i] = avz[i] = 0.0;
    }

/* loop through the conformations */

    fstep = 0.0;

    for(i = FirstF ; i <= LastF ; i += StepF) {  /* frames         */

        if(gomp_GetOneFrame(i , File_p , TRAJ_OLD)) {
            gomp_PrintMessage("?ERROR - can't read trajectory frame");
            fclose(File_p);
            free(avx); /* free the gomp_ratch space */
            free(avy);
            free(avz);
            return(1);
        }

        rewind(File_p);

        fstep += 1.0;

        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(Wstr) ; j++) {  /* atoms */

            work1 = x[j];
            work2 = y[j];
            work3 = z[j];

            if(j > 0) {
                temp1 = x[j] - avx[j - 1]/((float)j);
                work1 = x[j] - cella * nearbyint(cellaa * temp1);
                temp2 = y[j] - avy[j - 1]/((float)j);
                work2 = y[j] - cellb * nearbyint(cellbb * temp2);
                temp3 = z[j] - avz[j - 1]/((float)j);
                work3 = z[j] - cellc * nearbyint(cellcc * temp3);
            }

            avx[j] += work1;
            avy[j] += work2;
            avz[j] += work3;
        }
    }

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {   
        /* put average structure in x,y and z */
        x[i] = avx[i] / fstep /* - sumxyz[0]*/;
        y[i] = avy[i] / fstep /* - sumxyz[1]*/;
        z[i] = avz[i] / fstep /* - sumxyz[2]*/;
    }

    free(avx); /* free the gomp_ratch space */
    free(avy);
    free(avz);

    gomp_PrintMessage(".... Average structure calculated ....");

    fclose(File_p);

    return(0);
}

