/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2004 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <ctype.h>

#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "trajectory.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "selection.h"
#include "trace.h"

#include "stdafx.h"

/* structure to hold the trace of atoms       */

static struct {
    int trace_on;        /* switch to indicate that a trace is saved (=1)   */
    int trace_sets;      /* number of trace sets                            */
    int *trace_atoms;    /* number of traced atoms in each set              */
    int *trace_list;     /* list of atoms to be traced                      */
    float *trcx;         /* array to contain the x coordinates of the trace */
    float *trcy;
    float *trcz;
} trace_info;

static gom_Plotter *TracePlotterHandle;

/* trace the movement of the atoms in the selection list.
   The data come from the trajectory file */
/************************************************************************/
int gomp_TraceAtoms(const char *text1,const char *text2, const char *text3, int append )

/*
  const char *text1;  sgment name
  const char *text2;  residue name 
  const char *text3;  atom name    
  int   append; append     
*/
/************************************************************************/
{
    static int *sel_list;
    static int    slong;
    static float *tx,*ty,*tz;
    static int    i,j,jj,loop,new_len,new_len1,atom_max;
    static int    first_frame;
    static int    last_frame;
    static int    step_frame;
    static float *x;
    static float *y;
    static float *z;
    static const float *sumxyz;
    static int    Wstr;
    static char   OutText[BUFF_LEN];
    static FILE  *File_p;
    static float *SumX;
    static float *SumY;
    static float *SumZ;
    static int    Total;

    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintERROR("Number of frames is not defined ");
        return(1);
    }

    (void)gomp_GetTrajectoryDisplayParams(&first_frame ,
                                        &last_frame  ,
                                        &step_frame);

    atom_max = gomp_GetNumAtomsInMolecStruct(Wstr);
    sel_list = gomp_AllocateIntVector(atom_max);

    slong = gomp_MakeSelectionList(Wstr,text1,text2,text3,sel_list);

    if(slong > 0) { 

        Wstr     = 0;
        x        = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y        = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z        = gomp_GetModifiableAtomZCoordPointer(Wstr);
        sumxyz   = gomp_GetTranslateArray();
     
/* check that the trajectory file is known */
        if(*gomp_GetTrajectoryFileName() == '\0') {
            gomp_PrintMessage("?ERROR - trajectory file is not defined \n");
            return(1);
        }

/* do it ... */
        sprintf(OutText,"Selection is segment: '%s', residue: '%s', atom: '%s'",
                text1,text2,text3);
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
            sprintf(OutText,"Can't open input trajectory file: %s",
                    gomp_GetTrajectoryFileName());
            gomp_PrintMessage(OutText);
            return(1);
        }

        if(append) { 
/* check if there is something to be appended to */
            if(trace_info.trace_sets == 0) { /* no there is nothing */
                append = 0;
                gomp_PrintWARNING("Can't append to an empty list ");
            }
            gomp_PrintMessage("Selection list will be appended to old list");
        }

/* get some gomp_ratch space and save the x,y and z coordinates */
        tx   = gomp_AllocateFloatVector(atom_max);
        ty  = gomp_AllocateFloatVector(atom_max);
        tz = gomp_AllocateFloatVector(atom_max);

/* save coordinates */
        for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
            tx[i] = x[i];
            ty[i] = y[i];
            tz[i] = z[i];
        }

        switch(append) {

        case 0:  /* new list */

/* prepare the trace info structure */
            if(trace_info.trace_on > 0 || trace_info.trace_atoms != (const int *)NULL) {
                free(trace_info.trcx);
                free(trace_info.trcy);
                free(trace_info.trcz);
                free(trace_info.trace_list);
                free(trace_info.trace_atoms);
                trace_info.trace_atoms = (int *)NULL;
            }

            trace_info.trace_atoms = gomp_AllocateIntVector(1);

/* get space to save the trace                               */
            trace_info.trace_list = gomp_AllocateIntVector(slong);

            for(i = 0 ; i < slong ; i++)
                trace_info.trace_list[i] = sel_list[i];

            jj = slong * ((last_frame - first_frame + 1) / step_frame);

            gomp_FormatMessage("!INFO - will grab '%d' bytes of memory now",
                    (int)(3*jj*sizeof(int)));

            trace_info.trcx   = gomp_AllocateFloatVector(jj);
            trace_info.trcy  = gomp_AllocateFloatVector(jj);
            trace_info.trcz = gomp_AllocateFloatVector(jj);

            SumX = gomp_AllocateFloatVector(slong);
            SumY = gomp_AllocateFloatVector(slong);
            SumZ = gomp_AllocateFloatVector(slong);

            for(j = 0 ; j < slong ; j++)
                SumX[j] = SumY[j] = SumZ[j] = 0.0;

/* loop through the conformations */

            loop = 0;
            for(i = first_frame ; i <= last_frame ; i += step_frame) {  /* frames   */
        
                if(gomp_GetOneFrame(i , File_p , TRAJ_OLD)) {
                    gomp_PrintMessage("?ERROR - can't read trajectory frame");
                    fclose(File_p);
                    return(1);
                }

                rewind(File_p);

                for(j = 0 ; j < slong ; j++) {  /* atoms */
                    jj = sel_list[j];

                    SumX[j] += x[jj];
                    SumY[j] += y[jj];
                    SumZ[j] += z[jj];

                    trace_info.trcx[loop]   = x[jj];
                    trace_info.trcy[loop]  = y[jj];
                    trace_info.trcz[loop] = z[jj]; 
                    loop++;
                }
            }
            Total    = (last_frame - first_frame + 1);

            for(j = 0 ; j < slong ; j++) {  /* atoms */
                SumX[j] = SumX[j]/(float)Total;
                SumY[j] = SumY[j]/(float)Total;
                SumZ[j] = SumZ[j]/(float)Total;
                sprintf(OutText,"Center of coordinates: x: '%f', y: '%f', z: '%f'",
                        SumX[j]+sumxyz[0],SumY[j]+sumxyz[1],SumZ[j]+sumxyz[2]);
                gomp_PrintMessage(OutText);
            }

            trace_info.trace_on       = 1; /* trace is ready to be displayed */
            trace_info.trace_atoms[0] = slong;
            trace_info.trace_sets     = 1;

            if( !TracePlotterHandle )
                TracePlotterHandle = gomp_RegisterPlotter(
                    gomp_PlotAtomTrace,NULL,
                    PLOTTER_NAME_TRACE,PLOTTER_ORDER_TRACE);

            break; /* end of new list */

        case 1: /* append to old list */

/* realloc further space            */

            new_len = 0;

            for( i = 0 ; i < trace_info.trace_sets ; i++)
                new_len += trace_info.trace_atoms[i];

            new_len1 = (new_len + slong) * ((last_frame - first_frame + 1) / step_frame);

            trace_info.trace_list =
                realloc(trace_info.trace_list,
                (new_len + slong) * sizeof(*trace_info.trace_list));
            trace_info.trace_atoms =
                realloc(trace_info.trace_atoms, 
                (trace_info.trace_sets + 1) * sizeof(*trace_info.trace_atoms));
            trace_info.trcx = 
                realloc(trace_info.trcx,
                new_len1 * sizeof(*trace_info.trcx));
            trace_info.trcy =
                realloc(trace_info.trcy,
                new_len1 * sizeof(*trace_info.trcy));
            trace_info.trcz =
                realloc(trace_info.trcz,
                new_len1 * sizeof(*trace_info.trcz));

            for(i = 0 ; i < slong ; i++)
                trace_info.trace_list[i + new_len] = sel_list[i];

/* loop through the conformations */
            SumX = gomp_AllocateFloatVector(slong);
            SumY = gomp_AllocateFloatVector(slong);
            SumZ = gomp_AllocateFloatVector(slong);

            for(j = 0 ; j < slong ; j++) {
                SumX[j] = SumY[j] = SumZ[j] = 0.0;
            }

            loop = new_len * (last_frame - first_frame + 1) / step_frame;
            for(i = first_frame ; i <= last_frame ; i += step_frame) {  /* frames  */
            
                if(gomp_GetOneFrame(i , File_p , TRAJ_OLD)) {
                    gomp_PrintMessage("?ERROR - can't read trajectory frame");
                    fclose(File_p);
                    return(1);
                }

                rewind(File_p);

                for(j = 0 ; j < slong ; j++) {  /* atoms */
                    jj = sel_list[j];

                    SumX[j] += x[jj];
                    SumY[j] += y[jj];
                    SumZ[j] += z[jj];

                    trace_info.trcx[loop]   = x[jj];
                    trace_info.trcy[loop]  = y[jj];
                    trace_info.trcz[loop] = z[jj]; 
                    loop++;
                }
            }
            Total    = (last_frame - first_frame + 1);

            for(j = 0 ; j < slong ; j++) {  /* atoms */
                SumX[j] = SumX[j]/(float)Total;
                SumY[j] = SumY[j]/(float)Total;
                SumZ[j] = SumZ[j]/(float)Total;
                sprintf(OutText,"Center of coordinates: x: '%f', y: '%f', z: '%f'",
                        SumX[j]+sumxyz[0],SumY[j]+sumxyz[1],SumZ[j]+sumxyz[2]);
                gomp_PrintMessage(OutText);
            }

            trace_info.trace_on = 1; /* trace is ready to be displayed */

            trace_info.trace_atoms[trace_info.trace_sets] = slong;

            trace_info.trace_sets++;

            if( !TracePlotterHandle )
                TracePlotterHandle = gomp_RegisterPlotter(
                    gomp_PlotAtomTrace,NULL,
                    PLOTTER_NAME_TRACE,PLOTTER_ORDER_TRACE);

            break; /* end of append list */
        }

/* put values back */
        for(i = 0 ; i <  gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
            x[i] = tx[i];
            y[i] = ty[i];
            z[i] = tz[i];
        }

        free(tx); /* free the gomp_ratch space */
        free(ty);
        free(tz);

        free(sel_list);

        free(SumX);
        free(SumY);
        free(SumZ);

        fclose(File_p);

        return(0);
    }

/* done */

    else {
/* prepare the trace info structure */
        trace_info.trace_sets = slong;
        gomp_PrintMessage("?WARNING - no atoms in the selection list \n");
        return(0);
    }
}

/************************************************************************/
int    gomp_GetTraceState()
/************************************************************************/
{
    return(trace_info.trace_on);
}
/************************************************************************/
int    gomp_GetTraceSets()
/************************************************************************/
{
    return(trace_info.trace_sets);
}
#if 0
/************************************************************************/
int    gomp_GetTraceStep()
/************************************************************************/
{
    return(trace_info.trace_step);
}
#endif
/************************************************************************/
const int *gomp_GetTraceAtomsInSet()
/************************************************************************/
{
    return(trace_info.trace_atoms);
}
/************************************************************************/
const int *gomp_GetTraceAtomList()
/************************************************************************/
{
    return(trace_info.trace_list);
}
/************************************************************************/
const float *gomp_GetTraceAtomXCoord()
/************************************************************************/
{
    return(trace_info.trcx);
}

/************************************************************************/
const float *gomp_GetTraceAtomYCoord()
/************************************************************************/
{
    return(trace_info.trcy);
}

/************************************************************************/
const float *gomp_GetTraceAtomZCoord()
/************************************************************************/
{
    return(trace_info.trcz);
}

/************************************************************************/
int    gomp_DeleteTrace()
/************************************************************************/
{
    if(trace_info.trace_on == 1) {
        free(trace_info.trcx);
        free(trace_info.trcy);
        free(trace_info.trcz);
        free(trace_info.trace_list);
        free(trace_info.trace_atoms);
    }
    trace_info.trace_atoms = (int *)NULL;
    trace_info.trace_on    = 0;
    trace_info.trace_sets  = 0;

    gomp_UnregisterPlotter(TracePlotterHandle);
    TracePlotterHandle = NULL;
    
    return(0);
}

/************************************************************************/
int gomp_WriteAtomTrace(int Type , FILE *File_point)
/************************************************************************/
{

    static int i,j,k,jj,loopI,in_set,in_set1;
    static int    first_frame;
    static int    last_frame;
    static int    step_frame;
    static int    Wstr;
    static const float *sumxyz;
    static int    Total;
    static int    RunningAtoms;
    static float  Size;

    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintERROR("Number of frames is not defined ");
        return(1);
    }


    (void)gomp_GetTrajectoryDisplayParams(&first_frame ,
                                        &last_frame  ,
                                        &step_frame);

    Total    = gomp_AtomsInTrace();
    sumxyz   = gomp_GetTranslateArray();

    Wstr    = 0;
    in_set  = 0;
    in_set1 = 0;

    switch(Type) {  /* write type #1 */
    case PROBESURF_INPUT:
        fprintf(File_point,"Default title from Trace Atom facility\n");
        fprintf(File_point,"Default title from Trace Atom facility\n");
        fprintf(File_point,"%d\n",Total);
        break;
    }

    RunningAtoms = 0;

    for(k = 0 ; k < trace_info.trace_sets ; k++) {

        for(j = 0 ; j < trace_info.trace_atoms[k] ; j++) {

            loopI = 0;

            for(i = (first_frame - 1) ; 
                i <  last_frame       ; 
                i += step_frame) {

                jj = in_set + loopI * trace_info.trace_atoms[k];

                loopI++;

                switch(Type) {  /* write type #2 */
                case PROBESURF_INPUT:
                    fprintf(File_point,"%d %s %s %s %f %f %f %f\n",
                            (RunningAtoms + 1),
                            "seg",
                            "res",
                            "atm",
                            trace_info.trcx[j + jj] +   sumxyz[0],
                            trace_info.trcy[j + jj] +   sumxyz[1], 
                            trace_info.trcz[j + jj] +   sumxyz[2],
                            gomp_GetAtomVdwRad(Wstr , (trace_info.trace_list[j + in_set1]))); 
                    break;
                }

            }

        }
        in_set  += (last_frame -first_frame + 1) *
            trace_info.trace_atoms[k];
        in_set1 += trace_info.trace_atoms[k];
    }

    switch(Type) {  /* write type  #3 */
    case PROBESURF_INPUT:
        Size = gomp_GetSizeOfSystem();
        fprintf(File_point,"%f %f %f %f %f %f \n", -Size , Size , 
                -Size , Size ,
                -Size , Size);
        fprintf(File_point,"%d %d %d\n",60,60,60);
        fprintf(File_point,"%f\n",1.8);
        break;
    }
    return(0);
}

/************************************************************************/
int gomp_AtomsInTrace()
/************************************************************************/
{

    static int i,j,k,loop,in_set,in_set1;
    static int    first_frame;
    static int    last_frame;
    static int    step_frame;


    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintERROR("Number of frames is not defined ");
        return(1);
    }


    (void)gomp_GetTrajectoryDisplayParams(&first_frame ,
                                        &last_frame  ,
                                        &step_frame);


    in_set  = 0;
    in_set1 = 0;

    loop = 0;

    for(k = 0 ; k < trace_info.trace_sets ; k++) {

        for(j = 0 ; j < trace_info.trace_atoms[k] ; j++) {

            for(i = (first_frame - 1) ; 
                i <  last_frame       ; 
                i += step_frame) {

                loop++;
            }
        }
    }

    return(loop);
}
