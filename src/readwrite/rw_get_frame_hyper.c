/*

Copyright (c) 1997 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gomendian.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))

#define FIXED_SIZED_HEADER_LEN      43
#define VARIABLE_SIZED_HEADER_LEN   17
#define COORD_AND_VELOC_DATA_LEN     6
#define DYNA_INFO_DATA_LEN           8

static int ReadHyperChemTrajInfo(FILE *HyperChem_oc);

static int swap_bytes;

/***************************************************************************/
int gomp_GetFrameHyperChem(int alt,FILE *File_p, int iappend)  
    /* read one frame from a HyperChem trajectory   
       mode of operation (=0) first time in read
       (>0) trajectory number 
       if = 0 no append , if = 1 append */
/***************************************************************************/
{
    static float *x,*y,*z;

    static long record_len=0;         /* the "record" length in bytes of one
                                         record containing the x,y and z
                                         coordinates plus the information
                                         in between the coordinates */
    static long ret_fseek;

    static int    natom;             /* number of HyperChem atoms */
    static int    icount;
    static int    i;
    static int    FixedSizedHeader[FIXED_SIZED_HEADER_LEN];
    static float  CoordAndVelocData[COORD_AND_VELOC_DATA_LEN];
    static const float *sumxyz;
    static int    Wstr;
    static char   OutText[BUFF_LEN];


    if(!alt) {
        return(ReadHyperChemTrajInfo(File_p));
    }

/* always from first atom list */

/* determine the "record length"    */
    natom = gomp_GetNumAtomsInMolecStruct(0);
    record_len = natom * sizeof(float) * 6 + 43 * sizeof(float);

/*  start reading  */

    icount = fread(FixedSizedHeader,sizeof(int), 
                   FIXED_SIZED_HEADER_LEN ,File_p);
    if ( swap_bytes ) {
        gomp_Reverse_int_array( FixedSizedHeader , FIXED_SIZED_HEADER_LEN);
    }

    if(FixedSizedHeader[2] != natom) {
        gomp_PrintMessage("?ERROR - number of atoms in traj file does not match defined atoms");
        return(1);
    }

    ret_fseek = fseek(File_p,(FixedSizedHeader[1] + (alt - 1) * record_len),
                      SEEK_SET);
    if(ret_fseek) {
        sprintf(OutText,"?ERROR - can't read trajectory file : %s \n",
                gomp_GetTrajectoryFileName());
        gomp_PrintERROR(OutText);
        return(1);
    }

    sumxyz   = gomp_GetTranslateArray();

/* next records are specific for one frame */

    if(iappend == 0) {

/*  get pointer to coordinate vectors */
        Wstr  = 0;
        x     = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y     = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z     = gomp_GetModifiableAtomZCoordPointer(Wstr);
      
/* read controll record */
        for(i = 0 ; i < natom ; i++) {

            icount = fread(CoordAndVelocData,sizeof(float), 
                           COORD_AND_VELOC_DATA_LEN ,File_p);
            if (swap_bytes ) {
                gomp_Reverse_float_array( CoordAndVelocData, COORD_AND_VELOC_DATA_LEN );
            }

            x[i] = CoordAndVelocData[0]     - sumxyz[0];
            y[i] = CoordAndVelocData[1]    - sumxyz[1];
            z[i] = CoordAndVelocData[2]   - sumxyz[2];
        }

    }
    else {

/*  get pointer to coordinate vectors */
        Wstr = gomp_CreateMolecStruct("HyperChem frame" , natom , APPEND);
        if ( Wstr < 0 )
            return(1);
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y    = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z    = gomp_GetModifiableAtomZCoordPointer(Wstr);

        for(i = 0 ; i < natom ; i++) {

            icount = fread(CoordAndVelocData,sizeof(float),
                           COORD_AND_VELOC_DATA_LEN ,File_p);
            if (swap_bytes ) {
                gomp_Reverse_float_array( CoordAndVelocData, COORD_AND_VELOC_DATA_LEN );
            }

            x[i] = CoordAndVelocData[0];
            y[i] = CoordAndVelocData[1];
            z[i] = CoordAndVelocData[2];
        }

    }

    return(0);
/*                                  */

}
/***************************************************************************/
int ReadHyperChemTrajInfo(FILE *HyperChem_oc)
/***************************************************************************/
{
    static int natom;       /* number of HyperChem atoms */
    static int icount;
    static int record;
    static int nstep;
    static int FileOffset;
    static int SnapshotPeriod;
    static int OptimizerType;
    static int UsingPeriod;
    static int ConstantTemp;
    static float TimeStep;
    static float HeatTime;
    static float RunTime;
    static float CoolTime;
    static float TempStep;
    static int   DataCollPeriod;
    static char OutText[BUFF_LEN];


/*  start reading  */

    icount = fread(&nstep,sizeof(int), 1 ,HyperChem_oc);
    icount = fread(&FileOffset,sizeof(int), 1 ,HyperChem_oc);
    icount = fread(&natom,sizeof(int), 1 ,HyperChem_oc);

    if(natom != gomp_GetNumAtomsInMolecStruct(0)) {
        gomp_Reverse_int( & natom );
        if(natom != gomp_GetNumAtomsInMolecStruct(0)) {
            gomp_PrintMessage("?ERROR - number of atoms in traj file does not match defined atoms");
            return(1);
        }
        gomp_Reverse_int( & nstep );
        gomp_Reverse_int( & FileOffset );
        gomp_PrintMessage("Enabling automatic byte_swapping...");
        swap_bytes = 1;
    }

    icount = fread(&SnapshotPeriod,sizeof(int), 1 ,HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_int( &SnapshotPeriod );
    }
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&OptimizerType,sizeof(int), 1 , HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_int( &OptimizerType );
    }
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&UsingPeriod,sizeof(int), 1 , HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_int( &UsingPeriod );
    }
    icount = fread(&ConstantTemp,sizeof(int), 1 , HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_int( &ConstantTemp );
    }
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&HeatTime,sizeof(float), 1 , HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_float( &HeatTime );
    }
    icount = fread(&RunTime,sizeof(float), 1 , HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_float( &RunTime );
    }
    icount = fread(&CoolTime,sizeof(float), 1 , HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_float( &CoolTime );
    }
    icount = fread(&TimeStep,sizeof(float), 1 , HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_float( &TimeStep );
    }
    icount = fread(&TempStep,sizeof(float), 1 , HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_float( &TempStep );
    }
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);
    icount = fread(&DataCollPeriod,sizeof(int), 1 , HyperChem_oc);
    if ( swap_bytes ) {
        gomp_Reverse_int( &DataCollPeriod );
    }
    icount = fread(&record,sizeof(int), 1 , HyperChem_oc);


/*
  trajectory_info.time_bw_steps    = 
  DataCollPeriod * (int)(TimeStep * 1.e+3);
*/
/* update trajectory info ... */
    (void)gomp_SetNumberOfTrajectoryAtoms(natom);
    (void)gomp_SetNumberOfFrames(nstep);
    (void)gomp_SetTrajectoryTimeInfo(DataCollPeriod * (int)(TimeStep * 1.e+3) , 0);
    (void)gomp_SetNumberOfFreeAtoms(0 , natom);
    (void)gomp_SetTrajectoryDisplayParams(1 , nstep , 1);
    (void)gomp_PutDisplayFrameNumber(1);

    sprintf(OutText," Info for trajectory file     : %s",gomp_GetTrajectoryFileName());
    gomp_PrintMessage(OutText);
    sprintf(OutText," Atoms found                  : %d",natom);
    gomp_PrintMessage(OutText);
    sprintf(OutText," Dynamics steps               : %d",nstep);
    gomp_PrintMessage(OutText);
    sprintf(OutText," Heat time (ps)               : %f",HeatTime);
    gomp_PrintMessage(OutText);
    sprintf(OutText," Run time (ps)                : %f",RunTime);
    gomp_PrintMessage(OutText);
    sprintf(OutText," Cool time (ps)               : %f",CoolTime);
    gomp_PrintMessage(OutText);
    sprintf(OutText," Time step (ps)               : %f",TimeStep);
    gomp_PrintMessage(OutText);

    fclose(HyperChem_oc);
    return(0);

/*                                  */

}


