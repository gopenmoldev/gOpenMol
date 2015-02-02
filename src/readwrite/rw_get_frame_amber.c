/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gomendian.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))

static int CheckIfBoxIncluded(FILE * , int);
static int BoxInformationIncluded;

/***************************************************************************/
int gomp_GetFrameAmber(int alt,FILE *File_p , int iappend)  
    /* read one frame from a amber trajectory   
       mode of operation (=0) first time in read
       (>0) trajectory number 
       if = 0 no append , if = 1 append */
/***************************************************************************/
{
    static long record_len=0;         /* the "record" length in bytes of one
                                         record containing the x,y and z
                                         coordinates plus the information
                                         in between the coordinates */
    static long ret_fseek;

    static char   title[BUFF_LEN];  /* amber title */
    static int    natom;       /* number of amber atoms */
    static int    nstep;
    static int    icount;
    static int    record;
    static int    i;
    static int    idx1;
    static float *Master;
    static float  SimBox[3];
    static float *x,*y,*z;
    static const float *sumxyz;
    static int    StartRecord;
    static int    Wstr;
    static int    swap_bytes;

/* just to be sure ... */
    rewind(File_p);
/* determine the "record length"    */
    natom = gomp_GetNumAtomsInMolecStruct(0);

/*  start reading  */

    if(!alt) {
/* read controll record */
        icount = fread(&record,sizeof(int), 1 , File_p);
/* read 80 characters   */
        swap_bytes = 0;
        if(record != (80 * sizeof(char))) {
            gomp_Reverse_int( & record );
            if(record != (80 * sizeof(char))) {
                gomp_PrintERROR("wrong internal structure of trajectory file");
                return(1);
            } else {
                gomp_PrintMessage("Enabling automatic byte_swapping...");
                swap_bytes = 1;
            }
        }

        icount = fread(title,sizeof(char), 80 ,File_p);
        icount = fread(&record,sizeof(int), 1 , File_p);
        gomp_PrintMessage("Title:");
        title[79] = '\0';
        gomp_PrintMessage(title);
/* mark the starting point for reading a frame */
        StartRecord = ftell(File_p);

        BoxInformationIncluded = CheckIfBoxIncluded(File_p , swap_bytes);

/* update trajectory info ... */
        icount = fseek(File_p,0L, SEEK_END);
        if(icount != 0) {
            gomp_PrintMessage("?ERROR - can't find end of trajectory file");
            return(0);
        }

        record = ftell(File_p);
        icount = fseek(File_p, 88L ,SEEK_SET);
        icount = ftell(File_p);

        if(BoxInformationIncluded)
            record_len = 4 * sizeof(int) +  natom * sizeof(float) * 3 +
                sizeof(float) * 3;
        else
            record_len = 2 * sizeof(int) +  natom * sizeof(float) * 3;

        nstep = (record-icount)/record_len;

        (void)gomp_SetNumberOfTrajectoryAtoms(natom);
        (void)gomp_SetNumberOfFrames(nstep);
        (void)gomp_SetTrajectoryTimeInfo(0 , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , 0);
        (void)gomp_SetTrajectoryDisplayParams(1 , nstep , 1);
        (void)gomp_PutDisplayFrameNumber(1);

        sprintf(title," Atoms found                : %d   ",natom);
        gomp_PrintMessage(title);
        sprintf(title," Free atoms                 : %d   ",natom);
        gomp_PrintMessage(title);
        sprintf(title," Dynamics steps             : %d   ",nstep);
        gomp_PrintMessage(title);
      
        return(0);
    }
    else 
        icount = fseek(File_p, StartRecord ,SEEK_SET);

    if(alt > 0) {
        ret_fseek = fseek(File_p,((alt - 1) * record_len),SEEK_CUR);
        if(ret_fseek) {
            sprintf(title,"?ERROR - can't read trajectory file");
            gomp_PrintMessage(title);
            return(1);
        }
    }

    sumxyz   = gomp_GetTranslateArray();

/* next records are specific for one frame */
/* read controll record */

    Master = gomp_AllocateFloatVector(3 * natom);

    if(iappend == 0) {

/*  get pointer to coordinate vectors */ 
        Wstr = 0;
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y   = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z  = gomp_GetModifiableAtomZCoordPointer(Wstr);

/* read controll record */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(Master,sizeof(float), 3 * natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_float_array( Master, 3 * natom );
        }


        if(BoxInformationIncluded) {
            icount = fread(&record,sizeof(int),   1 ,File_p);
            icount = fread(SimBox ,sizeof(float), 3 ,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);
            if ( swap_bytes ) {
                gomp_Reverse_float_array( SimBox, 3 );
            }

        }

        icount = 0;
        for(i = 0 ; i < natom ; i++) {
            x[i] = Master[icount]       - sumxyz[0];
            y[i] = Master[icount + 1]  - sumxyz[1];
            z[i] = Master[icount + 2] - sumxyz[2];
            icount += 3;
        }
    }
    else {

        idx1 = 0;
        sprintf(title,"Amber frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(title , natom , APPEND);
        if ( Wstr < 0 )
            goto end;
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y   = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z  = gomp_GetModifiableAtomZCoordPointer(Wstr);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(Master,sizeof(float),3 * natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_float_array( Master, 3 * natom );
        }

        icount = 0;
        for(i = 0 ; i < natom ; i++) {
            x[idx1 + i] = Master[icount]       - sumxyz[0];
            y[idx1 + i] = Master[icount + 1]  - sumxyz[1];
            z[idx1 + i] = Master[icount + 2] - sumxyz[2];
            icount += 3;
        }

    }

end:
    free(Master);
    return(Wstr>=0 ? 0 : 1);
/*                                  */
}

/*
  this function returns:

  (1) on ERROR it returns a negative value
  (2) if box information is included it returns a positive value
  (3) if box information is not included it returns zero (0)
*/
/***************************************************************************/
int CheckIfBoxIncluded(FILE *file_pointer , int swap_bytes)
/***************************************************************************/
{
    int icount1;
    int icount2;
    int record1;
    int record2;

/* place pointer to record after after the title */
    if(fseek(file_pointer , (long)88 , SEEK_SET)) return(-1);

/* number of bytes in next record */
    icount1 = fread(&record1, (size_t)(sizeof(int)), (size_t)1 , file_pointer);
    if ( swap_bytes ) {
        gomp_Reverse_int( & record1 );
    }

/* jump to next record */
    if(fseek(file_pointer , (long)record1 , SEEK_CUR)) return(1);

/* check end of record */
    icount2 = fread(&record2, (size_t)(sizeof(int)), (size_t)1 , file_pointer);
    if ( swap_bytes ) {
        gomp_Reverse_int( & record2 );
    }
    
    if(record1 != record2) return(-1);

/* read new record block */
    icount2 = fread(&record2, (size_t)(sizeof(int)), (size_t)1 , file_pointer);
    if ( swap_bytes ) {
        gomp_Reverse_int( & record2 );
    }

    if(record1 == record2) 
        return(0);   /* no box information included */
    else
        return(1);   /* box information is included */
}
