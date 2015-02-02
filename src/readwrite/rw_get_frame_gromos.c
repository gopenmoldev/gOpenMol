/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "memalloc.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))

/***************************************************************************/
int gomp_GetFrameGromos(int alt, FILE *File_p , int iappend)  
    /* read one frame from a gromos trajectory   
       mode of operation (=0) first time in read
       (>0) trajectory number 
       if = 0 no append , if = 1 append */
/***************************************************************************/
{
    static float *x,*y,*z;
    static long   record_len=0;      /* the "record" length in bytes of one
                                        record containing the x,y and z
                                        coordinates plus the information
                                        in between the coordinates */
    static long   ret_fseek;

    static char   title[BUFF_LEN];  /* gromos title */
    static int    natom;       /* number of gromos atoms */
    static int    nstep;
    static int    icount;
    static int    record;
    static int    i;
    static int    idx1;
    static int    LabelRecord;
    static float *Master;
    static const float *sumxyz;
    static int    StartRecord;
    static int    Wstr;

    /* always from first atom list */

/* determine the "record length"    */
    natom = gomp_GetNumAtomsInMolecStruct(0);


/*  start reading  */

    if(!alt) {

        record_len = 2 * sizeof(int) +  natom * sizeof(float) * 3;

/* read controll record */
        icount = fread(&LabelRecord,sizeof(int), 1 , File_p);
        if(LabelRecord >= BUFF_LEN) {
            gomp_PrintMessage("?ERROR - buffer for label too small");
            return(-1);
        }
/* read LabelRecord characters   */
        icount = fread(title,sizeof(char), LabelRecord ,File_p);
        icount = fread(&record,sizeof(int), 1 , File_p);

        StartRecord = ftell(File_p);

        title[LabelRecord - 1] = '\0';
        gomp_PrintMessage("File title:");
        gomp_PrintMessage(title);
/* go first to the end of file */
        icount = fseek(File_p,0L,2);
        if(icount != 0) {
            gomp_PrintMessage("?ERROR - can't find end of trajectory file");
            return(-1);
        }

        record = ftell(File_p);
        icount = fseek(File_p, (long)(LabelRecord + 8) ,0);
        icount = ftell(File_p);

        nstep = (record-icount)/record_len;

        sprintf(title," Atoms found                : %d",natom);
        gomp_PrintMessage(title);
        sprintf(title," Dynamics steps             : %d",nstep);
        gomp_PrintMessage(title);

        (void)gomp_SetNumberOfTrajectoryAtoms(natom);
        (void)gomp_SetNumberOfFrames(nstep);
        (void)gomp_SetTrajectoryTimeInfo(0 , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , 0);
        (void)gomp_SetTrajectoryDisplayParams(1 , nstep , 1);
        (void)gomp_PutDisplayFrameNumber(1);
        return(0);
    }
    else
        icount = fseek(File_p , StartRecord , SEEK_SET);

    if(alt > 0) {
        ret_fseek = fseek(File_p,((alt - 1) * record_len), SEEK_CUR);
        if(ret_fseek) {
            sprintf(title,"?ERROR - can't read trajectory file");
            gomp_PrintMessage(title);
            return(-1);
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
        y    = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z    = gomp_GetModifiableAtomZCoordPointer(Wstr);

/* read controll record */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(Master,sizeof(float), 3 * natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = 0;
        for(i = 0 ; i < natom ; i++) {
            x[i] = 10. * Master[icount]       - sumxyz[0];
            y[i] = 10. * Master[icount + 1]  - sumxyz[1];
            z[i] = 10. * Master[icount + 2] - sumxyz[2];
            icount += 3;
        }
    }
    else {

        idx1 = 0;
/*  get pointer to coordinate vectors */
        sprintf(title,"Gromos frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(title , natom , APPEND);
        if ( Wstr < 0 )
            goto end;
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y    = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z    = gomp_GetModifiableAtomZCoordPointer(Wstr);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(Master,sizeof(float),3 * natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = 0;
        for(i = 0 ; i < natom ; i++) {
            x[idx1 + i] = 10. * Master[icount]       - sumxyz[0];
            y[idx1 + i] = 10. * Master[icount + 1]  - sumxyz[1];
            z[idx1 + i] = 10. * Master[icount + 2] - sumxyz[2];
            icount += 3;
        }

    }

end:
    free(Master);
    return(Wstr >= 0 ? 0 : 1);
/*                                  */

}

