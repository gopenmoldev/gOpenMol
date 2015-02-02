/*
  Function to read the MUMOD trajectory file

  Copyright (c) 1991 - 2005 by:
  Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
  Confidential unpublished property of 
  Leif Laaksonen
  All rights reserved

  Enhancements 2003, 2005 by:
  Eero Häkkinen

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))

/***************************************************************************/
int gomp_GetFrameMumod(int alt, FILE *File_p, int iappend)  
    /* read one frame from a mumod trajectory 
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

    static int    icount;
    static int    record;
    static float  tstep;       /* time step in 10**(-12) sec */
    static int    ifre;        /* sampling frequency         */
    static int    natom;       /* number of mumod atoms      */
    static int    npts7;       /* number of frames in file   */
    static float  x2;          /* box size in x, y and z     */
    static float  y2;
    static float  z2;
    static int    i , idx1;
    static float *x,*y,*z;
    static const float *sumxyz;
    static char   OutText[BUFF_LEN];
    static int    StartRecord;
    static int    nstep;
    static int    Wstr;

/*  start reading  */

    rewind(File_p);

    if(!alt) {

/* read controll record */
        icount  = fread(&record,sizeof(int), 1 , File_p) * sizeof(int);
/* read ...             */
        icount = fread(&tstep,sizeof(float), 1 ,File_p) * sizeof(float);
        icount = fread(&ifre ,sizeof(int)  , 1 ,File_p) * sizeof(int);
        icount = fread(&natom,sizeof(int)  , 1 ,File_p) * sizeof(int);
        icount = fread(&npts7,sizeof(int)  , 1 ,File_p) * sizeof(int);
#ifdef W64BITS
        icount = fread(&x2   ,sizeof(double), 1 ,File_p) * sizeof(double);
        icount = fread(&y2   ,sizeof(double), 1 ,File_p) * sizeof(double);
        icount = fread(&z2   ,sizeof(double), 1 ,File_p) * sizeof(double);
#else
        icount = fread(&x2   ,sizeof(float), 1 ,File_p) * sizeof(float);
        icount = fread(&y2   ,sizeof(float), 1 ,File_p) * sizeof(float);
        icount = fread(&z2   ,sizeof(float), 1 ,File_p) * sizeof(float);
#endif

/* read controll record (icount gives now number of bytes from start) */
        icount = fread(&record,sizeof(int), 1 , File_p) * sizeof(int);

        record = ftell(File_p);;

/* header information done   */
/* determine the "record length"    */
        if(natom != gomp_GetNumAtomsInMolecStruct(0)) {
            gomp_PrintMessage(
                "?ERROR - number of atoms in file != current number of atoms");
            return(1);
        }
        record_len = 
            natom * sizeof(float) * 3 + 2 * sizeof(int);

        StartRecord = ftell(File_p);

/* go first to the end of file */
        icount = fseek(File_p,0L,SEEK_END);
        if(icount != 0) {
            gomp_PrintMessage("?ERROR - can't find end of trajectory file");
            return(1);
        }

        icount = ftell(File_p) - StartRecord;

        nstep = icount / record_len;
    
        sprintf(OutText," Atoms found                : %d ",natom);
        gomp_PrintMessage(OutText);
        sprintf(OutText," Dynamics steps             : %d ",nstep);
        gomp_PrintMessage(OutText);
        sprintf(OutText," Time between data sets     : %d ",
                ((int)(tstep * 1000.)*ifre)/10);
        gomp_PrintMessage(OutText);

/* update trajectory info ... */
        (void)gomp_SetNumberOfTrajectoryAtoms(natom);
        (void)gomp_SetNumberOfFrames(nstep);
        (void)gomp_SetTrajectoryTimeInfo(
            ((int)(tstep * 1000.)*ifre)/10 , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , natom);
        (void)gomp_SetTrajectoryDisplayParams(1 , nstep , 1);
        (void)gomp_PutDisplayFrameNumber(1);
    
        return(0);
    }
/* check things and go ahead */
/* read controll record */

/* so far icount bytes ... */
    if(alt > 0) {
        ret_fseek = fseek(File_p,((alt - 1) * record_len + StartRecord),1);
        if(ret_fseek) {
            printf("?ERROR - can't read trajectory file");
            return(1);
        }
    }

    sumxyz   = gomp_GetTranslateArray();
   
/* next records are specific for one frame */

    if(!iappend) {

/*  get pointer to coordinate vectors */
        Wstr = 0;
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y    = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z    = gomp_GetModifiableAtomZCoordPointer(Wstr);

/* read controll record */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        for(i = 0 ; i < natom ; i++) {
            icount = fread(&x[i],sizeof(float), 1 ,File_p);
            icount = fread(&y[i],sizeof(float), 1 ,File_p);
            icount = fread(&z[i],sizeof(float), 1 ,File_p);
        }
        icount = fread(&record,sizeof(int), 1 , File_p);

    }
    else {

        idx1 = 0;
/*  get pointer to coordinate vectors */
        sprintf(OutText,"Mumod frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(OutText , natom , APPEND);
        if ( Wstr < 0 )
            return(1);
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y    = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z    = gomp_GetModifiableAtomZCoordPointer(Wstr);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        for(i = idx1 ; i < (idx1 + natom) ;i++) {
            icount = fread(&x[i],sizeof(float), 1 ,File_p);
            icount = fread(&y[i],sizeof(float), 1 ,File_p);
            icount = fread(&z[i],sizeof(float), 1 ,File_p);
        }
        icount = fread(&record,sizeof(int), 1 ,File_p);

    }

/*    update_mlist(natom); */

    return(0);
/*                                  */

}



