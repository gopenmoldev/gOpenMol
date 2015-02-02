/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved

Enhancements 2003, 2005 by:
Eero Häkkinen
  
Additions:
* Byte swap added by Robert Best (rbest@hydrogen.cem.uct.ac.za)
March 2000

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

#define CHARMM_TEXT_FIELD 80

/* Read a frame from the trajectory file. This function opens the file,
   reads the header information 'fseeks' the place and read the frame.
   after that the file is closed again. */

/*  Symmetric shape index data for a Crystal/Constant pressure calculation */
struct sym_shape {
    float xtlabc[6];
};       /* lower triangle is used */

static struct SimCryst {
    double xtlabc[6];
} ShapeIndex;

void gomp_Reverse(void *p,size_t size)
{
    char *low  = p;
    char *high = low + size - 1;
    const char * const mid  = low + size / 2;
    char  temp;

    while (low != mid) {
        temp  = *low;
        *low  = *high;
        *high = temp;
        low++;
        high--;
    }
}

void gomp_Reverse_array(void *p,size_t size,size_t count)
{
    size_t i;
    char *q = p;
    for ( i = 0 ; i < count ; i++ ) {
        gomp_Reverse(q,size);
        q += size;
    }
}

#define DECLARE_REVERSE(type) \
void gomp_Reverse_##type(type *p) \
{ \
    gomp_Reverse(p,sizeof(type)); \
} \
void gomp_Reverse_##type##_array(type array[], size_t size) \
{ \
    gomp_Reverse_array(array,sizeof(type),size); \
}

DECLARE_REVERSE(int)
DECLARE_REVERSE(float)
DECLARE_REVERSE(double)

/***************************************************************************/
int gomp_GetFrameCharmm(int alt, FILE *File_p, int iappend)  
    /* read frame number 'alt' from charmm trajectory   */
    /* mode of operation:
       ( = 0) Check file and print trajectory information
       ( > 0) get frame number trajectory number        */
    /* *File_p is file pointer to the trajectory file   */
    /* append, if = 0 no append , if = 1 append         */
/***************************************************************************/
{
/*  
    Program to read CHARMm binary files

    leif laaksonen  1989 modified for gOpenMol 1995

    1992-11-17. Added constant pressure/crystal information record. LUL

*/

/* 1 HDR , ICONTR  */
    static char hdr[4];      /* unknown header information */
    static int icntrl[20];/* contains information about the datasets held in file

                          (1) number of data sets in the file. Not necessary correct
                          (2) time of the first data set. Usually in femtoseconds
                          (3) time steps between data sets
                          (9) number of fixed atoms (NFIXED)   
                          (11) 1 for crystal/constant pressure calculation
                          0 otherwise
                          (20) version number (22 for CHARMm 22, 0 for
                          previous version */

/* 2 TITLE         */
    static int ntitl;      /* (title(i,j), i=1,10),j=1,ntitl)
                              double precision title(10,10) */

/* 3 NATOM          (NFREAT = NATOM - NFIXED) */
    static int natom;


    static int icount,i;
    static int swap_bytes = 0;
    static int record;
    static int nstep,ifbeg,ifstep;
    static long record_len1=0;        /* the "record" length in bytes of one
                                         record containing the x,y and z
                                         coordinates plus the information
                                         in between the coordinates (first 
                                         record)*/
    static long record_len2;          /* rest of the records if icntrl[8] > 0 */

    static long ret_fseek;
    static int idx1;
    static int nfreat; /* number of free atoms */

    static float *tax; /* temp vectors in case of fixed atoms */
    static float *tay;
    static float *taz;
    static const int *Free_p;
    static int    qcrys;
    static int    CrystLen;

    static char label[BUFF_LEN];
    static int  numset;  /* number of data sets */
    static int  StartRecord;
    static int  Wstr;

/* pointers to coordinates */
    float *Xcoord;
    float *Ycoord;
    float *Zcoord;
    const float *sumxyz;

    CrystLen = 0;

    sumxyz   = gomp_GetTranslateArray();

/*  start reading  */
  
    if(!alt) {

/* check for right type of trajectory
   record has to be = (4 chars + 20 * sizeof(int))!
*/
        icount = fread(&record,sizeof(int), 1 ,File_p);

        swap_bytes = 0;
        if(record != (4 + 20 * sizeof(int))) {
            gomp_Reverse_int( & record );
            if(record != (4 + 20 * sizeof(int))) {
                gomp_PrintERROR("wrong internal structure of trajectory file");
                return(1);
            } else {
                gomp_PrintMessage("Enabling automatic byte_swapping...");
                swap_bytes = 1;
            }
        }

        icount = fread(hdr,sizeof(char),    4 ,File_p);
        icount = fread(icntrl,sizeof(int), 20 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int_array( icntrl, 20 );
        }

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&ntitl,sizeof(int),  1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( & ntitl );
        }

        for(i = 0 ; i < ntitl ; i++) {
            icount = fread(label, CHARMM_TEXT_FIELD ,1,File_p);
            label[CHARMM_TEXT_FIELD - 1] = '\0';
            gomp_PrintMessage(label);
        }

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&natom,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( & natom );
        }

        icount = fread(&record,sizeof(int), 1 ,File_p);

        nfreat = natom;

        qcrys = 0;
        if(icntrl[10]) {
            qcrys = 1;
            CrystLen = 2 * sizeof(int) + 6 * sizeof(double);
        }

        if(icntrl[8] > 0) {  /* there are fixed atoms */
            nfreat = natom - icntrl[8];
            (void)gomp_SetFreeAtomListPointer( 0 , nfreat);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(gomp_GetModifiableFreeAtomListPointer(0),
                           sizeof(int), nfreat, File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);
            if ( swap_bytes ) {
                gomp_Reverse_int_array( gomp_GetModifiableFreeAtomListPointer(0), nfreat );
            }
        }

/* determine the "record length"    */
        record_len1 = 6 * sizeof(int) + /* the record contribution */
            3 * natom * sizeof(float); /* x-,y- and z-contribution */
/* record length if there are fixed atoms */
        record_len2 = 6 * sizeof(int) + /* the record contribution */
            3 * nfreat * sizeof(float); /* x-,y- and z-contribution */

        if(qcrys) {
            record_len1 += CrystLen;
            record_len2 += CrystLen;
        }

/*                                  */
        StartRecord = ftell(File_p);
        ret_fseek = fseek(File_p,0L,SEEK_END);
        if(ret_fseek) {
            sprintf(label,"?ERROR - can't read trajectory file : %s ",
                    gomp_GetTrajectoryFileName());
            gomp_PrintMessage(label);
            return(1);
        }

        numset  = (ftell(File_p) - StartRecord - record_len1)/record_len2 + 1;
        if(numset != icntrl[0]) {
            gomp_PrintERROR("Number of frames in file header != actual number of frames");
            gomp_PrintMessage("Will change number of frames to real number");
            icntrl[0] = numset;
        }
        nstep   = icntrl[0];
        numset  = nstep;
        ifbeg   = icntrl[1];
        ifstep  = icntrl[2];

        sprintf(label," Info for trajectory file   : %s   ",
                gomp_GetTrajectoryFileName());
        gomp_PrintMessage(label);
        if(icntrl[19]) {
            sprintf(label," CHARMM version             : %d   ",icntrl[19]);
            gomp_PrintMessage(label);
        }
        sprintf(label," Atoms found                : %d   ",natom);
        gomp_PrintMessage(label);
        sprintf(label," Free atoms                 : %d   ",nfreat);
        gomp_PrintMessage(label);
        sprintf(label," Dynamics steps             : %d   ",nstep);
        gomp_PrintMessage(label);
        sprintf(label," Time between data sets     : %d   ",icntrl[2]);
        gomp_PrintMessage(label);
        sprintf(label," Time of the first data set : %d   ",icntrl[1]);
        gomp_PrintMessage(label);

/* test that number of atoms in the file is the same as for the displayed
   molecule  */

        if(natom != gomp_GetNumAtomsInMolecStruct(0)) {
            sprintf(label,"**** ERROR number of atoms in the file does not match");
            gomp_PrintERROR(label);
            return(1);
        }

/* update trajectory info ... */
        (void)gomp_SetNumberOfTrajectoryAtoms(natom);
        (void)gomp_SetNumberOfFrames(nstep);
        (void)gomp_SetTrajectoryTimeInfo(icntrl[2] , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , nfreat);
        (void)gomp_SetTrajectoryDisplayParams(1 , nstep , 1);
        (void)gomp_PutDisplayFrameNumber(1);

        return(0);
    }
    else
        ret_fseek = fseek(File_p , StartRecord , SEEK_SET);

/* there are fixed atoms */
    if(icntrl[8]) {

/* get space for the temporary coordinate array */
        tax   = gomp_AllocateFloatVector(nfreat);
        tay  = gomp_AllocateFloatVector(nfreat);
        taz = gomp_AllocateFloatVector(nfreat);

        if(alt > 1) {
            ret_fseek = fseek(File_p,(record_len1 + (alt - 2) * record_len2),1);
            if(ret_fseek) {
                sprintf(label,"?ERROR - can't read trajectory file : %s ",
                        gomp_GetTrajectoryFileName());
                gomp_PrintMessage(label);
                return(1);
            }
        }
    }
/* there are no fixed atoms */
    else {
        if(alt) {
            ret_fseek = fseek(File_p,((alt - 1) * record_len1),1);
            if(ret_fseek) {
                sprintf(label,"?ERROR - can't read trajectory file : %s ",
                        gomp_GetTrajectoryFileName());
                gomp_PrintMessage(label);
                return(1);
            }
        }
    }

#ifdef DEBUG
    printf(" Retrieving frame number : %d \n",alt+1);
#endif

    if(!iappend) {  /* start NO append */

/*  get pointer to coordinate vectors */
        Wstr     = 0;

        Xcoord   = gomp_GetModifiableAtomXCoordPointer(Wstr);
        Ycoord  = gomp_GetModifiableAtomYCoordPointer(Wstr);
        Zcoord = gomp_GetModifiableAtomZCoordPointer(Wstr);

        if(icntrl[8] > 0 ) {

            if(qcrys) {
                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(ShapeIndex.xtlabc,sizeof(double),6,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
            }
            if (swap_bytes ) {
                gomp_Reverse_double_array( ShapeIndex.xtlabc, 6 );
            }

/* special case first frame contains ALL ATOMS */
            if(alt == 1) {
                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(Xcoord,sizeof(float),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(Ycoord,sizeof(float),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(Zcoord,sizeof(float),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
                if (swap_bytes ) {
                    gomp_Reverse_float_array( Xcoord, natom );
                    gomp_Reverse_float_array( Ycoord, natom );
                    gomp_Reverse_float_array( Zcoord, natom );
                }
            }
            else {
                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(float),nfreat,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(float),nfreat,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(float),nfreat,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
                if (swap_bytes ) {
                    gomp_Reverse_float_array( tax, nfreat );
                    gomp_Reverse_float_array( tay, nfreat );
                    gomp_Reverse_float_array( taz, nfreat );
                }
    
                for(i = 0 ; i < natom ; i++) {
                    Xcoord[i]   += sumxyz[0];
                    Ycoord[i]  += sumxyz[1];
                    Zcoord[i] += sumxyz[2];
                }

/* put new values in */
                Free_p = gomp_GetFreeAtomListPointer(0);

                for(i = 0 ; i < nfreat; i++) {
                    Xcoord[Free_p[i] - 1]   = tax[i];
                    Ycoord[Free_p[i] - 1]  = tay[i];
                    Zcoord[Free_p[i] - 1] = taz[i];
                }
            }
        }
        else {

            if(qcrys) {
                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(ShapeIndex.xtlabc,sizeof(double),6,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
                if (swap_bytes ) {
                    gomp_Reverse_double_array( ShapeIndex.xtlabc, 6 );
                }
            }

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(Xcoord,sizeof(float),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(Ycoord,sizeof(float),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(Zcoord,sizeof(float),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            if (swap_bytes ) {
                gomp_Reverse_float_array( Xcoord, natom );
                gomp_Reverse_float_array( Ycoord, natom );
                gomp_Reverse_float_array( Zcoord, natom );
            }
        }
    }  /*   end of NO append */
    else {   /*   start of append  */

        idx1   = 0;
/*  get pointer to coordinate vectors */
        Wstr   = gomp_CreateMolecStruct("CHARMm frame" , natom , APPEND);
        if ( Wstr < 0 )
            goto end;
        Xcoord = gomp_GetModifiableAtomXCoordPointer(Wstr);
        Ycoord = gomp_GetModifiableAtomYCoordPointer(Wstr);
        Zcoord = gomp_GetModifiableAtomZCoordPointer(Wstr);

        if(icntrl[8] > 0 ) {

            if(qcrys) {
                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(ShapeIndex.xtlabc,sizeof(double),6,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
                if ( swap_bytes ) {
                    gomp_Reverse_double_array( ShapeIndex.xtlabc, 6 );
                }
            }

/* special case first frame contains ALL ATOMS */
            if(alt == 1) {

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(Xcoord,sizeof(float),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(Ycoord,sizeof(float),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(Zcoord,sizeof(float),natom,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
                if (swap_bytes ) {
                    gomp_Reverse_float_array( Xcoord, natom );
                    gomp_Reverse_float_array( Ycoord, natom );
                    gomp_Reverse_float_array( Zcoord, natom );
                }
            }
            else {

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tax,sizeof(float),nfreat,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(tay,sizeof(float),nfreat,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);

                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(taz,sizeof(float),nfreat,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
                if (swap_bytes ) {
                    gomp_Reverse_float_array( tax, nfreat );
                    gomp_Reverse_float_array( tay, nfreat );
                    gomp_Reverse_float_array( taz, nfreat );
                }

                for(i = 0 ; i < natom ; i++) {
                    Xcoord[i + idx1]   = Xcoord[i] + sumxyz[0];
                    Ycoord[i + idx1]  = Ycoord[i] + sumxyz[1];
                    Zcoord[i + idx1] = Zcoord[i] + sumxyz[2] ;
                }

/* put new values in */
                Free_p = gomp_GetFreeAtomListPointer(Wstr);

                for(i = 0 ; i < nfreat; i++) {
                    Xcoord[idx1 + Free_p[i] - 1]   = tax[i];
                    Ycoord[idx1 + Free_p[i] - 1]  = tay[i];
                    Zcoord[idx1 + Free_p[i] - 1] = taz[i];
                }
            }
        }
        else {

            if(qcrys) {
                icount = fread(&record,sizeof(int), 1 ,File_p);
                icount = fread(ShapeIndex.xtlabc,sizeof(double),6,File_p);
                icount = fread(&record,sizeof(int), 1 ,File_p);
            }
            if ( swap_bytes ) {
                gomp_Reverse_double_array( ShapeIndex.xtlabc, 6);
            }

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(&Xcoord[idx1],sizeof(float),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(&Ycoord[idx1],sizeof(float),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(&Zcoord[idx1],sizeof(float),natom,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);

            if (swap_bytes ) {
                gomp_Reverse_float_array( Xcoord, natom );
                gomp_Reverse_float_array( Ycoord, natom );
                gomp_Reverse_float_array( Zcoord, natom );
            }
        }
    }   /* end of append */


/* shift back using the translation info */
    for(i = 0 ; i < natom ; i++) {
        Xcoord[i]   -= sumxyz[0];
        Ycoord[i]  -= sumxyz[1];
        Zcoord[i] -= sumxyz[2];
    }
/* end of translation                    */

end:
    if(icntrl[8]) {
        free(tax);
        free(tay);
        free(taz);
    }

    return(Wstr >=0 ? 0 : 1);
/*                                  */

}


