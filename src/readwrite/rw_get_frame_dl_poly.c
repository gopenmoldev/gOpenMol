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
#include <tcl.h>

#include "gomendian.h"
#include "gomstring.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))

#define DL_POLY_TEXT_FIELD 80

/* Read a frame from the trajectory file. This function opens the file,
   reads the header information 'fseeks' the place and read the frame.
   after that the file is closed again. */

static int PlaceDL_POLYFrame(FILE * , int , int);
static int Look4DL_POLYFrames(FILE * , int * , int * , int *);

static int *DL_PolyFrameStartPoint = NULL;

#if 0
/***************************************************************************/
int gomp_GetFrameDL_Poly(int alt, FILE *File_p, int iappend)  
    /* read frame number 'alt' from charmm trajectory   */
    /* mode of operation:
       ( = 0) Check file and print trajectory information
       ( > 0) get frame number trajectory number        */
    /* *File_p is file pointer to the trajectory file   */
    /* append, if = 0 no append , if = 1 append         */
/***************************************************************************/
{
/*  
    Program to read DL_Poly unformatted trajectory files

    leif laaksonen  2000

*/

    return(gomp_GetFrameDL_PolyUNFORMATTED(alt , File_p , iappend));
}
#endif
/***************************************************************************/
int gomp_GetFrameDL_PolyUNFORMATTED(int alt, FILE *File_p, int iappend)  
/***************************************************************************/
{
/* 1 Header  */
    static char hdr[80];      /* header information*/

/* 2 NATMS   */
    static double natms;

/*
  Jump over these records:

  record 3
  atname(1,...,natms)        atom names or symbols (character*8)
  record 4
  weight(1,...,natms)        atomic masses (real*8)
  record 5
  charge(1,...,natms)        atomic charges (real*8)

*/
    static double nstep;
    static double keytrj;
    static double imcon;
    static double tstep;
    static double cell[9];
    static long   record_len=0;

    static int natom;
    static double  *tax; /* temp vectors  */
    static double  *tay;
    static double  *taz;
    static float *tax1; /* temp vectors  */
    static float *tay1;
    static float *taz1;
    static int icount,i;
    static int swap_bytes = 0;
    static int record;
    static long ret_fseek;
    static int  StartRecord;
    static int  Wstr;
    static char label[BUFF_LEN];
    static int  numset;

/* pointers to coordinates */
    static float *Xcoord;
    static float *Ycoord;
    static float *Zcoord;
    static const float *sumxyz;

    sumxyz   = gomp_GetTranslateArray();

/*  start reading  */
  
    if(!alt) {

/* check for right type of trajectory
   record has to be = (80 chars)!
*/
/*#1*/
        icount = fread(&record,sizeof(int), 1 ,File_p);

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

        icount = fread(&hdr,sizeof(char),  80 ,File_p);
        hdr[DL_POLY_TEXT_FIELD - 1] = '\0';
        gomp_PrintMessage(hdr);
        icount = fread(&record,sizeof(int), 1 ,File_p);
/*#2*/
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&natms,sizeof(double),   1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_double( &natms );
        }
        natom = (int)natms;
        icount = fread(&record,sizeof(int), 1 ,File_p);
/*#3*/
        icount = fread(&record,sizeof(int), 1 ,File_p);
        ret_fseek = fseek(File_p,(8 * natom * sizeof(char)),SEEK_CUR);
        if(ret_fseek) {
            sprintf(label,"?ERROR - can't read trajectory file : %s (#1)",
                    gomp_GetTrajectoryFileName());
            gomp_PrintMessage(label);
            return(1);
        }
        icount = fread(&record,sizeof(int), 1 ,File_p);
/*#4*/
        icount = fread(&record,sizeof(int), 1 ,File_p);
        ret_fseek = fseek(File_p,(natom * sizeof(double)),SEEK_CUR);
        if(ret_fseek) {
            sprintf(label,"?ERROR - can't read trajectory file : %s (#2)",
                    gomp_GetTrajectoryFileName());
            gomp_PrintMessage(label);
            return(1);
        }
        icount = fread(&record,sizeof(int), 1 ,File_p);
/*#5*/
        icount = fread(&record,sizeof(int), 1 ,File_p);
        ret_fseek = fseek(File_p,(natom * sizeof(double)),SEEK_CUR);
        if(ret_fseek) {
            sprintf(label,"?ERROR - can't read trajectory file : %s (#3)",
                    gomp_GetTrajectoryFileName());
            gomp_PrintMessage(label);
            return(1);
        }
        icount = fread(&record,sizeof(int), 1 ,File_p);

        StartRecord = ftell(File_p);

/*#i*/
/*
  record i
  nstep          the current time-step (real*8)
  natms          number of atoms in configuration (real*8)
  keytrj         trajectory key (real*8)
  imcon          image convention key (real*8)
  tstep          integration timestep (real*8) 
*/
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&nstep,sizeof(double),   1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_double( &nstep );
        }
        icount = fread(&natms,sizeof(double),   1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_double( &natms );
        }
        icount = fread(&keytrj,sizeof(double),   1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_double( &keytrj );
        }
        icount = fread(&imcon,sizeof(double),   1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_double( &imcon );
        }
        icount = fread(&tstep,sizeof(double),   1 ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_double( &tstep );
        }
        icount = fread(&record,sizeof(int), 1 ,File_p);
/* determine the "record length"    */
        record_len = 6 * sizeof(int) + /* the record contribution */
            3 * natom * sizeof(double) + /* x-,y- and z-contribution */
            5 * sizeof(double) + 2 * sizeof(int);

        if(imcon > 0.0)  record_len += 9 * sizeof(double) + 2 * sizeof(int); /* cell contribution */

        if(keytrj > 0.0) {
            record_len += natom * sizeof(double) + 2 * sizeof(int); /* v x */
            record_len += natom * sizeof(double) + 2 * sizeof(int); /* v y */
            record_len += natom * sizeof(double) + 2 * sizeof(int); /* v z */
        }

        if(keytrj > 1.0) {
            record_len += natom * sizeof(double) + 2 * sizeof(int); /* f x */
            record_len += natom * sizeof(double) + 2 * sizeof(int); /* f y */
            record_len += natom * sizeof(double) + 2 * sizeof(int); /* f z */
        }
/*                                  */
        ret_fseek = fseek(File_p,0L,SEEK_END);
        if(ret_fseek) {
            sprintf(label,"?ERROR - can't read trajectory file : %s (#4)",
                    gomp_GetTrajectoryFileName());
            gomp_PrintMessage(label);
            return(1);
        }

        numset  = (ftell(File_p) - StartRecord)/record_len;

        sprintf(label," Info for trajectory file   : %s   ",
                gomp_GetTrajectoryFileName());
        sprintf(label," Atoms found                : %d   ",natom);
        gomp_PrintMessage(label);
        sprintf(label," Dynamics steps             : %d   ",numset);
        gomp_PrintMessage(label);
        if(keytrj > 0.0) {
            gomp_PrintMessage(" Velocities available       : YES");
        } else {
            gomp_PrintMessage(" Velocities available       : NO");
        }
        if(keytrj > 1.0) {
            gomp_PrintMessage(" Force components available : YES");
        } else {
            gomp_PrintMessage(" Force components available : NO");
        }      

/* test that number of atoms in the file is the same as for the displayed
   molecule  */

        if(natom != gomp_GetNumAtomsInMolecStruct(0)) {
            sprintf(label,"**** ERROR number of atoms in the file does not match");
            gomp_PrintERROR(label);
            return(1);
        }

/* update trajectory info ... */
        (void)gomp_SetNumberOfTrajectoryAtoms(natom);
        (void)gomp_SetNumberOfFrames(numset);
        (void)gomp_SetTrajectoryTimeInfo(0 , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , 0);
        (void)gomp_SetTrajectoryDisplayParams(1 , numset , 1);
        (void)gomp_PutDisplayFrameNumber(1);

        return(0);
    } else { 
        ret_fseek = fseek(File_p , StartRecord , SEEK_SET);
    }
    ret_fseek = fseek(File_p,((alt - 1) * record_len + 
                              5 * sizeof(double) + 2 * sizeof(int)),SEEK_CUR);
    if(ret_fseek) {
        sprintf(label,"?ERROR - can't read trajectory file : %s ",
                gomp_GetTrajectoryFileName());
        gomp_PrintMessage(label);
        return(1);
    }

#ifdef DEBUG
    printf(" Retrieving frame number : %d \n",alt+1);
#endif

    if(!iappend) {  /* start NO append */

/* get space for the temporary coordinate array */
        tax   = gomp_AllocateDoubleVector(natom);
        tay  = gomp_AllocateDoubleVector(natom);
        taz = gomp_AllocateDoubleVector(natom);

/*  get pointer to coordinate vectors */

        Xcoord   = gomp_GetModifiableAtomXCoordPointer(0);
        Ycoord  = gomp_GetModifiableAtomYCoordPointer(0);
        Zcoord = gomp_GetModifiableAtomZCoordPointer(0);

        if(imcon > 0.0 ) {
            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(cell,sizeof(double),9,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);
            if (swap_bytes ) {
                gomp_Reverse_double_array( cell, 9 );
            }
        }

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(tax,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(tay,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(taz,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        if (swap_bytes ) {
            gomp_Reverse_double_array( tax, natom );
            gomp_Reverse_double_array( tay, natom );
            gomp_Reverse_double_array( taz, natom );
        }

        for(i = 0 ; i < natom; i++) {
            Xcoord[i]   = (float)tax[i];
            Ycoord[i]  = (float)tay[i];
            Zcoord[i] = (float)taz[i];
        }
    }  /*   end of NO append */
    else {   /*   start of append  */

/*  get pointer to coordinate vectors */
        Wstr     = gomp_CreateMolecStruct("DL_Poly frame" , natom , APPEND);
        if ( Wstr < 0 )
            return(1);
        Xcoord   = gomp_GetModifiableAtomXCoordPointer(Wstr);
        Ycoord  = gomp_GetModifiableAtomYCoordPointer(Wstr);
        Zcoord = gomp_GetModifiableAtomZCoordPointer(Wstr);

/* get space for the temporary coordinate array */
        tax   = gomp_AllocateDoubleVector(natom);
        tay  = gomp_AllocateDoubleVector(natom);
        taz = gomp_AllocateDoubleVector(natom);

        if(imcon > 0.0 ) {
            icount = fread(&record,sizeof(int), 1 ,File_p);
            icount = fread(cell,sizeof(double),9,File_p);
            icount = fread(&record,sizeof(int), 1 ,File_p);
            if (swap_bytes ) {
                gomp_Reverse_double_array( cell, 9 );
            }
        }

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(tax,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(tay,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(taz,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        if (swap_bytes ) {
            gomp_Reverse_double_array( tax, natom );
            gomp_Reverse_double_array( tay, natom );
            gomp_Reverse_double_array( taz, natom );
        }

        for(i = 0 ; i < natom; i++) {
            Xcoord[i]   = (float)tax[i];
            Ycoord[i]  = (float)tay[i];
            Zcoord[i] = (float)taz[i];
        }
    }   /* end of append */


/* shift back using the translation info */
    for(i = 0 ; i < natom ; i++) {
        Xcoord[i]   -= sumxyz[0];
        Ycoord[i]  -= sumxyz[1];
        Zcoord[i] -= sumxyz[2];
    }
/* end of translation                    */

/* atom velocities */
/* retrieve velocities and forces if requested */
    if((keytrj > 0.0) && gomp_GetVelocityRetrieveState()) {

        if(gomp_GetVelocitySpace(natom))
            return(1);

        tax1 = gomp_GetModifiableVelocityXComponentPointer();
        tay1 = gomp_GetModifiableVelocityYComponentPointer();
        taz1 = gomp_GetModifiableVelocityZComponentPointer();

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(tax,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(tay,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(taz,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        if (swap_bytes ) {
            gomp_Reverse_double_array( tax, natom );
            gomp_Reverse_double_array( tay, natom );
            gomp_Reverse_double_array( taz, natom );
        }

        for(i = 0 ; i < natom; i++) {
            tax1[i] = (float)tax[i];
            tay1[i] = (float)tay[i];
            taz1[i] = (float)taz[i];
        }

        if(gomp_CalculateVelocityMinMax())
            gomp_PrintERROR("can't calculate velocity min/max values");
    }
/* atom forces     */
    if((keytrj > 1.0) && gomp_GetForceRetrieveState()) {

        if(gomp_GetForceSpace(natom))
            return(1);

        tax1 = gomp_GetModifiableForceXComponentPointer();
        tay1 = gomp_GetModifiableForceYComponentPointer();
        taz1 = gomp_GetModifiableForceZComponentPointer();

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(tax,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(tay,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(taz,sizeof(double),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        if (swap_bytes ) {
            gomp_Reverse_double_array( tax, natom );
            gomp_Reverse_double_array( tay, natom );
            gomp_Reverse_double_array( taz, natom );
        }
/* put new values in */

        for(i = 0 ; i < natom; i++) {
            tax1[i] = (float)tax[i];
            tay1[i] = (float)tay[i];
            taz1[i] = (float)taz[i];
        }

        if(gomp_CalculateForceMinMax())
            gomp_PrintERROR("can't calculate force min/max values");
    }

    free(tax);
    free(tay);
    free(taz);

    return(0);
/*                                  */

}
/***************************************************************************/
int gomp_GetFrameDL_PolyFORMATTED(int alt, FILE *File_p, int iappend)  
/***************************************************************************/
{
    static int    natom;
    static int    natomx;
    static int    nstep;
    static int    Wstr;
    static char   title[BUFF_LEN];
    static int    keytrj;

/* just to be sure ... */
    rewind(File_p);
/* determine the "record length"    */
    natomx = gomp_GetNumAtomsInMolecStruct(0);

/*  start reading  */

    if(!alt) {

        (void)Look4DL_POLYFrames(File_p , &natom , &keytrj , &nstep);

        if(natomx != natom) {
            gomp_PrintERROR("number of atoms in trajectory file do no match current system");
            return(1);
        }
        (void)gomp_SetNumberOfTrajectoryAtoms(natom);
        (void)gomp_SetNumberOfFrames(nstep);
        (void)gomp_SetTrajectoryTimeInfo(0 , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , 0);
        (void)gomp_SetTrajectoryDisplayParams(1 , nstep , 1);
        (void)gomp_PutDisplayFrameNumber(1);

        sprintf(title," Atoms found                : %d   ",natom);
        gomp_PrintMessage(title);
        sprintf(title," Dynamics steps             : %d   ",nstep);
        gomp_PrintMessage(title);
        if(keytrj > 0) {
            gomp_PrintMessage(" Velocities available       : YES");
        } else {
            gomp_PrintMessage(" Velocities available       : NO");
        }
        if(keytrj > 1) {
            gomp_PrintMessage(" Force components available : YES");
        } else {
            gomp_PrintMessage(" Force components available : NO");
        }      
      
        return(0);
    }

    if(iappend == 0)
/*  get pointer to coordinate vectors */ 
        Wstr = 0;
    else {
        sprintf(title,"DL_POLY frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(title , natom , APPEND);
        if ( Wstr < 0 )
            return(1);
    }
    (void)PlaceDL_POLYFrame(File_p , alt , Wstr);
    return(0);
}
#define DL_POLY_LINE_LEN   120   /*DL_POLY coord file line length */

/*********************************************************************/
int Look4DL_POLYFrames(FILE *chm_in , int *Atoms , int *keytrj , int *Frames)
/*********************************************************************/
{
    static int i;
    static int DL_POLYAtoms;
    static int DL_POLYFrames;

    char  inputl[DL_POLY_LINE_LEN];
    char  OutText[BUFF_LEN];
    int   items;
    int   imcon;
    int   natms;
    int   nstep;
    float tstep;


    rewind(chm_in);

    *Atoms  = 0;
    *Frames = 0;

    if(DL_PolyFrameStartPoint != NULL) {
        free(DL_PolyFrameStartPoint);
        DL_PolyFrameStartPoint  = NULL;
    }
/*
  Start reading file

*/
/*#1*/
    fgets(inputl,DL_POLY_LINE_LEN,chm_in);
/*#2*/
    fgets(inputl,DL_POLY_LINE_LEN,chm_in);
    items = sscanf(inputl,"%10d%10d%10d",keytrj,&imcon,&natms);

    DL_POLYFrames = 0;
    DL_POLYAtoms  = 0;

    while(!feof(chm_in)) {

/* save the pointers in the DL_POLY file */
        if(DL_POLYFrames < 1) {
            DL_PolyFrameStartPoint     = gomp_AllocateIntVector(1);
            DL_PolyFrameStartPoint[0]  = ftell(chm_in);
        } else {
            DL_PolyFrameStartPoint = 
                gomp_ReallocateIntVector(DL_PolyFrameStartPoint , DL_POLYFrames + 1);
            DL_PolyFrameStartPoint[DL_POLYFrames]  = ftell(chm_in);
        }

        if(fgets(inputl,DL_POLY_LINE_LEN,chm_in) == NULL) break;

        (void)gomp_StringTrim(inputl);

        if(!sscanf(inputl,"%8s%10d%10d%10d%10d%12f",
                   OutText,&nstep,&natms,keytrj,&imcon,&tstep)) break;

        if(imcon) {
            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
        }
/*
  record ii (3g12.4) for  imcon  0
  cell(1)        real        x component of  cell vector
  cell(2)        real        y component of  cell vector
  cell(3)        real        z component of  cell vector
  record iii (3g12.4) for  imcon  0
  cell(4)        real        x component of  cell vector
  cell(5)        real        y component of  cell vector
  cell(6)        real        z component of  cell vector
  record iv  (3g12.4) for  imcon  0
  cell(7)        real        x component of  cell vector
  cell(8)        real        y component of  cell vector
  cell(9)        real        z component of  cell vector
*/

        for(i = 0 ; i < natms ; i++) {
            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            if(*keytrj > 0) {
                fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            }
            if(*keytrj > 1) {
                fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            }
        }
        DL_POLYFrames++;
    }

    *Frames  = DL_POLYFrames;
    *Atoms   = natms;

    rewind(chm_in);

    return(0);
}

/*********************************************************************/
int PlaceDL_POLYFrame(FILE *chm_in , int Alt , int Wstr)
/*********************************************************************/
{
    static int i,j,k;
    static int DL_POLYFrames;

    float   TXc,TYc,TZc;
    const float *sumxyz;

    char inputl[DL_POLY_LINE_LEN];
    char OutText[BUFF_LEN];
    int  items;
    int  method;
    int   keytrj;
    int   imcon;
    int   natms;
    int   nstep;
    float tstep;
    float *tax1 = NULL;
    float *tay1 = NULL;
    float *taz1 = NULL;
    float *tax2 = NULL;
    float *tay2 = NULL;
    float *taz2 = NULL;

/*
  Start reading file

*/
    rewind(chm_in);
    sumxyz         = gomp_GetTranslateArray();
    DL_POLYFrames  = 0;

    method = gomp_GetFormattedTrajectoryReader();

    if(method) {
/* position -1 from the right one */
/*
  Start reading file

*/
/*#1*/
        fgets(inputl,DL_POLY_LINE_LEN,chm_in);
/*#2*/
        fgets(inputl,DL_POLY_LINE_LEN,chm_in);
        items = sscanf(inputl,"%10d%10d%10d",&keytrj,&imcon,&natms);

/* atom velocities */
        if((keytrj > 0) && gomp_GetVelocityRetrieveState()) {

            if(gomp_GetVelocitySpace(natms))
                return(1);

            tax1 = gomp_GetModifiableVelocityXComponentPointer();
            tay1 = gomp_GetModifiableVelocityYComponentPointer();
            taz1 = gomp_GetModifiableVelocityZComponentPointer();
        }
/* atom forces     */
        if((keytrj > 1) && gomp_GetForceRetrieveState()) {

            if(gomp_GetForceSpace(natms))
                return(1);

            tax2 = gomp_GetModifiableForceXComponentPointer();
            tay2 = gomp_GetModifiableForceYComponentPointer();
            taz2 = gomp_GetModifiableForceZComponentPointer();

        }

        for(i = 0 ; i < Alt ; i++) {

            if(fgets(inputl,DL_POLY_LINE_LEN,chm_in) == NULL) break;
            (void)gomp_StringTrim(inputl);
            if(!sscanf(inputl,"%8s%10d%10d%10d%10d%12f",
                       OutText,&nstep,&natms,&keytrj,&imcon,&tstep)) break;

            if(imcon) {
                fgets(inputl,DL_POLY_LINE_LEN,chm_in);
                fgets(inputl,DL_POLY_LINE_LEN,chm_in);
                fgets(inputl,DL_POLY_LINE_LEN,chm_in);

            }

            for(j = 0 ; j < natms ; j++) {

                fgets(inputl,DL_POLY_LINE_LEN,chm_in);
                fgets(inputl,DL_POLY_LINE_LEN,chm_in);
                if(i == (Alt - 1)) {
                    items = sscanf(inputl,"%f %f %f",&TXc,&TYc,&TZc);

                    k = gomp_PutAtomXCoord(Wstr , TXc 
                                  - sumxyz[0] , j);
                    k = gomp_PutAtomYCoord(Wstr , TYc 
                                  - sumxyz[1] , j);
                    k = gomp_PutAtomZCoord(Wstr , TZc 
                                  - sumxyz[2] , j);

                    if(keytrj > 0) {
                        fgets(inputl,DL_POLY_LINE_LEN,chm_in);
                        if(gomp_GetVelocityRetrieveState()) {
                            items = sscanf(inputl,"%f %f %f",&TXc,&TYc,&TZc);
                            tax1[j] = TXc; 
                            tay1[j] = TYc; 
                            taz1[j] = TZc;
                        }
                    }
                    if(keytrj > 1) {
                        fgets(inputl,DL_POLY_LINE_LEN,chm_in);
                        if(gomp_GetForceRetrieveState()) {
                            items = sscanf(inputl,"%f %f %f",&TXc,&TYc,&TZc);
                            tax2[j] = TXc; 
                            tay2[j] = TYc; 
                            taz2[j] = TZc;
                        }
                    }

                }
            }

        }
    } else {
 
        if(fseek(chm_in , (long)DL_PolyFrameStartPoint[Alt - 1] , SEEK_SET)) {
            gomp_PrintERROR("can't position at frame for DL_Poly trajectory");
            rewind(chm_in);
            return(1);
        }

        fgets(inputl,DL_POLY_LINE_LEN,chm_in);
        items = sscanf(inputl,"%8s%10d%10d%10d%10d%12f",
                       OutText,&nstep,&natms,&keytrj,&imcon,&tstep);

        if(imcon) {
            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
        }

/* atom velocities */
        if((keytrj > 0) && gomp_GetVelocityRetrieveState()) {

            if(gomp_GetVelocitySpace(natms))
                return(1);

            tax1 = gomp_GetModifiableVelocityXComponentPointer();
            tay1 = gomp_GetModifiableVelocityYComponentPointer();
            taz1 = gomp_GetModifiableVelocityZComponentPointer();
        }
/* atom forces     */
        if((keytrj > 1) && gomp_GetForceRetrieveState()) {

            if(gomp_GetForceSpace(natms))
                return(1);

            tax2 = gomp_GetModifiableForceXComponentPointer();
            tay2 = gomp_GetModifiableForceYComponentPointer();
            taz2 = gomp_GetModifiableForceZComponentPointer();

        }

        for(i = 0 ; i < natms ; i++) {

            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            fgets(inputl,DL_POLY_LINE_LEN,chm_in);
            items = sscanf(inputl,"%f %f %f",&TXc,&TYc,&TZc);
            j = gomp_PutAtomXCoord(Wstr , TXc 
                          - sumxyz[0] , i);
            j = gomp_PutAtomYCoord(Wstr , TYc 
                          - sumxyz[1] , i);
            j = gomp_PutAtomZCoord(Wstr , TZc 
                          - sumxyz[2] , i);

            if(keytrj > 0) {
                fgets(inputl,DL_POLY_LINE_LEN,chm_in);
                if(gomp_GetVelocityRetrieveState()) {
                    items = sscanf(inputl,"%f %f %f",&TXc,&TYc,&TZc);
                    tax1[i] = TXc; 
                    tay1[i] = TYc; 
                    taz1[i] = TZc;
                }
            }
            if(keytrj > 1) {
                fgets(inputl,DL_POLY_LINE_LEN,chm_in);
                if(gomp_GetForceRetrieveState()) {
                    items = sscanf(inputl,"%f %f %f",&TXc,&TYc,&TZc);
                    tax2[i] = TXc; 
                    tay2[i] = TYc; 
                    taz2[i] = TZc; 
                }
            }
        }
    }


    rewind(chm_in);
    return(0);

}
