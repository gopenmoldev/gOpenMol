/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved
  

This code is created from information on code from the GROMACS package.
Thanks a lot to Berk Hess (B.Hess@chem.rug.nl) for providing
code and information about the file formats.

Currently there is also some code included take from the CROMACS
package. Details and information about the authors are included.

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#include "cell.h"
#include "gomendian.h"
#include "gommain.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

#define MOD(a,b)  (((a+1)/b)*b - (a+1))

#define GROMACS_TEXT_FIELD      256
#define GROMACS_MAGIC_TRR      1993
#define GROMACS_MAGIC_XTC      1995
#define u_int  unsigned int

/* Read a frame from the trajectory file. This function opens the file,
   reads the header information 'fseeks' the place and read the frame.
   after that the file is closed again. */

static int GetFrameGromacsTrr(int, FILE *, int);
static int GetFrameGromacsXtc(int, FILE *, int);
static int GetGromacsPrecisionTrr(FILE *File_p);
   
static int xdr3dfcoord(FILE *, float *, int *, float *);
static int xdr3dfcoordDummy(FILE *, float *, int *, float *);
static int xdr_opaque(FILE *, void *, u_int );
static int xdr_float(FILE *, float *);
static int xdr_vector(FILE *, void *, u_int , u_int);
#if 0
static int xdr_double(FILE *, double  *);
#endif
static int xdr_float(FILE *, float *);
static int xdr_int(FILE *, int *);
static int TrajectorySwapBytes = 0;
static int TrajectoryPrecision = 0;

static int FileType;

typedef struct      /* This struct describes the order and the  */
/* sizes of the structs in a trjfile, sizes are given in bytes. */
{
    int ir_size;    /* Backward compatibility               */
    int e_size;     /* Backward compatibility               */
    int box_size;   /* Non zero if a box is present         */
    int   vir_size;   /* Backward compatibility             */
    int   pres_size;  /* Backward compatibility             */
    int top_size;   /* Backward compatibility               */
    int sym_size;   /* Backward compatibility               */
    int x_size;     /* Non zero if coordinates are present  */
    int v_size;     /* Non zero if velocities are present   */
    int f_size;     /* Non zero if forces are present       */

    int natoms;     /* The total number of atoms            */
    int step;       /* Current step number                  */
    int nre;        /* Backward compatibility               */
    float   t;          /* Current time                         */
    float   lambda;     /* Current value of lambda              */
} t_trnheader;

#define u_int unsigned int

/***************************************************************************/
int GetFrameGromacsTrr(int alt, FILE *File_p, int iappend)  
    /* read frame number 'alt' from charmm trajectory   */
    /* mode of operation:
       ( = 0) Check file and print trajectory information
       ( > 0) get frame number trajectory number        */
    /* *File_p is file pointer to the trajectory file   */
    /* append, if = 0 no append , if = 1 append         */
/***************************************************************************/
{
/*  
    Program to read a GROMACS binary files
*/

    static int  Wstr;
    static int  i;
    static int  icount;
    static int  record;
    static int  natoms;

/* pointers to coordinates */
    static float *Xcoord;
    static float *Ycoord;
    static float *Zcoord;
    static const float *sumxyz;
    static char   GromacsVersion[256];
    static t_trnheader   GromacsHeader;
    static float  Box[9];
    static double DBox[9];
    static int nsteps;
    static float *cvec  = NULL;
    static double  *dcvec = NULL;
    static char   label[BUFF_LEN];
    static int *DiskPointer = NULL;
    static int    DiskPointerLen = 0;
    static double DTemp1;
    static double DTemp2;

    sumxyz   = gomp_GetTranslateArray();

/*  start reading  */
  
    if(!alt) {

        if(DiskPointerLen) {
            free(DiskPointer);
            DiskPointer = NULL;
            DiskPointerLen = 0;
        }

        TrajectorySwapBytes  = 0;
        TrajectoryPrecision  = 0;
        nsteps = 0;

        if(GetGromacsPrecisionTrr(File_p)) return(1);

/* looping from here ... */

        while(!feof(File_p)) {
            icount = fread(&record,sizeof(int), 1 ,File_p);
            if(!icount) /* end of file */ 
                break;

            if(record != GROMACS_MAGIC_TRR) {
                gomp_Reverse_int( & record );
                if(record != GROMACS_MAGIC_TRR) {
                    if(nsteps) {
                        record = ftell(File_p);
                        fseek(File_p , 0L , SEEK_END);
                        record = record - ftell(File_p);
                        if(record) 
                            gomp_PrintMessage("incomplete frame");
                        break;
                    }
                    gomp_PrintERROR("can't figure out GROMACS internal file format");
                    return(1);
                }
                TrajectorySwapBytes = 1;
            }

            icount = fread(&record,sizeof(int), 1 ,File_p);

            icount = fread(&record,sizeof(int), 1 ,File_p);
            if ( TrajectorySwapBytes ) {
                gomp_Reverse_int( & record );
            }

            if(record > 255) {
                gomp_PrintERROR("version record too long > 255");
                return(1);
            }
            icount = fread(GromacsVersion,sizeof(char), record ,File_p);
            GromacsVersion[record] = (char)NULL;
            if(TrajectoryPrecision) {
                icount = fread(&GromacsHeader.ir_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.e_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.box_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.vir_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.pres_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.top_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.sym_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.x_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.v_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.f_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.natoms,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.step,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.nre,sizeof(int), 1 ,File_p);
                icount = fread(&DTemp1,sizeof(double), 1 ,File_p);
                icount = fread(&DTemp2,sizeof(double), 1 ,File_p);
                if ( TrajectorySwapBytes ) {
                    gomp_Reverse_int( & GromacsHeader.natoms );
                    gomp_Reverse_int( & GromacsHeader.box_size );
                    gomp_Reverse_int( & GromacsHeader.vir_size );
                    gomp_Reverse_int( & GromacsHeader.pres_size );
                    gomp_Reverse_int( & GromacsHeader.x_size );
                    gomp_Reverse_int( & GromacsHeader.v_size );
                    gomp_Reverse_int( & GromacsHeader.f_size );
                    gomp_Reverse_double( & DTemp1 );
                    gomp_Reverse_double( & DTemp2 );
                }
                GromacsHeader.t      = (float)DTemp1;
                GromacsHeader.lambda = (float)DTemp2;
            } else {
                icount = fread(&GromacsHeader.ir_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.e_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.box_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.vir_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.pres_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.top_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.sym_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.x_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.v_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.f_size,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.natoms,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.step,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.nre,sizeof(int), 1 ,File_p);
                icount = fread(&GromacsHeader.t,sizeof(float), 1 ,File_p);
                icount = fread(&GromacsHeader.lambda,sizeof(float), 1 ,File_p);
                if ( TrajectorySwapBytes ) {
                    gomp_Reverse_int( & GromacsHeader.natoms );
                    gomp_Reverse_int( & GromacsHeader.box_size );
                    gomp_Reverse_int( & GromacsHeader.vir_size );
                    gomp_Reverse_int( & GromacsHeader.pres_size );
                    gomp_Reverse_int( & GromacsHeader.x_size );
                    gomp_Reverse_int( & GromacsHeader.v_size );
                    gomp_Reverse_int( & GromacsHeader.f_size );
                    gomp_Reverse_float( & GromacsHeader.t );
                    gomp_Reverse_float( & GromacsHeader.lambda );
                }
            }

            natoms = GromacsHeader.natoms;

            if(natoms != gomp_GetNumAtomsInMolecStruct(0)) {
                sprintf(label,"**** ERROR number of atoms in the trajectory file does not match current number of atoms");
                gomp_PrintERROR(label);
                return(1);
            }

            if(DiskPointerLen) {
                DiskPointerLen++;
                DiskPointer = realloc(DiskPointer , (DiskPointerLen * sizeof(int)));
            } else {
                DiskPointerLen++;
                DiskPointer = malloc(sizeof(int));
            }

            if(DiskPointer == (const int *)NULL) {
                gomp_PrintERROR("can't assign memory for disk pointers");
                return(1);
            }
            DiskPointer[DiskPointerLen - 1] = ftell(File_p);

            if(GromacsHeader.box_size) {
                record = fseek(File_p , GromacsHeader.box_size , SEEK_CUR);
                if(record) {
                    sprintf(label,"?ERROR - can't read trajectory file : %s ",
                            gomp_GetTrajectoryFileName());
                    gomp_PrintMessage(label);
                    return(1);
                }
/*           icount = fread(Box,sizeof(float),GromacsHeader.box_size/sizeof(float) ,File_p);
             if (TrajectorySwapBytes ) {
             gomp_Reverse_float_array( Box, GromacsHeader.box_size/sizeof(float) );
             }
*/
            }

            if(GromacsHeader.vir_size) {
                record = fseek(File_p , GromacsHeader.vir_size , SEEK_CUR);
                if(record) {
                    sprintf(label,"?ERROR - can't read trajectory file : %s ",
                            gomp_GetTrajectoryFileName());
                    gomp_PrintMessage(label);
                    return(1);
                }
/*           icount = fread(Box,sizeof(float),GromacsHeader.vir_size/sizeof(float) ,File_p);
 */
            }
      
            if(GromacsHeader.pres_size) {
                record = fseek(File_p , GromacsHeader.pres_size , SEEK_CUR);
                if(record) {
                    sprintf(label,"?ERROR - can't read trajectory file : %s ",
                            gomp_GetTrajectoryFileName());
                    gomp_PrintMessage(label);
                    return(1);
                }
/*           icount = fread(Box,sizeof(float),GromacsHeader.pres_size/sizeof(float) ,File_p);
 */
            }

            if(GromacsHeader.x_size) {
                record = fseek(File_p , GromacsHeader.x_size , SEEK_CUR);
                if(record) {
                    sprintf(label,"?ERROR - can't read trajectory file : %s ",
                            gomp_GetTrajectoryFileName());
                    gomp_PrintMessage(label);
                    return(1);
                }
/*           cvec = malloc(GromacsHeader.x_size);
             icount = fread(cvec,sizeof(float), GromacsHeader.x_size/sizeof(float) ,File_p);
             free(cvec);
*/
            }
            if(GromacsHeader.v_size) {
                record = fseek(File_p , GromacsHeader.v_size , SEEK_CUR);
                if(record) {
                    sprintf(label,"?ERROR - can't read trajectory file : %s ",
                            gomp_GetTrajectoryFileName());
                    gomp_PrintMessage(label);
                    return(1);
                }
/*           cvec = malloc(GromacsHeader.v_size);
             icount = fread(cvec,sizeof(float), GromacsHeader.v_size/sizeof(float) ,File_p);
             free(cvec);
*/
            }
            if(GromacsHeader.f_size) {
                record = fseek(File_p , GromacsHeader.f_size , SEEK_CUR);
                if(record) {
                    sprintf(label,"?ERROR - can't read trajectory file : %s ",
                            gomp_GetTrajectoryFileName());
                    gomp_PrintMessage(label);
                    return(1);
                }
/*           cvec = malloc(GromacsHeader.f_size);
             icount = fread(cvec,sizeof(float), GromacsHeader.f_size/sizeof(float) ,File_p);
             free(cvec);
*/
            }

            nsteps++;
        }

/* update trajectory info ... */
        (void)gomp_SetNumberOfTrajectoryAtoms(natoms);
        (void)gomp_SetNumberOfFrames(nsteps);
        (void)gomp_SetNumberOfFreeAtoms(0 , 0);
        (void)gomp_SetTrajectoryDisplayParams(1 , nsteps , 1);
        (void)gomp_PutDisplayFrameNumber(1);

        sprintf(label," Info for trajectory file   : %s",gomp_GetTrajectoryFileName());
        gomp_PrintMessage(label);
        sprintf(label," GROMACS version            : %s",GromacsVersion);
        gomp_PrintMessage(label);
        sprintf(label," Atoms found                : %d",natoms);
        gomp_PrintMessage(label);
        sprintf(label," Dynamics steps             : %d",nsteps);
        gomp_PrintMessage(label);
        if(GromacsHeader.v_size) {
            gomp_PrintMessage(" Velocities available       : YES");
        } else {
            gomp_PrintMessage(" Velocities available       : NO");
        }
        if(GromacsHeader.f_size) {
            gomp_PrintMessage(" Force components available : YES");
        } else {
            gomp_PrintMessage(" Force components available : NO");
        }      
    } else {

        record = fseek(File_p , DiskPointer[alt - 1] , SEEK_SET);
        if(record) {
            sprintf(label,"?ERROR - can't read trajectory file : %s ",
                    gomp_GetTrajectoryFileName());
            gomp_PrintMessage(label);
            return(1);
        }

        if(GromacsHeader.box_size) {

            if(TrajectoryPrecision) {
                icount = fread(Box,sizeof(double),GromacsHeader.box_size/sizeof(double) ,File_p);
                if (TrajectorySwapBytes ) {
                    gomp_Reverse_double_array( DBox, GromacsHeader.box_size/sizeof(double) );
                }
                gomp_SetCellA((float)(10.0 * DBox[0]));
                gomp_SetCellB((float)(10.0 * DBox[4]));
                gomp_SetCellC((float)(10.0 * DBox[8]));
            } else {
                icount = fread(Box,sizeof(float),GromacsHeader.box_size/sizeof(float) ,File_p);
                if (TrajectorySwapBytes ) {
                    gomp_Reverse_float_array( Box, GromacsHeader.box_size/sizeof(float) );
                }
                gomp_SetCellA(10.0 * Box[0]);
                gomp_SetCellB(10.0 * Box[4]);
                gomp_SetCellC(10.0 * Box[8]);
            }
        }
/*  get pointer to coordinate vectors */

        if(!iappend) {
            Wstr   = 0;
            Xcoord = gomp_GetModifiableAtomXCoordPointer(Wstr);
            Ycoord = gomp_GetModifiableAtomYCoordPointer(Wstr);
            Zcoord = gomp_GetModifiableAtomZCoordPointer(Wstr);
            if(GromacsHeader.x_size) {
                if(TrajectoryPrecision) {
                    dcvec = malloc(GromacsHeader.x_size);
                    if(dcvec == (const double  *)NULL) {
                        gomp_PrintERROR("can't assign memory for temporary coordinates");
                        return(1);
                    }
                    icount = fread(dcvec,sizeof(double), GromacsHeader.x_size/sizeof(double) ,File_p);
                    if (TrajectorySwapBytes ) {
                        gomp_Reverse_double_array( dcvec , GromacsHeader.x_size/sizeof(double) );
                    }

                    for(i = 0 ; i < natoms ; i++) {
                        Xcoord[i] = (float)(10.0 * dcvec[i * 3    ] - sumxyz[0]); 
                        Ycoord[i] = (float)(10.0 * dcvec[i * 3 + 1] - sumxyz[1]);
                        Zcoord[i] = (float)(10.0 * dcvec[i * 3 + 2] - sumxyz[2]);
                    }
                    free(dcvec);
                } else {
                    cvec = malloc(GromacsHeader.x_size);
                    if(cvec == (const float *)NULL) {
                        gomp_PrintERROR("can't assign memory for temporary coordinates in 'rw_get_frame_gromacs'");
                        return(1);
                    }
                    icount = fread(cvec,sizeof(float), GromacsHeader.x_size/sizeof(float) ,File_p);
                    if (TrajectorySwapBytes ) {
                        gomp_Reverse_float_array( cvec , GromacsHeader.x_size/sizeof(float) );
                    }

                    for(i = 0 ; i < natoms ; i++) {
                        Xcoord[i] = 10.0 * cvec[i * 3    ] - sumxyz[0]; 
                        Ycoord[i] = 10.0 * cvec[i * 3 + 1] - sumxyz[1];
                        Zcoord[i] = 10.0 * cvec[i * 3 + 2] - sumxyz[2];
                    }
                    free(cvec);
                }
            }
        } else {
/*  get pointer to coordinate vectors */
            Wstr   = gomp_CreateMolecStruct("GROMACS frame" , natoms , APPEND);
            if ( Wstr < 0 )
                return(1);
            Xcoord = gomp_GetModifiableAtomXCoordPointer(Wstr);
            Ycoord = gomp_GetModifiableAtomYCoordPointer(Wstr);
            Zcoord = gomp_GetModifiableAtomZCoordPointer(Wstr);

            if(GromacsHeader.x_size) {
                if(TrajectoryPrecision) {
                    dcvec = malloc(GromacsHeader.x_size);
                    if(dcvec == (const double  *)NULL) {
                        gomp_PrintERROR("can't assign memory for temporary coordinates");
                        return(1);
                    }
                    icount = fread(dcvec,sizeof(double), GromacsHeader.x_size/sizeof(double) ,File_p);
                    if (TrajectorySwapBytes ) {
                        gomp_Reverse_double_array( dcvec , GromacsHeader.x_size/sizeof(double) );
                    }

                    for(i = 0 ; i < natoms ; i++) {
                        Xcoord[i] = (float)(10.0 * dcvec[i * 3    ] - sumxyz[0]); 
                        Ycoord[i] = (float)(10.0 * dcvec[i * 3 + 1] - sumxyz[1]);
                        Zcoord[i] = (float)(10.0 * dcvec[i * 3 + 2] - sumxyz[2]);
                    }
                    free(dcvec);
                } else {
                    cvec = malloc(GromacsHeader.x_size);
                    if(cvec == (const float *)NULL) {
                        gomp_PrintERROR("can't assign memory for temporary coordinates in 'rw_get_frame_gromacs'");
                        return(1);
                    }
                    icount = fread(cvec,sizeof(float), GromacsHeader.x_size/sizeof(float) ,File_p);
                    if (TrajectorySwapBytes ) {
                        gomp_Reverse_float_array( cvec , GromacsHeader.x_size/sizeof(float) );
                    }

                    for(i = 0 ; i < natoms ; i++) {
                        Xcoord[i] = 10.0 * cvec[i * 3    ] - sumxyz[0]; 
                        Ycoord[i] = 10.0 * cvec[i * 3 + 1] - sumxyz[1];
                        Zcoord[i] = 10.0 * cvec[i * 3 + 2] - sumxyz[2];
                    }
                    free(cvec);
                }
            }
        }

/* retrieve velocities and forces if requested */
        if(GromacsHeader.v_size && gomp_GetVelocityRetrieveState()) {

            if(gomp_GetVelocitySpace(natoms))
                return(1);

            Xcoord = gomp_GetModifiableVelocityXComponentPointer();
            Ycoord = gomp_GetModifiableVelocityYComponentPointer();
            Zcoord = gomp_GetModifiableVelocityZComponentPointer();
            if(GromacsHeader.v_size) {
                if(TrajectoryPrecision) {
                    dcvec = malloc(GromacsHeader.v_size);
                    if(dcvec == (const double  *)NULL) {
                        gomp_PrintERROR("can't assign memory for temporary coordinates");
                        return(1);
                    }
                    icount = fread(dcvec,sizeof(double), GromacsHeader.v_size/sizeof(double) ,File_p);
                    if (TrajectorySwapBytes ) {
                        gomp_Reverse_double_array( dcvec , GromacsHeader.v_size/sizeof(double) );
                    }

                    for(i = 0 ; i < natoms ; i++) {
                        Xcoord[i] = (float)(10.0 * dcvec[i * 3    ]); 
                        Ycoord[i] = (float)(10.0 * dcvec[i * 3 + 1]);
                        Zcoord[i] = (float)(10.0 * dcvec[i * 3 + 2]);

                    }
                    if(gomp_CalculateVelocityMinMax())
                        gomp_PrintERROR("can't calculate velocity min/max values");
                 
                    free(dcvec);
                } else {
                    cvec = malloc(GromacsHeader.v_size);
                    if(cvec == (const float *)NULL) {
                        gomp_PrintERROR("can't assign memory for temporary coordinates in 'rw_get_frame_gromacs'");
                        return(1);
                    }
                    icount = fread(cvec,sizeof(float), GromacsHeader.v_size/sizeof(float) ,File_p);
                    if (TrajectorySwapBytes ) {
                        gomp_Reverse_float_array( cvec , GromacsHeader.v_size/sizeof(float) );
                    }

                    for(i = 0 ; i < natoms ; i++) {
                        Xcoord[i] = 10.0 * cvec[i * 3    ]; 
                        Ycoord[i] = 10.0 * cvec[i * 3 + 1];
                        Zcoord[i] = 10.0 * cvec[i * 3 + 2];

                    }
                    if(gomp_CalculateVelocityMinMax())
                        gomp_PrintERROR("can't calculate velocity min/max values");
                    free(cvec);
                }
            }
        }

        if(GromacsHeader.f_size && gomp_GetForceRetrieveState()) {

            if(gomp_GetForceSpace(natoms))
                return(1);

            Xcoord = gomp_GetModifiableForceXComponentPointer();
            Ycoord = gomp_GetModifiableForceYComponentPointer();
            Zcoord = gomp_GetModifiableForceZComponentPointer();
            if(GromacsHeader.f_size) {
                if(TrajectoryPrecision) {
                    dcvec = malloc(GromacsHeader.f_size);
                    if(dcvec == (const double  *)NULL) {
                        gomp_PrintERROR("can't assign memory for temporary coordinates");
                        return(1);
                    }
                    icount = fread(dcvec,sizeof(double), GromacsHeader.f_size/sizeof(double) ,File_p);
                    if (TrajectorySwapBytes ) {
                        gomp_Reverse_double_array( dcvec , GromacsHeader.f_size/sizeof(double) );
                    }

                    for(i = 0 ; i < natoms ; i++) {
                        Xcoord[i] = (float)(10.0 * dcvec[i * 3    ]); 
                        Ycoord[i] = (float)(10.0 * dcvec[i * 3 + 1]);
                        Zcoord[i] = (float)(10.0 * dcvec[i * 3 + 2]);

                    }
                    if(gomp_CalculateForceMinMax())
                        gomp_PrintERROR("can't calculate force min/max values");
                    free(dcvec);
                } else {
                    cvec = malloc(GromacsHeader.f_size);
                    if(cvec == (const float *)NULL) {
                        gomp_PrintERROR("can't assign memory for temporary coordinates in 'rw_get_frame_gromacs'");
                        return(1);
                    }
                    icount = fread(cvec,sizeof(float), GromacsHeader.f_size/sizeof(float) ,File_p);
                    if (TrajectorySwapBytes ) {
                        gomp_Reverse_float_array( cvec , GromacsHeader.f_size/sizeof(float) );
                    }

                    for(i = 0 ; i < natoms ; i++) {
                        Xcoord[i] = 10.0 * cvec[i * 3    ]; 
                        Ycoord[i] = 10.0 * cvec[i * 3 + 1];
                        Zcoord[i] = 10.0 * cvec[i * 3 + 2];
                    }
                    if(gomp_CalculateForceMinMax())
                        gomp_PrintERROR("can't calculate force min/max values");
                    free(cvec);
                }
            }
        }
/* .......... */
    }

    return(0);
/*                                 */

}

/***************************************************************************/
int gomp_GetFrameGromacs(int alt, FILE *File_p, int iappend)  
/***************************************************************************/
{
    char TempText[BUFF_LEN];

    if(!alt) {
        FileType = 0;
        sprintf(TempText,"%s",gomp_GetTrajectoryFileName());
        if(Tcl_StringCaseMatch(TempText, "*.trr", 1) || 
           Tcl_StringCaseMatch(TempText, "*.trj", 1)) {
            FileType = 1;
            if(GetFrameGromacsTrr(alt , File_p , iappend)) 
                return(1);
        } else if(Tcl_StringCaseMatch(TempText, "*.xtc", 1)) {
            FileType = 2;
            if(GetFrameGromacsXtc(alt , File_p , iappend)) 
                return(1);
        } else {
            gomp_PrintERROR("unknow GROMACS file type extension (trr/trn/xtc)");
            FileType = 0;
            return(1);
        }
    } else {
        switch(FileType) {

        case 0: /* unknow file type */
            gomp_PrintERROR("GROMACS file type not defined");
            return(1);
            break;
        case 1: /* trr/trn file type */
            if(GetFrameGromacsTrr(alt , File_p , iappend)) 
                return(1);
            break;
        case 2: /* xtc file type */
            if(GetFrameGromacsXtc(alt , File_p , iappend)) 
                return(1);
            break;
        }

    }

    return(0);
}

/*
  The next XDR routines are not real ones. They are just some dummies
  for being able to READ values from the disk.
*/
/***************************************************************************/
static int xdr_int(FILE *xdrs, int *value)
/***************************************************************************/
{
    int icount;

    icount = fread(value,sizeof(int), 1 , xdrs);
    if ( TrajectorySwapBytes ) {
        gomp_Reverse_int( value );
    }

    return(1);
}
#if 0
/***************************************************************************/
static int xdr_double(FILE *xdrs, double  *dp)
/***************************************************************************/
{
    int icount;

    icount = fread(dp,sizeof(double), 1 , xdrs);
    if ( TrajectorySwapBytes ) {
        gomp_Reverse_double( dp );
    }

    return(1);
}
#endif
/***************************************************************************/
static int xdr_float(FILE *xdrs, float *fp)
/***************************************************************************/
{
    int icount;

    icount = fread(fp,sizeof(float), 1 , xdrs);
    if ( TrajectorySwapBytes ) {
        gomp_Reverse_float( fp );
    }

    return(1);
}


/***************************************************************************/
static int xdr_vector(FILE *xdrs, void *arrp, u_int size, u_int elsize)
/***************************************************************************/
{
    int  icount;

    icount = fread(arrp, elsize , size , xdrs);
    if (TrajectorySwapBytes) {
        gomp_Reverse_array(arrp, elsize, size);
    }

    return(1);
}

/***************************************************************************/
static int xdr_opaque(FILE *xdrs, void *cp, u_int cnt)
/***************************************************************************/
{
    int  icount;
    int  i,j;
    char dummy;

    icount = fread(cp, sizeof(char), cnt , xdrs);

    icount = cnt % sizeof(int);
    if(icount) {
        j = sizeof(int) - icount;
        for(i = 0 ; i < j ; i++) {
            icount = fread(&dummy , sizeof(char) , 1 , xdrs);
        }
    }
    return(1);
}
/*
  A filter primitive that translates between fixed-length arrays and
  their corresponding external representations.  The parameter arrp is
  the address of the pointer to the array, while size is the element
  count of the array.  The parameter elsize is the sizeof each of the
  array's elements, and elproc is an XDR filter that translates
  between the array elements' C form, and their external
  representation.  This routine returns 1 if it succeeds, 0 otherwise.



  A filter primitive that translates between fixed size opaque data
  and its external representation.  The parameter cp is the address of
  the opaque object, and cnt is its size in bytes.  This routine
  returns 1 if it succeeds, 0 otherwise.

*/

/*____________________________________________________________________________
  |
  | libxdrf - portable fortran interface to xdr. some xdr routines
  |      are C routines for compressed coordinates
  |
  | version 1.1
  |
  | This collection of routines is intended to write and read
  | data in a portable way to a file, so data written on one type
  | of machine can be read back on a different type.
  |
  | all fortran routines use an integer 'xdrid', which is an id to the
  | current xdr file, and is set by xdrfopen.
  | most routines have in integer 'ret' which is the return value.
  | The value of 'ret' is zero on failure, and most of the time one
  | on succes.
  |
  | There are three routines useful for C users:
  |  xdropen(), xdrclose(), xdr3dfcoord().
  | The first two replace xdrstdio_create and xdr_destroy, and *must* be
  | used when you plan to use xdr3dfcoord(). (they are also a bit
  | easier to interface). For writing data other than compressed coordinates 
  | you should use the standard C xdr routines (see xdr man page)
  |
  | xdrfopen(xdrid, filename, mode, ret)
  | character *(*) filename
  | character *(*) mode
  |
  | this will open the file with the given filename (string)
  | and the given mode, it returns an id in xdrid, which is
  | to be used in all other calls to xdrf routines.
  | mode is 'w' to create, or update an file, for all other
  | values of mode the file is opened for reading
  |
  | you need to call xdrfclose to flush the output and close
  | the file.
  | Note that you should not use xdrstdio_create, which comes with the
  | standard xdr library
  |
  | xdrfclose(xdrid, ret)
  | flush the data to the file, and closes the file;
  | You should not use xdr_destroy (which comes standard with
  | the xdr libraries.
  |
  | xdrfbool(xdrid, bp, ret)
  | integer pb
  |
  |     This filter produces values of either 1 or 0    
  |
  | xdrfchar(xdrid, cp, ret)
  | character cp
  |
  | filter that translate between characters and their xdr representation
  | Note that the characters in not compressed and occupies 4 bytes.
  |
  | xdrfdouble(xdrid, dp, ret)
  | double dp
  |
  | read/write a double.
  |
  | xdrffloat(xdrid, fp, ret)
  | float fp
  |
  | read/write a float.
  |
  | xdrfint(xdrid, ip, ret)
  | integer ip
  |
  | read/write integer.
  |
  | xdrflong(xdrid, lp, ret)
  | integer lp
  |
  | this routine has a possible portablility problem due to 64 bits longs.
  |
  | xdrfshort(xdrid, sp, ret)
  | integer *2 sp
  |
  | xdrfstring(xdrid, sp, maxsize, ret)
  | character *(*)
  | integer maxsize
  |
  | read/write a string, with maximum length given by maxsize
  |
  | xdrfwrapstring(xdris, sp, ret)
  | character *(*)
  |
  | read/write a string (it is the same as xdrfstring accept that it finds
  | the stringlength itself.
  |
  | xdrfvector(xdrid, cp, size, xdrfproc, ret)
  | character *(*)
  | integer size
  | external xdrfproc
  |
  | read/write an array pointed to by cp, with number of elements
  | defined by 'size'. the routine 'xdrfproc' is the name
  | of one of the above routines to read/write data (like xdrfdouble)
  | In contrast with the c-version you don't need to specify the
  | byte size of an element.
  | xdrfstring is not allowed here (it is in the c version)
  | 
  | xdrf3dfcoord(xdrid, fp, size, precision, ret)
  | real (*) fp
  | real precision
  | integer size
  |
  | this is *NOT* a standard xdr routine. I named it this way, because
  | it invites people to use the other xdr routines.
  |     It is introduced to store specifically 3d coordinates of molecules
  | (as found in molecular dynamics) and it writes it in a compressed way.
  | It starts by multiplying all numbers by precision and
  | rounding the result to integer. effectively converting
  | all floating point numbers to fixed point.
  | it uses an algorithm for compression that is optimized for
  | molecular data, but could be used for other 3d coordinates
  | as well. There is subtantial overhead involved, so call this
  | routine only if you have a large number of coordinates to read/write
  |
  | ________________________________________________________________________
  |
  | Below are the routines to be used by C programmers. Use the 'normal'
  | xdr routines to write integers, floats, etc (see man xdr)   
  |
  | int xdropen(XDR *xdrs, const char *filename, const char *type)
  | This will open the file with the given filename and the 
  | given mode. You should pass it an allocated XDR struct
  | in xdrs, to be used in all other calls to xdr routines.
  | Mode is 'w' to create, or update an file, and for all 
  | other values of mode the file is opened for reading. 
  | You need to call xdrclose to flush the output and close
  | the file.
  |
  | Note that you should not use xdrstdio_create, which
  | comes with the standard xdr library.
  |
  | int xdrclose(XDR *xdrs)
  | Flush the data to the file, and close the file;
  | You should not use xdr_destroy (which comes standard
  | with the xdr libraries).
  |  
  | int xdr3dfcoord(XDR *xdrs, const float *fp, const int *size, const float *precision)
  | This is \fInot\fR a standard xdr routine. I named it this 
  | way, because it invites people to use the other xdr 
  | routines.
  |
  | frans van hoesel hoesel@chem.rug.nl
*/  

/*___________________________________________________________________________
  |
  | what follows are the C routines for opening, closing xdr streams
  | and the routine to read/write compressed coordinates together
  | with some routines to assist in this task (those are marked
  | static and cannot be called from user programs)
*/
#define MAXABS INT_MAX-2

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x):(y))
#endif
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
static int magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645,
    812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501,
    8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536,
    82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042,
    8388607, 10568983, 13316085, 16777216 };

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))


/*____________________________________________________________________________
  |
  | sendbits - encode num into buf using the specified number of bits
  |
  | This routines appends the value of num to the bits already present in
  | the array buf. You need to give it the number of bits to use and you
  | better make sure that this number of bits is enough to hold the value
  | Also num must be positive.
  |
*/
#if 0
static void sendbits(int buf[], int num_of_bits, int num) {
    
    unsigned int cnt, lastbyte;
    int lastbits;
    unsigned char * cbuf;
    
    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = (unsigned int) buf[0];
    lastbits = buf[1];
    lastbyte =(unsigned int) buf[2];
    while (num_of_bits >= 8) {
        lastbyte = (lastbyte << 8) | ((num >> (num_of_bits -8)) /* & 0xff*/);
        cbuf[cnt++] = lastbyte >> lastbits;
        num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
        lastbyte = (lastbyte << num_of_bits) | num;
        lastbits += num_of_bits;
        if (lastbits >= 8) {
            lastbits -= 8;
            cbuf[cnt++] = lastbyte >> lastbits;
        }
    }
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    if (lastbits>0) {
        cbuf[cnt] = lastbyte << (8 - lastbits);
    }
}
#endif
/*_________________________________________________________________________
  |
  | sizeofint - calculate bitsize of an integer
  |
  | return the number of bits needed to store an integer with given max size
  |
*/

static int sizeofint(const int size) {
    int num = 1;
    int num_of_bits = 0;
    
    while (size >= num && num_of_bits < 32) {
        num_of_bits++;
        num <<= 1;
    }
    return num_of_bits;
}

/*___________________________________________________________________________
  |
  | sizeofints - calculate 'bitsize' of compressed ints
  |
  | given the number of small unsigned integers and the maximum value
  | return the number of bits needed to read or write them with the
  | routines receiveints and sendints. You need this parameter when
  | calling these routines. Note that for many calls I can use
  | the variable 'smallidx' which is exactly the number of bits, and
  | So I don't need to call 'sizeofints for those calls.
*/

static int sizeofints( const int num_of_ints, unsigned int sizes[]) {
    int i;
    unsigned int num,num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i=0; i < num_of_ints; i++) {   
        tmp = 0;
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
            tmp = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
        }
        while (tmp != 0) {
            bytes[bytecnt++] = tmp & 0xff;
            tmp >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) {
        num_of_bits++;
        num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;

}
    
#if 0
/*____________________________________________________________________________
  |
  | sendints - send a small set of small integers in compressed format
  |
  | this routine is used internally by xdr3dfcoord, to send a set of
  | small integers to the buffer. 
  | Multiplication with fixed (specified maximum ) sizes is used to get
  | to one big, multibyte integer. Allthough the routine could be
  | modified to handle sizes bigger than 16777216, or more than just
  | a few integers, this is not done, because the gain in compression
  | isn't worth the effort. Note that overflowing the multiplication
  | or the byte buffer (32 bytes) is unchecked and causes bad results.
  |
*/
 
static void sendints(int buf[], const int num_of_ints, const int num_of_bits,
                     unsigned int sizes[], unsigned int nums[]) {

    int i,num_of_bytes, bytecnt;
    unsigned int bytes[32], tmp;

    tmp = nums[0];
    num_of_bytes = 0;
    do {
        bytes[num_of_bytes++] = tmp & 0xff;
        tmp >>= 8;
    } while (tmp != 0);

    for (i = 1; i < num_of_ints; i++) {
        if (nums[i] >= sizes[i]) {
            fprintf(stderr,"major breakdown in sendints num %u doesn't "
                    "match size %u\n", nums[i], sizes[i]);
            gomp_Exit(1);
        }
        /* use one step multiply */    
        tmp = nums[i];
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
            tmp = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
        }
        while (tmp != 0) {
            bytes[bytecnt++] = tmp & 0xff;
            tmp >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    if (num_of_bits >= num_of_bytes * 8) {
        for (i = 0; i < num_of_bytes; i++) {
            sendbits(buf, 8, bytes[i]);
        }
        sendbits(buf, num_of_bits - num_of_bytes * 8, 0);
    } else {
        for (i = 0; i < num_of_bytes-1; i++) {
            sendbits(buf, 8, bytes[i]);
        }
        sendbits(buf, num_of_bits- (num_of_bytes -1) * 8, bytes[i]);
    }
}
#endif

/*___________________________________________________________________________
  |
  | receivebits - decode number from buf using specified number of bits
  | 
  | extract the number of bits from the array buf and construct an integer
  | from it. Return that value.
  |
*/

static int receivebits(int buf[], int num_of_bits) {

    int cnt, num, lastbits; 
    unsigned int lastbyte;
    unsigned const char * cbuf;
    int mask = (1 << num_of_bits) -1;

    cbuf = ((unsigned const char *)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = (unsigned int) buf[1];
    lastbyte = (unsigned int) buf[2];
    
    num = 0;
    while (num_of_bits >= 8) {
        lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
        num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
        num_of_bits -=8;
    }
    if (num_of_bits > 0) {
        if (lastbits < num_of_bits) {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | cbuf[cnt++];
        }
        lastbits -= num_of_bits;
        num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num; 
}

/*____________________________________________________________________________
  |
  | receiveints - decode 'small' integers from the buf array
  |
  | this routine is the inverse from sendints() and decodes the small integers
  | written to buf by calculating the remainder and doing divisions with
  | the given sizes[]. You need to specify the total number of bits to be
  | used from buf in num_of_bits.
  |
*/

static void receiveints(int buf[], const int num_of_ints, int num_of_bits,
                        unsigned int sizes[], int nums[]) {
    int bytes[32];
    int i, j, num_of_bytes, p, num;
    
    bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes = 0;
    while (num_of_bits > 8) {
        bytes[num_of_bytes++] = receivebits(buf, 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
        bytes[num_of_bytes++] = receivebits(buf, num_of_bits);
    }
    for (i = num_of_ints-1; i > 0; i--) {
        num = 0;
        for (j = num_of_bytes-1; j >=0; j--) {
            num = (num << 8) | bytes[j];
            p = num / sizes[i];
            bytes[j] = p;
            num = num - p * sizes[i];
        }
        nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}
    
/*____________________________________________________________________________
  |
  | xdr3dfcoord - read or write compressed 3d coordinates to xdr file.
  |
  | this routine reads or writes (depending on how you opened the file with
  | xdropen() ) a large number of 3d coordinates (stored in *fp).
  | The number of coordinates triplets to write is given by *size. On
  | read this number may be zero, in which case it reads as many as were written
  | or it may specify the number if triplets to read (which should match the
  | number written).
  | Compression is achieved by first converting all floating numbers to integer
  | using multiplication by *precision and rounding to the nearest integer.
  | Then the minimum and maximum value are calculated to determine the range.
  | The limited range of integers so found, is used to compress the coordinates.
  | In addition the differences between succesive coordinates is calculated.
  | If the difference happens to be 'small' then only the difference is saved,
  | compressing the data even more. The notion of 'small' is changed dynamically
  | and is enlarged or reduced whenever needed or possible.
  | Extra compression is achieved in the case of GROMOS and coordinates of
  | water molecules. GROMOS first writes out the Oxygen position, followed by
  | the two hydrogens. In order to make the differences smaller (and thereby
  | compression the data better) the order is changed into first one hydrogen
  | then the oxygen, followed by the other hydrogen. This is rather special, but
  | it shouldn't harm in the general case.
  |
*/
 


/*____________________________________________________________________________
  |
  | xdr3dfcoord - read or write compressed 3d coordinates to xdr file.
  |
  | this routine reads or writes (depending on how you opened the file with
  | xdropen() ) a large number of 3d coordinates (stored in *fp).
  | The number of coordinates triplets to write is given by *size. On
  | read this number may be zero, in which case it reads as many as were written
  | or it may specify the number if triplets to read (which should match the
  | number written).
  | Compression is achieved by first converting all floating numbers to integer
  | using multiplication by *precision and rounding to the nearest integer.
  | Then the minimum and maximum value are calculated to determine the range.
  | The limited range of integers so found, is used to compress the coordinates.
  | In addition the differences between succesive coordinates is calculated.
  | If the difference happens to be 'small' then only the difference is saved,
  | compressing the data even more. The notion of 'small' is changed dynamically
  | and is enlarged or reduced whenever needed or possible.
  | Extra compression is achieved in the case of GROMOS and coordinates of
  | water molecules. GROMOS first writes out the Oxygen position, followed by
  | the two hydrogens. In order to make the differences smaller (and thereby
  | compression the data better) the order is changed into first one hydrogen
  | then the oxygen, followed by the other hydrogen. This is rather special, but
  | it shouldn't harm in the general case.
  |
*/
 
static int xdr3dfcoord(FILE *xdrs, float *fp, int *size, float *precision) {
    

    static int *ip = NULL;
    static int oldsize;
    static int *buf = NULL;

    int minint[3], maxint[3], *lip;
    int smallidx;
    int minidx, maxidx;
    unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
    int flag, k;
    int small, smaller, larger, i, is_smaller, run;
    float *lfp;
    int tmp, *thiscoord,  prevcoord[3];

    int bufsize, lsize;
    unsigned int bitsize;
    float inv_precision;
    
    /* xdrs is open for reading */
    
    if (xdr_int(xdrs, &lsize) == 0) 
        return 0;
    if (*size != 0 && lsize != *size) {
        fprintf(stderr, "wrong number of coordinates in xdr3dfcoor; "
                "%d arg vs %d in file", *size, lsize);
    }
    *size = lsize;
    size3 = *size * 3;
    if (*size <= 9) {
        return (xdr_vector(xdrs, fp, (u_int)size3,(u_int)sizeof(*fp)));
    }
    xdr_float(xdrs, precision);
    if (ip == NULL) {
        ip = malloc((size_t)(size3 * sizeof(*ip)));
        if (ip == NULL) {
            gomp_PrintERROR("malloc failed in xdr3dfcoord (1)");
            return(1);
        }
        bufsize = (int)(size3 * 1.2);
        buf = malloc((size_t)(bufsize * sizeof(*buf)));
        if (buf == NULL) {
            gomp_PrintERROR("malloc failed in xdr3dfcoord (2)");
            return(1);
        }
        oldsize = *size;
    } else if (*size > oldsize) {
        ip = realloc(ip, (size_t)(size3 * sizeof(*ip)));
        if (ip == NULL) {
            gomp_PrintERROR("malloc failed in xdr3dfcoord (3)");
            return(1);
        }
        bufsize = (int)(size3 * 1.2);
        buf = realloc(buf, (size_t)(bufsize * sizeof(*buf)));
        if (buf == NULL) {
            gomp_PrintERROR("malloc failed in xdr3dfcoord (4)");
            return(1);
        }
        oldsize = *size;
    }
    buf[0] = buf[1] = buf[2] = 0;
    
    xdr_int(xdrs, &(minint[0]));
    xdr_int(xdrs, &(minint[1]));
    xdr_int(xdrs, &(minint[2]));

    xdr_int(xdrs, &(maxint[0]));
    xdr_int(xdrs, &(maxint[1]));
    xdr_int(xdrs, &(maxint[2]));
        
    sizeint[0] = maxint[0] - minint[0]+1;
    sizeint[1] = maxint[1] - minint[1]+1;
    sizeint[2] = maxint[2] - minint[2]+1;
    
    /* check if one of the sizes is to big to be multiplied */
    if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
        bitsizeint[0] = sizeofint(sizeint[0]);
        bitsizeint[1] = sizeofint(sizeint[1]);
        bitsizeint[2] = sizeofint(sizeint[2]);
        bitsize = 0; /* flag the use of large sizes */
    } else {
        bitsize = sizeofints(3, sizeint);
    }
    
    if (xdr_int(xdrs, &smallidx) == 0)  
        return 0;
    maxidx = MIN((int)LASTIDX, smallidx + 8) ;
    minidx = maxidx - 8; /* often this equal smallidx */
    smaller = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
    small = magicints[smallidx] / 2;
    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
    larger = magicints[maxidx];

    /* buf[0] holds the length in bytes */

    if (xdr_int(xdrs, &(buf[0])) == 0)
        return 0;
    if (xdr_opaque(xdrs, &(buf[3]), (u_int)buf[0]) == 0)
        return 0;
    buf[0] = buf[1] = buf[2] = 0;
    
    lfp = fp;
    inv_precision = 1.0 / * precision;
    run = 0;
    i = 0;
    lip = ip;
    while ( i < lsize ) {
        thiscoord = &lip[i * 3];

        if (bitsize == 0) {
            thiscoord[0] = receivebits(buf, bitsizeint[0]);
            thiscoord[1] = receivebits(buf, bitsizeint[1]);
            thiscoord[2] = receivebits(buf, bitsizeint[2]);
        } else {
            receiveints(buf, 3, bitsize, sizeint, thiscoord);
        }
        
        i++;
        thiscoord[0] += minint[0];
        thiscoord[1] += minint[1];
        thiscoord[2] += minint[2];
        
        prevcoord[0] = thiscoord[0];
        prevcoord[1] = thiscoord[1];
        prevcoord[2] = thiscoord[2];
        
       
        flag = receivebits(buf, 1);
        is_smaller = 0;
        if (flag == 1) {
            run = receivebits(buf, 5);
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller--;
        }
        if (run > 0) {
            thiscoord += 3;
            for (k = 0; k < run; k+=3) {
                receiveints(buf, 3, smallidx, sizesmall, thiscoord);
                i++;
                thiscoord[0] += prevcoord[0] - small;
                thiscoord[1] += prevcoord[1] - small;
                thiscoord[2] += prevcoord[2] - small;
                if (k == 0) {
                    /* interchange first with second atom for better
                     * compression of water molecules
                     */
                    tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
                    prevcoord[0] = tmp;
                    tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
                    prevcoord[1] = tmp;
                    tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
                    prevcoord[2] = tmp;
                    *lfp++ = prevcoord[0] * inv_precision;
                    *lfp++ = prevcoord[1] * inv_precision;
                    *lfp++ = prevcoord[2] * inv_precision;
                } else {
                    prevcoord[0] = thiscoord[0];
                    prevcoord[1] = thiscoord[1];
                    prevcoord[2] = thiscoord[2];
                }
                *lfp++ = thiscoord[0] * inv_precision;
                *lfp++ = thiscoord[1] * inv_precision;
                *lfp++ = thiscoord[2] * inv_precision;
            }
        } else {
            *lfp++ = thiscoord[0] * inv_precision;
            *lfp++ = thiscoord[1] * inv_precision;
            *lfp++ = thiscoord[2] * inv_precision;      
        }
        smallidx += is_smaller;
        if (is_smaller < 0) {
            small = smaller;
            if (smallidx > FIRSTIDX) {
                smaller = magicints[smallidx - 1] /2;
            } else {
                smaller = 0;
            }
        } else if (is_smaller > 0) {
            smaller = small;
            small = magicints[smallidx] / 2;
        }
        sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
    }

/* save the space and return */
    free(ip);
    ip = NULL;
    free(buf);
    buf = NULL;
    oldsize = 0;
/* ......................... */

    return 0;
}

static int xdr3dfcoordDummy(FILE *xdrs, float *fp, int *size, float *precision) {
    static int buf[3];

    int minint[3], maxint[3];
    int smallidx;
    unsigned size3;
    int k;

    int lsize;
    
    /* xdrs is open for reading */
    
    if (xdr_int(xdrs, &lsize) == 0) 
        return 0;
    if (*size != 0 && lsize != *size) {
        fprintf(stderr, "wrong number of coordinates in xdr3dfcoor; "
                "%d arg vs %d in file", *size, lsize);
    }
    *size = lsize;
    size3 = *size * 3;
    xdr_float(xdrs, precision);
    
    xdr_int(xdrs, &(minint[0]));
    xdr_int(xdrs, &(minint[1]));
    xdr_int(xdrs, &(minint[2]));

    xdr_int(xdrs, &(maxint[0]));
    xdr_int(xdrs, &(maxint[1]));
    xdr_int(xdrs, &(maxint[2]));
            
    if (xdr_int(xdrs, &smallidx) == 0)  
        return 0;

    /* buf[0] holds the length in bytes */

    if (xdr_int(xdrs, &(buf[0])) == 0)
        return 0;

/* leftover to full word boundary */
    k = buf[0] % sizeof(int);
    if(k) k = sizeof(int) - k;
    k = fseek(xdrs , (buf[0] + k) , SEEK_CUR);
    if(k) {
        gomp_PrintMessage("can't move pointer in trajectory file");
        return(1);
    }

    return 0;
}

/***************************************************************************/
int GetFrameGromacsXtc(int alt, FILE *File_p, int iappend)  
    /* read frame number 'alt' from charmm trajectory   */
    /* mode of operation:
       ( = 0) Check file and print trajectory information
       ( > 0) get frame number trajectory number        */
    /* *File_p is file pointer to the trajectory file   */
    /* append, if = 0 no append , if = 1 append         */
/***************************************************************************/
{
/*  
    Program to read a GROMACS binary files
*/

    static int  Wstr;
    static int  i;
    static int  icount;
    static int  record;
    static int  natoms;

/* pointers to coordinates */
    static float *Xcoord;
    static float *Ycoord;
    static float *Zcoord;
    static const float *sumxyz;
    static float Box[9];
    static int nsteps;
    static float *cvec = (float *)NULL;
    static char   label[BUFF_LEN];
    static int *DiskPointer = NULL;
    static int    DiskPointerLen = 0;
    static int    FrameNumber;
    static float  Time;
    static int    Size;
    static float  Precision;

    sumxyz   = gomp_GetTranslateArray();

/*  start reading  */
  
    if(!alt) {

        if(DiskPointerLen) {
            free(DiskPointer);
            DiskPointer = NULL;
            DiskPointerLen = 0;
        }

        TrajectorySwapBytes = 0;

/* magic number */
        icount = fread(&record,sizeof(int), 1 ,File_p);

        if(record != GROMACS_MAGIC_XTC) {
            gomp_Reverse_int( & record );
            if(record != GROMACS_MAGIC_XTC) {
                gomp_PrintERROR("can't figure out GROMACS internal file format");
                return(1);
            }
            TrajectorySwapBytes = 1;
        }
/* natoms */
        icount = fread(&natoms,sizeof(int), 1 ,File_p);
        if ( TrajectorySwapBytes ) {
            gomp_Reverse_int( & natoms );
        }

        if(natoms != gomp_GetNumAtomsInMolecStruct(0)) {
            sprintf(label,"**** ERROR number of atoms in the trajectory file does not match current number of atoms");
            gomp_PrintERROR(label);
            return(1);
        }

/* frame number */
        icount = fread(&FrameNumber,sizeof(int), 1 ,File_p);
        if ( TrajectorySwapBytes ) {
            gomp_Reverse_int( & FrameNumber );
        }
/* time */
        icount = fread(&Time,sizeof(float), 1 ,File_p);
        if ( TrajectorySwapBytes ) {
            gomp_Reverse_float( & Time );
        }

        DiskPointerLen = 1;
        DiskPointer = malloc(sizeof(int));
        if(DiskPointer == (const int *)NULL) {
            gomp_PrintERROR("can't assign memory for disk pointers");
            return(1);
        }
        DiskPointer[DiskPointerLen - 1] = ftell(File_p);

/* box */
        icount = fread(Box,sizeof(float), 9 ,File_p);
        if (TrajectorySwapBytes ) {
            gomp_Reverse_float_array( Box, 9 );
        }
        gomp_SetCellA(10.0 * Box[0]);
        gomp_SetCellB(10.0 * Box[4]);
        gomp_SetCellC(10.0 * Box[8]);

/* coordinate array */

        if(xdr3dfcoordDummy(File_p, NULL , &natoms, &Precision)) {
            gomp_PrintERROR("can't get coordinates");
            return(1);
        }

/* looping from here ... */
        nsteps = 1;
        while(!feof(File_p)) {

/* magic number */
            icount = fread(&record,sizeof(int), 1 ,File_p);
            if(!icount) break;

/* natoms */
            icount = fread(&natoms,sizeof(int), 1 ,File_p);
            if ( TrajectorySwapBytes ) {
                gomp_Reverse_int( & natoms );
            }
/* frame number */
            icount = fread(&FrameNumber,sizeof(int), 1 ,File_p);
            if ( TrajectorySwapBytes ) {
                gomp_Reverse_int( & FrameNumber );
            }
/* time */
            icount = fread(&Time,sizeof(float), 1 ,File_p);
            if ( TrajectorySwapBytes ) {
                gomp_Reverse_float( & Time );
            }

            DiskPointerLen++;
            DiskPointer = realloc(DiskPointer , (DiskPointerLen * sizeof(int)));
            if(DiskPointer == (const int *)NULL) {
                gomp_PrintERROR("can't assign memory for disk pointers");
                return(1);
            }
            DiskPointer[DiskPointerLen - 1] = ftell(File_p);

/* box */
            icount = fread(Box,sizeof(float), 9 ,File_p);
            if (TrajectorySwapBytes ) {
                gomp_Reverse_float_array( Box, 9 );
            }
            gomp_SetCellA(10.0 * Box[0]);
            gomp_SetCellB(10.0 * Box[4]);
            gomp_SetCellC(10.0 * Box[8]);

/* coordinate array */

            if(xdr3dfcoordDummy(File_p, NULL , &natoms, &Precision)) {
                gomp_PrintERROR("can't get coordinates");
                return(1);
            }

            nsteps++;
        }

/* update trajectory info ... */
        (void)gomp_SetNumberOfTrajectoryAtoms(natoms);
        (void)gomp_SetNumberOfFrames(nsteps);
        (void)gomp_SetNumberOfFreeAtoms(0 , 0);
        (void)gomp_SetTrajectoryDisplayParams(1 , nsteps , 1);
        (void)gomp_PutDisplayFrameNumber(1);

        sprintf(label," Info for trajectory file   : %s",gomp_GetTrajectoryFileName());
        gomp_PrintMessage(label);
        sprintf(label," Atoms found                : %d",natoms);
        gomp_PrintMessage(label);
        sprintf(label," Dynamics steps             : %d",nsteps);
        gomp_PrintMessage(label);      
    } else {

        record = fseek(File_p , DiskPointer[alt - 1] , SEEK_SET);
        if(record) {
            sprintf(label,"?ERROR - can't read trajectory file : %s ",
                    gomp_GetTrajectoryFileName());
            gomp_PrintMessage(label);
            return(1);
        }

/*  get pointer to coordinate vectors */

        if(!iappend) {
            Wstr   = 0;
            Xcoord = gomp_GetModifiableAtomXCoordPointer(Wstr);
            Ycoord = gomp_GetModifiableAtomYCoordPointer(Wstr);
            Zcoord = gomp_GetModifiableAtomZCoordPointer(Wstr);

/* box */
            icount = fread(Box,sizeof(float), 9 ,File_p);
            if (TrajectorySwapBytes ) {
                gomp_Reverse_float_array( Box, 9 );
            }
            gomp_SetCellA(10.0 * Box[0]);
            gomp_SetCellB(10.0 * Box[4]);
            gomp_SetCellC(10.0 * Box[8]);

/* coordinate array */

            Size = 3 * natoms * sizeof(float); 
            cvec = malloc(Size);
            if(cvec == (const float *)NULL) {
                gomp_PrintERROR("can't assign memory for temporary coordinate array");
                return(1);
            }
            if(xdr3dfcoord(File_p, cvec , &natoms, &Precision)) {
                gomp_PrintERROR("can't get coordinates");
                return(1);
            }

            for(i = 0 ; i < natoms ; i++) {
                Xcoord[i] = 10.0 * cvec[i * 3    ] - sumxyz[0]; 
                Ycoord[i] = 10.0 * cvec[i * 3 + 1] - sumxyz[1];
                Zcoord[i] = 10.0 * cvec[i * 3 + 2] - sumxyz[2];
            }
            free(cvec);
        } else {
/*  get pointer to coordinate vectors */
            Wstr   = gomp_CreateMolecStruct("GROMACS frame" , natoms , APPEND);
            if ( Wstr < 0 )
                return(1);
            Xcoord = gomp_GetModifiableAtomXCoordPointer(Wstr);
            Ycoord = gomp_GetModifiableAtomYCoordPointer(Wstr);
            Zcoord = gomp_GetModifiableAtomZCoordPointer(Wstr);

/* box */
            icount = fread(Box,sizeof(float), 9 ,File_p);
            if (TrajectorySwapBytes ) {
                gomp_Reverse_float_array( Box, 9 );
            }
            gomp_SetCellA(10.0 * Box[0]);
            gomp_SetCellB(10.0 * Box[4]);
            gomp_SetCellC(10.0 * Box[8]);

/* coordinate array */

            Size = 3 * natoms * sizeof(float); 
            cvec = malloc(Size);
        
            if(xdr3dfcoord(File_p, cvec , &natoms, &Precision)) {
                gomp_PrintERROR("can't get coordinates");
                return(1);
            }

            for(i = 0 ; i < natoms ; i++) {
                Xcoord[i] = 10.0 * cvec[i * 3    ] - sumxyz[0]; 
                Ycoord[i] = 10.0 * cvec[i * 3 + 1] - sumxyz[1];
                Zcoord[i] = 10.0 * cvec[i * 3 + 2] - sumxyz[2];
            }
            free(cvec);
        }
    
    }

    return(0);
/*                                 */

}

/***************************************************************************/
int GetGromacsPrecisionTrr(FILE *File_p)  
/***************************************************************************/
{
/*  
    Program to read a GROMACS binary files
*/

    static int  icount;
    static int  record;
    static int  natoms;

/* pointers to coordinates */
    static char   GromacsVersion[256];
    static t_trnheader   GromacsHeader;
    static char   label[BUFF_LEN];

  
    TrajectoryPrecision = 0;
    TrajectorySwapBytes = 0;

    rewind(File_p);

    icount = fread(&record,sizeof(int), 1 ,File_p);

    if(record != GROMACS_MAGIC_TRR) {
        gomp_Reverse_int( & record );
        if(record != GROMACS_MAGIC_TRR) {
            gomp_PrintERROR("can't figure out GROMACS internal file format");
            rewind(File_p);
            return(1);
        }
        TrajectorySwapBytes = 1;
    }

    icount = fread(&record,sizeof(int), 1 ,File_p);

    icount = fread(&record,sizeof(int), 1 ,File_p);
    if ( TrajectorySwapBytes ) {
        gomp_Reverse_int( & record );
    }

    if(record > 255) {
        gomp_PrintERROR("version record too long > 255");
        rewind(File_p);
        return(1);
    }
    icount = fread(GromacsVersion,sizeof(char), record ,File_p);
    GromacsVersion[record] = (char)NULL;

    icount = fread(&GromacsHeader.ir_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.e_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.box_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.vir_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.pres_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.top_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.sym_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.x_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.v_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.f_size,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.natoms,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.step,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.nre,sizeof(int), 1 ,File_p);
    icount = fread(&GromacsHeader.t,sizeof(float), 1 ,File_p);
    icount = fread(&GromacsHeader.lambda,sizeof(float), 1 ,File_p);
      
    if ( TrajectorySwapBytes ) {
        gomp_Reverse_int( & GromacsHeader.natoms );
        gomp_Reverse_int( & GromacsHeader.box_size );
        gomp_Reverse_int( & GromacsHeader.vir_size );
        gomp_Reverse_int( & GromacsHeader.pres_size );
        gomp_Reverse_int( & GromacsHeader.x_size );
        gomp_Reverse_int( & GromacsHeader.v_size );
        gomp_Reverse_int( & GromacsHeader.f_size );
    }

    natoms = GromacsHeader.natoms;

    if(natoms != gomp_GetNumAtomsInMolecStruct(0)) {
        sprintf(label,"**** ERROR number of atoms in the trajectory file does not match current number of atoms");
        gomp_PrintERROR(label);
        rewind(File_p);
        return(1);
    }

    if(GromacsHeader.box_size/sizeof(float) != 9) {
        if(GromacsHeader.box_size/sizeof(double) != 9) {
            gomp_PrintERROR("unknown internal format (neither real nor double) for cell dimensions");
            rewind(File_p);
            return(1);
        }
        TrajectoryPrecision = 1;
    }

        
    rewind(File_p);
    return(0);
}
