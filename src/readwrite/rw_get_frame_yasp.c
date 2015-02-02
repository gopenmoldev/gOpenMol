/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002, 2005 by:
Eero Häkkinen
*/

/*

Formal specification of the YASP (Version 2) trajectory file
============================================================

FMP 23 Aug 90

Abbreviations:
i*4 four-byte integer
r*4 four-byte floating point number
r*8 eight-byte floating point number

1) Header
=========

The Header  contains information  which is stored  only  once  on  the
trajectory  file. It precedes  all individual configurations (frames).

- character*80 title
write (io_unit) title
write (io_unit) number_of_records_per_frame (i*4)

2) Frames
=========

A "frame" is the information characterising the simulated  system at a
certain  time.    Each frame contains   as many  (Fortran)  records as
specified in the header.

- "Frame"
write (io_unit) frame_number, configuration_number,
$                number_of_atoms, time 
(i*4, i*4, i*4, r*8)

The configuration_number is the index number of the frame as it  ap-
pears in the simulation,  which is  in general not  the same as  the 
index number of  the frame on  the trajectory file.  The time is the 
simulated time (e.g. ps), not the cpu time.

- "Box"
write (io_unit) (a(i), i = 1, 3), (b(i), i = 1, 3), (c(i), i = 1, 3)
(all r*8)

The current size and shape of the periodic box is  described by  the
its cell unit vectors a, b and c.  Note that a, b and c refer to the
dimensions of the entire box.  The coordinate  origin for the atomic
positions is,  however,  in the centre of the box.  At present, YASP 
only uses rectangular boxes (all angles between cell vectors are  90
degrees). Hence the cell vectors look like
a = (xbox, 0,    0   )
b = (0,    ybox, 0   )
c = (0,    0,    zbox)
with xbox, ybox, zbox beind the box  lengths in the  three cartesian 
directions.

- "Pressure"   
write (io_unit) total_isotropic_pressure,
$         Pxx, 
$                Pyx, Pyy,
$                Pzx, Pyx, Pzz
(all r*8)
The total pressure is the isotropic pressure of the system which can
in principle be calculated from the components of the pressure ten-
sor that follow. This is just to avoid recalculation. All pressures
are in kPa.

- "Energies" 
write (io_unit) number_of_energies, (e(i), i = 1, number_of_energies)
(i*4, rest r*8)
As many energy-related quantities as there are. Presently there are
e(1) = total_energy
e(2) = potential_energy
e(3) = kinetic_energy
e(4) = temperature
YASP energies are in kJ/mole, temperatures in K.

- "X-coordinates" 
write (io_unit) (x(i), i = 1, natom) (all r*4)
Coordinates (atomic positions) are in nm.

- "Y-coordinates" 
write (io_unit) (y(i), i = 1, natom) (all r*4)
Coordinates (atomic positions) are in nm.

- "Z-coordinates"
write (io_unit) (z(i), i = 1, natom) (all r*4)
Coordinates (atomic positions) are in nm.

- "X-velocities"
write (io_unit) (vx(i), i = 1,natom) (all r*4)
Atomic velocities are in nm/ps.

- "Y-velocities"
write (io_unit) (vy(i), i = 1,natom) (all r*4)
Atomic velocities are in nm/ps.

- "Z-velocities"
write (io_unit) (vz(i), i = 1,natom) (all r*4)
Atomic velocities are in nm/ps.

- "Others"     
write (io_unit) (q(i), i = 1, natom) (all r*4)
As many of these records as necessary, all describing order(natom)
quantities (energy per atom and the like). These are presently
not used.

>From this description it can be seen that the information in a frame
naturally decomposes into three parts.
a) system description
"Frame", "Box", and "Pressure"
b) "Scalar" quantities. Their number does not depend on the number of
atoms. Let's call them O(natom**0)
"Energies"
c) "Vector" quantities. Their number depends linearly on the number of 
atoms. O(natom**1).
There is probably no need to include atom-pair-related quantities
of O(natom**2) or higher into the file.


For the current implementation of this format see subroutine xtrj
($YASP_HOME/mdmacro/src/xtrj.f).

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

/***************************************************************************/
int gomp_GetFrameYasp(int alt,FILE *File_p , int iappend)  
                          /* read one frame from a yasp trajectory   
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

    static char title[BUFF_LEN];  /* yasp title                */
    static int conf_num;          /* yasp configuration number */
    static int natom;             /* number of yasp atoms      */
    static int record;
    static int nrec;
    static int kstep;
    static int i;
    static int idx1;
    static int num_yasp_energ;
    static int icount;
    static double time_ps;              /* yasp simulation time in ps   */
    static double yasp_box[9];          /* yasp box                     */
    static double yasp_energ[20];       /* yasp energies                */
    static double yasp_ppress[7];       /* yasp pressure                */
    static double yasp_press;           /* total isotropic pressure     */
    static float *vx,*vy,*vz;           /* pointers to yasp velocities  */

    static const float *sumxyz;
    static int    StartRecord;
    static int    nstep;
    static int    Wstr;
    static int    swap_bytes;
/* ................................ */
    rewind(File_p);

/* determine the "record length"    */
    natom = gomp_GetNumAtomsInMolecStruct(0);

    if(!alt) {

/*  start reading  */

/* read controll record */
        icount = fread(&record,sizeof(int), 1 , File_p);
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
/* read 80 characters   */
        icount = fread(title,sizeof(char), 80 ,File_p);
        icount = fread(&record,sizeof(int), 1 , File_p);
        title[79] = '\0';
        gomp_PrintMessage("Title:");
        gomp_PrintMessage(title);

/* read controll record */
        icount = fread(&record,sizeof(int), 1 , File_p);
/* read 1 const int */
        icount = fread(&nrec,sizeof(int), 1  ,File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( & nrec );
        }
        icount = fread(&record,sizeof(int), 1 , File_p);

        StartRecord = ftell(File_p);

/* next records are specific for one frame */
/* read controll record */
        icount = fread(&record,sizeof(int), 1 , File_p);
        icount = fread(&conf_num,sizeof(int), 1 ,File_p);
        icount = fread(&kstep,sizeof(int), 1 ,File_p);
        icount = fread(&natom,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( & natom );
        }
        icount = fread(&time_ps,sizeof(double), 1 , File_p);
        icount = fread(&record,sizeof(int), 1 , File_p);

/* read controll record */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(yasp_box,sizeof(double), 9 ,File_p);
        icount = fread(&record,sizeof(int), 1 , File_p);

/* read controll record */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&yasp_press,sizeof(double), 1 , File_p);
        icount = fread(yasp_ppress,sizeof(double), 6 , File_p);
        icount = fread(&record,sizeof(int), 1 , File_p);


/* read controll record */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&num_yasp_energ,sizeof(int), 1 , File_p);
        if ( swap_bytes ) {
            gomp_Reverse_int( & num_yasp_energ );
        }
        icount = fread(yasp_energ,sizeof(double), num_yasp_energ , File_p);
        icount = fread(&record,sizeof(int), 1 , File_p);

        record_len = (ftell(File_p) - StartRecord) + 
            (2 * sizeof(int) + natom * sizeof(float)) * 6; 

/* check number of records */

/* go first to the end of file */
        icount = fseek(File_p,0L,SEEK_END);
        if(icount != 0) {
            gomp_PrintMessage("?ERROR - can't find end of trajectory file");
            return(1);
        }

        record = ftell(File_p);
        icount = fseek(File_p, StartRecord ,SEEK_SET);
        icount = ftell(File_p);

        nstep = (record-icount)/record_len;

        sprintf(title," Atoms found                : %d",natom);
        gomp_PrintMessage(title);
        sprintf(title," Dynamics steps             : %d",nstep);
        gomp_PrintMessage(title);
        gomp_PrintMessage(" Velocities available       : YES");

        (void)gomp_SetNumberOfTrajectoryAtoms(natom);
        (void)gomp_SetNumberOfFrames(nstep);
        (void)gomp_SetTrajectoryTimeInfo(0 , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , 0);
        (void)gomp_SetTrajectoryDisplayParams(1 , nstep , 1);
        (void)gomp_PutDisplayFrameNumber(1);
    
        return(0);
    }
/* so far 100 bytes ... */
    if(alt > 0) {
        ret_fseek = fseek(File_p,((alt - 1) * record_len + StartRecord),
                          SEEK_SET);
        if(ret_fseek) {
            printf("?ERROR - can't read trajectory file");
            return(1);
        }
    }

/* next records are specific for one frame */
/* read controll record */
    icount = fread(&record,sizeof(int), 1 , File_p);
    icount = fread(&conf_num,sizeof(int), 1 ,File_p);
    icount = fread(&kstep,sizeof(int), 1 ,File_p);
    icount = fread(&natom,sizeof(int), 1 , File_p);
    if ( swap_bytes ) {
        gomp_Reverse_int( & natom );
    }
    icount = fread(&time_ps,sizeof(double), 1 , File_p);
    icount = fread(&record,sizeof(int), 1 , File_p);

/* read controll record */
    icount = fread(&record,sizeof(int), 1 ,File_p);
    icount = fread(yasp_box,sizeof(double), 9 ,File_p);
    icount = fread(&record,sizeof(int), 1 , File_p);

/* read controll record */
    icount = fread(&record,sizeof(int), 1 ,File_p);
    icount = fread(&yasp_press,sizeof(double), 1 , File_p);
    icount = fread(yasp_ppress,sizeof(double), 6 , File_p);
    icount = fread(&record,sizeof(int), 1 , File_p);

/* read controll record */
    icount = fread(&record,sizeof(int), 1 ,File_p);
    icount = fread(&num_yasp_energ,sizeof(int), 1 , File_p);
    icount = fread(yasp_energ,sizeof(double), num_yasp_energ , File_p);
    icount = fread(&record,sizeof(int), 1 , File_p);

    sumxyz   = gomp_GetTranslateArray();

    if(iappend == 0) {

/*  get pointer to coordinate vectors */
        Wstr = 0;
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y    = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z    = gomp_GetModifiableAtomZCoordPointer(Wstr);

/* read controll record */
        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(x,sizeof(float),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(y,sizeof(float),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(z,sizeof(float),natom,File_p);
        icount = fread(&record,sizeof(int), 1 , File_p);

        if (swap_bytes ) {
            gomp_Reverse_float_array( x , natom );
            gomp_Reverse_float_array( y , natom );
            gomp_Reverse_float_array( z , natom );
        }

        for(i = 0 ; i < natom ; i++) {
            x[i] = 10.0 * x[i] - sumxyz[0];
            y[i] = 10.0 * y[i] - sumxyz[1];
            z[i] = 10.0 * z[i] - sumxyz[2];
        }
    }
    else {

        idx1 = 0;
/*  get pointer to coordinate vectors */
        sprintf(title,"Yasp frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(title , natom , APPEND);
        if ( Wstr < 0 )
            return(1);
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y    = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z    = gomp_GetModifiableAtomZCoordPointer(Wstr);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&x[idx1],sizeof(float),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&y[idx1],sizeof(float),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(&z[idx1],sizeof(float),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        if (swap_bytes ) {
            gomp_Reverse_float_array( &x[idx1] , natom );
            gomp_Reverse_float_array( &y[idx1] , natom );
            gomp_Reverse_float_array( &z[idx1] , natom );
        }

        for(i = 0 ; i < natom ; i++) {
            x[i + idx1] = 10.0 * x[i + idx1] - sumxyz[0];
            y[i + idx1] = 10.0 * y[i + idx1] - sumxyz[1];
            z[i + idx1] = 10.0 * z[i + idx1] - sumxyz[2];
        }
    }

/* atom velocities */
/* retrieve velocities and forces if requested */
    if(gomp_GetVelocityRetrieveState()) {

        if(gomp_GetVelocitySpace(natom))
            return(1);

        x = gomp_GetModifiableVelocityXComponentPointer();
        y = gomp_GetModifiableVelocityYComponentPointer();
        z = gomp_GetModifiableVelocityZComponentPointer();


        vx = gomp_AllocateFloatVector(natom);
        vy = gomp_AllocateFloatVector(natom);
        vz = gomp_AllocateFloatVector(natom);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(vx,sizeof(float),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(vy,sizeof(float),natom,File_p);
        icount = fread(&record,sizeof(int), 1 ,File_p);

        icount = fread(&record,sizeof(int), 1 ,File_p);
        icount = fread(vz,sizeof(float),natom,File_p);
        icount = fread(&record,sizeof(int), 1 , File_p);

        if (swap_bytes ) {
            gomp_Reverse_float_array( vx , natom );
            gomp_Reverse_float_array( vy , natom );
            gomp_Reverse_float_array( vz , natom );
        }

        for(i = 0 ; i < natom; i++) {
            x[i]   = vx[i];
            y[i]  = vy[i];
            z[i] = vz[i];
        }

        if(gomp_CalculateVelocityMinMax())
            gomp_PrintERROR("can't calculate velocity min/max values");

        free(vx);
        free(vy);
        free(vz);
    }

    return(0);
/*                                  */
}

