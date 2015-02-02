/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Enhancements 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include "gommath.h"
#include <ctype.h>
#include <sys/types.h>
#include <stdlib.h>

#include "cell.h"
#include "gomclipbrd.h"
#include "math_oper.h"
#include "measure.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "msd.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"
#include "selection.h"

#include "stdafx.h"

#define Rabs(a)        ( ( a ) > 0 ? (a) : -(a))
#define MAX_DIFF_RES   5000
#define MAXdeg         10

static int CalcCMassList(float *tx,float *ty,float *tz,const int *sel_list);

const float *GetDynamVec();
int    GetNumDynamVec();

/*static int MSDdataSet = 0;*/
static struct {
    float *MSDarray;
    int    Start;
    int    Stop;
    int    Step;
} MSDdata = {NULL , 0 , 0 , 0};

static struct {
    float *Value;
    int *AtomIndex;
    float *ValueX;
    float *ValueXY;
    float *ValueXZ;
    float *ValueY;
    float *ValueYZ;
    float *ValueZ;
    int    Long;
} MSFdata = {NULL , 0};

static int    MSDdataSet();
static int    MSDdataDelete();
static int    MSDdataDefine(int , int , int);
static int    MSDdataFill(int , float);

static struct ResSplit {
    int    Hits[MAX_DIFF_RES];
    int    Levels;
} ResSplitList;

static int GetNumResiduesList(void);
static int SplitIntoResiduesList(const int *List, int ListLong);


/************************************************************************/
int gomp_MeanSquareDisplacement(int MASSCoord,
                              const char *text1,const char *text2,const char *text3) 
    /*  calculate mean square displacement */
/*
  int   MASSCoord;   if  = 0 calculate from atom positions
  != 0 calculate from residue mass centre 
  const char *text1;  segment name 
  const char *text2;  residue name 
  const char *text3;  atom name    
*/
/************************************************************************/
{

    static float *tx,*ty,*tz;
    static float msfx,msfy,msfz;
    static float *Xcalc,*Ycalc;
/*  static float fhelp;*/
    static int i,j,jj;
    static int atom_max,loop;
    static int slong;
    static int *sel_list;
    static float vt2,vx2,vy2,vz2;
    static float  cellaa,cellbb,cellcc;
    static int    from_frame;
    static int    to_frame;
    static int    delta_frame;
    static float  cella,cellb,cellc;
     
    float  RegCoeff[MAXdeg];
    char   OutText[BUFF_LEN];
    int    Wstr;
    float *x;
    float *y;
    float *z;
    FILE  *fp;

/* check limits */
    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintMessage("?ERROR - number of frames is not defined ");
        return(1);
    }
   

    (void)gomp_GetTrajectoryDisplayParams(&from_frame,
                                        &to_frame,
                                        &delta_frame);

    sprintf(OutText,
            "> Mean square displacement, first frame: %d , last frame: %d , step: %d",
            from_frame,to_frame,delta_frame);
    gomp_PrintMessage(OutText);
    if(MASSCoord)
        gomp_PrintMessage("Calculating from residue mass centre");
    else
        gomp_PrintMessage("Calculating from atom positions");

/* done ... */

    Wstr     = 0;
    atom_max = gomp_GetNumAtomsInMolecStruct(0);
    x        = gomp_GetModifiableAtomXCoordPointer(Wstr);
    y        = gomp_GetModifiableAtomYCoordPointer(Wstr);
    z        = gomp_GetModifiableAtomZCoordPointer(Wstr);

    sel_list = gomp_AllocateIntVector(atom_max);

    slong = gomp_MakeSelectionList(0,text1,text2,text3,sel_list);

    if(slong > 0) { 

        sprintf(OutText,"Calculating rmd for a box: %g A x %g A x %g A",
                gomp_GetCellA(),gomp_GetCellB(),gomp_GetCellC());
        gomp_PrintMessage(OutText);

        (void)MSDdataDefine(from_frame , to_frame , delta_frame);

/* check that there is a file connected */

        if(gomp_GetTrajectoryFileName() == '\0') {
            gomp_PrintERROR("no trajectory file is defined ");
            return(1);
        }

        if(gomp_GetNumberOfTrajectoryAtoms() != gomp_GetNumAtomsInMolecStruct(Wstr)) {
            gomp_PrintERROR("number of atoms defined is not the same as that in md file");
            return(1);
        }

/* save coordinates (major savings here) */
        (void)gomp_SaveAtomCoords(Wstr);

/* open trajectory file */

        if(gomp_GetTrajectoryFileType() == XMOL_TRAJ      ||
           gomp_GetTrajectoryFileType() == GROMOS96A_TRAJ ||
           gomp_GetTrajectoryFileType() == TINKER_TRAJ    ||
           gomp_GetTrajectoryFileType() == FAMBER_TRAJ    ||
           gomp_GetTrajectoryFileType() == FDL_POLY_TRAJ) {
            fp = fopen(gomp_GetTrajectoryFileName(),"r");
        } else {
            fp = fopen(gomp_GetTrajectoryFileName(),"rb");
        }

        if(fp == NULL) {
            sprintf(OutText,">>> Can't open input file: '%s'",
                    gomp_GetTrajectoryFileName());
            gomp_PrintERROR(OutText);
            free(sel_list);
            (void)gomp_GetSavedAtomCoords(Wstr);
            return(1);
        }

/* read first frame */

        if(gomp_GetOneFrame(from_frame  , fp , TRAJ_OLD)) {
            free(sel_list);
            (void)gomp_GetSavedAtomCoords(Wstr);
            fclose(fp);
            return(1);
        }

/* time zero coordinates */
        tx = gomp_AllocateFloatVector(atom_max);
        ty = gomp_AllocateFloatVector(atom_max);
        tz = gomp_AllocateFloatVector(atom_max);

        for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
            tx[i] = x[i];
            ty[i] = y[i];
            tz[i] = z[i];
        }

        if(MASSCoord) {
            (void) SplitIntoResiduesList(sel_list,slong);
            if(CalcCMassList(tx,ty,tz,sel_list)) {
                free(sel_list);
                free(tx);
                free(ty);
                free(tz);
                (void)gomp_GetSavedAtomCoords(Wstr);
                fclose(fp);
                return(1);
            }
        }
        rewind(fp);

        i = (to_frame - from_frame + 1) / delta_frame;

        Xcalc = gomp_AllocateFloatVector(i);
        Ycalc = gomp_AllocateFloatVector(i);

        loop = 0;

/* loop through the conformations */
        cellaa = 1./gomp_GetCellA();
        cellbb = 1./gomp_GetCellB();
        cellcc = 1./gomp_GetCellC();

        for(i = from_frame ; i <= to_frame ; i += delta_frame) {  
            /* frames         */

            if(gomp_GetOneFrame(i , fp , TRAJ_OLD)) {
                free(sel_list);
                free(tx);
                free(ty);
                free(tz);
                (void)gomp_GetSavedAtomCoords(Wstr);
                fclose(fp);
                return(1);
            }

            rewind(fp);

            vt2 = 0.0; 

            cella = gomp_GetCellA();
            cellb = gomp_GetCellB();
            cellc = gomp_GetCellC();

            if(!MASSCoord) {     /* calculate msd from atoms */

                for(j = 0 ; j < slong ; j++) {  /* atoms */
                    jj = sel_list[j];
                    msfx     = (tx[jj] - x[jj]);
                    msfx     = msfx  - cella * nearbyint(cellaa * msfx);
                    vx2      = msfx*msfx;
                    msfy    = (ty[jj] - y[jj]);
                    msfy    = msfy  - cellb * nearbyint(cellbb * msfy);
                    vy2     = msfy*msfy;
                    msfz   = (tz[jj] - z[jj]);
                    msfz   = msfz  - cellc * nearbyint(cellcc * msfz);
                    vz2    = msfz*msfz; 

                    vt2 += vx2 + vy2 + vz2;

                }

                vt2 /= (float)slong;
            }
            else {          /* from residue mass centre */

                if(CalcCMassList(x,y,z,sel_list)) {
                    free(sel_list);
                    free(tx);
                    free(ty);
                    free(tz);
                    (void)gomp_GetSavedAtomCoords(Wstr);
                    return(1);
                }

                for(j = 0 ; j < ResSplitList.Levels ; j++) { 
                    msfx     = (tx[j] - x[j]);
                    msfx     = msfx  - cella * nearbyint(cellaa * msfx);
                    vx2      = msfx*msfx;
                    msfy    = (ty[j] - y[j]);
                    msfy    = msfy  - cellb * nearbyint(cellbb * msfy);
                    vy2     = msfy*msfy;
                    msfz   = (tz[j] - z[j]);
                    msfz   = msfz  - cellc * nearbyint(cellcc * msfz);
                    vz2    = msfz*msfz; 

                    vt2 += vx2 + vy2 + vz2;

                }

                vt2 /= (float)ResSplitList.Levels;
            }

            (void)MSDdataFill(loop , vt2);

            Ycalc[loop] = vt2;
/*
  if(trajectory_info.time_bw_steps) 
  Xcalc[loop] = (float)(loop * trajectory_info.time_bw_steps * delta_frame);
  else
*/
            Xcalc[loop] = (float)loop;

            loop++;

        }


        if(MASSCoord) {
            sprintf(OutText,"Number of center of mass points: %d",GetNumResiduesList());
        }
        else {
            sprintf(OutText,"Number of particles: %d",slong);
        }

        gomp_PrintMessage(OutText);

        i = 0;
        (void) gomp_nreg((Xcalc-1+i) , (Ycalc-1+i) , 
            (to_frame - from_frame + 1) / delta_frame , 1 , RegCoeff);


/*
  if(trajectory_info.time_bw_steps) {
  sprintf(OutText,"Fitted line: a0 = %f , a1 = %f",RegCoeff[0],RegCoeff[1]);
  gomp_PrintMessage(OutText); 

  fhelp = 0.1 * RegCoeff[1] / 6.0;
  sprintf(OutText,"Diffusion coefficient: %.4e cm2/s",fhelp);
  gomp_PrintMessage(OutText);
  }
  else {
  
        gomp_PrintMessage("?WARNING - can't calculate diffusion coefficient");
        gomp_PrintMessage("           trajectory timing is missing");
        gomp_PrintMessage(" ** X-axis values are just running numbers, not the right time **");
        sprintf(OutText,"Fitted line: a0 = %f , a1 = %f",RegCoeff[0],RegCoeff[1]);
        gomp_PrintMessage(OutText); 
*/ 

        gomp_PrintMessage("Mean square displacement is calculated and data is saved internally");

/* put coordinates back */
        (void)gomp_GetSavedAtomCoords(Wstr);
        free(sel_list);
        free(tx);
        free(ty);
        free(tz);

        fclose(fp);
    }
    else {
        gomp_PrintMessage("?ERROR - no atoms in the selection list ");
        return(0);
    }

    return(0);
}
/*
  on return
  1:    Structure is set
  0:    Structure is not set
*/
/************************************************************************/
static int MSDdataSet()
/************************************************************************/
{
    if(MSDdata.Stop > 0) return(1);

    return(0);
}
/************************************************************************/
static int MSDdataDelete()
/************************************************************************/
{
    if(MSDdataSet()) {
        free(MSDdata.MSDarray);
        MSDdata.Start = 0;
        MSDdata.Stop  = 0;
        MSDdata.Step  = 0;
    }

    return(0);
}
/************************************************************************/
int MSDdataDefine(int Start,int Stop,int Step)
/************************************************************************/
{
    int GetNum;

    if(MSDdataSet())
        (void)MSDdataDelete();

    MSDdata.Start = Start;
    MSDdata.Stop  = Stop;
    MSDdata.Step  = Step;

    GetNum = (Stop - Start + 1)/Step;

    if(GetNum < 1) {
        (void)MSDdataDelete();
        return(1);
    }

    MSDdata.MSDarray = gomp_AllocateFloatVector(GetNum);

    return(0);
}
/************************************************************************/
int gomp_MSDdataWrite(const char *FileName)
/************************************************************************/
{
    FILE *write_p;
    int i,j;
    char OutText[BUFF_LEN];

    if(!MSDdataSet()) {
        gomp_PrintERROR("?ERROR - no array for the mean square displacement defined");
        return(1);
    }

    if(FileName[0] == '?') {
        gomp_PrintMessage("?ERROR - file name can't contain '?' ");
        return(1);
    }

    write_p = fopen(FileName,"w");
    if(write_p == NULL) {
        sprintf(OutText,"?ERROR - unable to open file : '%s' ",FileName);
        gomp_PrintMessage(OutText);
        return(1);
    }

/* ready to write now ... */
    sprintf(OutText,"Writing mean square displacement to file '%s' on disk",FileName);
    gomp_PrintMessage(OutText);

    for(i = 0 ; i < (MSDdata.Stop - MSDdata.Start + 1)/MSDdata.Step ; i++) {
/*
  if(trajectory_info.time_bw_steps) 
  j = i * trajectory_info.time_bw_steps * MSDdata.Step;
  else
*/
        j = i;
        fprintf(write_p," %d %f \n",j,MSDdata.MSDarray[i]);
    }

    fclose(write_p);
    gomp_PrintMessage("Done...");

    return(0);
}

/*
  WARNING !!!!!!!!!!

  No index checking, be aware!!!
*/
/************************************************************************/
int MSDdataFill(int Indx , float Value)
/************************************************************************/
{
    MSDdata.MSDarray[Indx] = Value;

    return(0);
}
/************************************************************************/
int SplitIntoResiduesList(const int *List, int ListLong)
/************************************************************************/
{
    int  i,j;
    int  ListH;
    int  Level;
    int  Wstr;

    Wstr                    = 0;
    ListH                   = gomp_GetAtomResNum1(Wstr , 0);
    Level                   = 0;
    ResSplitList.Hits[0]    = 1;
    ResSplitList.Levels     = 1;

    if(ListLong < 2) return(0);

    for(i = 1 ; i < ListLong ; i++) {

        j =  gomp_GetAtomResNum1(Wstr , List[i]);

        if(ListH == j) {
            ResSplitList.Hits[Level]++;
        }
        else {
            Level++;
            if(Level == MAX_DIFF_RES) {
                gomp_PrintMessage("?ERROR - MAX_DIFF_RES list is too small");
                return(1);
            }
            ListH                    = j;
            ResSplitList.Hits[Level] = 1;
        }
    }

    ResSplitList.Levels = Level;

    if(Level < 2) 
        gomp_PrintMessage("?WARNING - only one residue type found");

    return(0);
}
/************************************************************************/
int CalcCMassList(float *tx,float *ty,float *tz,const int *sel_list)
/************************************************************************/
{

    int Istart = 0;
    int Loop1  = 0;
    int Loop2  = 0;
    int i,j,jj;
    float fhelp,Tmass;
    const float *x;
    const float *y;
    const float *z;
    int     Wstr;

    Wstr = 0;
    x    = gomp_GetAtomXCoordPointer(Wstr);
    y    = gomp_GetAtomYCoordPointer(Wstr);
    z    = gomp_GetAtomZCoordPointer(Wstr);

    for(j = 0 ; j < ResSplitList.Levels ; j++) {

        Tmass = 0.0;

        for(i = Istart ; i < ResSplitList.Hits[j] + Istart ; i++) {
            jj    = sel_list[Loop2];
            fhelp = gomp_GetAtomMass(Wstr , jj);

            if(fhelp < 0.001) {
                gomp_PrintERROR("the mass is undefined for an atom, can't continue");
                return(1);
            }

            tx[Loop1] += x[jj] * fhelp;
            ty[Loop1] += y[jj] * fhelp;
            tz[Loop1] += z[jj] * fhelp;

            Tmass += fhelp;
            Loop2++;
        }
        tx[Loop1] /= Tmass;
        ty[Loop1] /= Tmass;
        tz[Loop1] /= Tmass;

        if(Tmass < 1.e-05) {
            gomp_PrintMessage("?ERROR - atom masses not defined");
            return(1);
        }

        Istart += ResSplitList.Hits[j];
        Loop1++;
    }
    return(0);
}

#if 0
/************************************************************************/
int gomp_GetMSDObs()
/************************************************************************/
{
    return((MSDdata.Stop - MSDdata.Start + 1)/MSDdata.Step);
}
/************************************************************************/
const float *gomp_GetMSDVec()
/************************************************************************/
{
    return(MSDdata.MSDarray);
}
#endif
/************************************************************************/
int GetNumResiduesList()
/************************************************************************/
{
    return(ResSplitList.Levels);
}
/*

(rms fluct)i =  SQRT(<dRi**2>)

dRi**2 = (xi - <xi>)**2 + (yi - <yi>)**2 + (zi - <zi>)**2

where xi,yi,zi are cartesian coordinates of atom i,
x - <x>, ...etc. are the displacements from the average positions,
<...> denotes averaging over the trajectory

Leif Laaksonen, CSC, 1993, 1995, 1996

*/
/************************************************************************/
int gomp_RMS_Fluctuation(int What , const char *text1, const char *text2, const char *text3) 
/*  calculate mean square fluctuation 
    int   What;   == 0 , no fitting , != 0 , fitting 
    const char *text1;  segment name 
    const char *text2;  residue name 
    const char *text3;  atom name    
*/
/************************************************************************/
{
    static float *mx,*my,*mz; /* average coordinate set */
    static float msfx,*msfx2,msfy,*msfy2,msfz,*msfz2;
    static float msfxy,*msfxy2,msfxz,*msfxz2,msfyz,*msfyz2;
    static float fhelp,fhelp1,fhelp2,fhelp3;
    static float  cellaa,cellbb,cellcc;
    static float  cella,cellb,cellc;
    static int i,j,jj;
    static int atom_max;
    static int slong;
    static int *sel_list;
    static int FirstF,LastF,StepF;
    static const float *x;
    static const float *y;
    static const float *z;
    static int    Wstr;
    static const float *sumxyz;
    static FILE  *fp;
    char OutText[BUFF_LEN];
     

    Wstr = 0;

/* check that there is a file connected */

    if(gomp_GetTrajectoryFileName() == '\0') {
        gomp_PrintMessage("?ERROR - no trajectory file is defined ");
        return(1);
    }

    if(gomp_GetNumberOfTrajectoryAtoms() != gomp_GetNumAtomsInMolecStruct(Wstr)) {
        gomp_PrintMessage("?ERROR - number of atoms defined is not the same as that in md file");
        return(1);
    }

    (void)gomp_GetTrajectoryDisplayParams(&FirstF,
                                        &LastF,
                                        &StepF);
    sprintf(OutText,
            "Root mean square displacement, first frame: %d , last frame: %d , step: %d",
            FirstF,LastF,StepF);
    gomp_PrintMessage(OutText);

    atom_max = gomp_GetNumAtomsInMolecStruct(0);

    sel_list = gomp_AllocateIntVector(atom_max);

    slong = gomp_MakeSelectionList(0,text1,text2,text3,sel_list);

    if(slong > 0) { 

        x     = gomp_GetAtomXCoordPointer(Wstr);
        y     = gomp_GetAtomYCoordPointer(Wstr);
        z     = gomp_GetAtomZCoordPointer(Wstr);

/* save coordinates */
        (void)gomp_SaveAtomCoords(Wstr);

/* first calculate average structure */

        if(gomp_TrajAvStructure()) {
            free(sel_list);
            (void)gomp_GetSavedAtomCoords(Wstr);
            return(1);
        }

/* reserve space for it */
        mx = gomp_AllocateFloatVector(atom_max);
        my = gomp_AllocateFloatVector(atom_max);
        mz = gomp_AllocateFloatVector(atom_max);

/* save it           */
        for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
            mx[i] = x[i];
            my[i] = y[i];
            mz[i] = z[i];
        }

/* open trajectory file */

        if(gomp_GetTrajectoryFileType() == XMOL_TRAJ      ||
           gomp_GetTrajectoryFileType() == GROMOS96A_TRAJ ||
           gomp_GetTrajectoryFileType() == TINKER_TRAJ    ||
           gomp_GetTrajectoryFileType() == FAMBER_TRAJ    ||
           gomp_GetTrajectoryFileType() == FDL_POLY_TRAJ) {
            fp = fopen(gomp_GetTrajectoryFileName(),"r");
        } else {
            fp = fopen(gomp_GetTrajectoryFileName(),"rb");
        }
    
        if(fp == NULL) {
            sprintf(OutText,"Can't open input file: '%s'",
                    gomp_GetTrajectoryFileName());
            gomp_PrintERROR(OutText);
            free(sel_list);
            free(mx);free(my);free(mz);
            return(1); }


        msfx2  = gomp_AllocateFloatVector(slong);
        msfxy2 = gomp_AllocateFloatVector(slong);
        msfxz2 = gomp_AllocateFloatVector(slong);
        msfy2  = gomp_AllocateFloatVector(slong);
        msfyz2 = gomp_AllocateFloatVector(slong);
        msfz2  = gomp_AllocateFloatVector(slong);

        for(i = 0 ; i < slong ; i++) {
            msfx2[i]  = msfy2[i]  = msfz2[i]  = 0.0;
            msfxy2[i] = msfxz2[i] = msfyz2[i] = 0.0;
        }


        if(What) {
            gomp_PrintMessage("Will fit the structures before RMSD calculation");
        }
        else
            gomp_PrintMessage("No fit of the structures before RMSD calculation");

/* loop through the conformations */

        /* cell dimensions */
        cella = gomp_GetCellA();
        cellb = gomp_GetCellB();
        cellc = gomp_GetCellC();

        cellaa = 1./gomp_GetCellA();
        cellbb = 1./gomp_GetCellB();
        cellcc = 1./gomp_GetCellC();

        sumxyz = gomp_GetTranslateArray();
        fhelp  = 0.0;

        for(i = FirstF ; i <= LastF ; i += StepF) {  /* frames         */

            if(gomp_GetOneFrame(i , fp , TRAJ_OLD)) {
                fclose(fp);
                free(mx);free(my);free(mz);
                free(msfx2);free(msfy2);free(msfz2);
                free(msfxy2);free(msfxz2);free(msfyz2);
                free(sel_list);
                (void)gomp_GetSavedAtomCoords(Wstr);
                return(1);
            }
            rewind(fp);

/* put up system and superimpose the two conformations */

            if(What) {
                (void)gomp_QuatFit(mx ,
                              my ,
                              mz ,
                              gomp_GetNumAtomsInMolecStruct(Wstr) ,
                              gomp_GetModifiableAtomXCoordPointer(Wstr) ,
                              gomp_GetModifiableAtomYCoordPointer(Wstr) ,
                              gomp_GetModifiableAtomZCoordPointer(Wstr) ,
                              gomp_GetNumAtomsInMolecStruct(Wstr) ,
                              sel_list , slong ,
                              sel_list , slong ,
                              OFF);
            }

            fhelp += 1.0;

            for(j = 0 ; j < slong ; j++) {  /* atoms */
                jj = sel_list[j];
                msfx       = (x[jj] - mx[jj]);
                msfx       = msfx   - cella * nearbyint(cellaa * msfx);
                msfx2[j]  += msfx*msfx;
                msfy       = (y[jj] - my[jj]);
                msfy       = msfy   - cellb * nearbyint(cellbb * msfy);
                msfy2[j]  += msfy*msfy;
                msfz       = (z[jj] - mz[jj]);
                msfz       = msfz   - cellc * nearbyint(cellcc * msfz);
                msfz2[j]  += msfz*msfz; 

                msfxy2[j]   += msfx * msfy;
                msfxz2[j]   += msfx * msfz;
                msfyz2[j]   += msfy * msfz;
            }
        }

        (void)gomp_CreateMSFset(slong);

        msfx   = 0.0;
        msfy   = 0.0;
        msfz   = 0.0;
        msfxy  = 0.0;
        msfxz  = 0.0;
        msfyz  = 0.0;
        fhelp2 = 0.0;
        fhelp -= 1.0;

/* check that there is > 1 frames ... */
        if(fhelp < 0.5) {
            gomp_PrintERROR("Number of frames has to be > 1");
            fclose(fp);
            free(mx);free(my);free(mz);
            free(msfx2);free(msfy2);free(msfz2);
            free(msfxy2);free(msfxz2);free(msfyz2);
            free(sel_list);
            (void)gomp_GetSavedAtomCoords(Wstr);
            return(1);
        }
         
        fhelp3 = fhelp * fhelp;

        for(j = 0 ; j < slong ; j++) {
            jj = sel_list[j];

            fhelp1 = (msfx2[j] + msfy2[j] + msfz2[j])/fhelp;
            fhelp2                += fhelp1;
            MSFdata.Value[j]       = fhelp1;
            MSFdata.AtomIndex[j]   = jj;

/* construct the covariance matrix */
            MSFdata.ValueX[j]      = msfx2[j] / fhelp;
            msfx                  += MSFdata.ValueX[j];
            MSFdata.ValueY[j]      = msfy2[j] / fhelp;
            msfy                  += MSFdata.ValueY[j];
            MSFdata.ValueZ[j]      = msfz2[j] / fhelp;
            msfz                  += MSFdata.ValueZ[j];

            MSFdata.ValueXY[j]     = msfxy2[j] / fhelp;
            msfxy                 += MSFdata.ValueXY[j];
            MSFdata.ValueXZ[j]     = msfxz2[j] / fhelp;
            msfxz                 += MSFdata.ValueXZ[j];
            MSFdata.ValueYZ[j]     = msfyz2[j] / fhelp;
            msfyz                 += MSFdata.ValueYZ[j];

            sprintf(OutText,
                    "(%d) RMS2 for %s:%s(%d):%s %f (x: %f y: %f z: %f [xy: %f xz: %f yz: %f])",
                    (j+1),
                    gomp_GetAtomSegName(Wstr , jj),
                    gomp_GetAtomResName(Wstr , jj),
                    gomp_GetAtomResNum1(Wstr , jj),
                    gomp_GetAtomAtmName(Wstr , jj),
                    fhelp1,MSFdata.ValueX[j],MSFdata.ValueY[j],MSFdata.ValueZ[j],
                    MSFdata.ValueXY[j],MSFdata.ValueXZ[j],MSFdata.ValueYZ[j]);
            gomp_PrintMessage(OutText);
        }

        sprintf(OutText,"Total RMS2: %f (x: %f y: %f z: %f [xy: %f xz: %f yz: %f])",
                fhelp2/(float)slong,
                msfx/(float)slong,
                msfy/(float)slong,
                msfz/(float)slong,
                msfxy/(float)slong,
                msfxz/(float)slong,
                msfyz/(float)slong);
        gomp_PrintMessage(OutText);

        free(mx);
        free(my);
        free(mz);
        free(sel_list);
        free(msfx2);free(msfy2);free(msfz2);
        free(msfxy2);free(msfxz2);free(msfyz2);

        fclose(fp);
    }
    else {
        gomp_PrintMessage("?ERROR - no atoms in the selection list ");
        return(1);
    }

    (void)gomp_GetSavedAtomCoords(Wstr);

    return(0);
}
/************************************************************************/
int    gomp_CreateMSFset(int Long)
/************************************************************************/
{
    if(MSFdata.Long) 
        (void)gomp_DeleteMSFset();

    MSFdata.Value      = gomp_AllocateFloatVector(Long);
    MSFdata.AtomIndex  = gomp_AllocateIntVector(Long);
    MSFdata.ValueX     = gomp_AllocateFloatVector(Long);
    MSFdata.ValueXY    = gomp_AllocateFloatVector(Long);
    MSFdata.ValueXZ    = gomp_AllocateFloatVector(Long);
    MSFdata.ValueY     = gomp_AllocateFloatVector(Long);
    MSFdata.ValueYZ    = gomp_AllocateFloatVector(Long);
    MSFdata.ValueZ     = gomp_AllocateFloatVector(Long);
    MSFdata.Long       = Long;

    return(0);
}
/************************************************************************/
int    gomp_DeleteMSFset()
/************************************************************************/
{
    if(MSFdata.Long) {
        free(MSFdata.Value);
        free(MSFdata.AtomIndex);
        free(MSFdata.ValueX);
        free(MSFdata.ValueXY);
        free(MSFdata.ValueXZ);
        free(MSFdata.ValueY);
        free(MSFdata.ValueYZ);
        free(MSFdata.ValueZ);
        MSFdata.Long = 0;
    }

    return(0);
}
/************************************************************************/
int    gomp_ExportMSFset(const char *FileName)
/************************************************************************/
{
    FILE *File_p;
    int i;
    char OutText[BUFF_LEN];

    if(!MSFdata.Long) {
        gomp_PrintERROR("no msfluctuation data available");
        return(1);
    }

    File_p = fopen(FileName,"w");
    if(File_p == NULL) {
        sprintf(OutText,"can't open file '%s' for writing",FileName);
        gomp_PrintERROR(OutText);
        return(1);
    }

    for(i = 0 ; i < MSFdata.Long ; i++)
        fprintf(File_p,"%d  %f %f %f %f\n",(i+1),
                MSFdata.Value[i],
                MSFdata.ValueX[i],
                MSFdata.ValueY[i],
                MSFdata.ValueZ[i]);

    fclose(File_p);

    return(0);
}

/************************************************************************/
int    gomp_GetRMSlength()
/************************************************************************/
{
    return(MSFdata.Long);
}
/************************************************************************/
const float *gomp_GetRMSfluctuation(int Which)
/************************************************************************/
{
    static float Temp[4];

    if(!MSFdata.Long) {
        gomp_PrintERROR("no rms fluctuation list is set up");
        return((const float *)NULL);
    }

    if(Which < 1 || Which > MSFdata.Long) {
        gomp_PrintERROR("your index is out of range");
        return((const float *)NULL);
    }

    Temp[0] = MSFdata.Value[Which-1];
    Temp[1] = MSFdata.ValueX[Which-1];
    Temp[2] = MSFdata.ValueY[Which-1];
    Temp[3] = MSFdata.ValueZ[Which-1];

    return(Temp);
}
/************************************************************************/
int    gomp_GetRMSatomindex(int Which)
/************************************************************************/
{
    if(!MSFdata.Long) {
        gomp_PrintERROR("no rms fluctuation list is set up");
        return(0);
    }

    if(Which < 1 || Which > MSFdata.Long) {
        gomp_PrintERROR("your index is out of range");
        return(0);
    }

    return(MSFdata.AtomIndex[Which - 1] + 1);
}

/************************************************************************/
int   gomp_CopyMSDarray2Clipboard()
/************************************************************************/
{
    static int  i;
    static char Text[BUFF_LEN];

    if(!MSDdataSet()) {
        gomp_PrintERROR("no array for the mean square displacement defined");
        return(1);
    }

    {
        char *String = NULL;

        for(i = 0 ; 
            i < (MSDdata.Stop - MSDdata.Start + 1)/MSDdata.Step ; 
            i++) {
#if defined(WIN32)
            sprintf(Text,"%d  %f\r\n",(i+1),MSDdata.MSDarray[i]);
#else
            sprintf(Text,"%d  %f\n",(i+1),MSDdata.MSDarray[i]);
#endif

            if(!i) {
                String = gomp_AllocateCharVector(strlen(Text) + 1);
                strncpy(String,Text,strlen(Text));
                String[strlen(Text)] = (char)NULL;
            }
            else {
                String = gomp_ReallocateCharVector(String , 
                                                strlen(Text) + strlen(String) + 1);
                strncat(String,Text,strlen(Text));
                String[strlen(String)] = (char)NULL;
            }
        }
        (void)gomp_CopyText2Clipboard(String);
        free(String);
    }

    return(0);
}

