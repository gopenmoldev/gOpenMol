/*  

Copyright (c) 1995 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing , ESPOO , FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  


Enhancements 2002 - 2004 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <tcl.h>

#include "bond.h"
#include "coord_file.h"
#include "gomstring.h"
#include "memalloc.h"
#include "molecoord.h"
#include "printmsg.h"
#include "tclutils.h"
#include "trajectory.h"
#include "gomfile.h"

#include "stdafx.h"

/*  trajectory functions and handles */

static int DefaultTrajectoryFileType = CHARMM_TRAJ;
static struct {
    const char *tag;
    char file_ext[BUFF_LEN];
} TrajectoryFileTypes[] = {
#define  AMBER_TRAJ      1
    {"ambe$r",     ""},
#define  CERIUS2_TRAJ    2
    {"ceri$us2",   "trj"},
#define  CHARMM_TRAJ     3
    {"char$mm",    "dcd"},
#define  DISCOVER_TRAJ   4
    {"disc$over",  "his"},
#define  FDL_POLY_TRAJ   5
    {"fdl_$poly",  ""},
#define  UDL_POLY_TRAJ   6
    {"udl_$poly",  ""},
#define  FAMBER_TRAJ     7
    {"famb$er",    ""},
#define  GROMACS_TRAJ    8
    {"gromacs",    ""},
#define  GROMOS_TRAJ     9
    {"gromos",     ""},
#define  GROMOS96A_TRAJ  10
    {"gromos96a",  "dat"},
#define  GROMOS96B_TRAJ  11
    {"gromos96b",  "dat"},
#define  HYPERCHEM_TRAJ  12
    {"hype$rchem", "snp"},
#define  MUMOD_TRAJ      13
    {"mumo$d",     ""},
#define  OPENMOL_TRAJ    14
    {"open",       ""},
#define  TINKER_TRAJ     15
    {"tink$er",    "arc"},
#define  XMOL_TRAJ       16
    {"xmol",       "xmol"},
#define  XPLOR_TRAJ      17
    {"xplo$r",     "xplor"},
#define  YASP_TRAJ       18
    {"yasp",       "yasp"},
    {NULL,         ""}
};

static TrajectoryInfo_t TrajectoryInfo;

static struct {
    int  nfreat;         /* number of free atoms           */
    int *free_atom_list; /* list of free (not fixed) atoms */
} FreeAtomList = { 0 , (int *)NULL};

/*  end of trajectory functions and handles */

/*************************************************************************/
int gomp_ParseTrajType(const char *intext) /* determine trajectory file type */
/*************************************************************************/
{
    char text[BUFF_LEN];
    int Type;

    if(intext[0] == '\0') 
        return(DefaultTrajectoryFileType); /* default trajectory type (CHARMM)*/

    gomp_CopyString(text,intext,BUFF_LEN);
    text[BUFF_LEN-1] = '\0';
    (void)gomp_String2Lower(text);

/* do it ... */
    /* Type numbering start from 1. C indeces start from 0. */
    for ( Type = 1 ; TrajectoryFileTypes[Type - 1].tag ; Type++ ) {
        if ( gomp_StringMatch(text,TrajectoryFileTypes[Type - 1].tag))
            return(Type);
    }

    gomp_FormatERROR("can't parse trajectory type '%s'",text);
    return(0);
}

/*************************************************************************/
int gomp_SetTrajectoryFileName(const char *FileName)
/*************************************************************************/
{
    gomp_CopyString(TrajectoryInfo.traj_file,FileName,BUFF_LEN);

    return(0);
}
/*************************************************************************/
const char *gomp_GetTrajectoryFileName(void)
/*************************************************************************/
{
    return(TrajectoryInfo.traj_file);
}

/*************************************************************************/
int   gomp_SetTrajectoryFileType(int FileType)
/*************************************************************************/
{
    TrajectoryInfo.type = FileType;

    return(0);
}
/*************************************************************************/
int   gomp_GetTrajectoryFileType(void)
/*************************************************************************/
{
    return(TrajectoryInfo.type);
}
/*************************************************************************/
const char *gomp_GetTrajectoryFileTypeName(void)
/*************************************************************************/
{
    static char OutText[BUFF_LEN];
    int         Type = gomp_GetTrajectoryFileType();
    char       *p;
    const char *q;

    if ( Type == 0 )
        /* Type isn't set. */
/*        Type = CHARMM_TRAJ; */
/* Changed type to point to the default or currently defined type */
          Type = gomp_GetTrajectoryTypeDefault();

    /* Copy the type name but omit $ characters. */
    for ( p = OutText , q = TrajectoryFileTypes[Type - 1].tag ; *q ; q++ ) {
        if ( *q != '$' )
            *p++ = *q;
    }
    *p = '\0';

    return OutText;
}
/*************************************************************************/
int gomp_SetNumberOfFrames(int Frames)
/*************************************************************************/
{

    TrajectoryInfo.nstep = Frames;

    return(0);
}
/*************************************************************************/
int gomp_GetNumberOfFrames(void)
/*************************************************************************/
{

    return(TrajectoryInfo.nstep);
}

/*************************************************************************/
int   gomp_SetTrajectoryTimeInfo(int  tbf, int toff)
/*************************************************************************/
{
    TrajectoryInfo.time_bw_steps    = tbf;
    TrajectoryInfo.time_first_frame = toff;

    return(0);
}
/*************************************************************************/
int   gomp_GetTrajectoryTimeInfo(int *tbf, int *toff)
/*************************************************************************/
{
    *tbf  = TrajectoryInfo.time_bw_steps;
    *toff = TrajectoryInfo.time_first_frame;

    return(0);
}

/*************************************************************************/
int   gomp_SetTrajectoryDisplayParams(int  ff, int lf , int df)
/*************************************************************************/
{
    if(!TrajectoryInfo.nstep) {
        gomp_PrintERROR("no trajectory file is defined or can't resolve number of frames");
        return(1);
    }

    if(ff < 1 || ff > TrajectoryInfo.nstep) {
        gomp_PrintERROR("first index is < 0 or > number of frames");
        return(1);
    }
    if(lf < 1 || lf > TrajectoryInfo.nstep) {
        gomp_PrintERROR("last index is < 0 or > number of frames");
        return(1);
    }
    if(df < 1 || df > TrajectoryInfo.nstep) {
        gomp_PrintERROR("step index is < 0 or > number of frames");
        return(1);
    }
    TrajectoryInfo.first_frame = ff;
    TrajectoryInfo.last_frame  = lf;
    TrajectoryInfo.delta_frame = df;

    return(0);
}
/*************************************************************************/
int   gomp_GetTrajectoryDisplayParams(int *ff, int *lf, int *df)
/*************************************************************************/
{
    *ff = TrajectoryInfo.first_frame;
    *lf = TrajectoryInfo.last_frame;
    *df = TrajectoryInfo.delta_frame;

    return(0);
}


/*************************************************************************/
int   gomp_SetTrajectoryStructureInUse(void)
/*************************************************************************/
{
    TrajectoryInfo.inuse = 1;

    return(0);
}
/*************************************************************************/
int   gomp_GetTrajectoryStructureState(void)
/*************************************************************************/
{
    return(TrajectoryInfo.inuse);
}
/*************************************************************************/
int   gomp_SetTrajectoryStructureOff(void)
/*************************************************************************/
{
    const char *Value;

    Value = Tcl_GetVar(gomp_GetTclInterp(),"lulPostReadFrameTrigger",TCL_GLOBAL_ONLY);

    if(Value) {
        if(TCL_OK != Tcl_UnsetVar(gomp_GetTclInterp(),"lulPostReadFrameTrigger",TCL_GLOBAL_ONLY)) {
            gomp_PrintERROR("can't unset the 'lulPostReadFrameTrigger' variable");
        }
    }

    TrajectoryInfo.inuse        = 0;
    TrajectoryInfo.type         = 0;
    TrajectoryInfo.natom        = 0;
    TrajectoryInfo.nstep        = 0;
    TrajectoryInfo.traj_file[0] = '\0';

    return(0);
}
/*************************************************************************/
int   gomp_DeleteTrajectoryStructure(void)
/*************************************************************************/
{
    const char *Value;

    Value = Tcl_GetVar(gomp_GetTclInterp(),"lulPostReadFrameTrigger",TCL_GLOBAL_ONLY);

    if(Value) {
        if(TCL_OK != Tcl_UnsetVar(gomp_GetTclInterp(),"lulPostReadFrameTrigger",TCL_GLOBAL_ONLY)) {
            gomp_PrintERROR("can't unset the 'lulPostReadFrameTrigger' variable");
        }
    }

    TrajectoryInfo.inuse        = 0;
    TrajectoryInfo.type         = 0;
    TrajectoryInfo.natom        = 0;
    TrajectoryInfo.nstep        = 0;
    TrajectoryInfo.traj_file[0] = '\0';

    return(0);
}

/*************************************************************************/
int   gomp_GetTrajectoryFirstFrame(void)
/*************************************************************************/
{
    return(TrajectoryInfo.first_frame);
}

/*************************************************************************/
int   gomp_GetTrajectoryLastFrame(void)
/*************************************************************************/
{
    return(TrajectoryInfo.last_frame);
}
/*************************************************************************/
int   gomp_GetTrajectoryDeltaFrame(void)
/*************************************************************************/
{
    return(TrajectoryInfo.delta_frame);
}

/*************************************************************************/
int   gomp_GetNumberOfFreeAtoms(int Set)
/*************************************************************************/
{
    return(FreeAtomList.nfreat);
}
/*************************************************************************/
int   gomp_SetNumberOfFreeAtoms(int Set , int NFreat)
/*************************************************************************/
{
    FreeAtomList.nfreat = NFreat;

    return(0);
}
/*************************************************************************/
/* gomp_GetFreeAtomListPointer                                             */
/*************************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    int,FreeAtomListPointer,(int Set),
    FreeAtomList.free_atom_list,;)
/*************************************************************************/
int gomp_SetFreeAtomListPointer(int Set , int NFreat)
/*************************************************************************/
{
    if(FreeAtomList.nfreat) free(FreeAtomList.free_atom_list);

    FreeAtomList.free_atom_list = gomp_AllocateIntVector(NFreat);

    FreeAtomList.nfreat         = NFreat;

    return(0);
}
/*************************************************************************/
int   gomp_SetNumberOfTrajectoryAtoms(int Atoms)
/*************************************************************************/
{
    TrajectoryInfo.natom = Atoms;

    return(0);
}

/*************************************************************************/
int   gomp_GetNumberOfTrajectoryAtoms(void)
/*************************************************************************/
{
    return(TrajectoryInfo.natom);
}

/*************************************************************************/
int   gomp_GetDisplayFrameNumber(void)
/*************************************************************************/
{
    return(TrajectoryInfo.display);
}
/*************************************************************************/
int   gomp_PutDisplayFrameNumber(int Frame)
/*************************************************************************/
{
    TrajectoryInfo.display = Frame;

    return(0);
}
/*************************************************************************/
int gomp_GetOneFrame(int alt , FILE *File_p , int Append)
/*************************************************************************/
{
    int RetVal;

    switch(gomp_GetTrajectoryFileType()) {
    
    case AMBER_TRAJ:
        RetVal = gomp_GetFrameAmber(alt , File_p , Append);
        break;

    case CHARMM_TRAJ:
        RetVal = gomp_GetFrameCharmm(alt , File_p , Append);
        break;

    case CERIUS2_TRAJ:
        RetVal = gomp_GetFrameCerius2(alt , File_p , Append);
        break;

    case DISCOVER_TRAJ:
        RetVal = gomp_GetFrameDiscover(alt , File_p , Append);
        break;

    case FDL_POLY_TRAJ:
        RetVal = gomp_GetFrameDL_PolyFORMATTED(alt , File_p , Append);
        break;

    case UDL_POLY_TRAJ:
        RetVal = gomp_GetFrameDL_PolyUNFORMATTED(alt , File_p , Append);
        break;

    case FAMBER_TRAJ:
        RetVal = gomp_GetFrameFAmber(alt , File_p , Append);
        break;

    case GROMACS_TRAJ:
        RetVal = gomp_GetFrameGromacs(alt , File_p , Append);
        break;

    case GROMOS_TRAJ:
        RetVal = gomp_GetFrameGromos(alt , File_p , Append);
        break;

    case GROMOS96A_TRAJ:
        RetVal = gomp_GetFrameGROMOS96A(alt , File_p , Append);
        break;

    case HYPERCHEM_TRAJ:
        RetVal = gomp_GetFrameHyperChem(alt , File_p , Append);
        break;

    case MUMOD_TRAJ:
        RetVal = gomp_GetFrameMumod(alt , File_p , Append);
        break;

    case TINKER_TRAJ:
        RetVal = gomp_GetFrameTINKER(alt , File_p , Append);
        break;

    case XMOL_TRAJ:
        RetVal = gomp_GetFrameXmol(alt , File_p , Append);
        break;

    case XPLOR_TRAJ:
        RetVal = gomp_GetFrameXplor(alt , File_p , Append);
        break;

    case YASP_TRAJ:
        RetVal = gomp_GetFrameYasp(alt , File_p , Append);
        break;

    default:
        gomp_PrintMessage(
            "$ERROR - undefined file type when reading the trajectory");
        RetVal = 1;
        break;
    }
    if(alt) {
/* check if the atom connectivity should be recalculated */
        if(gomp_GetBondReconnectivityState()) {
            if(gomp_CalcAtomConn(0))
                return(1);
        }
        if(gomp_GetHBondReconnectivityState()) {
            if(gomp_RecalculateHbonds(0))
                return(1);
        }

        (void)gomp_PeekRunningFrameNumberProperty();
        (void)gomp_PutDisplayFrameNumber(alt);
        (void)gomp_PutRunningFrameNumber(alt);
    }

    return(RetVal);
}

/*************************************************************************/
int gomp_ReadCoordinatesFRAME(int FrameNumber, int Append)
/*************************************************************************/
{
    FILE *fp;
    int   retv;

    if(!gomp_GetTrajectoryStructureState() ||
       gomp_GetTrajectoryFileName()       == '\0') {
        gomp_PrintERROR("no trajectory information available");
        return(1);                                        /* not active */
    }

    if(FrameNumber < 1) {
        gomp_PrintERROR("Defined frame number is < 0 (has to be >= 1)");
        return(1);
    }
    if(FrameNumber > gomp_GetTrajectoryLastFrame()) {
        gomp_PrintERROR("Defined frame number is > max frames");
        return(1);
    }

    if(gomp_GetTrajectoryFileType() == XMOL_TRAJ      ||
       gomp_GetTrajectoryFileType() == GROMOS96A_TRAJ ||
       gomp_GetTrajectoryFileType() == TINKER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FAMBER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FDL_POLY_TRAJ) {
           if(gomp_CheckTextFileLineEnding(gomp_GetTrajectoryFileName())) {
             fp = fopen(gomp_GetTrajectoryFileName() , "rb");
             printf("Hello Leif\n");
             gomp_PrintMessage("it looks as if your trajectory is not a proper text file!\nWill try to solve the problem");
           } else {
             fp = fopen(gomp_GetTrajectoryFileName() , "r");
           }
    } else {
        fp = fopen(gomp_GetTrajectoryFileName() , "rb");
    }

    if(fp == NULL) {
        gomp_PrintERROR("Can't open trajectory file to read a frame");
        return(1);
    }

    retv = gomp_GetOneFrame(FrameNumber , fp , Append);  

    fclose(fp);

    return(retv);
}

/*************************************************************************/
int   gomp_GetTrajectoryDisplayFrames(void)
/*************************************************************************/
{
    return((TrajectoryInfo.last_frame - TrajectoryInfo.first_frame + 1) /
           TrajectoryInfo.delta_frame);
}

/*************************************************************************/
const char *gomp_GetTrajectoryTypeFileExtension(int Type)
/*************************************************************************/
{
    return(TrajectoryFileTypes[Type - 1].file_ext);
}
/*************************************************************************/
int         gomp_PutTrajectoryTypeFileExtension(int Type,
                                                const char *ext)
/*************************************************************************/
{
    gomp_CopyString(TrajectoryFileTypes[Type - 1].file_ext, ext,
        sizeof(TrajectoryFileTypes[0].file_ext));

    return(1);
}
/*************************************************************************/
int         gomp_GetTrajectoryTypeDefault(void)
/*************************************************************************/
{
    return(DefaultTrajectoryFileType);
}
/*************************************************************************/
int         gomp_PutTrajectoryTypeDefault(int Type)
/*************************************************************************/
{
    DefaultTrajectoryFileType = Type;

    return(1);
}
/*************************************************************************/
int   gomp_SetFormattedTrajectoryReader(int Value)
/*************************************************************************/
{
    TrajectoryInfo.retrieve_type = Value;

    return(0);
}
/*************************************************************************/
int   gomp_GetFormattedTrajectoryReader(void)
/*************************************************************************/
{
    return(TrajectoryInfo.retrieve_type);
}

/*************************************************************************/
int   gomp_SetVelocityRetrieveState(int Value)
/*************************************************************************/
{
    TrajectoryInfo.retrieve_velocity = Value;

    return(0);
}
/*************************************************************************/
int   gomp_GetVelocityRetrieveState(void)
/*************************************************************************/
{
    return(TrajectoryInfo.retrieve_velocity);
}
/*************************************************************************/
int   gomp_SetForceRetrieveState(int Value)
/*************************************************************************/
{
    TrajectoryInfo.retrieve_force = Value;

    return(0);
}
/*************************************************************************/
int   gomp_GetForceRetrieveState(void)
/*************************************************************************/
{
    return(TrajectoryInfo.retrieve_force);
}
/*************************************************************************/
int   gomp_DeleteVelocityForceData(void)
/*************************************************************************/
{
    if(TrajectoryInfo.retrieve_velocity_ready) {
        free(TrajectoryInfo.Xvelocities);
        free(TrajectoryInfo.Yvelocities);
        free(TrajectoryInfo.Zvelocities);
        TrajectoryInfo.retrieve_velocity_ready = 0;
    }

    if(TrajectoryInfo.retrieve_force_ready) {
        free(TrajectoryInfo.Xforces);
        free(TrajectoryInfo.Yforces);
        free(TrajectoryInfo.Zforces);
        TrajectoryInfo.retrieve_force_ready = 0;
    }

    return(0);
}
/*************************************************************************/
int   gomp_SetVelocityRetrieveReadyState(int Value)
/*************************************************************************/
{
    TrajectoryInfo.retrieve_velocity_ready = Value;

    return(0);
}
/*************************************************************************/
int   gomp_GetVelocityRetrieveReadyState(void)
/*************************************************************************/
{
    return(TrajectoryInfo.retrieve_velocity_ready);
}
/*************************************************************************/
int   gomp_SetForceRetrieveReadyState(int Value)
/*************************************************************************/
{
    TrajectoryInfo.retrieve_force_ready = Value;

    return(0);
}
/*************************************************************************/
int   gomp_GetForceRetrieveReadyState(void)
/*************************************************************************/
{
    return(TrajectoryInfo.retrieve_force_ready);
}

/*************************************************************************/
int    gomp_GetVelocitySpace(int Atoms)
/*************************************************************************/
{
    if(TrajectoryInfo.retrieve_velocity_ready) {
        free(TrajectoryInfo.Xvelocities);
        free(TrajectoryInfo.Yvelocities);
        free(TrajectoryInfo.Zvelocities);
        TrajectoryInfo.retrieve_velocity_ready = 0;
    }

    TrajectoryInfo.Xvelocities = gomp_AllocateFloatVector(Atoms);
    TrajectoryInfo.Yvelocities = gomp_AllocateFloatVector(Atoms);
    TrajectoryInfo.Zvelocities = gomp_AllocateFloatVector(Atoms);

    if((TrajectoryInfo.Xvelocities == (const float *)NULL) ||
       (TrajectoryInfo.Yvelocities == (const float *)NULL) ||
       (TrajectoryInfo.Zvelocities == (const float *)NULL)) 
        return(1);

    memset(TrajectoryInfo.Xvelocities , 0 , sizeof(float) * Atoms);
    memset(TrajectoryInfo.Yvelocities , 0 , sizeof(float) * Atoms);
    memset(TrajectoryInfo.Zvelocities , 0 , sizeof(float) * Atoms);

    TrajectoryInfo.retrieve_velocity_ready = 1;

    return(0);
}
/*************************************************************************/
/* gomp_GetVelocityXComponentPointer,                                    */
/* gomp_GetModifiableVelocityXComponentPointer                           */
/*************************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    float,VelocityXComponentPointer,(void),
    TrajectoryInfo.retrieve_velocity_ready ?
    TrajectoryInfo.Xvelocities : NULL,;)
/*************************************************************************/
/* gomp_GetVelocityYComponentPointer,                                    */
/* gomp_GetModifiableVelocityYComponentPointer                           */
/*************************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    float,VelocityYComponentPointer,(void),
    TrajectoryInfo.retrieve_velocity_ready ?
    TrajectoryInfo.Yvelocities : NULL,;)
/*************************************************************************/
/* gomp_GetVelocityXComponentPointer,                                    */
/* gomp_GetModifiableVelocityXComponentPointer                           */
/*************************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    float,VelocityZComponentPointer,(void),
    TrajectoryInfo.retrieve_velocity_ready ?
    TrajectoryInfo.Zvelocities : NULL,;)
/*************************************************************************/
int    gomp_GetForceSpace(int Atoms)
/*************************************************************************/
{

    if(TrajectoryInfo.retrieve_force_ready) {
        free(TrajectoryInfo.Xforces);
        free(TrajectoryInfo.Yforces);
        free(TrajectoryInfo.Zforces);
        TrajectoryInfo.retrieve_force_ready = 0;
    }

    TrajectoryInfo.Xforces = gomp_AllocateFloatVector(Atoms);
    TrajectoryInfo.Yforces = gomp_AllocateFloatVector(Atoms);
    TrajectoryInfo.Zforces = gomp_AllocateFloatVector(Atoms);

    if((TrajectoryInfo.Xforces == (const float *)NULL) ||
       (TrajectoryInfo.Yforces == (const float *)NULL) ||
       (TrajectoryInfo.Zforces == (const float *)NULL)) 
        return(1);

    memset(TrajectoryInfo.Xforces , 0 , sizeof(float) * Atoms);
    memset(TrajectoryInfo.Yforces , 0 , sizeof(float) * Atoms);
    memset(TrajectoryInfo.Zforces , 0 , sizeof(float) * Atoms);

    TrajectoryInfo.retrieve_force_ready = 1;

    return(0);
}
/*************************************************************************/
/* gomp_GetForceXComponenPointer,                                        */
/* gomp_GetModifiableForceXComponentPointer                              */
/*************************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    float,ForceXComponentPointer,(void),
    TrajectoryInfo.retrieve_force_ready ? TrajectoryInfo.Xforces : NULL,;)
/*************************************************************************/
/* gomp_GetForceYComponenPointer,                                        */
/* gomp_GetModifiableForceYComponentPointer                              */
/*************************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    float,ForceYComponentPointer,(void),
    TrajectoryInfo.retrieve_force_ready ? TrajectoryInfo.Yforces : NULL,;)
/*************************************************************************/
/* gomp_GetForceZComponenPointer,                                        */
/* gomp_GetModifiableForceZComponentPointer                              */
/*************************************************************************/
IMPLEMENT_GET_MODIFIABLE_POINTER(
    float,ForceZComponentPointer,(void),
    TrajectoryInfo.retrieve_force_ready ? TrajectoryInfo.Zforces : NULL,;)
/*************************************************************************/
int    gomp_SetVelocityMinMax(float Vmin , float Vmax)
/*************************************************************************/
{
    TrajectoryInfo.Vmin = Vmin;
    TrajectoryInfo.Vmax = Vmax;

    return(0);
}
/*************************************************************************/
int    gomp_GetVelocityMinMax(float *Vmin, float *Vmax)
/*************************************************************************/
{
    *Vmin = TrajectoryInfo.Vmin;
    *Vmax = TrajectoryInfo.Vmax;

    return(0);
}
/*************************************************************************/
int    gomp_SetForceMinMax(float  Fmin, float Fmax)
/*************************************************************************/
{
    TrajectoryInfo.Fmin = Fmin;
    TrajectoryInfo.Fmax = Fmax;

    return(0);
}
/*************************************************************************/
int    gomp_GetForceMinMax(float *Fmin, float *Fmax)
/*************************************************************************/
{
    *Fmin = TrajectoryInfo.Fmin;
    *Fmax = TrajectoryInfo.Fmax;

    return(0);
}
/*************************************************************************/
int    gomp_CalculateVelocityMinMax()
/*************************************************************************/
{
    float Ttemp;
    float VMax;
    float VMin;
    int   i;

    VMin = 1.0E+20f;
    VMax = 0.0f;

    if(!TrajectoryInfo.retrieve_velocity_ready) 
        return(1);

    for(i = 0 ; i < TrajectoryInfo.natom ; i++) {

        Ttemp = sqrt(
            TrajectoryInfo.Xvelocities[i] * TrajectoryInfo.Xvelocities[i] +
            TrajectoryInfo.Yvelocities[i] * TrajectoryInfo.Yvelocities[i] +
            TrajectoryInfo.Zvelocities[i] * TrajectoryInfo.Zvelocities[i]);

        if(Ttemp < VMin) {
            VMin = Ttemp;
        } else if(Ttemp > VMax) {
            VMax = Ttemp;
        }
    }

    if(gomp_SetVelocityMinMax(VMin, VMax))
        return(1);

    return(0);
}
/*************************************************************************/
int    gomp_CalculateForceMinMax()
/*************************************************************************/
{
    float Ttemp;
    float FMax;
    float FMin;
    int   i;

    FMin = 1.0E+20f;
    FMax = 0.0f;

    if(!TrajectoryInfo.retrieve_force_ready) 
        return(1);

    for(i = 0 ; i < TrajectoryInfo.natom ; i++) {

        Ttemp = sqrt(TrajectoryInfo.Xforces[i] * TrajectoryInfo.Xforces[i] +
                     TrajectoryInfo.Yforces[i] * TrajectoryInfo.Yforces[i] +
                     TrajectoryInfo.Zforces[i] * TrajectoryInfo.Zforces[i]);

        if(Ttemp < FMin) {
            FMin = Ttemp;
        } else if(Ttemp > FMax) {
            FMax = Ttemp;
        }
    }

    if(gomp_SetForceMinMax(FMin, FMax))
        return(1);

    return(0);
}

