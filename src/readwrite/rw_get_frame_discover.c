/*
  Copyright (c) 1991 - 2005 by:
  Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
  Confidential unpublished property of 
  Leif Laaksonen
  All rights reserved

  Enhancements 2003, 2005 by:
  Eero Häkkinen

  Based on code from BIOSYM Technologies Inc. 
  and hacked further to c by
  Leif Laaksonen Centre for Scientific Computing 1991 - 2003

*/

#include "maindefs.h"

#include "gomstdio.h"

#include "gomendian.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"

#include "stdafx.h"

/* define some useful things */

#define MAX_DISCM  10000      /* max Discover molecules */

#define RECORD()   { record = \
                   fread(&record, sizeof(int) , 1L , File_p);\
                   if(record < 1) {\
                     gomp_PrintMessage("?ERROR1 - in reading trajectory file");\
                     return(1);}}

#define FREAD(value_p , size)    { record = \
                                 fread(value_p, size , 1L , File_p);\
                                 if(record < 1) {\
                     gomp_PrintMessage("?ERROR2 - in reading trajectory file");\
                     return(1);}}

#define FREADN(value_p , num , size) { record = \
                                fread(value_p, size , num , File_p);\
                   if(record != num) {\
                     gomp_PrintMessage("?ERROR3 - in reading trajectory file");\
                     return(1);}}

#define NAME_LEN4                4   /* name is 4 bytes long version < 2.9*/
#define NAME_LEN5                5   /* name is 5 bytes long version >= 2.9*/

#define DISC_MAX_ATOMS 50000        /* max atoms in DISCOVER */

/* Common Block Declarations */

static struct {
    int nat;
    char at[800];
} atdat_;

#define atdat_1 atdat_

static struct {
    int *knd;
    char *jname;
} kindat_;

#define kindat_1 kindat_

static struct {
    double cmprsb, demax, dseed, dtemp, dtoten, pressb, temp, tempb, 
        temp0, timprs, timtmp, tmstpi, tmstpo, toten, vscale, vershn;
    int iprint, iread;
    char ititlr[80], iver[80], jtitlr[80];
    int katpfc, katscl, kavpr, kcnstr, kexact, kinit, kinitv, kpresb, 
        kprint, krsclv, krstfc, kshape, kstppr, ktempb, ktsscl, kvrscl, 
        mlhis, lhis, hiscnt, loprst, nclnan, nrlnan, nsteps, nstpav, 
        nwrstp;
} mdopt_;

#define mdopt_1 mdopt_

#ifdef TRAJ_ENERGY
static struct {
    double etotl, eb, et, ep, eop, ebb, ebt, ett, ebp, etph, ettp, ebb2, 
        eoo, elj, erp, edp, est, ehb, dobelj, doberp, dobedp, dobest;
} energ_;

#define energ_1 energ_

static struct {
    double etotml[MAX_DISCM];
    double ebmol[MAX_DISCM];
    double etmol[MAX_DISCM];
    double epmol[MAX_DISCM];
    double eopmol[MAX_DISCM];
    double ebbmol[MAX_DISCM];
    double ebtmol[MAX_DISCM];
    double ettmol[MAX_DISCM];
    double ettpml[MAX_DISCM];
    double eoomol[MAX_DISCM];
    double eljmol[MAX_DISCM];
    double erpmol[MAX_DISCM];
    double edpmol[MAX_DISCM];
    double estmol[MAX_DISCM];
    double ebpmol[MAX_DISCM];
    double ehbmol[MAX_DISCM];
    double etpmol[MAX_DISCM];
    double ebb2ml[MAX_DISCM];
} enmol_;

#define enmol_1 enmol_
#endif

static struct {
    double atemp, etk, etot, etp;
    int nenrg;
} mdeng_;

#define mdeng_1 mdeng_

static struct {
    int *ibw;   /* was [2][11400] */
} intdat_;

#define intdat_1 intdat_

#ifdef TRAJ_ENERGY
static struct {
    double erpmlr[MAX_DISCM];
    double edpmlr[MAX_DISCM];
    double eljmlr[MAX_DISCM];
    double estmlr[MAX_DISCM];
    double ehbmlr[MAX_DISCM];
} interm_;

#define interm_1 interm_
#endif

static struct {
    double ctemp, ctmpav, tpresm, tpress, tprmav, tprsav, uenavs[17], 
        uens[17], utk, utkav, utot, utotav, utp, utpav;
    int isteps, istpav, jstfor;
} mdavs_;

#define mdavs_1 mdavs_

#ifdef TRAJ_ENERGY
static struct {
    double presur, etkv[9]  /* was [3][3] */, virial[9] /* was [3][3] 
                                                         */, presst[9]  /* was [3][3] */, presrm, etkm[9]   /* 
                                                                                                               was [3][3] */, virilm[9] /* was [3][3] */, prestm[9] /* 
                                                                                                                                                                       was [3][3] */;
} mdpres_;

#define mdpres_1 mdpres_
#endif

#ifdef DISCOALL
static struct {
    double atmas[200], timcnv, fconv, tmpcnv, tmpacv, econv, ekconv, 
        ekcnvt, prscnv, avono, boltz, gaskcn, aclcnv, velcnv, ratmas[200],
        timstp, h2d2, rh2d2;
    int nlinep;
} mdkons_;

#define mdkons_1 mdkons_
#endif

static struct {
    int molxyz[MAX_DISCM];
    int natmol[MAX_DISCM];
    int molres[MAX_DISCM];
    int nrsmol[MAX_DISCM]; 
    int molbp[MAX_DISCM];
    int moltp[MAX_DISCM];
    int molpp[MAX_DISCM];
    int molop[MAX_DISCM];
    int molttp[MAX_DISCM];
    int molgrp[MAX_DISCM];
    int molbbp[MAX_DISCM];
} molptr_;

#define molptr_1 molptr_

static struct {
    double bufzon;
    int nmove, natmov, natmv3, nat32m, nbmv, nthmv, nphmv, nopmv, ntthmv, 
        nbbmv, imove[DISC_MAX_ATOMS], lmove[DISC_MAX_ATOMS];
} movinf_;

#define movinf_1 movinf_

static struct {
    int nmol, nres, natom, nbond, ntheta, nphi, nopln, nthth, nbb, nint, 
        nbnded, nat3, ihydnb, ngroup;
} numinf_;

#define numinf_1 numinf_

static struct {
    int nnmres;
    char namres[400];
} resdat_;

#define resdat_1 resdat_

static struct {
    int lisres[MAX_DISCM];
    int ifirst[MAX_DISCM];
    int last[MAX_DISCM];
    int mresml[MAX_DISCM];
    int irsgpf[MAX_DISCM];
    int irsgpl[MAX_DISCM];
} resinf_;

#define resinf_1 resinf_

static struct {
    double aboxa, aboxb, aboxc, boxlim, boxmnx, boxmny, boxmnz, boxmxx, 
        boxmxy, boxmxz, carcry[9]   /* was [3][3] */, crycar[9] /* 
                                                                   was [3][3] */, symcr[1764]   /* was [3][3][196] */, symops[1764] 
        /* was [3][3][196] */, tranop[588]  /* was [3][196] */, ucvcar[9]   
        /* was [3][3] */, ucvec[6], volum;
    int ncellx, ncelly, ncellz, ncelsz, ncelx2, ncely2, 
        ncelz2, nclszy, nsymcd, nsymm2, nsymop;
} symdat_;

#define symdat_1 symdat_

static struct {
    float *xcoor;
    float *ycoor;
    float *zcoor;
    float *xvlcty;
    float *yvlcty; 
    float *zvlcty;
} tempry_;

#define tempry_1 tempry_


/***************************************************************************/
int gomp_GetFrameDiscover(int alt, FILE *File_p , int iappend)  
    /* read one frame from a insight trajectory 
       mode of operation (=0) first time in read , 
       take just trajectory info
       (>0) trajectory number 
       if = 0 no append , if = 1 append */
/***************************************************************************/
{
    /* System generated locals */
    int i_1, i_2;
#ifdef TRAJ_ENERGY
    int i_3, i_4, i_5, i_6, i_7, i_8, i_9;
#endif

    /* Table of constant values */

    static int c__6 = 6;
    static int c__9 = 9;
    static int c__1764 = 1764;
    static int c__588 = 588;

    /* Local variables */
    static int ibond, nsymsm;
#define enrgys ((const double  *)&energ_1)
#define enmols ((const double  *)&enmol_1 + MAX_DISCM)
#ifdef TRAJ_ENERGY
    static int ien, iml;
#endif
    static int iat, ires;
    static int record,kntrl;
    static int FixedAtoms;

    static char OutText[BUFF_LEN];
    static struct {
        long WhereInFile;
        long JumpPoint;
        long RetSeek;
    } Discover;

    static long RecordL;
    static long RecordL1;
    static long RecordL2;
    static long RecordA;

    static float *xtemp;
    static float *ytemp;
    static float *ztemp;
    static float *x;
    static float *y;
    static float *z;
    static const float *sumxyz;
    static int    numat;
    static int    dynamics_frames;
    static int    Wstr;
    static int    EndOfFirstFrame;
    static int    JumpLength;
    static int    swap_bytes;

/* Function */
/*       Read one frame from an unformatted history file written by */
/*       DISCOVER version 2.2 or later. */

/*  Description of items in his file                      number of bytes */
/*------------------------------------------------------------------------
  kntrl  = 0 for first record, 1 for all others           4   
  iver   = program header                                 80   
  vershn = Discover version number                        8   
  jtitlr = run title                                      80   
  ititlr = run subtitle                                   80   
  nat    = number of atom types                           4   
  at     = names of atom types                            4*nat   
  atmas  = masses of atom types                           8*nat   
  nnmres = number of unique residues in system            4   
  namres = names of residues included in system           4*nnmres   
  natom  = number of atoms in system                      4   
  knd    = pointer to atom type for each atom             4*natom   
  jname  = unique name for each atom                      4*natom   
  nmove  =                                                4   
  natmov = number of moving atoms in system               4   
  nmol   = number of molecules in system                  4   
  natmol = number of atoms in each molecule               4*nmol   
  nrsmol = number of residues in each molecule            4*nmol   
  nres   = number of residues in system                   4   
  ifirst = pointer to first atom of each residue in knd   4*nres   
  last   = pointer to last atom of each residue in knd    4*nres   
  lisres = pointer to residue name in namres              4*nres   
  nbond  = number of bonds                                4   
  ibw    = pointers to bonded atoms in knd                8*nbond   
  ucvec  = unit cell vectors                              72   
  ucvcar = cartesian components of unit cell vectors      72   
  crycar = matrix to convert cryst. to cart. coordinates  72   
  carcry = matrix to convert cart. to cryst. coordinates  72   
  symops = symmetry operations                            72*mxsym   
  tranop = translational part of symmetry operations      24*mxsym   
  symcr  = symmetry operations for cartesian coordinates  72*mxsym   
  nsymop = number of symmetry operations                  4   
  aboxa, aboxb, aboxc                                     24   
  boxlim =   
  nsymsm =   
  boxmnx =   
  boxmny =   
  boxmnz =   
  ncellx = # of unit cells in one direction along x axis  4   
  ncelly = # of unit cells in one direction along y axis  4   
  ncellz = # of unit cells in one direction along z axis  4   
  ncelx2 = 2*ncellx + 1                                   4   
  ncely2 = 2*ncelly + 1                                   4   
  ncelz2 = 2*ncellz + 1                                   4   
  nenrg  = number of energies                             4   
  tmstpi = time per step in femtoseconds                  8   
  nwrstp = number of steps before writing history file    4   
  isteps =                                                4   
  etot   = total energy                                   8   
  etp    = total potential energy                         8   
  etk    = total kinetic energy                           8   
  enrgys = system energies (bond, angle, torsion, etc.)   8*nenrg   
  etotml = total energy in each molecule                  8*nmol   
  enmols = molecular energies                             8*nmol*nenrg   

  edpmlr = dipole energy in each molecule                 8*nmol   
  erpmlr = repulsion energy in each molecule              8*nmol   
  eljmlr = Lennard-Jones energy in each molecule          8*nmol   
  estmlr = electrostatic energy in each molecule          8*nmol   
  presur =                                                8   
  presrm =                                                8   
  presst =                                                72   
  prestm =                                                72   
  etkv   =                                                72   
  etkm   =                                                72   
  virial =                                                72   
  virilm =                                                72   
  xcoor  = x,y,z coordinates (in angstroms) of each atom  24*natom   
  velcty = x,y,z components of velocities for each atom   24*natmov   
*/
/* ====================================================================== 
 */
/*  To increase (or decrease) the capacity of DISCOVER to handle larger */
/*  or smaller systems, one need only change MXATM, MXRES, and MXMOL as */
/*  desired and the rest of the program will adjust as necessary. */
/*  Remember that any changes to this file will, in general, require */
/*  recompilation of the entire program. */
/*  How the dimensions of the various internal coordinate arrays scale */
/*  with system size is determined by experience with typical systems. */
/*  For example, the number of bonds in a polypeptide is roughly 1.02 */
/*  times the number of atoms.  To allow a safety margin, this observed */
/*  ratio is increased slightly and then used to scale MXBND (e.g., */
/*  MXBND = 1.2 * MXATM).  The empirical ratios found in large protein */
/*  systems along with the actual factors used to determine the values */
/*  found in param.ins are given below.  These gomp_aling factors are */
/*  used explicitly in the parameter statements below to automatically */
/*  scale the appropriate parameters as MXATM and MXRES are changed. */

/*                                             ave #/atom   #/atom ratio 
 */
/*                                               found in  used to scale 
 */
/*        parameter   meaning                    proteins     parameters 
 */
/*        ---------   ---------------------    ----------    ---------- */

/*            mxbnd   maximum # of bonds            1.02          1.20 */
/*           mxthet   maximum # of angles           1.83          2.00 */
/*            mxphi   maximum # of torsions         2.68          3.00 */
/*           mxopln   maximum # of out-of-planes    0.20          0.30 */
/*           mxthth   maximum # of angle-angles     3.03          3.50 */
/*            mxecl   maximum # of excluded atoms  10.00         13.00 */
/* ---MXSYM s the maximum number of symmetry operators permitted */
/* ---ipulmx is the maximum allowable number of 'pulled' atom pairs */
/*     CTLOPT contains primarily control flags. */
/*     added 6/86: PRMPTR - indicates if there are valid pointers to the 
 */
/*                   parameters in the CPV array.  This is used to avoid 
 */
/*                   redundant and potentially time consuming re */
/*                   generation of internals/pointer in MOLIN/GENTOP. */
/*     added 6/86: NONEIB - indicates if there is "NO NEIghBOR list" for 
 */
/*                   the current conformation.  Again, this is to save */
/*                   time (generation can be time consuming) but need */
/*                   not be done unless energy calculations or print */
/*                   of the list is needed. */
/*    added 9/86: PAXFIT,PAXREF - indicate if the FITted or REFerence */
/*                   molecules have been rotated to a Principal Axis */
/*                   coordinate system. */
/*    added 1/87: KASYM,USEGRP - kASYMmetric unit indicator for symmetry 
 */
/*                          USEGRP indicates if "charge groups" are used 
 */

/*     BIOVER is the version of Discover, set in BLKDAT. */
/*     FILVER is the version read in from topology/restart file. */
/*     CVFVER is the version read in from CVFF forcefield file. */
/*     SPACGP is the name of the space group (for PBC/symmetry) */

    rewind(File_p);

    if(!alt) {

/* #1 */
/*  START WITH THE CONTROL RECORD      */
        FREAD(&kntrl, sizeof(int));
/* check for right type of trajectory
   record has to be = (4 chars + 20 * sizeof(int))!
*/
        swap_bytes = 0;
        if(kntrl != sizeof(int)) {
            gomp_Reverse_int( & kntrl );
            if(kntrl != sizeof(int)) {
                gomp_PrintERROR("wrong internal structure of trajectory file");
                return(1);
            } else {
                gomp_PrintMessage("Enabling automatic byte_swapping...");
                swap_bytes = 1;
            }
        }

        FREAD(&kntrl, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( & kntrl );
        }
        RECORD();


/* #2 */
/*     DISCOVER HEADER AND VERSION NUMBER */
        RECORD();
        FREAD(mdopt_1.iver , (long)80 );

        sprintf(OutText,"Program header: '%s'",mdopt_1.iver);
        gomp_PrintMessage(OutText);

        FREAD(&mdopt_1.vershn, sizeof(double));
        if ( swap_bytes ) {
            gomp_Reverse_double( &mdopt_1.vershn );
        }

        sprintf(OutText,"Discover version number:        %f",mdopt_1.vershn);
        gomp_PrintMessage(OutText);

        RECORD();

        mdopt_1.lhis = 20;


/* #3 */
/*     READ INITIAL RECORD */
        RECORD();
        FREAD(mdopt_1.jtitlr, (long)80 );
        mdopt_1.jtitlr[79] = '\0';  /* null terminate string */

        gomp_PrintMessage("Run Title:");
        gomp_PrintMessage(mdopt_1.jtitlr);

        RECORD();

/* #4 */
        RECORD();
        FREAD(mdopt_1.ititlr, (long)80 );
        mdopt_1.ititlr[79] = '\0'; /* null terminate string */

        gomp_PrintMessage("Run Subtitle:");
        gomp_PrintMessage(mdopt_1.ititlr);

        RECORD();

/* #5 */
        RECORD();

        FREAD(&atdat_1.nat, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &atdat_1.nat );
        }

        i_1 = atdat_1.nat;

#ifdef DISCOALL

        FREADN(&atdat_1.at[0], i_1 , (long)4);
        if ( swap_bytes ) {
            gomp_Reverse_int( &atdat_1.at[0] );
        } 

        i_2 = atdat_1.nat;

        FREADN(&mdkons_1.atmas[0], i_2 , sizeof(double));
        if ( swap_bytes ) {
            gomp_Reverse_double( &mdkons_1.atmas[0] );
        } 
#else

        i_2 = i_1 * (sizeof(int) + sizeof(double));

        Discover.RetSeek = fseek(File_p , i_2 , SEEK_CUR);
        if(Discover.RetSeek) {
            sprintf(OutText,"?ERROR - can't read trajectory file");
            gomp_PrintMessage(OutText);
            return(1);
        }

#endif

        RECORD();

/* #6 */
        RECORD();
        FREAD(&resdat_1.nnmres, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &resdat_1.nnmres );
        } 
        i_1 = resdat_1.nnmres;

#ifdef DISCOALL

        FREADN(&resdat_1.namres[0] , i_1 , (long)4);
        if ( swap_bytes ) {
            gomp_Reverse_int( &resdat_1.namres[0] );
        } 

#else
        i_2 = i_1 * sizeof(int);

        Discover.RetSeek = fseek(File_p , i_2 , SEEK_CUR);
        if(Discover.RetSeek) {
            sprintf(OutText,"?ERROR - can't read trajectory file");
            gomp_PrintMessage(OutText);
            return(1);
        }

#endif

        RECORD();

/* #7 */
        RECORD();
        FREAD(&numinf_1.natom, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &numinf_1.natom );
        } 

        if(gomp_GetNumAtomsInMolecStruct(0) != numinf_1.natom) {
            gomp_PrintMessage("?ERROR - natom != DISCOVER trajectory natom");
            return(1);
        }

        if(numinf_1.natom >= DISC_MAX_ATOMS) {
            gomp_PrintMessage("?ERROR - too many atoms in DISCOVER trajectory ");
            return(1);
        }

        sprintf(OutText,"Number of atoms:               %d",numinf_1.natom);
        gomp_PrintMessage(OutText);

        i_1 = numinf_1.natom;

/* get space ... */
        kindat_1.knd = gomp_AllocateIntVector(i_1);

        FREADN(&kindat_1.knd[0] , i_1 , sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &kindat_1.knd[0] );
        } 

        i_2 = numinf_1.natom;

/* get space ... */
        if(mdopt_1.vershn < 2.9) {
            kindat_1.jname = gomp_AllocateCharVector(i_2 * NAME_LEN4);
            FREADN(&kindat_1.jname[0] , i_2 , (long)NAME_LEN4);
        }
        else {
            kindat_1.jname = gomp_AllocateCharVector(i_2 * NAME_LEN5);
            FREADN(&kindat_1.jname[0] , i_2 , (long)NAME_LEN5);
        }

        RECORD();

/* #8 */
        if (mdopt_1.vershn < 2.6) {
            RECORD();
            FREAD(&movinf_1.nmove  , sizeof(int));
            if ( swap_bytes ) {
                gomp_Reverse_int( &movinf_1.nmove );
            } 
            FREAD(&movinf_1.natmov , sizeof(int));
            RECORD();
        } else {
            RECORD();
            FREAD(&movinf_1.nmove  , sizeof(int));
            if ( swap_bytes ) {
                gomp_Reverse_int( &movinf_1.nmove );
            } 
            FREAD(&movinf_1.natmov , sizeof(int));
            if ( swap_bytes ) {
                gomp_Reverse_int( &movinf_1.natmov );
            } 
            i_1 = movinf_1.natmov;

            FREADN(&movinf_1.lmove[0], i_1 , sizeof(int));
            if ( swap_bytes ) {
                gomp_Reverse_int_array( &movinf_1.lmove[0] , i_1 );
            } 

            RECORD();
        }

        FixedAtoms = numinf_1.natom - movinf_1.natmov ;
        numat      = numinf_1.natom;

        sprintf(OutText,"Number of moving atoms:         %d",movinf_1.natmov);
        gomp_PrintMessage(OutText);
        sprintf(OutText,"Number of fixed atoms:          %d",FixedAtoms);
        gomp_PrintMessage(OutText);

/* #9 */
        RECORD();
        FREAD(&numinf_1.nmol, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &numinf_1.nmol );
        } 
        i_1 = numinf_1.nmol;

        if(i_1 >= MAX_DISCM) { /* too many molecules */
            gomp_PrintMessage("?ERROR - too many DISCOVER molecules ");
/* free gomp_ratch space ... */
            gomp_FreeVector(kindat_1.knd);
            gomp_FreeVector(kindat_1.jname);
            return(1);
        }

        FREADN(&molptr_1.natmol[0], i_1 , sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int_array( &molptr_1.natmol[0] , i_1 );
        } 

        i_2 = numinf_1.nmol;

        FREADN(&molptr_1.nrsmol[0], i_2 , sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int_array( &molptr_1.nrsmol[0] , i_2 );
        } 

        RECORD();

/* #10 */
        RECORD();
        FREAD(&numinf_1.nres, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &numinf_1.nres );
        } 
        i_1 = numinf_1.nres;
        for (ires = 0; ires < i_1; ++ires) {
            FREAD(&resinf_1.ifirst[ires], sizeof(int));
            if ( swap_bytes ) {
                gomp_Reverse_int( &resinf_1.ifirst[ires] );
            } 
            FREAD(&resinf_1.last[ires], sizeof(int));
            if ( swap_bytes ) {
                gomp_Reverse_int( &resinf_1.last[ires] );
            } 
        }
        i_2 = numinf_1.nres;

        FREADN(&resinf_1.lisres[0], i_2 , sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int_array( &resinf_1.lisres[0] , i_2 );
        } 

        RECORD();

/* #11 */
        RECORD();
        FREAD(&numinf_1.nbond, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &numinf_1.nbond );
        } 
        i_1 = numinf_1.nbond;
        intdat_1.ibw = gomp_AllocateIntVector(2 * i_1);

/* here are the latest changes from BIOSYM (1992-09-21) */

        if (mdopt_1.vershn > 2.7) {
            RECORD();
            RECORD();
        }

        for (ibond = 0; ibond < i_1; ++ibond) {
            for (iat = 0; iat < 2; ++iat) {
                FREAD(&intdat_1.ibw[iat + 2 * ibond], sizeof(int));
                if ( swap_bytes ) {
                    gomp_Reverse_int( &intdat_1.ibw[iat + 2 * ibond] );
                } 
            }
        }
        gomp_FreeVector(intdat_1.ibw);
        RECORD();

/* #12 */
        RECORD();
        FREADN(&symdat_1.ucvec[0] , c__6 , sizeof(double));
        FREADN(&symdat_1.ucvcar[0], c__9 , sizeof(double));
        FREADN(&symdat_1.crycar[0], c__9 , sizeof(double));
        FREADN(&symdat_1.carcry[0], c__9 , sizeof(double));
        FREADN(&symdat_1.symops[0], c__1764 , sizeof(double));
        FREADN(&symdat_1.tranop[0], c__588 , sizeof(double));
        FREADN(&symdat_1.symcr[0] ,  c__1764 ,sizeof(double));
        FREAD(&symdat_1.nsymop, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &symdat_1.nsymop );
        } 
        FREAD(&symdat_1.aboxa ,  sizeof(double));
        FREAD(&symdat_1.aboxb ,  sizeof(double));
        FREAD(&symdat_1.aboxc ,  sizeof(double));
        FREAD(&symdat_1.boxlim, sizeof(double));
        FREAD(&nsymsm         , sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &nsymsm );
        } 
        FREAD(&symdat_1.boxmnx, sizeof(double));
        FREAD(&symdat_1.boxmxx, sizeof(double));
        FREAD(&symdat_1.boxmny, sizeof(double));
        FREAD(&symdat_1.boxmxy, sizeof(double));
        FREAD(&symdat_1.boxmnz, sizeof(double));
        FREAD(&symdat_1.boxmxz, sizeof(double));
        FREAD(&symdat_1.ncellx, sizeof(int));
        FREAD(&symdat_1.ncelly, sizeof(int));
        FREAD(&symdat_1.ncellz, sizeof(int));
        FREAD(&symdat_1.ncelx2, sizeof(int));
        FREAD(&symdat_1.ncely2, sizeof(int));
        FREAD(&symdat_1.ncelz2, sizeof(int));
        RECORD();

/* #13 */
        RECORD();
        FREAD(&mdeng_1.nenrg , sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &mdeng_1.nenrg );
        } 
#ifdef DEBUG
        printf("Number of energies %d\n",mdeng_1.nenrg);
#endif
        FREAD(&mdopt_1.tmstpi, sizeof(double));
        if ( swap_bytes ) {
            gomp_Reverse_double( &mdopt_1.tmstpi );
        } 
        FREAD(&mdopt_1.nwrstp, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &mdopt_1.nwrstp );
        } 
        FREAD(&mdavs_1.isteps, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &mdavs_1.isteps );
        } 
        RECORD();

/* Ok the header field is now done ... */

/* calculate the length of record #14 */
        RecordL = 
            sizeof(int)                    +
            3 * sizeof(double)                 +
            mdeng_1.nenrg * sizeof(double)                 +
            numinf_1.nmol * sizeof(double)                 +
            mdeng_1.nenrg * numinf_1.nmol * sizeof(double) +
            4 * numinf_1.nmol * sizeof(double) +
            2 * sizeof(double)                 +
            54 * sizeof(double)                 +
            sizeof(int);

        RecordA =      sizeof(int)               +
            6  * sizeof(double)            +
            9  * sizeof(double)            +
            sizeof(int);

        RecordL1 = 3 * numinf_1.natom  * sizeof(float) + 2 * sizeof(int);

        RecordL2 = 3 * movinf_1.natmov * sizeof(float) + 2 * sizeof(int);

        EndOfFirstFrame      = ftell(File_p)   + RecordL + RecordL1 + RecordL1;
        Discover.WhereInFile = ftell(File_p)   + RecordL;
        JumpLength           = 3 * sizeof(int) + RecordL + RecordA  +
            RecordL2 + RecordL2; 

        /* no need to  do more, I was just interested in headers */
        /* just look for the number of frames                    */
        Discover.RetSeek = fseek(File_p,0L,SEEK_END);
        if(Discover.RetSeek) {
            gomp_PrintERROR("Can't read trajectory file");
            return(1);
        }

        dynamics_frames = (ftell(File_p) - EndOfFirstFrame) / JumpLength + 1;

        sprintf(OutText,"Number of frames:              %d",dynamics_frames);
        gomp_PrintMessage(OutText);

/* free gomp_ratch space ... */

        gomp_FreeVector(kindat_1.knd);
        gomp_FreeVector(kindat_1.jname);

/* update trajectory info ... */
        (void)gomp_SetNumberOfTrajectoryAtoms(numat);
        (void)gomp_SetNumberOfFrames(dynamics_frames);
        (void)gomp_SetTrajectoryTimeInfo((int) mdopt_1.tmstpi , 0);
        (void)gomp_SetNumberOfFreeAtoms(0 , movinf_1.natmov);
        (void)gomp_SetTrajectoryDisplayParams(1 , dynamics_frames , 1);
        (void)gomp_PutDisplayFrameNumber(1);

        return(0);
}


    if(alt == 1) { /* first frame is special */
        Discover.RetSeek = fseek(File_p,Discover.WhereInFile,SEEK_CUR);
        if(Discover.RetSeek) {
            gomp_PrintERROR("Can't read trajectory file");
            return(1);
        }
    }
    else {        
        /* for all other records there is first the control 
           block */

        Discover.RetSeek = fseek(File_p,EndOfFirstFrame,SEEK_CUR);
        if(Discover.RetSeek) {
            gomp_PrintERROR("Can't read trajectory file");
            return(1);
        }

        Discover.RetSeek = fseek(File_p,(alt - 2) * JumpLength,SEEK_CUR);
        if(Discover.RetSeek) {
            printf("?ERROR - can't read trajectory file");
            return(1);
        }

      
        RECORD();
        FREAD(&kntrl, sizeof(int));
        if ( swap_bytes ) {
            gomp_Reverse_int( &kntrl );
        } 
/*        if(mdopt_1.vershn > 2.81) { */
        if(kntrl < 1) {
            gomp_PrintMessage("?ERROR - in control record of Discover file ");
            return(1);
/*          } */
        }
        RECORD();

        Discover.RetSeek = fseek(File_p,(RecordL + RecordA),SEEK_CUR);
        if(Discover.RetSeek) {
            gomp_PrintERROR("Can't read trajectory file");
            return(1);
        }
    }
#ifdef TRAJ_ENERGY

    RECORD();
    FREAD(&mdeng_1.etot, sizeof(double));
    FREAD(&mdeng_1.etp , sizeof(double));
    FREAD(&mdeng_1.etk , sizeof(double));

    i_1 = mdeng_1.nenrg;
    FREADN(&enrgys[0], i_1 , sizeof(double));
    i_2 = numinf_1.nmol;
    FREADN(&enmol_1.etotml[0], i_2 , sizeof(double));

    i_3 = mdeng_1.nenrg;
    for (ien = 0; ien < i_3; ++ien) {
        i_4 = numinf_1.nmol;
        for (iml = 0; iml < i_4; ++iml) {
            FREAD(&enmols[iml + ien * MAX_DISCM - MAX_DISCM], sizeof(double));
        }
    }
    i_5 = numinf_1.nmol;
    FREADN(&interm_1.edpmlr[0], i_5 , sizeof(double));
    i_6 = numinf_1.nmol;
    FREADN(&interm_1.erpmlr[0], i_6 , sizeof(double));
    i_7 = numinf_1.nmol;
    FREADN(&interm_1.eljmlr[0], i_7 , sizeof(double));
    i_8 = numinf_1.nmol;
    FREADN(&interm_1.estmlr[0], i_8 , sizeof(double));

    FREAD(&mdpres_1.presur, sizeof(double));
    FREAD(&mdpres_1.presrm, sizeof(double));
    FREADN(&mdpres_1.presst[0], c__9 , sizeof(double));
    FREADN(&mdpres_1.prestm[0], c__9 , sizeof(double));
    FREADN(&mdpres_1.etkv[0]  , c__9 , sizeof(double));
    FREADN(&mdpres_1.etkm[0]  , c__9 , sizeof(double));
    FREADN(&mdpres_1.virial[0], c__9 , sizeof(double));
    FREADN(&mdpres_1.virilm[0], c__9 , sizeof(double));
    RECORD();

    if(alt) {
        RECORD();
        FREADN(&symdat_1.ucvec[0]  , c__6 , sizeof(double));
        FREADN(&symdat_1.ucvcar[0] , c__9 , sizeof(double));
        RECORD();
    }

#endif

    sumxyz   = gomp_GetTranslateArray();
         
    if(iappend) {   /* append atoms */
/*  get pointer to coordinate vectors */
        sprintf(OutText,"Discover frame (%d)",alt);
        Wstr = gomp_CreateMolecStruct(OutText , numat , APPEND);
        if ( Wstr < 0 )
            return(1);
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y   = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z  = gomp_GetModifiableAtomZCoordPointer(Wstr);

        i_1 = 0;

        RECORD();

        if( FixedAtoms ) {
            xtemp = gomp_AllocateFloatVector(movinf_1.natmov);
            ytemp = gomp_AllocateFloatVector(movinf_1.natmov);
            ztemp = gomp_AllocateFloatVector(movinf_1.natmov);

            for (iat = 0; iat < movinf_1.natmov; iat++) {
                FREAD(&xtemp[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &xtemp[iat] );
                } 
                FREAD(&ytemp[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &ytemp[iat] );
                } 
                FREAD(&ztemp[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &ztemp[iat] );
                } 
            }
            for(iat = 0 ; iat < numinf_1.natom ; iat++) {
                i_2 = iat + i_1;
                x[i_2] = x[iat] + sumxyz[0];
                y[i_2] = y[iat] + sumxyz[1];
                z[i_2] = z[iat] + sumxyz[2];
            }

            for(iat = 0 ; iat < movinf_1.natmov ; iat++) {
                i_2 = i_1 + movinf_1.lmove[iat] - 1;
                x[i_2] = xtemp[iat];
                y[i_2] = ytemp[iat];
                z[i_2] = ztemp[iat];
            }

            gomp_FreeVector(xtemp);
            gomp_FreeVector(ytemp);
            gomp_FreeVector(ztemp);
        }

        else {

            tempry_1.xcoor = &x[i_1]; 
            tempry_1.ycoor = &y[i_1]; 
            tempry_1.zcoor = &z[i_1];

            for (iat = 0; iat < numinf_1.natom; iat++) {
                FREAD(&tempry_1.xcoor[iat], sizeof(float)); 
                if ( swap_bytes ) {
                    gomp_Reverse_float( &tempry_1.xcoor[iat] );
                } 
                FREAD(&tempry_1.ycoor[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &tempry_1.ycoor[iat] );
                } 
                FREAD(&tempry_1.zcoor[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &tempry_1.zcoor[iat] );
                } 
            }
        }
    }

    else {  /* no append */

/*  get pointer to coordinate vectors */
        Wstr = 0;
        x    = gomp_GetModifiableAtomXCoordPointer(Wstr);
        y   = gomp_GetModifiableAtomYCoordPointer(Wstr);
        z  = gomp_GetModifiableAtomZCoordPointer(Wstr);
                     
        RECORD();

        if( FixedAtoms) {
            xtemp = gomp_AllocateFloatVector(movinf_1.natmov);
            ytemp = gomp_AllocateFloatVector(movinf_1.natmov);
            ztemp = gomp_AllocateFloatVector(movinf_1.natmov);

            for (iat = 0; iat < movinf_1.natmov; iat++) {
                FREAD(&xtemp[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &xtemp[iat] );
                } 
                FREAD(&ytemp[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &ytemp[iat] );
                } 
                FREAD(&ztemp[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &ztemp[iat] );
                } 
            }

            for(iat = 0 ; iat < numinf_1.natom ; iat++) {
                x[iat] -= sumxyz[0];
                y[iat] -= sumxyz[1];
                z[iat] -= sumxyz[2];
            }

            for(iat = 0 ; iat < movinf_1.natmov ; iat++) {
                i_2 = movinf_1.lmove[iat] - 1;
                x[i_2] = xtemp[iat]    - sumxyz[0];
                y[i_2] = ytemp[iat]   - sumxyz[1];
                z[i_2] = ztemp[iat]  - sumxyz[2];
            }

            gomp_FreeVector(xtemp);
            gomp_FreeVector(ytemp);
            gomp_FreeVector(ztemp);
        }

        else {

            tempry_1.xcoor = &x[0]; 
            tempry_1.ycoor = &y[0]; 
            tempry_1.zcoor = &z[0];

            for (iat = 0; iat < numinf_1.natom; iat++) {
                FREAD(&tempry_1.xcoor[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &tempry_1.xcoor[iat] );
                } 
                tempry_1.xcoor[iat] -=  sumxyz[0];
                FREAD(&tempry_1.ycoor[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &tempry_1.ycoor[iat] );
                } 
                tempry_1.ycoor[iat] -=  sumxyz[1];
                FREAD(&tempry_1.zcoor[iat], sizeof(float));
                if ( swap_bytes ) {
                    gomp_Reverse_float( &tempry_1.zcoor[iat] );
                } 
                tempry_1.zcoor[iat] -=  sumxyz[2];
            }
        }
    }

    RECORD();

#ifdef TRAJ_VELOCITY

/* get temp space for the velocities (they are not used now) */

    tempry_1.xvlcty = gomp_AllocateFloatVector(numinf_1.natom);
    tempry_1.yvlcty = gomp_AllocateFloatVector(numinf_1.natom);
    tempry_1.zvlcty = gomp_AllocateFloatVector(numinf_1.natom);

    RECORD();

    if(kntrl) 
        i_1 = movinf_1.natmov;
    else 
        i_1 = numinf_1.natom; 

    for (iat = 0; iat < i_1; ++iat) {
        FREAD(&tempry_1.xvlcty[iat], sizeof(float));
        FREAD(&tempry_1.yvlcty[iat], sizeof(float));
        FREAD(&tempry_1.zvlcty[iat], sizeof(float));
    }
    RECORD();

/* free the velocity array ...*/

    gomp_FreeVector(tempry_1.xvlcty);
    gomp_FreeVector(tempry_1.yvlcty);
    gomp_FreeVector(tempry_1.zvlcty);

#endif

    return(0);
#if 0
/*     PARTIAL RECORD HAS BEEN READ IN.  SET ERROR FLAG TO 1 AND RETURN */


/*
  L900:
*/

/* free gomp_ratch space ... */

    gomp_FreeVector(kindat_1.knd);
    gomp_FreeVector(kindat_1.jname);

    return (1);
#endif
} /* rdin22_ */

#undef enmols
#undef enrgys


