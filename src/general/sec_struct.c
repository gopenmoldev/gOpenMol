/*

Copyright (c) 1996 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include "gommath.h"
#include <tcl.h>

#include "gommonitor.h"
#include "measure.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "plumber.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

#define Rabs(a)    ( ( a ) > 0 ? (a) : -(a))

#define MAXppconf 10   /* maximum number of ppconf entries */

#define GAP_dist 25.          /*  Gap distance for different alpha-carbon */

struct secang {

    char phi[16];  /* space for 4 * 4 characters */
    char psi[16];
    char omega[16];
};

static struct secang ppept_ang  = {
    {"C   N   CA  C  "} ,
    /* phi:  C(i-1) N(i)  CA(i)  C(i)    */
    {"N   CA  C   N  "} ,
    /* psi:  N(i)   CA(i) C(i)   N(i+1)  */
    {"CA  C   N   CA "}
}; 
/* omega CA(i)  C(i)  N(i+1) CA(i+1) */


static int *phi_i;  /* index array to phi torsion angle array */
static int *psi_i;  /* index array to psi torsion angle array */
static int *omega_i;  /* index array to omega torsion angle array */
static int phi_e;     /* length of array */
static int psi_e;
static int omega_e;
static int *resi_num; /* residue numbers */
static float *phi_v;  /* the real values */
static float *psi_v;
static float *omega_v;

#if 0
struct pconf {

    char name[10];
    float phi;
    float psi;
    float omega;
    float rpt;
    float tpr;
    float delta;
    char color[BUFF_LEN];
}; 

static struct pconf ppept_conf[MAXppconf];
*/

/* structure to hold the trace of ramachandran       */

static struct {
    int trace_on;        /* switch to indicate that a trace is saved (=1)   */
    int trace_sets;      /* number of trace sets                            */
    int trace_step;      /* step length                                     */
    const int *trace_atoms;    /* number of traced atoms in each set              */
    const int *trace_list;     /* list of atoms to be traced                      */
    const float *phi;          /* array to contain the x coordinates of the trace */
    const float *psi;
} trace_info_rama;
#endif

static int   Num_KAA   = 20; /* number of known amino acids */
static const char *Name_KAA3 =
"ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL"; /* symbol as three letter names */
/*static const char *Name_KAA1 =  "ARNDCQEGHILKMFPSTWYV"; *//* symbol as one letter names */

static int rama_wind = 0; /* = 0 no display of Ramachandran plot , = 1 display */


static int IsItAminoAcid(const char *);
static int KillExtraSpaces(char *);

/***************************************************************************/
int gomp_SecondaryStructure()
/***************************************************************************/
{
    static int   i,j,chains;
    static int   first,first_aa;
    static int   num_ca;
    static char  chelp1[4];
    static char  chelp2[4];
    static float angle;
    static char  OutText[BUFF_LEN];
    static int   Wstr;
    static const float *x;
    static const float *y;
    static const float *z;

    if(!gomp_GetNumMolecStructs()) {
        rama_wind = 0;
        gomp_PrintMessage("?ERROR - no atoms for secondary structure analyzis");
        return(1);
    }

    Wstr   = 0;
    x       = gomp_GetAtomXCoordPointer(Wstr);
    y       = gomp_GetAtomYCoordPointer(Wstr);
    z       = gomp_GetAtomZCoordPointer(Wstr);

/* check that this is a protein */

    first_aa = -1;
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
        if((j = IsItAminoAcid(gomp_GetAtomResName(Wstr , i))) == 0) {
            first_aa = gomp_GetAtomResNum1(Wstr , i);
            break;
        }
    }

    if(j < 0) {
        rama_wind = 0;
        gomp_PrintMessage("?ERROR - this is most likely not a protein");
        return(1);
    }

/* check if there are several independent chains */
    chains = gomp_ProteinChainList(0, NULL, NULL, &num_ca);
    if(chains > 1) {
        gomp_PrintMessage("I'm sorry but this does not work for a protein with independent chains");
        sprintf(OutText,"There are now '%d' chains ",(chains+1));
        gomp_PrintMessage(OutText);
        rama_wind = 0;
        return(1);
    }
/*......................*/
    

/* check first if there is already memory reserved */
    gomp_FreeVector(phi_i);
    gomp_FreeVector(psi_i);
    gomp_FreeVector(omega_i);

    gomp_FreeVector(phi_v);
    gomp_FreeVector(psi_v);
    gomp_FreeVector(omega_v);

    phi_e = 0;
    psi_e = 0;
    omega_e = 0;

/* get space for the index arrays */
    phi_i     = gomp_AllocateIntVector(4 * (num_ca + 1));    
    psi_i     = gomp_AllocateIntVector(4 * (num_ca + 1));    
    omega_i   = gomp_AllocateIntVector(4 * (num_ca + 1));

    phi_v     = gomp_AllocateFloatVector(num_ca + 1);
    psi_v     = gomp_AllocateFloatVector(num_ca + 1);
    omega_v   = gomp_AllocateFloatVector(num_ca + 1);
    resi_num  = gomp_AllocateIntVector(num_ca + 1);

    phi_e     = 0;
    psi_e     = 0;
    omega_e   = 0;


    /* PHI ... */
    first = first_aa; /* first residue number */

/* find phi c */
    phi_e = 0;
    strncpy(chelp2,ppept_ang.phi, 4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                phi_i[4 * phi_e] = i;
                phi_e++;
                if(phi_e > num_ca) {
                    gomp_PrintERROR("reached limit for phi_e");
                    return(1);
                }
                first++;
            }
        }
    }         
/* find phi n */
    phi_e = 0;
    first = first_aa + 1; /* first residue number */
    strncpy(chelp2,ppept_ang.phi+ 4 , 4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1 , gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);
            if(gomp_StringMatch(chelp1,chelp2)) {
                phi_i[4 * phi_e + 1] = i;
                phi_e++;
                if(phi_e > num_ca) {
                    gomp_PrintERROR("reached limit for phi_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         
/* find phi ca*/
    phi_e = 0;
    first = first_aa + 1; /* first residue number */
    strncpy(chelp2,ppept_ang.phi+8, 4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                phi_i[4 * phi_e + 2] = i;
                phi_e++;
                if(phi_e > num_ca) {
                    gomp_PrintERROR("reached limit for phi_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         
/* find phi c */
    phi_e = 0;
    first = first_aa + 1; /* first residue number */
    strncpy(chelp2,ppept_ang.phi + 12, 4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                phi_i[4 * phi_e + 3] = i;
                phi_e++;
                if(phi_e > num_ca) {
                    gomp_PrintERROR("?ERROR - reached limit for phi_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         

    for(i = 0 ; i < phi_e ; i++)  {

        gomp_floDihedAngle(x[phi_i[4*i]]  ,y[phi_i[4*i]]  ,z[phi_i[4*i]],
                      x[phi_i[4*i+1]],y[phi_i[4*i+1]],z[phi_i[4*i+1]],
                      x[phi_i[4*i+2]],y[phi_i[4*i+2]],z[phi_i[4*i+2]],
                      x[phi_i[4*i+3]],y[phi_i[4*i+3]],z[phi_i[4*i+3]],
                      &angle);
#ifdef DEBUG
        print_names(phi_i[4*i]);
        print_names(phi_i[4*i+1]);
        print_names(phi_i[4*i+2]);
        print_names(phi_i[4*i+3]);
        phi_v[i] = 180.*angle/M_PI;
        printf(" = %f \n",180.*angle/M_PI);
#endif
        phi_v[i] = 180.*angle/M_PI;
    }

    /* PSI ... */

    first = gomp_GetAtomResNum1(Wstr , 0); /* first residue number */

/* find psi n */
    psi_e = 0;
    first = first_aa; /* first residue number */
    strncpy(chelp2,ppept_ang.psi, 4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                resi_num[psi_e]  = gomp_GetAtomResNum1(Wstr , i);
                psi_i[4 * psi_e] = i;
                psi_e++;
                if(psi_e > num_ca) {
                    gomp_PrintERROR("?ERROR - reached limit for psi_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         
/* find psi ca */
    psi_e = 0;
    first = first_aa; /* first residue number */
    strncpy(chelp2,ppept_ang.psi+4,4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                psi_i[4 * psi_e + 1] = i;
                psi_e++;
                if(psi_e > num_ca) {
                    gomp_PrintMessage("reached limit for psi_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         
/* find psi c */
    psi_e = 0;
    first = first_aa; /* first residue number */
    strncpy(chelp2,ppept_ang.psi+8, 4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                psi_i[4 * psi_e + 2] = i;
                psi_e++;
                if(psi_e > num_ca) {
                    gomp_PrintMessage("?ERROR - reached limit for psi_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         
/* find psi n */
    psi_e = 0;
    first = first_aa + 1; /* first residue number */
    strncpy(chelp2,ppept_ang.psi + 12, 4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                psi_i[4 * psi_e + 3] = i;
                psi_e++;
                if(psi_e > num_ca) {
                    gomp_PrintMessage("?ERROR - reached limit for psi_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         

    for(i = 0 ; i < psi_e ; i++)  {
        gomp_floDihedAngle(x[psi_i[4*i]]  ,y[psi_i[4*i]]  ,z[psi_i[4*i]],
                      x[psi_i[4*i+1]],y[psi_i[4*i+1]],z[psi_i[4*i+1]],
                      x[psi_i[4*i+2]],y[psi_i[4*i+2]],z[psi_i[4*i+2]],
                      x[psi_i[4*i+3]],y[psi_i[4*i+3]],z[psi_i[4*i+3]],
                      &angle);

#ifdef DEBUG
        print_names(psi_i[4*i]);
        print_names(psi_i[4*i+1]);
        print_names(psi_i[4*i+2]);
        print_names(psi_i[4*i+3]);
        psi_v[i] = 180.*angle/M_PI;
        printf(" = %f\n",180.*angle/M_PI);
#endif
        psi_v[i] = 180.*angle/M_PI;
    }

    /* OMEGA ... */

    first = first_aa; /* first residue number */

/* find omega ca */
    omega_e = 0;
    strncpy(chelp2,ppept_ang.omega,4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                omega_i[4 * omega_e] = i;
                omega_e++;
                if(omega_e > num_ca) {
                    gomp_PrintMessage("?ERROR - reached limit for omega_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         
/* find omega c */
    omega_e = 0;
    first = first_aa; /* first residue number */
    strncpy(chelp2,ppept_ang.omega+4,4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                omega_i[4 * omega_e + 1] = i;
                omega_e++;
                if(omega_e > num_ca) {
                    gomp_PrintMessage("?ERROR - reached limit for omega_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         
/* find omega n */
    omega_e = 0;
    first = first_aa + 1; /* first residue number */
    strncpy(chelp2,ppept_ang.omega+8,4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                omega_i[4 * omega_e + 2] = i;
                omega_e++;
                if(omega_e > num_ca) {
                    gomp_PrintMessage("?ERROR - reached limit for omega_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         
/* find omega ca */
    omega_e = 0;
    first = first_aa + 1; /* first residue number */
    strncpy(chelp2,ppept_ang.omega+12,4);
    (void)KillExtraSpaces(chelp2);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(gomp_GetAtomResNum1(Wstr , i) == first) {
            strncpy(chelp1,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);

            if(gomp_StringMatch(chelp1,chelp2)) {
                omega_i[4 * omega_e + 3] = i;
                omega_e++;
                if(omega_e > num_ca) {
                    gomp_PrintMessage("?ERROR - reached limit for omega_e");
                    rama_wind = 0;
                    return(1);
                }
                first++;
            }
        }
    }         

    for(i = 0 ; i < omega_e ; i++)  {
        gomp_floDihedAngle(x[omega_i[4*i]]  ,y[omega_i[4*i]]  ,z[omega_i[4*i]],
                      x[omega_i[4*i+1]],y[omega_i[4*i+1]],z[omega_i[4*i+1]],
                      x[omega_i[4*i+2]],y[omega_i[4*i+2]],z[omega_i[4*i+2]],
                      x[omega_i[4*i+3]],y[omega_i[4*i+3]],z[omega_i[4*i+3]],
                      &angle);

#ifdef DEBUG
        print_names(omega_i[4*i]);
        print_names(omega_i[4*i+1]);
        print_names(omega_i[4*i+2]);
        print_names(omega_i[4*i+3]);
        omega_v[i] = 180.*angle/M_PI;
        printf(" = %f\n",180.*angle/M_PI);
#endif
        omega_v[i] = 180.*angle/M_PI;
    }

    if(phi_e == 0 || psi_e == 0 || omega_e == 0) { /* probably not a protein*/
        gomp_PrintERROR(
            "can't find needed atoms (This is not a protein?)");
        rama_wind = 0;
        return(1);
    }

    return(0);
}

/***************************************************************************/
int IsItAminoAcid(const char *text)   
    /* checks to see if name in text is an amino acid */
    /* on return 0 = yes it is , < 0 it is not        */
/***************************************************************************/
{

    int i;
    char ctemp[MAX_RES_NAME_LEN];

    strncpy(ctemp,text,MAX_RES_NAME_LEN);

    for(i = 0 ; i < Num_KAA ; i++) {
        if(strncmp(ctemp,Name_KAA3+MAX_RES_NAME_LEN*i,3) == 0) return(0);
    }

    return(-1);
}
/***************************************************************************/
int gomp_ProteinChainList(int Wstr, int **chain_list, int **pAlpha_list, int *pAlpha_c)
/***************************************************************************/
{
    static int i,j,k;
    static int chains;
    static int alpha_c,*alpha_list;
    static float dist,dx,dy,dz;
    static char ctemp[MAX_ATM_NAME_LEN];
    static const float *x;
    static const float *y;
    static const float *z;

/* go and hunt for the alpha carbons */
    alpha_c = 0;
    x       = gomp_GetAtomXCoordPointer(Wstr);
    y       = gomp_GetAtomYCoordPointer(Wstr);
    z       = gomp_GetAtomZCoordPointer(Wstr);

    if( pAlpha_list && *pAlpha_list ) {
        alpha_list = *pAlpha_list;
        alpha_c    = *pAlpha_c;
    }
    else {
        alpha_list = gomp_AllocateIntVector(gomp_GetNumAtomsInMolecStruct(Wstr));
        
        for( i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
            strncpy(ctemp,gomp_GetAtomAtmName(Wstr , i),MAX_ATM_NAME_LEN);
            if(strncmp(ctemp,"CA",2) == 0) {
                alpha_list[alpha_c] = i;
                alpha_c++;
            }
        }

        if( pAlpha_list )
            *pAlpha_list = alpha_list;
        if( pAlpha_c)
            *pAlpha_c    = alpha_c;
    }

    if( alpha_c > 0 ) {
        chains = 1;
        if( chain_list ) {
            *chain_list    = gomp_AllocateIntVector(1);
            *chain_list[0] = 0;
        }
        
        for( i = 1 ; i < alpha_c ; i++) {
            j = alpha_list[i];
            k = alpha_list[i-1];
            dx = (x[j]-x[k]);
            dy = (y[j]-y[k]);
            dz = (z[j]-z[k]);
            dist=dx * dx + dy * dy + dz * dz;
            
            if(dist > GAP_dist ) {
                if( chain_list ) {
                    *chain_list           = gomp_ReallocateIntVector(*chain_list , chains + 1);
                    (*chain_list)[chains] = i;
                }
                chains++;
            }
        }
    }
    
    if( !pAlpha_list )
        gomp_FreeVector(alpha_list);
    return(chains);
}
/***************************************************************************/
int KillExtraSpaces(char *Text)
/***************************************************************************/
{
    int i;

    for(i = strlen(Text) - 1 ; i > 0 ; i--) {
        if(Text[i] == ' ') Text[i] = (char)NULL;
    }

    return(0);
}

/************************************************************************/
int  gomp_WriteBackboneDihedrals(const char *FileName)
/************************************************************************/
{
    FILE *File_p;
    int i,ii;
    char OutText[BUFF_LEN];
    float angle1,angle2,angle3;
    int   Wstr;
    const float *x;
    const float *y;
    const float *z;

    if(!psi_e || !phi_e || !omega_e) {
        gomp_PrintERROR("no backbone torsion angles available");
        return(1);
    }

    File_p = fopen(FileName,"w");
    if(File_p == NULL) {
        sprintf(OutText,"can't open file '%s' for writing",FileName);
        gomp_PrintERROR(OutText);
        return(1);
    }

    Wstr   = 0;
    x       = gomp_GetAtomXCoordPointer(Wstr);
    y       = gomp_GetAtomYCoordPointer(Wstr);
    z       = gomp_GetAtomZCoordPointer(Wstr);

    for(i = 0 ; i < phi_e ; i++) {

        if(!i) {
            gomp_floDihedAngle(x[psi_i[4*i]]  ,y[psi_i[4*i]]  ,z[psi_i[4*i]],
                          x[psi_i[4*i+1]],y[psi_i[4*i+1]],z[psi_i[4*i+1]],
                          x[psi_i[4*i+2]],y[psi_i[4*i+2]],z[psi_i[4*i+2]],
                          x[psi_i[4*i+3]],y[psi_i[4*i+3]],z[psi_i[4*i+3]],
                          &angle2);
            gomp_floDihedAngle(x[omega_i[4*i]]  ,y[omega_i[4*i]]  ,z[omega_i[4*i]],
                          x[omega_i[4*i+1]],y[omega_i[4*i+1]],z[omega_i[4*i+1]],
                          x[omega_i[4*i+2]],y[omega_i[4*i+2]],z[omega_i[4*i+2]],
                          x[omega_i[4*i+3]],y[omega_i[4*i+3]],z[omega_i[4*i+3]],
                          &angle3);
            fprintf(File_p,"%d               %f  %f\n",
                    (i+1),angle2 * 180./M_PI,angle3 * 180./M_PI);
        }
        else {
            ii = i - 1;
            gomp_floDihedAngle(x[phi_i[4*ii]]  ,y[phi_i[4*ii]]  ,z[phi_i[4*ii]],
                          x[phi_i[4*ii+1]],y[phi_i[4*ii+1]],z[phi_i[4*ii+1]],
                          x[phi_i[4*ii+2]],y[phi_i[4*ii+2]],z[phi_i[4*ii+2]],
                          x[phi_i[4*ii+3]],y[phi_i[4*ii+3]],z[phi_i[4*ii+3]],
                          &angle1);
            gomp_floDihedAngle(x[psi_i[4*i]]  ,y[psi_i[4*i]]  ,z[psi_i[4*i]],
                          x[psi_i[4*i+1]],y[psi_i[4*i+1]],z[psi_i[4*i+1]],
                          x[psi_i[4*i+2]],y[psi_i[4*i+2]],z[psi_i[4*i+2]],
                          x[psi_i[4*i+3]],y[psi_i[4*i+3]],z[psi_i[4*i+3]],
                          &angle2);
            gomp_floDihedAngle(x[omega_i[4*i]]  ,y[omega_i[4*i]]  ,z[omega_i[4*i]],
                          x[omega_i[4*i+1]],y[omega_i[4*i+1]],z[omega_i[4*i+1]],
                          x[omega_i[4*i+2]],y[omega_i[4*i+2]],z[omega_i[4*i+2]],
                          x[omega_i[4*i+3]],y[omega_i[4*i+3]],z[omega_i[4*i+3]],
                          &angle3);
            fprintf(File_p,"%d  %f  %f  %f\n",
                    (i+1),angle1 * 180./M_PI,
                    angle2 * 180./M_PI,
                    angle3 * 180./M_PI);
        }
    }
    i = phi_e - 1;
    gomp_floDihedAngle(x[phi_i[4*i]]  ,y[phi_i[4*i]]  ,z[phi_i[4*i]],
                  x[phi_i[4*i+1]],y[phi_i[4*i+1]],z[phi_i[4*i+1]],
                  x[phi_i[4*i+2]],y[phi_i[4*i+2]],z[phi_i[4*i+2]],
                  x[phi_i[4*i+3]],y[phi_i[4*i+3]],z[phi_i[4*i+3]],
                  &angle1);
    fprintf(File_p,"%d  %f\n",
            (i+1),angle1 * 180./M_PI);

    fclose(File_p);

    return(0);
}
