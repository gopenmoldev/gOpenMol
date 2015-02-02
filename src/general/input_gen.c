/*

Copyright (c) 1995 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "coord_file.h"
#include "gomstring.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"

#include "stdafx.h"

#define CHARGED_RESIDUES  16
#define NUM_ATOMS         19
#define CONV_ATOMS        57
#define DEFAULT_ICON8     "icon8.inp"
#define DEFAULT_MOPAC     "mopac6.inp"
#define DEFAULT_PROBESURF "probesurf.inp"
#define DEFAULT_VSS       "vss_scare.inp"
#define DEFAULT_DENSITY   "density.inp"

#define STRUCTURE         0

struct SURF_LIM {
    int set;             /* = 0 , not set , > 0 set */
    float Xmin;
    float Xmax;
    float Ymin;
    float Ymax;
    float Zmin;
    float Zmax;
    char  MeshFile[BUFF_LEN];
    char  WFFile[BUFF_LEN];
    float ProbeVal;
    int Orbital;
    int Dmethod;
    int Xpts;
    int Ypts;
    int Zpts;
} ;

void check_name_stack(int, const char *);
void build_name_stack(int);
static int which_element(char *);

/***********************************************************************/
/* generate icon8 input */
int gomp_ExternalInput4ICON8(int Which , const char *FileName)
/***********************************************************************/
{
   
    FILE *icon8_p;

    static const char *print = "FTTTTTTFFTTTTTTFTTTT";
    static const char *punch = "FFFFFFFFFTFFFFFFFFFF";
    static const char *def_title = "Default title for: ICON8";
    static char *element_vec;
    static char Icon8File[BUFF_LEN] = DEFAULT_ICON8;

    static int i,j,ret_val;
    static int nh; /* number of hydrogen atoms */
    static int na; /* number of heavy atoms    */
    static int charge=0; /* molecular charge   */
    static int meth=0;   /* calculational method desired */
    static int iprint=0;  
    static int ipunch=0;
    static char l1='F';
    static char l2='F';
    static char l3='F';
    static char l4='F';
    static char l5='F';
    static float con=0.0;  /* constant used in H(i,j) formula */
    static float peep=0.0; /* hydrogen orbital exponen (defaulr 1.30) */
    static float coulh=0.0;/* hydrogen H(i,i) (default -13.6) */

/*  static  int    ChargedResNum = CHARGED_RESIDUES;
    static  const char *ChargedRes =
        "ARGNARG ASP GLU HSC LYSNLYS GLYPCa  Mg  Fe  Zn  Li  Na  K   Cl  ";

    static float  ChargedResVal[CHARGED_RESIDUES] = {
        0.0f , 1.0f , -1.0f , -1.0f , 1.0f , 0.0f , 1.0f , 0.8f , 2.0f , 2.0f , 2.0f , 2.0f ,
        1.0f , 1.0f ,  1.0f , -1.0f};
*/
    static char   OutText[BUFF_LEN];
    static float  SCharge;
    static const char *DispListP;
/*  static int    Switch;*/
    static int    Wstr;
    static const float *x,*y,*z;
    static const float *sumxyz;

/* check first structures available */
    i = gomp_GetNumMolecStructs();
    if(!i) {
        gomp_PrintERROR("no structure available");
        return(1);
    }
    if(Which < 0 || Which >= i) {
        gomp_PrintERROR("structure index out of range");
        return(1);
    }

/* Check if an input name is supplied */
    if(*FileName != '\0') { /* ok there is a name coming ... */
        gomp_CopyString(Icon8File,FileName,BUFF_LEN);
    }
    else {
        gomp_CopyString(Icon8File,DEFAULT_ICON8,BUFF_LEN);
    }

    sprintf(OutText,"=> Writing ICON8 input to '%s'",Icon8File);
    gomp_PrintMessage(OutText);

    SCharge = 0.0;

/* open the output file and run ... */
    icon8_p = fopen(Icon8File,"w");
    if(icon8_p == NULL) {
        sprintf(OutText,"?ERROR - can't open output file '%s'",Icon8File);
        gomp_PrintMessage(OutText);
        return(1);
    }

/* calculate first number of hydrogens and heavy atoms */

    nh=0;
    na=0;

    Wstr = Which;
    DispListP = gomp_GetAtomDisplayStatePointer(Wstr);
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);
    sumxyz    = gomp_GetTranslateArray();

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(DispListP[i] == 0) continue;

        if(strncmp(gomp_GetAtomAtmName(Wstr , i),"H",1) == 0 ||
           strncmp(gomp_GetAtomAtmName(Wstr , i),"h",1) == 0) nh++;

        else na++;
    }


    gomp_PrintMessage("\n\n **** ICON8 input generator ");
    sprintf(OutText," Number of heavy atoms    : %d ",na);
    gomp_PrintMessage(OutText);
    sprintf(OutText," Number of hydrogen atoms : %d ",nh);
    gomp_PrintMessage(OutText);
    sprintf(OutText," System charge: %f",SCharge);
    gomp_PrintMessage(OutText);

    if((na+nh) != gomp_GetNumAtomsInMolecStruct(Wstr)) { 
        sprintf(OutText," Writing a subset of all atoms (na+nh) = %d",na+nh);
        gomp_PrintMessage(OutText);
    }

    element_vec = 
        gomp_AllocateCharVector(2 * gomp_GetNumAtomsInMolecStruct(Wstr));

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(STRUCTURE) ; i++) {
        if(DispListP[i] == 0) continue;
        strncpy(element_vec+2*i,gomp_GetAtomAtmName(Wstr , i),2);
    }

    ret_val = which_element(element_vec); 
    /* check the atoms if they are allowed */
    if(ret_val > 0) {
        free(element_vec);
        fclose(icon8_p);
        return(1);
    }

    fprintf(icon8_p,"%s\n",def_title); 

/* OBS! I'm killing nh now!  */

    na = nh + na;
    nh = 0;

    charge = SCharge > 0.0 ? (int) (SCharge + 0.5) : (int) (SCharge -0.5);

    fprintf(icon8_p,"%3d%3d%3d%3d%3d%3d%c%c%c%c%c%5.2f%6.3f%6.3f%s%s\n",
            nh,na,charge,meth,iprint,ipunch,l1,l2,l3,l4,l5,con,peep,coulh,
            print,punch);

/* write out the heavy atoms first      */
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(DispListP[i] == 0) continue;

        if(element_vec[2*i+1] == 'H' ||
           element_vec[2*i+1] == 'h') continue;

        fprintf(icon8_p,"%15.6f%15.6f%15.6f\n",(x[i]+sumxyz[0]),
                (y[i]+sumxyz[1]),
                (z[i]+sumxyz[2]));
    }

/* now the hydrogens                   */
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(DispListP[i] == 0) continue;

        if(element_vec[2*i+1] != 'H' &&
           element_vec[2*i+1] != 'h') continue;

        fprintf(icon8_p,"%15.6f%15.6f%15.6f\n",(x[i]+sumxyz[0]),
                (y[i]+sumxyz[1]),
                (z[i]+sumxyz[2]));
    }

/* now write the atom labels starting with the heavy atoms */
    j = 0;
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(DispListP[i] == 0) continue;
        if(element_vec[2*i+1] == 'H' ||
           element_vec[2*i+1] == 'h') continue;
        if(j == 40) {
            fprintf(icon8_p,"\n");
            j = 0;
        }

        fprintf(icon8_p,"%.2s",element_vec+2*i);
        j++;
    }
/* now the hydrogen atoms */

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        if(DispListP[i] == 0) continue;

        if(element_vec[2*i+1] != 'H' &&
           element_vec[2*i+1] != 'h') continue;

        if(j == 40) {
            fprintf(icon8_p,"\n");
            j = 0;
        }

        fprintf(icon8_p,"%.2s",element_vec+2*i);
        j++;
    }

    fprintf(icon8_p,"\n");
    fclose(icon8_p);
    free(element_vec);

    return(0);
}  


/**********************************************************************/
int which_element(char *element_vec)
/**********************************************************************/
{
/*  static int nsymbl= CONV_ATOMS;*/ /* atom symbol list */
/*
 *
 */
/*  static const char *pt =
        "?? HHeLiBe B C N O FNeNaMgAlSi P SClArCACBCGCDCECZNANBNGNDNENZOAOBOGODOEOZOHSASBSGSDHHHGHZNHCHOCHAHBHCHDHEHNHTHOCa";
    static const char *PT =
        " HHeLiBe B C N O FNeNaMgAlSi P SClArCa";

    static int velec[NUM_ATOMS]  = { 1 , 2 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 ,
                                     1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 2};
    static int vshell[NUM_ATOMS] = { 1 , 1 , 1 , 4 , 4 , 4 , 4 , 4 , 4 , 4 ,
                                     1 , 1 , 4 , 4 , 4 , 4 , 4 , 4 , 4};
                                       
    static int ihelpv[CONV_ATOMS] = {0,
                                     1,2,
                                     3,4,
                                     5,6,7,8,9,10,
                                     11,12,
                                     13,14,15,16,17,18,
                                     6,6,6,6,6,6,
                                     7,7,7,7,7,7,
                                     8,8,8,8,8,8,8,
                                     16,16,16,16,
                                     1,1,1,7,6,8,1,1,1,1,1,1,1,1,19};

    static int    morbit[3] = { 1 , 5 , 14 };
    static int    me[3]     = { 2 , 10 , 18 };
    static int    mshell    = 3;
*/
    static int    natom,norbit,ne;
    static int    i/*,j,swtch,nerror*/;
    static const char *disp_list;
    static int    Wstr;


/*  static char OutText[BUFF_LEN];*/

    natom  = 0;
    norbit = 0;
    ne     = 0;

    Wstr = 0;
    disp_list = gomp_GetAtomDisplayStatePointer(Wstr);

/* clear all numbers in the element vector */
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
        if(disp_list[i] == 0) continue;

        if(element_vec[2*i+1] == ' ' || 
           element_vec[2*i+1] == '\0') {
            element_vec[2*i+1] = element_vec[2*i];
            element_vec[2*i] = ' ';
           }

        if(isdigit(element_vec[2*i+1])) {
            element_vec[2*i+1] = element_vec[2*i];
            element_vec[2*i] = ' ';
        }
    }

/*
  for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(STRUCTURE) ; i++) {

  if(disp_list[i] == 0) continue;
  swtch = 0;
  for(j = 0 ; j < nsymbl ; j++) { 
  if(strncmp((element_vec+2*i),(pt+2*j),2) == 0)  {

  strncpy(element_vec+2*i,PT+2*(ihelpv[j]-1),2);
  ne     += velec[ihelpv[j] - 1];
  natom  += 1;
  norbit += vshell[ihelpv[j] - 1];
  swtch   = 1;
  break;

  break;
  }
  }
  if(swtch != 0) continue;

  sprintf(OutText," ?WARNING - undefined element  >%.2s<",
  element_vec+2*i);
  gomp_PrintMessage(OutText);
  gomp_PrintMessage("Allowed atoms are:\n");
  for(j = 0 ; j < strlen(PT)/ 30 ; j++) {
  sprintf(OutText,"%s",PT+j*30);
  gomp_PrintMessage(OutText);
  }
  if(strlen(PT) - 30 * ( strlen(PT) / 30)) {
  sprintf(OutText,"%s",PT+(strlen(PT) / 30)*30);
  gomp_PrintMessage(OutText);
  }
  nerror++;
  }
  if(nerror > 0) return(1);

  sprintf(OutText," Number of atoms : %d ",natom);
  gomp_PrintMessage(OutText);
  sprintf(OutText," Number of orbitals: %d ",norbit);
  gomp_PrintMessage(OutText);
  sprintf(OutText," Number of electrons: %d ",ne);
  gomp_PrintMessage(OutText);
*/

    return(0);
}
/***********************************************************************/
int gomp_ExternalInput4PROBESURF(int Which ,
                               const char *FileName , float MinX , float MaxX ,
                               float MinY , float MaxY ,
                               float MinZ , float MaxZ ,
                               int Xpts, int Ypts, int Zpts ,
                               float ProbeRad)
    /* generate probsurf input */
/***********************************************************************/
{
   
    FILE *probsurf_p;

    static const char *def_title = "Default title for: PROBESURF";
    static char probsurf_input[BUFF_LEN] = DEFAULT_PROBESURF;

    static int i;
    static int na; /* number of heavy atoms    */

    char OutText[BUFF_LEN];
    int FileNameDef = 0;  /* 0 not defined, 1 defined */
    const char *DispListP;
    const float *Move;

/* check first structures available */
    i = gomp_GetNumMolecStructs();
    if(!i) {
        gomp_PrintERROR("no structure available");
        return(1);
    }
    if(Which < 0 || Which >= i) {
        gomp_PrintERROR("structure index out of range");
        return(1);
    }

    FileNameDef = 0;

/* Check if an input name is supplied */
    if(*FileName != '\0' && FileName[strlen(FileName) - 1]  != '/') { 
        /* ok there is a name coming ... */
        gomp_CopyString(probsurf_input,FileName,BUFF_LEN);
        FileNameDef = 1;
    }
    else {
        gomp_CopyString(probsurf_input,DEFAULT_PROBESURF,BUFF_LEN);
    }

    sprintf(OutText,"=> Writing PROBESURF input to '%s'",probsurf_input);
    gomp_PrintMessage(OutText);

    probsurf_p = fopen(probsurf_input,"w");
    if(probsurf_p == NULL) {
        sprintf(OutText,"?ERROR - can't open PROBESURF input file '%s' \n",probsurf_input);
        gomp_PrintMessage(OutText);
        return(1);
    }

    DispListP = gomp_GetAtomDisplayStatePointer(Which);

    na = 0;
    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {
        if(DispListP[i]) na++;
    }
 
    gomp_PrintMessage("\n\n **** PROBESURF input generator \n");
    sprintf(OutText," Number of atoms    : %d ",na);
    gomp_PrintMessage(OutText);

/* start writing input file */
    fprintf(probsurf_p,"%s\n",def_title);
    fprintf(probsurf_p,"%s\n",def_title);
    fprintf(probsurf_p,"%d \n",na);

    Move  = gomp_GetTranslateArray();

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Which) ; i++) {

        if(!DispListP[i]) continue; /* look into display list */

        fprintf(probsurf_p,"%d %s %s %s %f %f %f %f \n",
                (i+1),
                gomp_GetAtomSegName(Which , i),
                gomp_GetAtomResName(Which , i),
                gomp_GetAtomAtmName(Which , i),
                (gomp_GetAtomXCoord(Which , i) + Move[0]),
                (gomp_GetAtomYCoord(Which , i) + Move[1]),
                (gomp_GetAtomZCoord(Which , i) + Move[2]),
                gomp_GetAtomVdwRad(Which , i));
    }

    fprintf(probsurf_p,"%f %f %f %f %f %f \n",MinX ,
            MaxX,
            MinY,
            MaxY,
            MinZ,
            MaxZ);
    fprintf(probsurf_p,"%d %d %d\n",Xpts,
            Ypts,
            Zpts);
    fprintf(probsurf_p,"%f \n",ProbeRad);

    fclose(probsurf_p);

    return(0);
}
