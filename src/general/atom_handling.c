/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <tcl.h>

#include "atom_param.h"
#include "bond.h"
#include "colouring.h"
#include "coord_file.h"
#include "g_status.h"
#include "gomenv.h"
#include "gomfile.h"
#include "gommain.h"
#include "gomstring.h"
#include "model_file.h"
#include "molecule.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "memalloc.h"
#include "plot_molec.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#include "stdafx.h"

#define DEFAULT_ATOM_TYPE    499
#define DEFAULT_ATOM_COLOUR "cyan"

#define MAX_FLOAT  1.e+35    /* This max float is just a guess */
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define Rabs(a)    ( ( a ) > 0   ? (a) : -(a))
#define Fabs(a)    ( ( a ) > 0.0 ? (a) : -(a))

/* 
   list of standard covalent radii as they are used in CSD 3.1 
   Cambridge Structural Data Base System. User's Manual
   Part II, page 84.
*/ 


/*
  const char *AtomSymbols={"\
  Ac  Ag  Al  Am  As  Au  B   Ba  Be  Bi  Br  C   Ca  Cd  \
  Ce  Cl  Co  Cr  Cs  Cu  D   Dy  Er  Eu  F   Fe  Ga  Gd  \
  Ge  H   Hf  Hg  Ho  I   In  Ir  K   La  Li  Lu  Mg  Mn  \
  Mo  N   Na  Nb  Nd  Ni  Np  O   Os  P   Pa  Pb  Pd  Pm  \
  Po  Pr  Pt  Pu  Ra  Rb  Re  Rh  Ru  S   Sb  Sc  Se  Si  \
  Sm  Sn  Sr  Ta  Tb  Tc  Te  Th  Ti  Tl  Tm  U   V   W   \
  Y   Yb  Zn  Zr  He  Ne  Ar  Kr  Xe  Rn  Fr  Cm  Bk  Cf  \
  Es  Fm  XX  YY  ZZ  "};
*/
/* ..................................................................*/
/*#define ATOM_NAMES_TOTAL 104
  const char *AtomNameTable[ATOM_NAMES_TOTAL] = {
  "*Unknown*",
  "Hydrogen","Helium",
  "Lithium","Beryllium","Boron","Carbon","Nitrogen","Oxygen","Fluorine","Neon",
  "Sodium","Magnesium","Aluminium","Silicon","Phosphorus","Sulfur","Chlorine",
  "Argon",
  "Potassium","Calcium","Scandium","Titanum","Vanadium","Chromium","Manganese",
  "Iron","Cobalt","Nickel","Copper","Zinc","Gallium","Germanium","Arsenic",
  "Selenium","Bromine","Krypton",
  "Rubidium","Strontium","Yttrium","Zirconium","Niobium","Molybdenum",
  "Technetium","Ruthenium","Rhodium","Palladium","Silver","Cadmium","Indium",
  "Tin","Antimony","Tellurium","Iodine","Xenon",
  "Caesium","Barium",
  "Lanthanum","Cerium","Praseodymium","Neodymium","Promethium","Samarium",
  "Europium","Gadolinium","Terbium","Dysprosium","Holmium","Erbium",
  "Thulium","Ytterbium","Lutetium",
  "Hafnium","Tantalum","Tungsten","Rhenium","Osmium","Iridium","Platinum",
  "Gold","Mercury","Thallium","Lead","Bismuth","Polonium","Astatine",
  "Radon",
  "Francium","Radium",
  "Actinium","Thorium","Protactinium","Uranium","Neptunium","Plutonium",
  "Americium","Curium","Bercelium","Californium","Einsteinum","Fermium" ,
  "PseudoAtomX" , "PseudoAtomY" , "PseudoAtomZ"};
*/
/* ..................................................................*/

/* this is a sort of temporary fix ... */
/*   int NumAtomSymbols = 104;

int AtomSymbol_p[] =
{ 89 , 47 , 13 , 95 , 33 , 79 , 5  , 56 , 4  , 83 , 35 , 6  , 20 , 48 , 
58 , 17 , 27 , 24 , 55 , 29 , 0  , 66 , 68 , 63 , 9  , 26 , 31 , 64 ,
32 , 1  , 72 , 80 , 67 , 53 , 49 , 77 , 19 , 57 , 3  , 71 , 12 , 25 ,
42 , 7  , 11 , 41 , 60 , 28 , 93 , 8  , 76 , 15 , 91 , 82 , 46 , 61 ,
84 , 59 , 78 , 94 , 88 , 37 , 75 , 45 , 44 , 16 , 51 , 21 , 34 , 14 , 
62 , 50 , 38 , 73 , 65 , 43 , 52 , 90 , 22 , 81 , 69 , 92 , 23 , 74 , 
39 , 70 , 30 , 40 , 2  , 10 , 18 , 36 , 54 , 86 , 87 , 96 , 97 , 98 , 
99 ,100 ,101 ,102 ,103};
 
float AtomCovar[] =
{1.88,1.59,1.35,1.51,1.21,1.50,0.83,1.34,0.35,1.54,1.21,0.68,0.99,1.69,
 1.83,0.99,1.33,1.35,1.67,1.52,0.23,1.75,1.73,1.99,0.64,1.34,1.22,1.79,
 1.17,0.23,1.57,1.70,1.74,1.40,1.63,1.32,1.33,1.87,0.68,1.72,1.10,1.35,
 1.47,0.68,0.97,1.48,1.81,1.50,1.55,0.68,1.37,1.05,1.61,1.54,1.50,1.80,
 1.68,1.82,1.50,1.53,1.90,1.47,1.35,1.45,1.40,1.02,1.46,1.44,1.22,1.20,
 1.80,1.46,1.12,1.43,1.76,1.35,1.47,1.79,1.47,1.55,1.72,1.58,1.33,1.37,
 1.78,1.94,1.45,1.56,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
 0.00,0.00,5.00,10.00,20.00};

float Atom_max_corad2 = 16.0 ; *//* trigger value in connection search */

static struct {
    int ElementCount;
    struct {
        char Symbol[MAX_ATM_NAME_LEN];
        float Covar;
    } *Elements;
    float AtomMaxCorad2;
} ElementData = { 0 , NULL , 16.0 };

static int ConnSearch_window = 500;
/* search window for the atom connection analyzis */
/* the bigger this value is the longer it takes */


static float MinimumSizeOfSystem  = 10.0;
static float GetMinimumSizeOfSystem(void);

#if 0
static int  ShowAllowedAtomSymbols();
#endif
static int  IdentifyAtomColours(int);
static int  AtomCompare(const char *, const char *);
static int  AtomCompareWithNumbers(const char *, const char *);
static int  GetAtomIndexInTable(int);

/*  */
static int CalculateMoleculeDimensions(void);
static struct {
    float TransDamping;
    float SizeOfSystem;
    float MinX,MaxX;
    float MinY,MaxY;
    float MinZ,MaxZ;
    float Min,Max;
    int   TranslateState;
} MoleculeDimensions;

/*  */

typedef struct {
    char  AtomName[BUFF_LEN];
    int   AtomType;
    int   AtomTypeIndex;
    int   AtomIndex;
    char  AtomColour[BUFF_LEN];
    char  AtomBasisSet[BUFF_LEN];
} AtomHit;

static const AtomHit *LookIntoAtomConversionTable(const char *);
static const AtomHit *LookIntoAtomConversionTableDirectHit(const char *);
static int     AssignAtomInfo(int , int , const AtomHit *);
static int     AssignColour2Atoms(int , int , const AtomHit *);
static int     ExpandAtomConversionTable(void);
static int     PutIntoAtomConversionTable(int , const char * , const char * , int , const char * ,const char *);
static int     AtomsInAtomConversionTable(void);

static struct {
    int Atoms;                   /* atoms in conversion table */
    int Reserved;                /* space reserved for number of atoms */
    int Step;                    /* increment in which space is reserved */
    char *Atom;                  /* real atom name */
    char *Equiv;                 /* equivalence name */
    int *ContInt;               /* the equivalence name contains an integer/integers */
    int *AtomType;              /* pointer to any atom type table */
    int *AtomTypeIndex;         /* pointer to any atom type table index */
    int *LabelIndex;            /* pointer to label index table   */
    char *Colour;                /* colour to be used with this atom */
    char *BasisSet;              /* atom basis set tag               */
} AtomConversionTable = {0 , 0 , 50 , 
                         NULL , NULL , NULL ,
                         NULL , NULL , NULL , NULL ,NULL };

static struct {
    int BondState;
    int HBondState;
} AtomReconnectivityState;

static int BondDisplayStyle  = 1; /* default "half bond display" */

#if 0
/************************************************************************/
int ShowAllowedAtomSymbols()
/************************************************************************/
{
  char OutText[BUFF_LEN];
  const char *p;

  (void)gomp_PrintMessage("Allowed atoms are:");

  p = AtomSymbols;

  while(*p != '\0') {
      sprintf(OutText,"%.40s",p);
      gomp_PrintMessage(OutText);
      p += 40;
  }

  return(0);
}
#endif
/***********************************************************************/
int gomp_ParseGetAtomInfo(const char *AtonName)
/***********************************************************************/
{
    int i,j,Hit;
    static char    OutText[5*BUFF_LEN];
    static char    TempAtom[BUFF_LEN],atype[MAX_ATM_NAME_LEN+1];
    float          Red,Green,Blue,Covar;
    const AtomHit *Atom_t;

    Atom_t = NULL;
    Hit    = 0;

#ifndef ELEM_PARAM
    NumAtomSymbols = strlen(AtomSymbols)/4;
#endif

    if(isdigit(AtonName[0])) {
        i = GetAtomIndexInTable(atoi(AtonName));

        if( i < 0  ) {
            gomp_PrintERROR("Atom index out of allowed range");
            return(1);
        }

        strncpy(TempAtom, gomp_GetAtomSymbol(i), MAX_ATM_NAME_LEN);
    }
    else
        strncpy(TempAtom, AtonName, MAX_ATM_NAME_LEN);

    if(*TempAtom == '\0')
        return(1);

/*   first level hit by looking into the atom conversion table        */    
    Atom_t = LookIntoAtomConversionTable(TempAtom);
    if(Atom_t->AtomName[0] != '\0')
        Hit = Atom_t->AtomIndex;
    else {
        Hit = gomp_MatchAtom(TempAtom);
        if(Hit) {
            Atom_t = LookIntoAtomConversionTableDirectHit(TempAtom);
            if(Atom_t->AtomName[0] == '\0')
                Atom_t = NULL;
        }
    }

    if(!Hit) {
/*   direct hit by decreasing the number of characters */
        i = strlen(TempAtom);
    
        for(j = 1 ; j < i ; j++) {
            if(isdigit(TempAtom[j])) {
                TempAtom[j] = '\0';
                break;
            }
        }

        for(j = strlen(TempAtom); j > 0 ; j--) {
            TempAtom[j] = '\0';

/* new try ... */
            Atom_t = LookIntoAtomConversionTable(TempAtom);
            if(Atom_t->AtomName[0] != '\0')
                Hit = gomp_MatchAtom(Atom_t->AtomName);
            else {
                Hit = gomp_MatchAtom(TempAtom);
                if(Hit) {
                    Atom_t = LookIntoAtomConversionTableDirectHit(TempAtom);
                    if(Atom_t->AtomName[0] == '\0')
                        Atom_t = NULL;
                }
            }

            if(Hit)
                break;

            Atom_t = NULL;
        }
    }

    if(Hit && Hit<=ElementData.ElementCount)
        Covar = ElementData.Elements[Hit-1].Covar;
    else
        Covar = 0.0f;

    if(Atom_t) {
        strncpy(TempAtom,Atom_t->AtomName,MAX_ATM_NAME_LEN);
        gomp_ColourName2RGB(Atom_t->AtomColour,&Red,&Green,&Blue);

        if(Hit) {
            strncpy(atype,gomp_GetAtom_atype(Atom_t->AtomTypeIndex),
                    MAX_ATM_NAME_LEN);

            sprintf(OutText,"{%s} %d {%f %f %f} {%s} %f {"
                    "type %d bndrad %f vdwrad %f plurad %f global %c "
                    "emin %f rmin %f patom %f mass %f cnct %d hbond %c atype %s"
                    "}",
                    TempAtom,Hit-1,
                    Red,Green,Blue,
                    Atom_t->AtomBasisSet,Covar,
                    gomp_GetAtom_type(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_bndrad(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_vdwrad(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_plurad(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_global(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_emin(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_rmin(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_patom(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_mass(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_cnct(Atom_t->AtomTypeIndex),
                    gomp_GetAtom_hbond(Atom_t->AtomTypeIndex),
                    atype);
        }
        else
            sprintf(OutText,"{%s} %d {%f %f %f} {%s} %f {}",
                    TempAtom,Hit-1,
                    Red,Green,Blue,
                    Atom_t->AtomBasisSet,Covar);
    }
    else
        sprintf(OutText,"{} -1 {1.0 1.0 1.0} {} %f {}",Covar);

    gomp_SendTclReturn(OutText);

    return(0);
}

/***********************************************************************/
int gomp_IdentifyAtom(int AtomStructure, int AtomIndex, int *pHits)
/***********************************************************************/
{
    int i,j,k;
    int Hit;
    const atom_name_t *Atom_p;
    const AtomHit *Atom_t;
    const AtomHit *Atom_tt;
    int      NumAtoms;
    int      Hits[10];
    static char  TempAtom[BUFF_LEN];
    static char  TempAtom1[BUFF_LEN];
    static char  OutText[BUFF_LEN];

#ifndef ELEM_PARAM
    NumAtomSymbols = strlen(AtomSymbols)/4;
#endif

    if(AtomStructure < 0 ||
       AtomStructure > gomp_GetNumMolecStructs()) {
        gomp_PrintMessage("?gOpenMol - ERROR");
        gomp_PrintMessage("Error in atom structure number");
        return(1);
    }

    NumAtoms = gomp_GetNumAtomsInMolecStruct(AtomStructure);

    if(AtomIndex < 0 ||
       AtomIndex >= NumAtoms ) {
        gomp_PrintMessage("?gOpenMol - ERROR");
        gomp_PrintMessage("Error in atom index");
        return(1);
    }

    Atom_p = gomp_GetAtomAtmNamePointer(AtomStructure);

    for(i = 0 ; i < 10 ; i++) Hits[i] = 0;

/* Put atom, residue and segment names into stacks */
    (void)gomp_UpdateNameStack(1, gomp_GetAtomAtmName(AtomStructure , AtomIndex));
    (void)gomp_UpdateNameStack(2, gomp_GetAtomResName(AtomStructure , AtomIndex));
    (void)gomp_UpdateNameStack(3, gomp_GetAtomSegName(AtomStructure , AtomIndex));

    Hit = 0;
    
    strncpy(TempAtom  , Atom_p[AtomIndex] , MAX_ATM_NAME_LEN);
    strncpy(TempAtom1 , Atom_p[AtomIndex] , MAX_ATM_NAME_LEN);

    if(*TempAtom == '\0') {
        sprintf(OutText,"atom name can not be empty string (NULL)\nUnknown atom #%d",AtomIndex+1);
        gomp_PrintEXIT(OutText);
    }

/*   first level hit by looking into the atom conversion table        */
    
    Atom_t = LookIntoAtomConversionTable(TempAtom);
    if(Atom_t->AtomName[0] != '\0') {
        Hit  = Atom_t->AtomIndex;
        if(Hit) 
            (void)AssignAtomInfo( AtomStructure , AtomIndex , Atom_t); 
        Hits[0] += 1;
        Hits[6] += 1;
    }
    else {
        Hit = gomp_MatchAtom(TempAtom);
        if(Hit) {
            Atom_tt = LookIntoAtomConversionTableDirectHit(TempAtom);
            if(Atom_tt->AtomName[0] != '\0') {
                (void)AssignAtomInfo( AtomStructure , AtomIndex , Atom_tt);
                Hits[1] += 1;
                Hits[6] += 1;
            } else {
                sprintf(OutText,"Can't assign type to atom '%s'",TempAtom);
                gomp_PrintMessage(OutText);
                Hits[2] += 1;
            }
        }
    }

/*   direct hit by decreasing the number of characters */
    if(!Hit) {
        j = strlen(TempAtom);
        
        for(k = 0 ; k < j ; k++) {
            if(isdigit(TempAtom[k])) {
                if(k == 0) {
                    /* PDB files contains hydrogens starting with a digit */
                    if(TempAtom[1]!='H')
                        gomp_PrintMessage(
                            "?WARNING - first character in atom name is a digit");
                    continue;
                }
                else {
                    TempAtom[k] = '\0';
                    break;
                }
            }
        }

        for(k = strlen(TempAtom); k > 0 ; k--) {
            TempAtom[k] = '\0';

/* new try ... */
            Atom_t = LookIntoAtomConversionTable(TempAtom);
            if(Atom_t->AtomName[0] != '\0') {
                Hit  = gomp_MatchAtom(Atom_t->AtomName);
                if(Hit) 
                    (void)AssignAtomInfo( AtomStructure , AtomIndex , Atom_t); 
                Hits[3] += 1;
                Hits[6] += 1;
            }
            else {
                Hit = gomp_MatchAtom(TempAtom);
                if(Hit) {
                    Atom_tt = LookIntoAtomConversionTableDirectHit(TempAtom);
                    if(Atom_tt->AtomName[0] != '\0') {
                        (void)AssignAtomInfo( AtomStructure , AtomIndex , Atom_tt);
                        Hits[4] += 1;
                        Hits[6] += 1;
                    } else {
                        sprintf(OutText,"Can't assign type to atom '%s'",TempAtom);
                        gomp_PrintMessage(OutText);
                        Hits[5] += 1;
                    }
                }
            }

            if(Hit)
                break;
        }   
    }

#ifndef ELEM_PARAM
    if(Hit) {
        (void)gomp_PutAtomCovar(AtomStructure , AtomCovar[Hit-1] , AtomIndex);
    }
#else
    if(Hit && Hit<=ElementData.ElementCount) {
        (void)gomp_PutAtomCovar(AtomStructure,
                            ElementData.Elements[Hit-1].Covar, AtomIndex);
    }
#endif
    else {
        gomp_PutAtomColour(AtomStructure, 1.0, 1.0, 1.0, AtomIndex);
        
        gomp_PutAtomAtype(AtomStructure , "???" , AtomIndex);
        sprintf(OutText,"Unknown atom #%d:'%.4s' can't calculate connectivity for this atom!",AtomIndex+1,TempAtom1);
        gomp_PrintWARNING(OutText);
    }

    if(pHits) {
        for( i=0; i < 10; i++ ) pHits[i] += Hits[i];
    }
    return(0);
}

/***********************************************************************/
int gomp_IdentifyAtoms(int AtomStructure)
/***********************************************************************/
{
    int i;
    int      NumAtoms;
    int      Hits[10];
    char  OutText[BUFF_LEN];

#ifndef ELEM_PARAM
    NumAtomSymbols = strlen(AtomSymbols)/4;
#endif

    if(AtomStructure < 0 ||
       AtomStructure > gomp_GetNumMolecStructs()) {
        gomp_PrintMessage("?gOpenMol - ERROR");
        gomp_PrintMessage("Error in atom structure number");
        return(1);
    }

    NumAtoms = gomp_GetNumAtomsInMolecStruct(AtomStructure);

    for(i = 0 ; i < 10 ; i++) Hits[i] = 0;

    for(i = 0 ; i < NumAtoms ; i++) {
/* link to the display system !!!!!!!!!!!! */
        (void)gomp_PutText2StatusLine2((10 * (i + 1))/ NumAtoms);
    
        gomp_IdentifyAtom(AtomStructure,i,Hits);
    }


/* link to the display system !!!!!!!!!!!! */
    (void)gomp_PutText2StatusLine2(0);

    gomp_PrintMessage("Hits in bins:");
    sprintf(OutText,"Fast: #1(direct):%d #2(indirect):%d #3(lost):%d",Hits[0],Hits[1],Hits[2]);
    gomp_PrintMessage(OutText);
    sprintf(OutText,"Slow: #4(direct):%d #5(indirect):%d #6(lost):%d",Hits[3],Hits[4],Hits[5]);
    gomp_PrintMessage(OutText);
    sprintf(OutText,"#1+#2+#4+#5 = %d (Atoms:%d, lost:%d)",Hits[6],NumAtoms,(NumAtoms - Hits[6]));
    gomp_PrintMessage(OutText);
    if(Hits[3]) {
        if(Hits[0]/Hits[3] < 1) {
            gomp_PrintMessage("please add more of your atoms to the 'atom_conversion.data' file");
        }
    }
    return(0);
}


/***********************************************************************/
int IdentifyAtomColours(int AtomStructure)
/***********************************************************************/
{
    int i,j,k;
    int Hit;
    const atom_name_t *Atom_p;
    const AtomHit *Atom_t;
    const AtomHit *Atom_tt;
    int      NumAtoms;
    char  TempAtom[BUFF_LEN];
    char  TempAtom1[BUFF_LEN];
    char  OutText[BUFF_LEN];

#ifndef ELEM_PARAM
    NumAtomSymbols = strlen(AtomSymbols)/4;
#endif

    if(AtomStructure < 0 ||
       AtomStructure > gomp_GetNumMolecStructs()) {
        gomp_PrintMessage("?gOpenMol - ERROR");
        gomp_PrintMessage("Error in atom structure number");
        return(1);
    }

    Atom_p = gomp_GetAtomAtmNamePointer(AtomStructure);

    NumAtoms = gomp_GetNumAtomsInMolecStruct(AtomStructure);

    for(i = 0 ; i < NumAtoms ; i++) {

/* link to the display system !!!!!!!!!!!! */
        (void)gomp_PutText2StatusLine2((10 * (i + 1))/ NumAtoms);

        Hit = 0;

        strncpy(TempAtom  , Atom_p[i] , MAX_ATM_NAME_LEN);
        strncpy(TempAtom1 , Atom_p[i] , MAX_ATM_NAME_LEN);

        if(*TempAtom == '\0') {
            sprintf(OutText,"atom name can not be empty string (NULL)\nUnknown atom #%d",i+1);
            gomp_PrintEXIT(OutText);
        }

/*   first level hit by looking into the atom conversion table        */

        Atom_t = LookIntoAtomConversionTable(TempAtom);
        if(Atom_t->AtomName[0] != '\0') {
            Hit  = Atom_t->AtomIndex;
            if(Hit) 
                (void)AssignColour2Atoms( AtomStructure , i , Atom_t); 
        }
        else {
            Hit = gomp_MatchAtom(TempAtom);
            if(Hit) {
                Atom_tt = LookIntoAtomConversionTableDirectHit(TempAtom);
                if(Atom_tt->AtomName[0] != '\0') {
                    (void)AssignColour2Atoms( AtomStructure , i , Atom_tt);
                } else {
                    sprintf(OutText,"Can't assign type to atom '%s'",TempAtom);
                    gomp_PrintMessage(OutText);
                }
            }
        }

/*   direct hit by decreasing the number of characters */
        if(!Hit) {
            j = strlen(TempAtom);

            for(k = 0 ; k < j ; k++) {
                if(isdigit(TempAtom[k])) {
                    if(k == 0) {
                        gomp_PrintMessage(
                            "?WARNING - first character in atom name is a digit");
                        continue;
                    }
                    else {
                        TempAtom[k] = '\0';
                        break;
                    }
                }
            }

            for(k = strlen(TempAtom) ; k > 0 ; k--) {
                TempAtom[k] = '\0';

/* new try ... */
                Atom_t = LookIntoAtomConversionTable(TempAtom);
                if(Atom_t->AtomName[0] != '\0') {
                    Hit  = gomp_MatchAtom(Atom_t->AtomName);
                    if(Hit) 
                        (void)AssignColour2Atoms( AtomStructure , i , Atom_t); 
                }
                else {
                    Hit = gomp_MatchAtom(TempAtom);
                    if(Hit) {
                        Atom_tt = LookIntoAtomConversionTableDirectHit(TempAtom);
                        if(Atom_tt->AtomName[0] != '\0') {
                            (void)AssignColour2Atoms( AtomStructure , i , Atom_tt);
                        } else {
                            sprintf(OutText,"Can't assign type to atom '%s'",TempAtom);
                            gomp_PrintMessage(OutText);
                        }
                    }
                }

                if(Hit)
                    break;
            }
        }

        if(!Hit) {
            gomp_PutAtomColour(AtomStructure, 1.0, 1.0, 1.0, i);
            
            sprintf(OutText,"Unknown atom #%d:'%.4s' can't calculate connectivity for this atom!",i+1,TempAtom1);
            gomp_PrintWARNING(OutText);
        }
    }

/* link to the display system !!!!!!!!!!!! */
    (void)gomp_PutText2StatusLine2(0);
    return(0);
}


/***********************************************************************/
int gomp_CalcAtomConn(int Wstr)
    /* calculates the atom connection  
       a window of 'serach_window' atoms backward and forward is used
       to find the connection. if you need more change
       the "grid" parameter.
       the algorithm is very simple. all atoms have a defined
       covalent radius rc(i). a distance r(i,j) is calculated
       between atoms i and j. there is a covalent bond between
       atom i and j if the following is true
       rc(i)+rc(j) - delta < r(i,j) < rc(i)+rc(j) + delta
       where the delta is a small value (look for COVdelta).

       Leif Laaksonen 1994       */
/***********************************************************************/
{
    int i,j,k;
    int grid, last, maxConn;
    float ipj, Atom_max_corad, Atom_max_corad2;
    int from,to;
    register float w1,w2,w3,r;
    register const float *X_p,*Y_p,*Z_p,*Covar_p;
    char   OutText[BUFF_LEN];
    int *AtmConn1;
    int *AtmConn2;
    const char *Value;
    float  CovalentDelta;

    grid = ConnSearch_window;

    from = 0;
    to   = gomp_GetNumAtomsInMolecStruct(Wstr);

    if(to - from < 1 ) return(0); /* no atoms */

    Value = Tcl_GetVar(gomp_GetTclInterp(),"gomCovalentDelta",TCL_GLOBAL_ONLY);

    if(Value) {
        CovalentDelta = atof(Value);
    } else {
        CovalentDelta = COVdelt;
    }

    sprintf(OutText,"Calculating atom connection matrix for structure %d ...",Wstr+1);
    gomp_PrintMessage(OutText);

    X_p = gomp_GetAtomXCoordPointer(Wstr);
    Y_p = gomp_GetAtomYCoordPointer(Wstr);
    Z_p = gomp_GetAtomZCoordPointer(Wstr);
    Covar_p = gomp_GetAtomCovarPointer(Wstr);

    maxConn = gomp_GetMaxAtomConnections();

/* kill old connection table */
    for( i = from ; i < to ; i++) {
        AtmConn1 = gomp_GetModifiableAtomConnection(Wstr, i);
        for(j = 0 ; j < maxConn ; j++) 
            AtmConn1[j]=0;
    }

    for( i = from ; i < to ; i++) {

        AtmConn1 = gomp_GetModifiableAtomConnection(Wstr, i);

/* Do only the search in a grid   */

        last  = ( (i + grid) < to   ? (i + grid) : to  );

        for( j = i + 1 ; j < last ; j++) {

            AtmConn2 = gomp_GetModifiableAtomConnection(Wstr, j);

            Atom_max_corad  = Covar_p[i] + Covar_p[j] + CovalentDelta;
            Atom_max_corad2 = Atom_max_corad * Atom_max_corad;

            w1 = X_p[i] - X_p[j];
            r = w1 * w1;
            if(r > Atom_max_corad2) continue;
            w2 = Y_p[i] - Y_p[j];
            r += w2 * w2;
            if(r > Atom_max_corad2) continue;
            w3 = Z_p[i] - Z_p[j];
            r += w3 * w3;
            if(r > Atom_max_corad2) continue;

            ipj = Covar_p[i] + Covar_p[j];

            r = (float)sqrt(r);

            if((r < ipj + CovalentDelta) && (r > ipj - CovalentDelta) ) {

/*  Special case nr 1:  */
/*  Hydrogens can be close to each other without a bond */

                if(AtmConn1[0]+1 >= maxConn || AtmConn2[0]+1 >= maxConn) {

                    gomp_PrintMessage("*** ERROR . Max atom connections exceeded");

                    sprintf(OutText,"Atom nr: %d, '%.4s:%.4s:%.4s'",(i+1),
                            gomp_GetAtomSegName(Wstr , i),
                            gomp_GetAtomResName(Wstr , i),
                            gomp_GetAtomAtmName(Wstr , i));
                    gomp_PrintMessage(OutText);

                    gomp_PrintMessage("Connecting to atoms:");
                    for(k = 1; k <= AtmConn1[0]; k++) {
                        sprintf(OutText,"Atom nr: %d, '%.4s:%.4s:%.4s'",(AtmConn1[k]+1),
                                gomp_GetAtomSegName(Wstr , AtmConn1[k]),
                                gomp_GetAtomResName(Wstr , AtmConn1[k]),
                                gomp_GetAtomAtmName(Wstr , AtmConn1[k]));
                        gomp_PrintMessage(OutText);
                    }

                    sprintf(OutText,"Atom nr: %d, '%.4s:%.4s:%.4s'",(j+1),
                            gomp_GetAtomSegName(Wstr , j),
                            gomp_GetAtomResName(Wstr , j),
                            gomp_GetAtomAtmName(Wstr , j));
                    gomp_PrintMessage(OutText);

                    gomp_PrintMessage("Connecting to atoms:");
                    for(k = 1; k <= AtmConn2[0]; k++) {
                        sprintf(OutText,"Atom nr: %d, '%.4s:%.4s:%.4s'",(AtmConn2[k]+1),
                                gomp_GetAtomSegName(Wstr , AtmConn2[k]),
                                gomp_GetAtomResName(Wstr , AtmConn2[k]),
                                gomp_GetAtomAtmName(Wstr , AtmConn2[k]));
                        gomp_PrintMessage(OutText);
                    }

                    gomp_PrintERROR("Problems with your coordinates.\nYou have reached max number of currently allowed bonds.\nPlease check internal format.\nNew value can be defined with the command:\n'define atom maxco Ivalue'.\nWill now restart gOpenMol!");
                    (void)gomp_ResetgOpenMol();
                    return(1);
                }

                AtmConn1[++AtmConn1[0]]=j;
                AtmConn2[++AtmConn2[0]]=i;
            }
        }
    }

    gomp_PrintMessage("Done!");

    return(0);
}

/*                                                                      */
/***********************************************************************/
int gomp_CalcAtomConnI(int Wstr , int Atom)
    /* calculates the atom connection  
       a window of 'serach_window' atoms backward and forward is used
       to find the connection. if you need more change
       the "grid" parameter.
       the algorithm is very simple. all atoms have a defined
       covalent radius rc(i). a distance r(i,j) is calculated
       between atoms i and j. there is a covalent bond between
       atom i and j if the following is true
       rc(i)+rc(j) - delta < r(i,j) < rc(i)+rc(j) + delta
       where the delta is a small value (look for COVdelta).

       Leif Laaksonen 1994

       adds also connections from other atoms to the atom (Wstr,Atom).

       Eero Häkkinen 2002    */
/***********************************************************************/
{
    int i,j,k;
    int grid, first , last, maxConn;
    float ipj, Atom_max_corad, Atom_max_corad2;
    register float w1,w2,w3,r;
    register const float *X_p,*Y_p,*Z_p,*Covar_p;
    char OutText[BUFF_LEN];
    int *AtmConn1;
    int *AtmConn2;
    const char *Value;
    float  CovalentDelta;

    grid = ConnSearch_window;

    if(gomp_GetNumAtomsInMolecStruct(Wstr) < 1 ) return(0); /* no atoms */

    Value = Tcl_GetVar(gomp_GetTclInterp(),"gomCovalentDelta",TCL_GLOBAL_ONLY);

    if(Value) {
        CovalentDelta = atof(Value);
    } else {
        CovalentDelta = (float)COVdelt;
    }

    X_p = gomp_GetAtomXCoordPointer(Wstr);
    Y_p = gomp_GetAtomYCoordPointer(Wstr);
    Z_p = gomp_GetAtomZCoordPointer(Wstr);
    Covar_p = gomp_GetAtomCovarPointer(Wstr);

    maxConn    = gomp_GetMaxAtomConnections();

    i    = Atom;

    sprintf(OutText,"Calculating atom connection matrix for structure nr: %d and atom nr: %d ...",Wstr+1,Atom+1);
    gomp_PrintMessage(OutText);

/* kill old connection table */
    k=0;
    AtmConn1 = gomp_GetModifiableAtomConnection(Wstr, i);
    for(j = AtmConn1[0]; j > 0; j--)
        gomp_BreakBond(0, Wstr, i, AtmConn1[j]);
    for(j = 0 ; j < maxConn ; j++)
        AtmConn1[j]=0;

/* Do only the search in a grid   */

    first = ( (i - grid) > 0
              ? (i - grid) : 0);
    last  = ( (i + grid) < gomp_GetNumAtomsInMolecStruct(Wstr)
              ? (i + grid) : gomp_GetNumAtomsInMolecStruct(Wstr));

    for( j = first ; j < last ; j++) {

        if( i == j )
            continue;

        Atom_max_corad  = Covar_p[i] + Covar_p[j] + CovalentDelta;
        Atom_max_corad2 = Atom_max_corad * Atom_max_corad;

        w1 = X_p[i] - X_p[j];
        r = w1 * w1;
        if(r > Atom_max_corad2)
            continue;
        w2 = Y_p[i] - Y_p[j];
        r += w2 * w2;
        if(r > Atom_max_corad2)
            continue;
        w3 = Z_p[i] - Z_p[j];
        r += w3 * w3;
        if(r > Atom_max_corad2)
            continue;

        ipj = Covar_p[i] + Covar_p[j];

        r = (float)sqrt(r);

        if((r < ipj + CovalentDelta) && (r > ipj - CovalentDelta) ) {

/*  Special case nr 1:  */
/*  Hydrogens can be close to each other without a bond */

/*
  printf("k: %d , l: %d %f %f %f %f %f %f %f\n",k,l,r,X_p[j],Y_p[j],Z_p[j],
  X_p[i],Y_p[i],Z_p[i]);
*/
            if( AtmConn1[0]+1 >= maxConn) {
                gomp_PrintMessage("*** ERROR . Max atom connections exceeded");
                sprintf(OutText,"Atom nr: %d, '%.4s:%.4s:%.4s'",(i+1),
                        gomp_GetAtomSegName(Wstr , i),
                        gomp_GetAtomResName(Wstr , i),
                        gomp_GetAtomAtmName(Wstr , i));
                gomp_PrintMessage(OutText);

                gomp_PrintMessage("Connecting to atoms:");
                for(k = 1; k <= AtmConn1[0]; k++) {
                    sprintf(OutText,"Atom nr: %d, '%.4s:%.4s:%.4s'",(AtmConn1[k]+1),
                            gomp_GetAtomSegName(Wstr , AtmConn1[k]),
                            gomp_GetAtomResName(Wstr , AtmConn1[k]),
                            gomp_GetAtomAtmName(Wstr , AtmConn1[k]));
                    gomp_PrintMessage(OutText);
                }

                gomp_PrintERROR("Problems with your coordinates.\nYou have reached max number of currently allowed bonds.\nPlease check internal format.\nNew value can be defined with the command:\n'define atom maxco Ivalue'.\nWill now restart gOpenMol!");
                (void)gomp_ResetgOpenMol();
                return(1);
            }

            if( gomp_CheckBond(0, Wstr, j, i) ) {
                AtmConn2 = gomp_GetModifiableAtomConnection(Wstr, j);

                if( AtmConn2[0]+1 >= maxConn) {
                    gomp_PrintMessage("*** ERROR . Max atom connections exceeded");
                    sprintf(OutText,"Atom nr: %d, '%.4s:%.4s:%.4s'",(j+1),
                            gomp_GetAtomSegName(Wstr , j),
                            gomp_GetAtomResName(Wstr , j),
                            gomp_GetAtomAtmName(Wstr , j));
                    gomp_PrintMessage(OutText);

                    gomp_PrintMessage("Connecting to atoms:");
                    for(k = 1; k <= AtmConn1[0]; k++) {
                        sprintf(OutText,"Atom nr: %d, '%.4s:%.4s:%.4s'",(AtmConn2[k]+1),
                                gomp_GetAtomSegName(Wstr , AtmConn2[k]),
                                gomp_GetAtomResName(Wstr , AtmConn2[k]),
                                gomp_GetAtomAtmName(Wstr , AtmConn2[k]));
                        gomp_PrintMessage(OutText);
                    }

                    gomp_PrintERROR("Problems with your coordinates.\nYou have reached max number of currently allowed bonds.\nPlease check internal format.\nNew value can be defined with the command:\n'define atom maxco Ivalue'.\nWill now restart gOpenMol!");
                    (void)gomp_ResetgOpenMol();
                    return(1);
                }

                AtmConn2[++AtmConn2[0]]=i;
            }

            AtmConn1[++AtmConn1[0]]=j;
        }
    }

    gomp_PrintMessage("Done!");
    return(0);
}


/************************************************************************/
int gomp_PostReadAtoms(int What,int Wstr)
/************************************************************************/
{
    const char *value;

    if(Wstr < 0) {
        gomp_PrintMessage("Molecule structure index number out of allowed range");
        return(1);
    }

    value  = Tcl_GetVar(gomp_GetTclInterp() , 
                                "lulIdentifyAtoms", 
                                TCL_GLOBAL_ONLY);

    if(value && gomp_StringMatch(value,"no"))
        gomp_PrintMessage("Atoms are not identified");
    else {
        if(value && !gomp_StringMatch(value,"yes"))
            gomp_PrintERROR("unknown option for 'lulIdentifyAtoms' found. Must be 'yes/no'\n"
                         "Will identify atoms");
        if(gomp_IdentifyAtoms(Wstr)) {
            gomp_PrintMessage("$gOpenMol - ERROR");
            gomp_PrintERROR("Can't identify atoms");
            return(1);
        }
    }


    value  = Tcl_GetVar(gomp_GetTclInterp() , 
                                "lulCalculateAtomConnectivity", 
                                TCL_GLOBAL_ONLY);

    if(value && gomp_StringMatch(value,"no"))
        gomp_PrintMessage("Atom connectivity is not calculated");
    else {
        if(value && !gomp_StringMatch(value,"yes"))
            gomp_PrintERROR("unknown option for 'lulCalculateAtomConnectivity' found. Must be 'yes/no'\n"
                         "Will calculate connectivity");
        if(gomp_CalcAtomConn(Wstr)) {
            gomp_PrintMessage("$gOpenMol - ERROR");
            gomp_PrintERROR("Can't calculate connectivity");
            return(1);
        }
    }

    if(What) {  /* this has to be done only when a new file is read */

        (void)gomp_CalculateGeometricalCenter(Wstr);

        if(CalculateMoleculeDimensions()) {
            gomp_PrintMessage("$gOpenMol - ERROR");
            gomp_PrintERROR("Can't calculate molecular system dimensions");
            return(1);
        }
    }

    return(0);
}
/***********************************************************************/
int gomp_ReadAtomConversionTable(const char *File)
/***********************************************************************/
{
    int  Loop;
    FILE *File_p;
    char FileName[BUFF_LEN];
    char InLine[BUFF_LEN];
    char TempAtom[BUFF_LEN];
    char TempEq[BUFF_LEN];
    char TempColour[BUFF_LEN];
    char TempBasisSet[BUFF_LEN];
    int  TempType;

    if((File_p = fopen(File,"r")) == NULL) {

#if defined(WIN32)
        sprintf(FileName,"%s\\%s",gomp_ShowDataDir(),File);
#else
        sprintf(FileName,"%s/%s",gomp_ShowDataDir(),File);
#endif

        File_p = fopen(FileName,"r");
        if(File_p == NULL) {
            sprintf(InLine,"**** ERROR. Can't open input file : '%s'",FileName);
            (void)gomp_PrintMessage(InLine);
            return(1);
        }
    }

/* ready to continue ... */

    if(!AtomConversionTable.Reserved) { /* no atoms so far */
        if(ExpandAtomConversionTable()) {
            gomp_PrintMessage("?gOpenMol - ERROR");
            gomp_PrintEXIT("Can't build atom conversion list");
        }
    }

    Loop = 0;

#ifdef DEBUG
    sprintf(InLine,"Reading atom comversion table: %s ...",FileName);
    gomp_PrintMessage(InLine);
#endif

    while(fgets(InLine,BUFF_LEN,File_p) != NULL) {

        if((InLine[0] == '#') || (InLine[1] == '#') ||
           (InLine[0] == '\0'))
            continue;  /* comments */

        if(!sscanf(InLine,"%s %s %d %s %s",
                   TempAtom,TempEq,&TempType,TempColour,TempBasisSet))
            continue;

        if(Loop == AtomConversionTable.Reserved) {
            if(ExpandAtomConversionTable()) {
                gomp_PrintMessage("?gOpenMol - ERROR");
                gomp_PrintEXIT("Can't build atom conversion list");
            }
        }

        if(TempType < 1)
            TempType = DEFAULT_ATOM_TYPE;
        if(TempColour[0] == '\0')
            gomp_CopyString(TempColour,DEFAULT_ATOM_COLOUR,BUFF_LEN);

        if(PutIntoAtomConversionTable(Loop , TempAtom , TempEq , TempType,
                                      TempColour , TempBasisSet)) {
            gomp_PrintMessage("?gOpenMol - ERROR");
            gomp_PrintEXIT("Can't set symbols into conversion table");
        }

        Loop++;

        AtomConversionTable.Atoms = Loop;
    }

    rewind(File_p);
    fclose(File_p);

#ifdef DEBUG
    gomp_PrintMessage(" done\n");
#endif 

    return(0);
}
/***********************************************************************/
const AtomHit *LookIntoAtomConversionTable(const char *Atom)
/***********************************************************************/
{
    register int i,j;
    static   AtomHit OutText;
    static   char Text[MAX_ATM_NAME_LEN + 1];
     

    OutText.AtomName[0] = '\0';
    strncpy(OutText.AtomColour, "cyan" , BUFF_LEN-1);
    OutText.AtomType = 499;

    j = 0;

    for(i = 0 ; i < AtomsInAtomConversionTable() ; i++) {

        strncpy(Text,AtomConversionTable.Equiv + j,MAX_ATM_NAME_LEN);

        if(AtomConversionTable.ContInt[i]) {
            if(!AtomCompareWithNumbers(Atom , Text)) {
                strncpy(OutText.AtomName,AtomConversionTable.Atom + 
                        (MAX_ATM_NAME_LEN * i) , MAX_ATM_NAME_LEN);
                strncpy(OutText.AtomColour,AtomConversionTable.Colour +
                        (BUFF_LEN * i) , BUFF_LEN-1);
                strncpy(OutText.AtomBasisSet,AtomConversionTable.BasisSet +
                        (BUFF_LEN * i) , BUFF_LEN-1);
                OutText.AtomType       = AtomConversionTable.AtomType[i];
                OutText.AtomTypeIndex  = AtomConversionTable.AtomTypeIndex[i];
                OutText.AtomIndex      = AtomConversionTable.LabelIndex[i];

                return(&OutText);
            }
        } else {
            if(!AtomCompare(Atom , Text)) {
                strncpy(OutText.AtomName,AtomConversionTable.Atom + 
                        (MAX_ATM_NAME_LEN * i) , MAX_ATM_NAME_LEN);
                strncpy(OutText.AtomColour,AtomConversionTable.Colour +
                        (BUFF_LEN * i) , BUFF_LEN-1);
                strncpy(OutText.AtomBasisSet,AtomConversionTable.BasisSet +
                        (BUFF_LEN * i) , BUFF_LEN-1);
                OutText.AtomType       = AtomConversionTable.AtomType[i];
                OutText.AtomTypeIndex  = AtomConversionTable.AtomTypeIndex[i];
                OutText.AtomIndex      = AtomConversionTable.LabelIndex[i];

                return(&OutText);
            }
        }

        j += MAX_ATM_NAME_LEN;
    }

    return(&OutText);
}

/***********************************************************************/
const AtomHit *LookIntoAtomConversionTableDirectHit(const char *Atom)
/***********************************************************************/
{
    register int i,j;
    static   AtomHit OutText;
    static   char Text[MAX_ATM_NAME_LEN + 1];
     

    OutText.AtomName[0] = '\0';
    strcpy(OutText.AtomColour, "cyan");
    OutText.AtomType = 499;

    j = 0;

    for(i = 0 ; i < AtomsInAtomConversionTable() ; i++) {

        strncpy(Text,AtomConversionTable.Atom + j,MAX_ATM_NAME_LEN);

        if(!AtomCompare(Atom , Text)) {
            strncpy(OutText.AtomName, Text, MAX_ATM_NAME_LEN);
            strncpy(OutText.AtomColour,AtomConversionTable.Colour +
                    (BUFF_LEN * i) , BUFF_LEN-1);
            strncpy(OutText.AtomBasisSet,AtomConversionTable.BasisSet +
                    (BUFF_LEN * i) , BUFF_LEN-1);
            OutText.AtomType       = AtomConversionTable.AtomType[i];
            OutText.AtomTypeIndex  = AtomConversionTable.AtomTypeIndex[i];
            OutText.AtomIndex      = AtomConversionTable.LabelIndex[i];

            return(&OutText);
        }

        j += MAX_ATM_NAME_LEN;
    }

    return(&OutText);
}

/***********************************************************************/
int AtomsInAtomConversionTable()
/***********************************************************************/
{
    return(AtomConversionTable.Atoms);
}
/***********************************************************************/
int PutIntoAtomConversionTable(int Place , const char *Atom , 
                               const char *Equiv ,
                               int         AtomType ,
                               const char *AtomColour,
                               const char *AtomBasisSet)
/***********************************************************************/
{
    static int i,j;

    strncpy(AtomConversionTable.Atom+(MAX_ATM_NAME_LEN*Place) ,
            Atom , MAX_ATM_NAME_LEN);

    strncpy(AtomConversionTable.Equiv+(MAX_ATM_NAME_LEN*Place) ,
            Equiv , MAX_ATM_NAME_LEN);

    AtomConversionTable.ContInt[Place]      = 0;

    j = strlen(Equiv) > MAX_ATM_NAME_LEN ? MAX_ATM_NAME_LEN : strlen(Equiv);

    for( i = 0 ; i < j ; i++) {

        if(isdigit(Equiv[i])) {
            AtomConversionTable.ContInt[Place]      = 1;
            break;
        }
    }

    AtomConversionTable.LabelIndex[Place]      = gomp_MatchAtom(Atom);

    AtomConversionTable.AtomType[Place]        = AtomType;

    AtomConversionTable.AtomTypeIndex[Place]   = gomp_GetAtomPointerForIndex(AtomType);

    strncpy(AtomConversionTable.Colour+(BUFF_LEN*Place) ,
            AtomColour , BUFF_LEN-1);

    strncpy(AtomConversionTable.BasisSet+(BUFF_LEN*Place) ,
            AtomBasisSet , BUFF_LEN-1);


    return(0);
}

/***********************************************************************/
int ExpandAtomConversionTable()
/***********************************************************************/
{
    if(!AtomConversionTable.Reserved) {
        AtomConversionTable.Reserved = AtomConversionTable.Step;
        AtomConversionTable.Atom  = 
            gomp_AllocateCharVector(MAX_ATM_NAME_LEN * 
                              AtomConversionTable.Step);
        AtomConversionTable.Equiv = 
            gomp_AllocateCharVector(MAX_ATM_NAME_LEN *
                              AtomConversionTable.Step);

        AtomConversionTable.ContInt  = 
            gomp_AllocateIntVector(AtomConversionTable.Step);

        AtomConversionTable.AtomType = 
            gomp_AllocateIntVector(AtomConversionTable.Step);

        AtomConversionTable.AtomTypeIndex = 
            gomp_AllocateIntVector(AtomConversionTable.Step);

        AtomConversionTable.LabelIndex = 
            gomp_AllocateIntVector(AtomConversionTable.Step);

        AtomConversionTable.Colour =
            gomp_AllocateCharVector(BUFF_LEN *
                              AtomConversionTable.Step);
        AtomConversionTable.BasisSet =
            gomp_AllocateCharVector(BUFF_LEN *
                              AtomConversionTable.Step);
    }
    else {
        AtomConversionTable.Reserved += AtomConversionTable.Step;
        AtomConversionTable.Atom  = 
            gomp_ReallocateCharVector(
                AtomConversionTable.Atom,
                MAX_ATM_NAME_LEN * AtomConversionTable.Reserved);
        AtomConversionTable.Equiv = 
            gomp_ReallocateCharVector(
                AtomConversionTable.Equiv,
                MAX_ATM_NAME_LEN * AtomConversionTable.Reserved);

        AtomConversionTable.ContInt  =
            gomp_ReallocateIntVector(
                AtomConversionTable.ContInt , 
                AtomConversionTable.Reserved);

        AtomConversionTable.AtomType =
            gomp_ReallocateIntVector(
                AtomConversionTable.AtomType , 
                AtomConversionTable.Reserved);

        AtomConversionTable.AtomTypeIndex =
            gomp_ReallocateIntVector(
                AtomConversionTable.AtomTypeIndex , 
                AtomConversionTable.Reserved);

        AtomConversionTable.LabelIndex =
            gomp_ReallocateIntVector(
                AtomConversionTable.LabelIndex , 
                AtomConversionTable.Reserved);

        AtomConversionTable.Colour =
            gomp_ReallocateCharVector(
                AtomConversionTable.Colour ,
                BUFF_LEN * AtomConversionTable.Reserved);

        AtomConversionTable.BasisSet =
            gomp_ReallocateCharVector(
                AtomConversionTable.BasisSet ,
                BUFF_LEN * AtomConversionTable.Reserved);
    }

    if(AtomConversionTable.Atom            == NULL ||
       AtomConversionTable.Equiv           == NULL ||
       AtomConversionTable.ContInt         == NULL ||
       AtomConversionTable.AtomType        == NULL ||
       AtomConversionTable.AtomTypeIndex   == NULL ||
       AtomConversionTable.LabelIndex      == NULL ||
       AtomConversionTable.BasisSet        == NULL ||
       AtomConversionTable.Colour          == NULL) {
        return(1);
    }

    return(0);
}
/***********************************************************************/
int gomp_MatchAtom(const char *Text)
/***********************************************************************/
{
    int i;

#ifndef ELEM_PARAM
    for(i = 0 ; i < (NumAtomSymbols = strlen(AtomSymbols)/MAX_ATM_NAME_LEN) ; i++) {

        if(!AtomCompare(Text , AtomSymbols+(MAX_ATM_NAME_LEN * i))) {
            return(i+1);
        }
    }
#else
    for( i=0; i<ElementData.ElementCount; i++ ) {
        if(!AtomCompare(Text,ElementData.Elements[i].Symbol))
            return(i+1);
    }
#endif

    return(0);
}

/* On return:

= 0 ; strings are equal
!= 0 ; strings are not equal

*/
/***********************************************************************/
int AtomCompare(const char *a, const char *b)
/***********************************************************************/
{
    register const char *p1;
    register const char *p2;
    register int   iMax;

    p1   = a;
    p2   = b;
    iMax = 0; 


    while(*p1 != '\0') {

        if(iMax >= MAX_ATM_NAME_LEN) return(0); 

        if(isdigit(*p1) ||
           isspace(*p1)) {
            p1++;
            continue;
        }

        if((*p2 == '\0') ||
           (*p2 == ' ')) return(1);
        if(*p1 != *p2)   return(1);
        p1++;
        p2++;
    }

    if(*p2 == '\0') 
        return(0);
    else if(*p2 == ' ')
        return(0);
    else
        return(1);
}
/***********************************************************************/
int AtomCompareWithNumbers(const char *a, const char *b)
/***********************************************************************/
{
    register const char *p1;
    register const char *p2;
    register int   iMax;

    p1   = a;
    p2   = b;
    iMax = 0; 


    while(*p1 != '\0') {

        if(iMax >= MAX_ATM_NAME_LEN) return(0); 

        if(isspace(*p1)) {
            p1++;
            continue;
        }

        if((*p2 == '\0') ||
           (*p2 == ' ')) return(1);
        if(*p1 != *p2)   return(1);
        p1++;
        p2++;
    }

    if(*p2 == '\0') 
        return(0);
    else if(*p2 == ' ')
        return(0);
    else
        return(1);
}

/***********************************************************************/
const char *gomp_GetAtomSymbol(int Elem)
/***********************************************************************/
{
#ifndef ELEM_PARAM
    return(AtomSymbols+MAX_ATM_NAME_LEN*Elem);
#else
    return(ElementData.Elements[Elem].Symbol);
#endif
}
/***********************************************************************/
int   GetAtomIndexInTable(int Elem)
/***********************************************************************/
{
#ifndef ELEM_PARAM
    register int i;

    for(i = 0 ; i < gomp_GetNumberOfAtomSymbols() ; i++) {

        if(Elem == AtomSymbol_p[i]) return(i);

    }

    return(-1);
#else
    if(Elem>=0 && Elem<ElementData.ElementCount)
        return(Elem);
    return(-1);
#endif
}
/***********************************************************************/
int gomp_GetNumberOfAtomSymbols()
/***********************************************************************/
{
#ifndef ELEM_PARAM
    return(strlen(AtomSymbols)/4);
#else
    return(ElementData.ElementCount);
#endif
}
/***********************************************************************/
int AssignAtomInfo(int Wstr , int Hit , const AtomHit *AtomInfo)
/***********************************************************************/
{

    static int         SmashHit;
    static const char *BSTag;

/*  if AtomInfo == NULL the info is from a dictionary file */
    if(AtomInfo == NULL) {
        SmashHit = gomp_GetAtomPointerForIndex(gomp_GetAtomType(Wstr , Hit));
        (void)gomp_AssignAtomProperties(  Wstr , Hit , SmashHit);
        return(0);
    }

    SmashHit = AtomInfo->AtomTypeIndex;

/* assign properties and colour to atom number 'Hit'   */

    (void)gomp_AssignAtomProperties(  Wstr , Hit , SmashHit);
    (void)gomp_AssignAtomColourProperties(Wstr , Hit , AtomInfo->AtomColour);
    BSTag = gomp_GetAtomBasisSetTag(Wstr , Hit);
/* if the basis set tag has a value don't put a new one */
    if(BSTag[0] == '\0')
        (void)gomp_AssignAtomBasisSet(Wstr     , Hit , AtomInfo->AtomBasisSet);

    return(0);
}

/***********************************************************************/
int AssignColour2Atoms(int Wstr , int Hit , const AtomHit *AtomInfo)
/***********************************************************************/
{

/*  if AtomInfo == NULL the info is from a dictionary file */
    if(AtomInfo == NULL) {
        return(0);
    }

/* assign colour to atom number 'Hit'   */

    (void)gomp_AssignAtomColourProperties(Wstr , Hit , AtomInfo->AtomColour);

    return(0);
}
#if 0
/***********************************************************************/
int      gomp_AtomConversionTableAtoms() /* atoms in the atom conversion
                                          tabel */
/***********************************************************************/
{
    return(AtomConversionTable.Atoms);
}
/***********************************************************************/
const char *gomp_AtomConversionTableAtom(int Place)
/***********************************************************************/
{
    static char Temp[MAX_ATM_NAME_LEN];

    strncpy(Temp,&AtomConversionTable.Atom[Place * MAX_ATM_NAME_LEN] ,
            MAX_ATM_NAME_LEN);

    return(Temp);
}  
/***********************************************************************/
const char *gomp_AtomConversionTableEquiv(int Place)
/***********************************************************************/
{
    static char Temp[MAX_ATM_NAME_LEN];

    strncpy(Temp,&AtomConversionTable.Equiv[Place * MAX_ATM_NAME_LEN] ,
            MAX_ATM_NAME_LEN);

    return(Temp);
}  

/***********************************************************************/
int     gomp_AtomConversionTableType(int Place)
/***********************************************************************/
{
    return(AtomConversionTable.AtomType[Place]);
}
/***********************************************************************/
const char *gomp_AtomConversionTableColour(int Place)
/***********************************************************************/
{
    static char Temp[BUFF_LEN];

    strncpy(Temp,&AtomConversionTable.Colour[Place * BUFF_LEN] , BUFF_LEN-1);

    return(Temp);
}  

/***********************************************************************/
const char *gomp_AtomConversionTableBS(int Place)
/***********************************************************************/
{
    static char Temp[BUFF_LEN];

    strncpy(Temp,&AtomConversionTable.BasisSet[Place * BUFF_LEN] ,
            BUFF_LEN-1);

    return(Temp);
}  
#endif
/*
  Calculate the geometrical center of the molecule or the
  total system.

  If 'Wstr' is >= 0 calculate the geometrical center for that
  structure and if 'Wstr' < 0 calculate the geometrical center
  for the total system.

  Leif Laaksonen 1994

*/
/***********************************************************************/
const float *gomp_CalculateGeometricalCenter(int Wstr)
/***********************************************************************/
{
    register int i,j;
    register float *XC,*YC,*ZC;
    static float  TotalXYZ[3],TotalSUM;
    static const float *TempXYZ;
    static int from,to;
    static char OutText[BUFF_LEN];

    TotalXYZ[0] = TotalXYZ[1] = TotalXYZ[2] = 0.0;
    TotalSUM = 0.0;

    if(Wstr > gomp_GetNumMolecStructs()) {
        gomp_PrintERROR(
            "$Structure number specified is > than number of structures)");
        return(TotalXYZ);
    }

/* system physical translation is off */
    if(!gomp_GetSystemTranslateState()) { 
        TempXYZ = gomp_CalculateGeometricalCenterMT( -1);
        (void)gomp_SaveTranslateArray(TotalXYZ[0] , 
                                    TotalXYZ[1] ,
                                    TotalXYZ[2]);
        return(TempXYZ);
    }

/* just one structure ... */
    if(gomp_GetNumMolecStructs() == 1) {

        from = 0;
        to   = gomp_GetNumAtomsInMolecStruct(0);

        XC = gomp_GetModifiableAtomXCoordPointer(0);
        YC = gomp_GetModifiableAtomYCoordPointer(0);
        ZC = gomp_GetModifiableAtomZCoordPointer(0);

        for(i = from ; i < to ; i++) {
            TotalXYZ[0] += XC[i];
            TotalXYZ[1] += YC[i];
            TotalXYZ[2] += ZC[i];

            TotalSUM += 1.0;
        }
    }
    else {  /* there are several structures ... */
        TempXYZ = gomp_GetTranslateArray();
        for(j = 0 ; j < gomp_GetNumMolecStructs() - 1 ; j++) {
            from = 0;
            to   = gomp_GetNumAtomsInMolecStruct(j);

            XC = gomp_GetModifiableAtomXCoordPointer(j);
            YC = gomp_GetModifiableAtomYCoordPointer(j);
            ZC = gomp_GetModifiableAtomZCoordPointer(j);

            /* put the translation back to structures  */

            for(i = from ; i < to ; i++) {
                XC[i]       += TempXYZ[0];
                YC[i]       += TempXYZ[1];
                ZC[i]       += TempXYZ[2];

                TotalXYZ[0] += XC[i];
                TotalXYZ[1] += YC[i];
                TotalXYZ[2] += ZC[i];

                TotalSUM += 1.0;
            }
/* ... done ...  */
        }

/* now calculate the leftower ... */

        j    = gomp_GetNumMolecStructs() - 1;
        from = 0;
        to   = gomp_GetNumAtomsInMolecStruct(j);

        XC = gomp_GetModifiableAtomXCoordPointer(j);
        YC = gomp_GetModifiableAtomYCoordPointer(j);
        ZC = gomp_GetModifiableAtomZCoordPointer(j);

        for(i = from ; i < to ; i++) {
            TotalXYZ[0] += XC[i];
            TotalXYZ[1] += YC[i];
            TotalXYZ[2] += ZC[i];

            TotalSUM += 1.0;
        }
    }
/* handle the 'shift' */

    TotalXYZ[0]   = TotalXYZ[0]/TotalSUM;
    TotalXYZ[1]  = TotalXYZ[1]/TotalSUM;
    TotalXYZ[2] = TotalXYZ[2]/TotalSUM;
      
    for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {
        from = 0;
        to   = gomp_GetNumAtomsInMolecStruct(j);

        XC = gomp_GetModifiableAtomXCoordPointer(j);
        YC = gomp_GetModifiableAtomYCoordPointer(j);
        ZC = gomp_GetModifiableAtomZCoordPointer(j);

        for(i = from ; i < to ; i++) {
            XC[i]   -= TotalXYZ[0];
            YC[i]  -= TotalXYZ[1];
            ZC[i] -= TotalXYZ[2];
        }
    }


    sprintf(OutText,"Will apply a physical translation (x,y,z): %f %f %f",
            TotalXYZ[0],TotalXYZ[1],TotalXYZ[2]);
    gomp_PrintMessage(OutText);

    (void)gomp_SaveTranslateArray(TotalXYZ[0] , 
                                TotalXYZ[1] ,
                                TotalXYZ[2]);  

    TempXYZ = gomp_CalculateGeometricalCenterMT(Wstr);

    return(TotalXYZ);
}
/***********************************************************************/
const float *gomp_CalculateGeometricalCenterMT(int Wstr)
/***********************************************************************/
{
    register int i,j;
    register const float *XC,*YC,*ZC;
    static float  TotalXYZ[3],TotalSUM;
    static int from,to;
    static char OutText[BUFF_LEN];

    TotalXYZ[0] = TotalXYZ[1] = TotalXYZ[2] = 0.0;

    if(Wstr > gomp_GetNumMolecStructs()) {
        gomp_PrintERROR(
            "structure number specified is > than number of structures)");
        return(TotalXYZ);
    }

    if(gomp_GetNumMolecStructs() < 1) {
        return(TotalXYZ);
    }


    if(Wstr < 0) { /* Total system */

        TotalXYZ[0] = TotalXYZ[1] = TotalXYZ[2] = 0.0;
        TotalSUM = 0.0;
/* just one structure ... */
        if(gomp_GetNumMolecStructs() == 1) {

            from = 0;
            to   = gomp_GetNumAtomsInMolecStruct(0);

            XC = gomp_GetAtomXCoordPointer(0);
            YC = gomp_GetAtomYCoordPointer(0);
            ZC = gomp_GetAtomZCoordPointer(0);

            for(i = from ; i < to ; i++) {
                TotalXYZ[0] += XC[i];
                TotalXYZ[1] += YC[i];
                TotalXYZ[2] += ZC[i];

                TotalSUM += 1.0;
            }

/* handle the 'shift' */

            TotalXYZ[0]   = TotalXYZ[0]/TotalSUM;
            TotalXYZ[1]  = TotalXYZ[1]/TotalSUM;
            TotalXYZ[2] = TotalXYZ[2]/TotalSUM;

            sprintf(OutText,"Will apply a MT translation (x,y,z): %f %f %f",
                    TotalXYZ[0],TotalXYZ[1],TotalXYZ[2]);
            gomp_PrintMessage(OutText);

            (void)gomp_SaveTranslateArrayMT(   0 , TotalXYZ[0] , 
                                             TotalXYZ[1] ,
                                             TotalXYZ[2]);  
        }
        else {  /* there are several structures ... */

            if(gomp_GetObjectCenterType()) { /* local */
                for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {
                    from = 0;
                    to   = gomp_GetNumAtomsInMolecStruct(j);

                    XC = gomp_GetAtomXCoordPointer(j);
                    YC = gomp_GetAtomYCoordPointer(j);
                    ZC = gomp_GetAtomZCoordPointer(j);

                    TotalXYZ[0] = TotalXYZ[1] = TotalXYZ[2] = 0.0;
                    TotalSUM = 0.0;

                    for(i = from ; i < to ; i++) {
                        TotalXYZ[0] += XC[i];
                        TotalXYZ[1] += YC[i];
                        TotalXYZ[2] += ZC[i];

                        TotalSUM += 1.0;
                    }
/* handle the 'shift' */

                    TotalXYZ[0]   = TotalXYZ[0]/TotalSUM;
                    TotalXYZ[1]  = TotalXYZ[1]/TotalSUM;
                    TotalXYZ[2] = TotalXYZ[2]/TotalSUM;

                    sprintf(OutText,"Will apply a MT translation #%d (x,y,z): %f %f %f",
                            j+1,TotalXYZ[0],TotalXYZ[1],TotalXYZ[2]);
                    gomp_PrintMessage(OutText);

                    (void)gomp_SaveTranslateArrayMT(   j , TotalXYZ[0] , 
                                                     TotalXYZ[1] ,
                                                     TotalXYZ[2]);  
                }
            } else { /* global */
                TotalXYZ[0] = TotalXYZ[1] = TotalXYZ[2] = 0.0;
                TotalSUM = 0.0;

                for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {
                    from = 0;
                    to   = gomp_GetNumAtomsInMolecStruct(j);

                    XC = gomp_GetAtomXCoordPointer(j);
                    YC = gomp_GetAtomYCoordPointer(j);
                    ZC = gomp_GetAtomZCoordPointer(j);

                    for(i = from ; i < to ; i++) {
                        TotalXYZ[0] += XC[i];
                        TotalXYZ[1] += YC[i];
                        TotalXYZ[2] += ZC[i];

                        TotalSUM += 1.0;
                    }
                }

/* handle the 'shift' */
                TotalXYZ[0]   = TotalXYZ[0]/TotalSUM;
                TotalXYZ[1]  = TotalXYZ[1]/TotalSUM;
                TotalXYZ[2] = TotalXYZ[2]/TotalSUM;

                sprintf(OutText,"Will apply a MT translation (x,y,z): %f %f %f",
                        TotalXYZ[0],TotalXYZ[1],TotalXYZ[2]);
                gomp_PrintMessage(OutText);

                for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {
                    (void)gomp_SaveTranslateArrayMT(   j , TotalXYZ[0] , 
                                                     TotalXYZ[1] ,
                                                     TotalXYZ[2]);  
                }
            }
        }

    } /* end of total system */
    else {

        if(gomp_GetObjectCenterType()) { /* local */
            from = 0;
            to   = gomp_GetNumAtomsInMolecStruct(Wstr);

            XC = gomp_GetAtomXCoordPointer(Wstr);
            YC = gomp_GetAtomYCoordPointer(Wstr);
            ZC = gomp_GetAtomZCoordPointer(Wstr);

            TotalXYZ[0] = TotalXYZ[1] = TotalXYZ[2] = 0.0;
            TotalSUM = 0.0;

            for(i = from ; i < to ; i++) {
                TotalXYZ[0] += XC[i];
                TotalXYZ[1] += YC[i];
                TotalXYZ[2] += ZC[i];

                TotalSUM += 1.0;
            }

/* handle the 'shift' */

            TotalXYZ[0]   = TotalXYZ[0]/TotalSUM;
            TotalXYZ[1]  = TotalXYZ[1]/TotalSUM;
            TotalXYZ[2] = TotalXYZ[2]/TotalSUM;

            sprintf(OutText,"Will apply a MT translation #%d (x,y,z): %f %f %f",
                    Wstr+1,TotalXYZ[0],TotalXYZ[1],TotalXYZ[2]);
            gomp_PrintMessage(OutText);

            (void)gomp_SaveTranslateArrayMT(   Wstr , TotalXYZ[0] , 
                                             TotalXYZ[1] ,
                                             TotalXYZ[2]);

        } else { /* global */

            TotalXYZ[0] = TotalXYZ[1] = TotalXYZ[2] = 0.0;
            TotalSUM = 0.0;

            for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {
                from = 0;
                to   = gomp_GetNumAtomsInMolecStruct(j);

                XC = gomp_GetAtomXCoordPointer(j);
                YC = gomp_GetAtomYCoordPointer(j);
                ZC = gomp_GetAtomZCoordPointer(j);

                for(i = from ; i < to ; i++) {
                    TotalXYZ[0] += XC[i];
                    TotalXYZ[1] += YC[i];
                    TotalXYZ[2] += ZC[i];

                    TotalSUM += 1.0;
                }
            }

/* handle the 'shift' */
            TotalXYZ[0]   = TotalXYZ[0]/TotalSUM;
            TotalXYZ[1]  = TotalXYZ[1]/TotalSUM;
            TotalXYZ[2] = TotalXYZ[2]/TotalSUM;

            for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {
                sprintf(OutText,"Will apply a MT translation (x,y,z): %f %f %f",
                        TotalXYZ[0],TotalXYZ[1],TotalXYZ[2]);
                gomp_PrintMessage(OutText);

                (void)gomp_SaveTranslateArrayMT(   j , TotalXYZ[0] , 
                                                 TotalXYZ[1] ,
                                                 TotalXYZ[2]);  
            }

        }

    }

    return(TotalXYZ);

}


/***********************************************************************/
int CalculateMoleculeDimensions()
/***********************************************************************/
{

    int i,j;
    float minx,maxx;
    float miny,maxy;
    float minz,maxz;
    const float *XC;
    const float *YC;
    const float *ZC;

    /* Calculate min and max of the x,y and z coordinates  */

    minx=1.e+30f; maxx= -1.e+30f;
    miny=1.e+30f; maxy= -1.e+30f;
    minz=1.e+30f; maxz= -1.e+30f;

    for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) { 
        /* handle the structures */

        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) {

            XC   = gomp_GetAtomXCoordPointer(i);
            YC  = gomp_GetAtomYCoordPointer(i);
            ZC = gomp_GetAtomZCoordPointer(i);

            if(XC[j] < minx)  minx=XC[j];
            if(XC[j] > maxx)  maxx=XC[j];

            if(YC[j] < miny)  miny=YC[j];
            if(YC[j] > maxy)  maxy=YC[j];

            if(ZC[j] < minz)  minz=ZC[j];
            if(ZC[j] > maxz)  maxz=ZC[j]; 
        }
    }


/* To prevent clipping in case of rotating a box (periodic boundary conditions)
   the clipping plane is calculated to the edge of the box ( sqrt(2)*r ) */


    MoleculeDimensions.MinX = (minx > 0.0 ? minx*0.67 : minx*1.41);
    MoleculeDimensions.MaxX = (maxx < 0.0 ? maxx*0.67 : maxx*1.41);
    MoleculeDimensions.MinY = (miny > 0.0 ? miny*0.67 : miny*1.41);
    MoleculeDimensions.MaxY = (maxy < 0.0 ? maxy*0.67 : maxy*1.41);
    MoleculeDimensions.MinZ = (minz > 0.0 ? minz*0.67 : minz*1.41);
    MoleculeDimensions.MaxZ = (maxz < 0.0 ? maxz*0.67 : maxz*1.41);


    MoleculeDimensions.Min  = MIN( minz , MIN(minx , miny));
    MoleculeDimensions.Max  = MAX( maxz , MAX(maxx , maxy));

    MoleculeDimensions.SizeOfSystem = 
        MoleculeDimensions.Max - MoleculeDimensions.Min;

    if(MoleculeDimensions.SizeOfSystem < GetMinimumSizeOfSystem())
        MoleculeDimensions.SizeOfSystem = GetMinimumSizeOfSystem();

    MoleculeDimensions.TransDamping = 
        MoleculeDimensions.SizeOfSystem / 200.;
/*  ok done ...                                                     */

    return(0);
}

/***********************************************************************/
float  gomp_GetSizeOfSystem()
/***********************************************************************/
{
    if(gomp_GetNumMolecStructs())
        return(MoleculeDimensions.SizeOfSystem);
    else
        return((float)1.0);
}
/***********************************************************************/
int  gomp_SetSizeOfSystem(float Value)
/***********************************************************************/
{
    MoleculeDimensions.SizeOfSystem = Value;

    return(0);
}
#if 0
/***********************************************************************/
float gomp_TranslationDamping()
/***********************************************************************/
{
    return( MoleculeDimensions.TransDamping);
}
/***********************************************************************/
int gomp_PrintAtomInfo(int Wstr, int Place)
/***********************************************************************/
{
    char OutString[BUFF_LEN];

/*  if there are no atoms return NOW */
    if(gomp_GetNumAtomsInMolecStruct(Wstr) < 1) return(1);

    gomp_PrintMessage("** Atom info for atom **");

    sprintf(OutString,"Structure nr: %d , atom nr: %d",
            (Wstr+1),(Place+1));
    gomp_PrintMessage(OutString);

    sprintf(OutString,"Segment name: '%s', residue name: '%s', atom name: '%s'",
            gomp_GetAtomSegName(Wstr , Place),
            gomp_GetAtomResName(Wstr , Place),
            gomp_GetAtomAtmName(Wstr , Place));

    gomp_PrintMessage(OutString);

    sprintf(OutString,"Basis set tag: '%s'",gomp_GetAtomBasisSetTag(Wstr,Place));
    gomp_PrintMessage(OutString);

    return(0);
}
#endif
/***********************************************************************/
float  gomp_GetTranslationDamping()
/***********************************************************************/
{
    return(MoleculeDimensions.TransDamping);
}
/***********************************************************************/
int   gomp_SetTranslationDamping(float Value)
/***********************************************************************/
{
    MoleculeDimensions.TransDamping = Value;

    return(0);
}
/***********************************************************************/
int   gomp_ParseSetAtomCovar(const char *Atom,const char *Value)
/***********************************************************************/
{
    int   index;
    float MaxCovar;
    float Covar;

    if(isdigit(Atom[0]))
        index = atoi(Atom);
    else
        index = gomp_MatchAtom(Atom) - 1;

    if( index < 0 || index >= ElementData.ElementCount ) {
        gomp_PrintERROR("Atom index out of allowed range or unknown atom symbol name");
        return(1);
    }

    MaxCovar = sqrt(ElementData.AtomMaxCorad2);
    Covar    = (float)atof(Value);
    if(Covar > MaxCovar)
        Covar = MaxCovar;

    ElementData.Elements[index].Covar = Covar;

    return(0);
}
/***********************************************************************/
int gomp_ReadElementParams(const char *File)
/***********************************************************************/
{
    char input[BUFF_LEN];
    char temp1[BUFF_LEN];
    FILE *File_p;

    char chelp[BUFF_LEN];

    if((File_p = fopen(File,"r")) == NULL) {

#if defined(WIN32)
        sprintf(chelp,"%s\\%s",gomp_ShowDataDir(),File);
#else
        sprintf(chelp,"%s/%s",gomp_ShowDataDir(),File);
#endif

        File_p = fopen(chelp,"r");
        if(File_p == NULL) {
            sprintf(temp1,"**** ERROR. Can't open input file : '%s'",chelp);
            (void)gomp_PrintMessage(temp1);
            return(1);
        }
    }

/*   Check if there is already a structure */
    if( ElementData.Elements ) {
        free(ElementData.Elements);
        ElementData.ElementCount = 0;
        ElementData.Elements     = NULL;
    }

    while(fgets(input,BUFF_LEN,File_p) != NULL) {

/* allow for a comment starting with '#' */
        if(input[0] == '#')
            continue;

        if(!ElementData.Elements)
            ElementData.Elements = gomp_AllocateVoidVector(
                sizeof(*ElementData.Elements));
        else
            ElementData.Elements = gomp_ReallocateVoidVector(
                ElementData.Elements,
                sizeof(*ElementData.Elements)*(ElementData.ElementCount+1));

        if(!sscanf(input,"%4s %f",
                   temp1,&ElementData.Elements[ElementData.ElementCount].Covar))
            continue;

        strncpy(
            ElementData.Elements[ElementData.ElementCount].Symbol,
            temp1, MAX_ATM_NAME_LEN );

        ElementData.ElementCount++;
    }

    fclose(File_p);

    return (0);
}

/*
  In fact this also looks for an input tcl gomp_ipt ...
*/
/***********************************************************************/
int gomp_StartWithCoordinateFile(const char *InputText0)
/***********************************************************************/
{
    char  InputText[BUFF_LEN];
    char  Temp1[BUFF_LEN];
    char  Temp2[2*BUFF_LEN];
    const char *result;

/* check for URL */
    gomp_CopyString(InputText,InputText0,BUFF_LEN);
    InputText[BUFF_LEN-1] = '\0';
    if(gomp_FileNameIsURL(InputText)) {
        sprintf(Temp1,"?Can't handle the given URL '%s'",InputText);
        gomp_PrintERROR(Temp1);
        return(1);
    }
/* determine first the type of input file */

/* gom coordinates */
    if(Tcl_StringCaseMatch(InputText,"*.gom",1)) {
        return(gomp_StartWithModelFile(InputText));
    }
/* tcl script name */
    else if(Tcl_StringCaseMatch(InputText,"*.tcl",1)) {
        gomp_FormatMessage("Reading tcl file '%s' at startup", InputText);
        if ( gomp_SendFile2TclParser(InputText) != TCL_OK ) {
            gomp_PrintMessage(Tcl_GetStringResult(gomp_GetTclInterp()));
            gomp_PrintERROR("can't evaluate supplied input tcl script");
            return(1);
        }
        return(0);
    } else {
        Tcl_Obj* list;
        list = gomp_CreateTclList(
            "%s %s %s",
            "gom::GetAtomCoordinateFileFormatByExtension",
            "import",
            InputText);
        if( Tcl_GlobalEvalObj(gomp_GetTclInterp(), list) != TCL_OK ) {
            gomp_PrintMessage(Tcl_GetStringResult(gomp_GetTclInterp()));
            gomp_PrintERROR("can't execute the coordinates script 'gom::GetAtomCoordinateFileFormatByExtension'");
            return(1);
        }
        result = Tcl_GetStringResult( gomp_GetTclInterp() );
        if( result && *result ) {
            gomp_CopyString(Temp2, result, BUFF_LEN);
            Temp2[BUFF_LEN-1] = '\0';
            if( gomp_ReadCoordinates( Temp2 , InputText , "append" ) == 0 )
                return(0);
        }
    }

    printf("$ERROR - gOpenMol does not support the file type in '%s'\n",InputText);

    return(1);
}
/***********************************************************************/
int gomp_StartWithModelFile(const char *InputText0)
/***********************************************************************/
{
    char InputText[BUFF_LEN];
    char Text[BUFF_LEN];

    gomp_CopyString(InputText,InputText0,BUFF_LEN);
    InputText[BUFF_LEN-1] = '\0';
    if(gomp_FileNameIsURL(InputText)) {
        sprintf(Text,"?Can't handle the given URL '%s'",InputText);
        gomp_PrintERROR(Text);
        return(1);
    }

    return(gomp_ReadOldModel(InputText));
}

/****************************************************************************/
const char *gomp_Number2Name(int AtomIndex)
/****************************************************************************/
{
    static int      i;
    static char     AtomName[3];

    i = GetAtomIndexInTable(AtomIndex);

    if(i < 0) perror("Error in atom index table");

    strncpy(AtomName,gomp_GetAtomSymbol(i),2);
    AtomName[2] = '\0';    /* NULL terminate always */

    return(AtomName);
}
/*
  Recalculated also the connection matrix for the selected structures
  1998-07-04 ... 2002-06-24
*/
/****************************************************************************/
int   gomp_SetSearchWindow(int SearchWindow)
/****************************************************************************/
{
    ConnSearch_window = SearchWindow;

    return(0);
}
/****************************************************************************/
int   gomp_GetConnectionSearchWindow()
/****************************************************************************/
{
    return(ConnSearch_window);
}
/****************************************************************************/
int   gomp_AssignAtomInfoFromDictionary( int Wstr , int Hit)
/****************************************************************************/
{
    AtomHit *Dummy;

    Dummy = NULL;

    (void)AssignAtomInfo( Wstr , Hit , Dummy);

    return(0);
}
#if 0
/************************************************************************/
int   gomp_ResetAtoms()
/************************************************************************/
{
    int i;

    if(!gomp_GetNumMolecStructs()) return(0);

    for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {
        if(gomp_IdentifyAtoms(i)) {
            gomp_PrintERROR("can't reset atoms");
            return(1);
        }
    }

    return(0);
}
/************************************************************************/
int   gomp_SetMinimumSizeOfSystem(float Size)
/************************************************************************/
{
    MinimumSizeOfSystem = Size;

    return(0);
}
#endif
/************************************************************************/
float  GetMinimumSizeOfSystem()
/************************************************************************/
{
    return(MinimumSizeOfSystem);
}
/************************************************************************/
int  gomp_SetBondReconnectivityState(int State)
/************************************************************************/
{

    AtomReconnectivityState.BondState = State;  

    return(0);  
}
/************************************************************************/
int  gomp_GetBondReconnectivityState()
/************************************************************************/
{
    return(AtomReconnectivityState.BondState);  
}
/************************************************************************/
int  gomp_SetHBondReconnectivityState(int State)
/************************************************************************/
{

    AtomReconnectivityState.HBondState = State;  

    return(0);  
}
/************************************************************************/
int  gomp_GetHBondReconnectivityState()
/************************************************************************/
{
    return(AtomReconnectivityState.HBondState);  
}

/************************************************************************/
int    gomp_GetBondDisplayStyle()
/************************************************************************/
{
    return(BondDisplayStyle);
}
/************************************************************************/
int    gomp_SetBondDisplayStyle(int Value)
/************************************************************************/
{
    BondDisplayStyle  = Value;

    return(0);
}


/************************************************************************/
int   gomp_IdentifyAtomColoursAllStructures()
/************************************************************************/
{
    int i;

    for( i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {

        if(IdentifyAtomColours(i))
            return(1);
    }

    return(0);
}
/***********************************************************************/
int      gomp_GetSystemTranslateState()
/***********************************************************************/
{

    return(MoleculeDimensions.TranslateState);
}
/***********************************************************************/
int      gomp_SetSystemTranslateState(int Value)
/***********************************************************************/
{
    register int i,j;
    register float *XC,*YC,*ZC;
    static float  TotalXYZ[3],TotalSUM;
    static const float *TempXYZ;
    static int from,to;
    static char OutText[BUFF_LEN];
    static float Identity[16] = {
        1.0 , 0.0 , 0.0 , 0.0,
        0.0 , 1.0 , 0.0 , 0.0,
        0.0 , 0.0 , 1.0 , 0.0,
        0.0 , 0.0 , 0.0 , 1.0};

    if(MoleculeDimensions.TranslateState == Value)
        return(0);

    if(gomp_GetNumMolecStructs() < 1) {
        MoleculeDimensions.TranslateState = Value;
        return(0);
    }

    if(!Value) {  /* translate on!  */

        TotalSUM    = 0.0;
        TotalXYZ[0] = TotalXYZ[1] = TotalXYZ[2] = 0.0;

        for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {

            from = 0;
            to   = gomp_GetNumAtomsInMolecStruct(j);

            XC = gomp_GetModifiableAtomXCoordPointer(j);
            YC = gomp_GetModifiableAtomYCoordPointer(j);
            ZC = gomp_GetModifiableAtomZCoordPointer(j);

            for(i = from ; i < to ; i++) {

                TotalXYZ[0] += XC[i];
                TotalXYZ[1] += YC[i];
                TotalXYZ[2] += ZC[i];

                TotalSUM += 1.0;
            }
        }

/* handle the 'shift' */

        TotalXYZ[0]   = TotalXYZ[0]/TotalSUM;
        TotalXYZ[1]  = TotalXYZ[1]/TotalSUM;
        TotalXYZ[2] = TotalXYZ[2]/TotalSUM;

/* do the shift */
        for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {

            from = 0;
            to   = gomp_GetNumAtomsInMolecStruct(j);

            XC = gomp_GetModifiableAtomXCoordPointer(j);
            YC = gomp_GetModifiableAtomYCoordPointer(j);
            ZC = gomp_GetModifiableAtomZCoordPointer(j);

            for(i = from ; i < to ; i++) {

                XC[i] -= TotalXYZ[0];
                YC[i] -= TotalXYZ[1];
                ZC[i] -= TotalXYZ[2];

            }
        }

        sprintf(OutText,"Will apply a translation (x,y,z): %f %f %f",
                TotalXYZ[0],TotalXYZ[1],TotalXYZ[2]);
        gomp_PrintMessage(OutText);

        (void)gomp_SaveTranslateArray(TotalXYZ[0] , 
                                    TotalXYZ[1] ,
                                    TotalXYZ[2]);  
        (void)gomp_ResetView();
        (void)gomp_SaveModelViewMatrix(Identity);

        MoleculeDimensions.TranslateState = 0;

    } else {      /* translate off! */

        TempXYZ = gomp_GetTranslateArray();

        for(j = 0 ; j < gomp_GetNumMolecStructs() ; j++) {

            from = 0;
            to   = gomp_GetNumAtomsInMolecStruct(j);

            XC = gomp_GetModifiableAtomXCoordPointer(j);
            YC = gomp_GetModifiableAtomYCoordPointer(j);
            ZC = gomp_GetModifiableAtomZCoordPointer(j);

            for(i = from ; i < to ; i++) {

                XC[i] += TempXYZ[0];
                YC[i] += TempXYZ[1];
                ZC[i] += TempXYZ[2];

            }
        }

        (void)gomp_SaveTranslateArray(0.0 , 
                                    0.0 ,
                                    0.0);

        TempXYZ = gomp_CalculateGeometricalCenterMT( - 1);

        (void)gomp_ResetView();
        (void)gomp_SaveModelViewMatrix(Identity);
        MoleculeDimensions.TranslateState = 1;
    }

    return(0);
}

