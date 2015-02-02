/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2003, 2005 by:
Eero Häkkinen

*/


/*
  This program reads a PDB coordinate file

  Leif Laaksonen 1989, rewritten and bugs fixed 1998
  2000-07-25: Added CRYST1 record (LUL)
  2001-05-24: Added the HELIX, SHEET and TURN record (LUL)
  2003-03-15: Added the MODEL/ENDMDL structure optio into the reader
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <tcl.h>

#include "cell.h"
#include "coord_file.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

#define PDB_LINE_LEN       120                  /* pdb file line length */
#define PDB_ATM_NAME_LEN     4                  /* atom name length     */
#define PDB_RES_NAME_LEN     3                  /* residue name length  */
#define PDB_SEG_NAME_LEN     1                  /* segment name length  */
#define PDB_RES_NUM_LEN      4                  /* residue number length */
#define PDB_ATM_COORD_LEN    8                  /* coord string length   */
#define PDB_TEMP_FACTOR      6

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

static int AtomsInPdbFile(FILE *);
static int ModelsInPdbFile(FILE *, int **);

/*************************************************************************/
int gomp_ReadCoordinatesPDB(const char *Text1, int Append)  
    /* Brookhaven format file reader */
/*************************************************************************/
{

    char inputl[PDB_LINE_LEN];   /* input line */

    char tmp_atm[BUFF_LEN];
    char tmp_res[BUFF_LEN];
    char tmp_seg[BUFF_LEN];
    char tmp_res_num[BUFF_LEN];
    char tmp_coord[BUFF_LEN];
    char tmp_temp[BUFF_LEN];

    char tmp_a[BUFF_LEN];
    char tmp_b[BUFF_LEN];
    char tmp_c[BUFF_LEN];
    char tmp_alpha[BUFF_LEN];
    char tmp_beta[BUFF_LEN];
    char tmp_gamma[BUFF_LEN];
    char tmp_sGroup[BUFF_LEN];
    char tmp_z[BUFF_LEN];

    char OutText[BUFF_LEN];

    float TXc,TYc,TZc,TBv;
    int   TRs1;
    int   Helix, Sheet, Turn;
    int   ITemp;

    static int i,j,loop;
    static int type_warning;
    static int PDBatoms;
    static int Wstr;
    static int RWstr;

    int  PDBmodels       = 0;
    int *PDBatomsInModel = NULL;


    FILE *pdb_in;

    type_warning = 0;

    Helix = 0;
    Sheet = 0;
    Turn  = 0;

    pdb_in=fopen(Text1,"r");
    if(pdb_in == NULL) {
        sprintf(OutText,"$Can't open input file : %s",Text1);
        gomp_PrintERROR(OutText);
        return(1);
    }

    sprintf(OutText,"********** Reading : %s **********",Text1);
    gomp_PrintMessage(OutText);

/* calculate number of models in the PDB file              */
    if((PDBmodels=ModelsInPdbFile(pdb_in,&PDBatomsInModel)) > 0) {

        RWstr = 0;
        for(i = 0 ; i < PDBmodels ; i++) {

            sprintf(OutText,"%s_%d",Text1,(i+1));
            if(!i) {

                if(Append)
                    Wstr = gomp_CreateMolecStruct(OutText , PDBatomsInModel[i] , APPEND);
                else
                    Wstr = gomp_CreateMolecStruct(OutText , PDBatomsInModel[i] , NEW);
                if ( Wstr < 0 )
                    goto end;

            } else {
                if(gomp_CreateMolecStruct(OutText , PDBatomsInModel[i] , APPEND)<0) {
                    while ( --i >= 0 )
                        gomp_DeleteMolecStruct(Wstr+i);
                    Wstr = -1;
                    goto end;
                }
            }

        }
    } else {
/* calculate number of atoms in the PDB file              */
        PDBatoms = AtomsInPdbFile(pdb_in);

        if(Append)
            Wstr = gomp_CreateMolecStruct(Text1 , PDBatoms , APPEND);
        else
            Wstr = gomp_CreateMolecStruct(Text1 , PDBatoms , NEW);
        if ( Wstr < 0 )
            goto end;
    }

/*
  Start reading file
  Tags recogniced by this routine are:

  ATOM:   Atom coordinate records for "standard groups"
  HETATM: Atom coordinate records for "non-standard" groups
  END:    End-of-entry record
*/

    loop = 0;

    while(fgets(inputl,PDB_LINE_LEN,pdb_in) != NULL) { 

/* MODEL */
        if((strncmp(inputl,"MODEL",5)   == 0)) { /*starts model*/
            loop = 0;
        }
/* ATOM and HETATM */
        if((strncmp(inputl,"ATOM",4)   == 0) ||
           (strncmp(inputl,"HETATM",6) == 0) ) { /*start atom or hetatm*/

/* check the the atom index */
            if(PDBmodels) {
                if(loop >= PDBatomsInModel[RWstr]) {
                    gomp_PrintERROR("atom indexout of allowed range");
                    free(PDBatomsInModel);
                    return(1);
                }
            }

            memset(tmp_atm , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_atm , (inputl+12) , PDB_ATM_NAME_LEN);
            memset(tmp_res , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_res , (inputl+17) , PDB_RES_NAME_LEN);
            memset(tmp_seg , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_seg , (inputl+21) , PDB_SEG_NAME_LEN);

            sscanf(tmp_atm,"%s",OutText);
            j = gomp_PutAtomAtmName(Wstr , OutText , loop);
            sscanf(tmp_res,"%s",OutText);
            j = gomp_PutAtomResName(Wstr , OutText , loop);
            sscanf(tmp_seg,"%s",OutText);
            j = gomp_PutAtomSegName(Wstr , OutText , loop);

/* check first postion 27 to see if it is present If a number length = 5 else 4 */
            if(isdigit(*(inputl+26))) {
                memset(tmp_res_num , 0 , (size_t)BUFF_LEN);
                strncpy(tmp_res_num , (inputl+22) , PDB_RES_NUM_LEN + 1);
                sscanf(tmp_res_num,"%d",&TRs1);
            } else {
                memset(tmp_res_num , 0 , (size_t)BUFF_LEN);
                strncpy(tmp_res_num , (inputl+22) , PDB_RES_NUM_LEN);
                sscanf(tmp_res_num,"%d",&TRs1);
            }

            j = gomp_PutAtomResNum1(Wstr , TRs1 , loop);
            j = gomp_PutAtomResNum2(Wstr , TRs1 , loop);

            memset(tmp_coord , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_coord , (inputl+30) , PDB_ATM_COORD_LEN);
            sscanf(tmp_coord,"%f",&TXc);
            j = gomp_PutAtomXCoord(Wstr , TXc , loop);
            memset(tmp_coord , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_coord , (inputl+38) , PDB_ATM_COORD_LEN);
            sscanf(tmp_coord,"%f",&TYc);
            j = gomp_PutAtomYCoord(Wstr , TYc , loop);
            memset(tmp_coord , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_coord , (inputl+46) , PDB_ATM_COORD_LEN);
            sscanf(tmp_coord,"%f",&TZc);
            j = gomp_PutAtomZCoord(Wstr , TZc , loop);

            memset(tmp_temp , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_temp , (inputl+60) , PDB_TEMP_FACTOR);
            sscanf(tmp_temp,"%f",&TBv);
            j = gomp_PutAtomBValue(Wstr , TBv , loop);
            j = gomp_PutAtomCharge(Wstr , 0.0 , loop);

            if(inputl[21] == ' ')
                gomp_PutAtomSegName(Wstr , "S1",loop);

            loop++;

        } /*end atom*/
        else  if(strncmp(inputl,"CRYST1",6) == 0) { /*start CRYST1 record*/
/*
  COLUMNS       DATA TYPE      FIELD         DEFINITION
  -------------------------------------------------------------
  1 -  6       Record name    "CRYST1"

  7 - 15       Real(9.3)      a             a (Angstroms).

  16 - 24       Real(9.3)      b             b (Angstroms).

  25 - 33       Real(9.3)      c             c (Angstroms).

  34 - 40       Real(7.2)      alpha         alpha (degrees).

  41 - 47       Real(7.2)      beta          beta (degrees).

  48 - 54       Real(7.2)      gamma         gamma (degrees).

  56 - 66       LString        sGroup        Space group.

  67 - 70       Integer        z             Z value.

*/
            memset(tmp_a , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_a , (inputl+6)  , 9);
            memset(tmp_b , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_b , (inputl+15) , 9);
            memset(tmp_c , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_c , (inputl+24) , 9);

            sscanf(tmp_a,"%f",&TXc);
            (void)gomp_SetCellA(TXc);
            sscanf(tmp_b,"%f",&TYc);
            (void)gomp_SetCellB(TYc);
            sscanf(tmp_c,"%f",&TZc);
            (void)gomp_SetCellC(TZc);

            memset(tmp_alpha , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_alpha , (inputl+33) , 7);
            sscanf(tmp_alpha,"%f",&TXc);
            (void)gomp_SetCellAlpha(TXc);
            memset(tmp_beta , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_beta , (inputl+40)  , 7);
            sscanf(tmp_beta,"%f",&TYc);
            (void)gomp_SetCellBeta(TYc);
            memset(tmp_gamma , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_gamma , (inputl+47) , 7);
            sscanf(tmp_gamma,"%f",&TZc);
            (void)gomp_SetCellGamma(TZc);

            memset(tmp_sGroup , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_sGroup , (inputl+55) , 11);

            memset(tmp_z , 0 , (size_t)BUFF_LEN);
            strncpy(tmp_z , (inputl+66) , 4);
            sscanf(tmp_z,"%f",&TBv);

        }
/* HELIX */
        else if((strncmp(inputl,"HELIX",5)   == 0)) {
            sprintf(OutText,"lulPlumber::SecondaryStructureSaver %d %d {%s}",Wstr+1,(Sheet+Helix+Turn+1),inputl);
            ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), OutText );
            if(ITemp != TCL_OK) {
                gomp_PrintERROR("can't execute the script 'lulPlumber::SecondaryStructureSaver'");
                continue;
            }

            Helix += 1;
        }
/* SHEET */
        else if((strncmp(inputl,"SHEET",5)   == 0)) {
            sprintf(OutText,"lulPlumber::SecondaryStructureSaver %d %d {%s}",Wstr+1,(Sheet+Helix+Turn+1),inputl);
            ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), OutText );
            if(ITemp != TCL_OK) {
                gomp_PrintERROR("can't execute the script 'lulPlumber::SecondaryStructureSaver'");
                continue;
            }

            Sheet += 1;
        }
/* TURN */
        else if((strncmp(inputl,"TURN",4)   == 0)) {
            sprintf(OutText,"lulPlumber::SecondaryStructureSaver %d %d {%s}",Wstr+1,(Sheet+Helix+Turn+1),inputl);
            ITemp = Tcl_GlobalEval(gomp_GetTclInterp(), OutText );
            if(ITemp != TCL_OK) {
                gomp_PrintERROR("can't execute the script 'lulPlumber::SecondaryStructureSaver'");
                continue;
            }

            Turn += 1;
        }

/* ENDMDL */
        if((strncmp(inputl,"ENDMDL",6)   == 0)) { /*starts model*/
            Wstr++;
            RWstr++;
            continue;
        }


/* END */
        if(strncmp(inputl,"END",3)   == 0) break; /*the end*/

    }

    if(!PDBmodels) {
        if(loop != PDBatoms) {
            gomp_PrintMessage("?ERROR - can't read correct number of atoms from PDB file");
            return(1);
        }
    }

    gomp_PrintMessage("**********   Done   **********");

end:
    fclose(pdb_in);

    free(PDBatomsInModel);

    return(Wstr >= 0 ? 0 : 1);

}
/*************************************************************************/
int AtomsInPdbFile(FILE *file_p)
/*************************************************************************/
{
    char inputl[PDB_LINE_LEN];   /* input line */
    int PDBatoms;

/* just for sure rewind the file */
    rewind(file_p);

/*
  Start reading file
  Tags recogniced by this routine are:

  ATOM:   Atom coordinate records for "standard groups"
  HETATM: Atom coordinate records for "non-standard" groups
  END:    End-of-entry record
*/
    PDBatoms = 0 ; /* start counting number of atoms */

    while(fgets(inputl,PDB_LINE_LEN,file_p) != NULL) { 

/* ATOM and HETATM */
        if((strncmp(inputl,"ATOM",4)   == 0) ||
           (strncmp(inputl,"HETATM",6) == 0) ) { /*start atom or hetatm*/

            PDBatoms++;

        } /*end atom*/

/* END */
        if(strncmp(inputl,"END",3)   == 0) break; /*the end*/

    }

/* rewind the file to be ready for the "real" processing */
    rewind(file_p);

    return(PDBatoms);
}
/*************************************************************************/
static int GetAtomsInModelSection(FILE *file_p)
/*************************************************************************/
{



    char inputl[PDB_LINE_LEN];   /* input line */
    int PDBatoms;

/*
  Start reading file
  Tags recogniced by this routine are:

  ATOM:   Atom coordinate records for "standard groups"
  HETATM: Atom coordinate records for "non-standard" groups
  END:    End-of-entry record
*/
    PDBatoms = 0 ; /* start counting number of atoms */

    while(fgets(inputl,PDB_LINE_LEN,file_p) != NULL) { 

/* ATOM and HETATM */
        if((strncmp(inputl,"ATOM",4)   == 0) ||
           (strncmp(inputl,"HETATM",6) == 0) ) { /*start atom or hetatm*/

            PDBatoms++;

        } /*end atom*/

/* ENDMDL */
        if((strncmp(inputl,"ENDMDL",6) == 0)) break; /*the end*/

    }

    return(PDBatoms);
}


/*************************************************************************/
int ModelsInPdbFile(FILE *file_p, int **PDBatomsInModel)
/*************************************************************************/
{



    char inputl[PDB_LINE_LEN];   /* input line */
    int PDBmodel;
    int PDBendmdl;

/* just for sure rewind the file */
    rewind(file_p);

/*
  Start reading file
  Tags recogniced by this routine are:

  MODEL:   Atom coordinate records for "standard groups"
  ENDMDL:  Atom coordinate records for "non-standard" groups
*/
    PDBmodel  = 0 ; /* start counting number of models */
    PDBendmdl = 0;

    while(fgets(inputl,PDB_LINE_LEN,file_p) != NULL) { 

/* MODEL */
        if((strncmp(inputl,"MODEL",5)   == 0)) { /*start MODEL */

            *PDBatomsInModel = gomp_ReallocateIntVector(*PDBatomsInModel , PDBmodel + 1);

            (*PDBatomsInModel)[PDBmodel] = GetAtomsInModelSection(file_p);

            PDBmodel++;
            PDBendmdl++;

        } /*end */
    }

/* rewind the file to be ready for the "real" processing */
    rewind(file_p);

    if(PDBmodel != PDBendmdl) {
        gomp_PrintERROR("number of 'MODEL' and 'ENDMDL' records do not match!");
        return(MIN(PDBmodel,PDBendmdl));
    }

    return(PDBmodel);
}
