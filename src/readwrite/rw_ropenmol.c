/*

Copyright (c) 1993 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002, 2005 by:
Eero Häkkinen
*/


/*
  This program reads a center ascii file

  Leif Laaksonen 1997
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <tcl.h>

#include "coord_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

#define OPENMOL_LINE_LEN   120   /* OpenMol file line length */

/* ... Stuff needed for the interface between OpenMol and gOpenMol */

static int AtomsInOpenMolFile(FILE *);

/*************************************************************************/
int gomp_ReadCoordinatesOPENMOL(const char *Text1 , int Append)  
    /* center binary file reader */
/*************************************************************************/
{

    char OutText[BUFF_LEN];
    char AtomName[BUFF_LEN];
    char inputl[BUFF_LEN];
    float NucCharge;
    float Xc,Yc,Zc;
    char BasisSetTag[BUFF_LEN];
    int  OpenMolAtoms;
    int  atoms;
    float conv = 0.52917715f;
    int Wstr;
    int  j;


    FILE *openmol_in;

    openmol_in=fopen(Text1,"r");
    if(openmol_in == NULL) {
        sprintf(OutText,"$Can't open input file : %s",Text1);
        gomp_PrintERROR(OutText);
        return(1);
    }

    sprintf(OutText,"********** Reading : %s **********",Text1);
    gomp_PrintMessage(OutText);

/* calculate number of atoms in the PDB file              */
    OpenMolAtoms = AtomsInOpenMolFile(openmol_in);

    if(OpenMolAtoms < 1) {
        gomp_PrintERROR("?ERROR - no atoms in OPENMOL file");
        return(1);
    }

/*  update atom structure */

    if(Append)
        Wstr = gomp_CreateMolecStruct(Text1 , OpenMolAtoms , APPEND);
    else
        Wstr = gomp_CreateMolecStruct(Text1 , OpenMolAtoms , NEW);
    if ( Wstr < 0 )
        goto end;

    atoms = 0;

    while(fgets(inputl,OPENMOL_LINE_LEN,openmol_in) != NULL) {

        if(Tcl_StringMatch(inputl , ": gaussian_basis_set_charge_centers*")) {

            while(fgets(inputl,OPENMOL_LINE_LEN,openmol_in) != NULL) {

                if(inputl[0] == ':') break;

                sscanf(inputl,"%s %f %f %f %f %s",AtomName,&NucCharge,&Xc,&Yc,&Zc,
                       BasisSetTag);

                j = gomp_PutAtomAtmName(Wstr , AtomName             , atoms);
                j = gomp_PutAtomResName(Wstr , DEFAULT_RESIDUE_NAME , atoms);
                j = gomp_PutAtomSegName(Wstr , DEFAULT_SEGMENT_NAME , atoms);

                j = gomp_PutAtomResNum1(Wstr ,  1 , atoms);
                j = gomp_PutAtomResNum2(Wstr ,  1 , atoms);

                j = gomp_PutAtomXCoord(Wstr , Xc * conv   , atoms);
                j = gomp_PutAtomYCoord(Wstr , Yc * conv  , atoms);
                j = gomp_PutAtomZCoord(Wstr , Zc * conv , atoms);

                j = gomp_PutAtomBValue(Wstr , 0.0 , atoms);
                j = gomp_PutAtomCharge(Wstr , 0.0 , atoms);

                j = gomp_PutAtomBasisSetTag(Wstr , BasisSetTag , atoms);
                j = gomp_PutAtomNucCharge(Wstr , NucCharge  , atoms);

                atoms++;

            }
        }
    }


    gomp_PrintMessage("**********   Done   **********");

end:
    fclose(openmol_in);

    return(Wstr >= 0 ? 0 : 1);

}
/*************************************************************************/
int AtomsInOpenMolFile(FILE *file_p)
/*************************************************************************/
{



    char inputl[OPENMOL_LINE_LEN];   /* input line */
    int OPENMOLatoms;

/* just for sure rewind the file */
    rewind(file_p);

/*
  Start reading file
  Tags recogniced by this routine are:

  : label => next record is a title 
  : gaussian_basis_set_charge_centers => atom records 
*/
    OPENMOLatoms = 0 ; /* start counting number of atoms */

    while(fgets(inputl,OPENMOL_LINE_LEN,file_p) != NULL) {

        if(Tcl_StringMatch(inputl , ": gaussian_basis_set_charge_centers*")) {

            while(fgets(inputl,OPENMOL_LINE_LEN,file_p) != NULL) {
                if(inputl[0] == ':') break;
                OPENMOLatoms++;
            }                /* end of atom loop    */
        }                  /* end of centers loop */
    }                    /* and of file loop    */

/* rewind the file to be ready for the "real" processing */
    rewind(file_p);

    return(OPENMOLatoms);
}

