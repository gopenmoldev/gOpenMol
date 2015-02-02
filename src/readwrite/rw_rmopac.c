/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Enhancements 2003, 2005 by:
Eero Häkkinen
*/

/*

This is the piece of code from MOPAC 5 which writes the data to disk

IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
* write to disk the following data for graphics calculation, in order:
*
*      number of atoms, orbital, electrons
*      all atomic coordinates
*      orbital counters
*      orbital exponents, s, p, and d, and atomic numbers
*      eigenvectors (m.o.s not re-normalized)
*      inverse-square root of the overlap matrix.
*
write(13)numat,norbs,nelecs,((xyz(i,j),j=1,numat),i=1,3)
write(13)(nlast(i),nfirst(i),i=1,numat)
write(13)(zs(nat(i)),i=1,numat),(zp(nat(i)),i=1,numat),
1         (zd(nat(i)),i=1,numat),(nat(i),i=1,numat)
linear=norbs*norbs
write(13)(c(i),i=1,linear)
write(13)(f(i),i=1,linear)
if(index(keywrd,'mullik').eq.0)return
endif
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "coord_file.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

/***************************************************************************/
int gomp_ReadCoordinatesMOPACgraph(const char *filen, int Append)
    /* mopac graphics output file */
/***************************************************************************/
{

    FILE *mopac_p;
    int record;
    int i,j,k,l,found_el;
    int numatm;  /* number of atoms */
    int norbs;   /*           orbitals */
    int nelecs;  /*           electrons */
    double  **xyz; /* coordinates pointer */
    int *nlast,*nfirst; /* orbital counters */
    double  *zs;   /* orbital exponents s,p,d */
    double  *zp;
    double  *zd;
    int *nat;
    char OutText[BUFF_LEN];
    int Wstr;

    mopac_p = fopen(filen,"r");
    if(mopac_p == NULL) {
        sprintf(OutText,"can't open MOPAC output file '%s'",filen);
        gomp_PrintERROR(OutText);
        return(1);
    }

    fread(&record,sizeof(int),1,mopac_p); /* controll record */
    fread(&numatm,sizeof(int),1,mopac_p);
    fread(&norbs,sizeof(int),1,mopac_p);
    fread(&nelecs,sizeof(int),1,mopac_p);

    sprintf(OutText,"Reading MOPAC graphics output file '%s'",filen);
    gomp_PrintMessage(OutText);
    sprintf(OutText,"numat: %d norbs: %d nelecs: %d",numatm,norbs,nelecs);
    gomp_PrintMessage(OutText);

    if(Append)
        Wstr = gomp_CreateMolecStruct(filen , numatm , APPEND);
    else
        Wstr = gomp_CreateMolecStruct(filen , numatm , NEW);
    if ( Wstr < 0 )
        goto end;

/* get space for the rest */

    nlast  = malloc( numatm * sizeof(int));
    nfirst = malloc( numatm * sizeof(int));

    if(nlast == NULL || nfirst == NULL) {
        gomp_PrintERROR("?ERROR - can't allocate space for arrays 'rmopac'");
        return(1);
    }

    xyz = malloc( numatm * sizeof(const double  *));
       
    if(xyz == NULL ) {
        gomp_PrintERROR("?ERROR - can't allocate space for arrays 'rmopac'");
        return(1);
    }

 
    for(i = 0 ; i < numatm ; i++) {
        xyz[i] = malloc( 3 * sizeof(double));
        if(xyz[i] == NULL ) {
            gomp_PrintERROR("?ERROR - can't allocate space for arrays 'rmopac'");
            return(1);
        }
    }

    zs = malloc(numatm * sizeof(double));
    zp = malloc(numatm * sizeof(double));
    zd = malloc(numatm * sizeof(double));

    nat = malloc(numatm * sizeof(int));

    for(j = 0 ; j < 3 ; j++) {
        for(i = 0 ; i < numatm ; i++) {
            fread(&xyz[i][j],sizeof(double), 1 ,mopac_p);
        }
    }

    fread(&record,sizeof(int),1,mopac_p); /* controll record */

    fread(&record,sizeof(int),1,mopac_p); /* controll record */
    for(i = 0 ; i < numatm ; i++) {
        fread(&nlast[i],sizeof(int), 1 ,mopac_p);
        fread(&nfirst[i],sizeof(int), 1 ,mopac_p);
    }
    fread(&record,sizeof(int),1,mopac_p); /* controll record */

    fread(&record,sizeof(int),1,mopac_p); /* controll record */
    fread(zs,sizeof(double),numatm,mopac_p);
    fread(zp,sizeof(double),numatm,mopac_p);
    fread(zd,sizeof(double),numatm,mopac_p);
    fread(nat,sizeof(int),numatm,mopac_p);
    fread(&record,sizeof(int),1,mopac_p); /* controll record */

/* flush coordinates to the gOpenMol vectors */
    k = 0;
    for(i = 0 ; i < numatm ; i++) {
        l = gomp_PutAtomXCoord(Wstr ,        xyz[i][0] , i); 
        l = gomp_PutAtomYCoord(Wstr ,       xyz[i][1] , i); 
        l = gomp_PutAtomZCoord(Wstr ,      xyz[i][2] , i);
        l = gomp_PutAtomCharge( Wstr , 0.0 , i);
/* segment name    */
        l = gomp_PutAtomSegName(Wstr , "mopc",i);
/* residue name    */
        l = gomp_PutAtomResName(Wstr , "mopc",i);
/* residue number  */
        l = gomp_PutAtomResNum1( Wstr ,  1 , i);
        l = gomp_PutAtomResNum2( Wstr , 1 , i);
/* set atom symbol */
        l = nat[k];
        found_el = 0;
#ifndef ELEM_PARAM
        for(j = 0 ; j < NumAtomSymbols ; j++) {
            if(l == AtomSymbol_p[j]) {
                found_el = 1;
                if(AtomSymbols[4*j+1] == ' ') 
                    sprintf(OutText,"%.1s%.2d",AtomSymbols+4*j,(i+1));
                else 
                    sprintf(OutText,"%.2s%.2d",AtomSymbols+4*j,(i+1));
                break;
            }
        }
#else
        if(l>0 && l<=gomp_GetNumberOfAtomSymbols()) {
            found_el = 1;
            if(gomp_GetAtomSymbol(l)[1] == ' ')
                sprintf(OutText,"%.1s%.2d",gomp_GetAtomSymbol(l),(i+1));
            else 
                sprintf(OutText,"%.2s%.2d",gomp_GetAtomSymbol(l),(i+1));
        }
#endif
        if(found_el == 0) {
            gomp_PrintWARNING("can't find element");
            sprintf(OutText,"Element number: %d is unknown",nat[i]);
            gomp_PrintMessage(OutText);
        }
        else {
            l = gomp_PutAtomAtmName(Wstr , OutText , i);
        }
        k++; /* increment mopac list pointer */
    }

/* free memory */
    free(nlast);
    free(nfirst);
    free(xyz);
    free(zs);
    free(zp);
    free(zd);
    free(nat);
end:
    fclose(mopac_p);
    return(Wstr >= 0 ? 0 : 1);
}
