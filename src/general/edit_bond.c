/*

Copyright (c) 1990 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2001, 2002 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "bond.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "selection.h"

#include "stdafx.h"

/* ........ */

static int * (*GetConnection[2])(int, int) = {
    gomp_GetModifiableAtomConnection,
    gomp_GetModifiableAtomHydrogenBond
};

#if 0
/************************************************************************/
int gomp_EditBond(int   bond_type,
                  int   alt,
                  const char *text1,
                  const char *text2,
                  const char *text3,
                  const char *text4,
                  const char *text5,
                  const char *text6)
/************************************************************************/
{
    static int i,j,k;
    static char chelp[BUFF_LEN];
    static int ihelp;
    static int *AtomConnectionTable1;
    static int *AtomConnectionTable2;
    static int *sel_list1,*sel_list2;
    static int slong1,slong2,atom_max;

    atom_max  = gomp_GetTotalNumberOfAtoms();

    if(!atom_max) {
        gomp_PrintWARNING("no molecular system is defined");
        return(1);
    }
    sel_list1 = gomp_AllocateIntVector(atom_max);
    sel_list2 = gomp_AllocateIntVector(atom_max);

    for(k = 0 ; k < gomp_GetNumMolecStructs() ; k++) {

        slong1 = gomp_MakeSelectionList( k , text1,text2,text3,sel_list1);
        slong2 = gomp_MakeSelectionList( k , text4,text5,text6,sel_list2);

        if((slong1 < 0) || (slong2 < 0))
            continue;

        if(slong1 > 1 || slong2 > 1) {
            sprintf(chelp,
                    "Structure (%d): more than one atom in selection list 1 and/or 2",(k+1));
            gomp_PrintERROR(chelp);
            free(sel_list1);
            free(sel_list2);
            return(1);
        }

        switch(alt) {

        case 1: /* add bond */

            if(slong1 == 1 && slong2 == 1) {

                i = sel_list1[0];
                j = sel_list2[0];

/* check first if there is a bond between atoms i and j */

                if(!gomp_CheckBond(bond_type,k,i,j)) {
                    gomp_PrintWARNING("there is already a bond between specified atoms ");
                    free(sel_list1);
                    free(sel_list2);
                    return(0);
                }

                AtomConnectionTable1 = GetConnection[bond_type](k, i);
                ihelp = AtomConnectionTable1[0] + 1;

                if(ihelp >= gomp_GetMaxAtomConnections()) {
                    gomp_PrintERROR("max connections reached");
                    free(sel_list1);
                    free(sel_list2);
                    return(2);
                }

                AtomConnectionTable1[0]     = ihelp;
                AtomConnectionTable1[ihelp] = j;

                AtomConnectionTable2 = GetConnection[bond_type](k, j);
                ihelp = AtomConnectionTable2[0] + 1;

                if(ihelp >= gomp_GetMaxAtomConnections()) {
                    gomp_PrintERROR("?ERROR - max connections reached");
                    free(sel_list1);
                    free(sel_list2);
                    return(2);
                }

                AtomConnectionTable2[0]     = ihelp;
                AtomConnectionTable2[ihelp] = i;
            }
            else {
                sprintf(chelp,
                        "Structure (%d): no atoms in the selection list 1 and/or 2 ",(k+1));
                gomp_PrintWARNING(chelp);
            }

            break;

        case 2:  /*break bond */

            if(slong1 == 1 && slong2 == 1) {

                i = sel_list1[0];
                j = sel_list2[0];

/* check first if there is a bond between atoms i and j */
                if(gomp_CheckBond(bond_type,k,i,j)) {
                    gomp_PrintMessage("?WARNING - there is no bond between specified atoms ");
                    free(sel_list1);
                    free(sel_list2);
                    return(1);
                }
/* remove bond i - j from the connectivity matrix */
                (void)gomp_BreakBond(bond_type,k,i,j);

            }
            else {
                sprintf(chelp,
                        "Structure (%d): no atoms in the selection list 1 and/or 2 ",(k+1));
                gomp_PrintWARNING(chelp);
            }

            break;

        default:
            gomp_PrintMessage("?ERROR - you should not be here in edit_bond \n");
            return(3);
        }
    }
    free(sel_list1);
    free(sel_list2);

    return(0);
}
#endif
/************************************************************************/
int gomp_EditBondI(int   bond_type,
                 int   alt,
                 int  Which,
                 const char *text1,
                 const char *text2,
                 const char *text3,
                 const char *text4,
                 const char *text5,
                 const char *text6)
/************************************************************************/
{
    static int i,j,k;
    static char chelp[BUFF_LEN];
    static int ihelp;
    static int *AtomConnectionTable1;
    static int *AtomConnectionTable2;
    static int *sel_list1,*sel_list2;
    static int slong1,slong2,atom_max;

    atom_max  = gomp_GetTotalNumberOfAtoms();

    if(!atom_max) {
        gomp_PrintWARNING("no molecular system is defined");
        return(1);
    }
    sel_list1 = gomp_AllocateIntVector(atom_max);
    sel_list2 = gomp_AllocateIntVector(atom_max);

    if(Which < 0 || Which > (gomp_GetNumMolecStructs() - 1)) {
        gomp_PrintERROR("structure index out of range");
        return(1);
    }

    k = Which;

    slong1 = gomp_MakeSelectionList( k , text1,text2,text3,sel_list1);
    slong2 = gomp_MakeSelectionList( k , text4,text5,text6,sel_list2);

    if(slong1 > 1 || slong2 > 1) {
        sprintf(chelp,
                "Structure (%d): more than one atom in selection list 1 and/or 2",(k+1));
        gomp_PrintERROR(chelp);
        free(sel_list1);
        free(sel_list2);
        return(1);
    }

    switch(alt) {

    case 1: /* add bond */

        if(slong1 == 1 && slong2 == 1) {

            i = sel_list1[0];
            j = sel_list2[0];

/* check first if there is a bond between atoms i and j */

            if(!gomp_CheckBond(bond_type,k,i,j)) {
                gomp_PrintMessage("?WARNING - there is already a bond between specified atoms ");
                free(sel_list1);
                free(sel_list2);
                return(0);
            }

            AtomConnectionTable1 = GetConnection[bond_type](Which, i);
            ihelp = AtomConnectionTable1[0] + 1;

            if(ihelp >= gomp_GetMaxAtomConnections()) {
                gomp_PrintMessage("?ERROR - max connections reached");
                free(sel_list1);
                free(sel_list2);
                return(2);
            }

            AtomConnectionTable1[0]     = ihelp;
            AtomConnectionTable1[ihelp] = j;

            AtomConnectionTable2 = GetConnection[bond_type](Which, j);
            ihelp = AtomConnectionTable2[0] + 1;

            if(ihelp >= gomp_GetMaxAtomConnections()) {
                gomp_PrintMessage("?ERROR - max connections reached");
                free(sel_list1);
                free(sel_list2);
                return(2);
            }

            AtomConnectionTable2[0]     = ihelp;
            AtomConnectionTable2[ihelp] = i;

        }
        else {
            sprintf(chelp,
                    "Structure (%d): no atoms in the selection list 1 and/or 2 ",(k+1));
            gomp_PrintWARNING(chelp);
        }

        break;

    case 2:  /*break bond */

        if(slong1 == 1 && slong2 == 1) {

            i = sel_list1[0];
            j = sel_list2[0];

/* check first if there is a bond between atoms i and j */
            if(gomp_CheckBond(bond_type,k,i,j)) {
                gomp_PrintMessage("?WARNING - there is no bond between specified atoms ");
                free(sel_list1);
                free(sel_list2);
                return(1);
            }
/* remove bond i - j from the connectivity matrix */
            (void)gomp_BreakBond(bond_type,k,i,j);
        }
        else {
            sprintf(chelp,
                    "Structure (%d): no atoms in the selection list 1 and/or 2 ",(k+1));
            gomp_PrintWARNING(chelp);
        }

        break;

    default:
        gomp_PrintMessage("?ERROR - you should not be here in edit_bond \n");
        return(3);
    }

    free(sel_list1);
    free(sel_list2);

    return(0);
}


/************************************************************************/
int gomp_CheckBond(int bond_type, int Wstr, int i, int j)           /* check bond i - j */
/************************************************************************/
{
    static int k,l;
    static const int *AtomConnectionTable;

    AtomConnectionTable = GetConnection[bond_type](Wstr, i);
    l = AtomConnectionTable[0]; /* number of connections for atom i */
    for(k = 1 ; k <= l ; k++) 
        if(AtomConnectionTable[k] == j) return(0); /* yes there is a bond */
     
    return(1); /* no bond found */
}
/************************************************************************/
int gomp_BreakBond(int bond_type, int Wstr , int i, int j)          /* break bond i - j */
/************************************************************************/
{
    static int k,l,n;
    static int *AtomConnectionTable;
    static int found;

    found = 0;

    AtomConnectionTable = GetConnection[bond_type](Wstr, i);
    l = AtomConnectionTable[0];

    for(k = 1 ; k <= l ; k++) {
        if(AtomConnectionTable[k] == j) { /* delete it */

            AtomConnectionTable[0]--;

            for(n = k ; n < l ; n++) { /* pop stack */
                AtomConnectionTable[n] = AtomConnectionTable[n+1];
            }
            found = 1;
        }
    }

    AtomConnectionTable = GetConnection[bond_type](Wstr, j);
    l = AtomConnectionTable[0];

    for(k = 1 ; k <= l ; k++) {
        if(AtomConnectionTable[k] == i) { /* delete it */

            AtomConnectionTable[0]--;

            for(n = k ; n < l ; n++) { /* pop stack */
                AtomConnectionTable[n] = AtomConnectionTable[n+1];
            }
            found = 1;
        }
    }

    return(!found); 
}

/************************************************************************/
int gomp_EditBondM(int   bond_type,
                 int   alt,
                 const char *text1,
                 const char *text2,
                 const char *text3,
                 const char *text4,
                 const char *text5,
                 const char *text6)
/************************************************************************/
{

    static int i,j,k;
    static int ii,jj;
    static int ihelp;
    static int *AtomConnectionTable1;
    static int *AtomConnectionTable2;
    static int *sel_list1,*sel_list2;
    static int slong1,slong2,atom_max;

    atom_max  = gomp_GetTotalNumberOfAtoms();

    if(!atom_max) {
        gomp_PrintWARNING("no molecular sustem is defined");
        return(1);
    }

    sel_list1 = gomp_AllocateIntVector(atom_max);
    sel_list2 = gomp_AllocateIntVector(atom_max);

    for(k = 0 ; k < gomp_GetNumMolecStructs() ; k++) {

        slong1 = gomp_MakeSelectionList( k , text1,text2,text3,sel_list1);
        slong2 = gomp_MakeSelectionList( k , text4,text5,text6,sel_list2);

        if((slong1 < 0) || (slong2 < 0)) continue;

        switch(alt) {

        case 1: /* add bond */

            if(slong1 > 0 && slong2 > 0) {

                for(ii = 0 ; ii < slong1 ; ii++) {
                    i = sel_list1[ii];
                    for(jj = 0 ; jj < slong2 ; jj++) {
                        j = sel_list2[jj];

                        if(j == i) continue;

/* check first if there is a bond between atoms i and j */

                        if(!gomp_CheckBond(bond_type, k , i, j)) {
                            gomp_PrintMessage("?WARNING - there is already a bond between specified atoms ");
                            return(0);
                        }

                        AtomConnectionTable1 = GetConnection[bond_type](k, i);
                        ihelp = AtomConnectionTable1[0] + 1;

                        if(ihelp >= gomp_GetMaxAtomConnections()) {
                            gomp_PrintMessage("?ERROR - max connections reached");
                            free(sel_list1);
                            free(sel_list2);
                            return(2);
                        }

                        AtomConnectionTable1[0]     = ihelp;
                        AtomConnectionTable1[ihelp] = j;

                        AtomConnectionTable2 = GetConnection[bond_type](k, j);
                        ihelp = AtomConnectionTable2[0] + 1;

                        if(ihelp >= gomp_GetMaxAtomConnections()) {
                            gomp_PrintMessage("?ERROR - max connections reached");
                            free(sel_list1);
                            free(sel_list2);
                            return(2);
                        }

                        AtomConnectionTable2[0]     = ihelp;
                        AtomConnectionTable2[ihelp] = i;

                    }
                }
            }
            else
                gomp_PrintMessage("?ERROR - no atoms in the selection list 1 or 2 ");

            break;

        case 2:  /*break bond */

            if(slong1 > 0 && slong2 > 0) {

                for(ii = 0 ; ii < slong1 ; ii++) {
                    i = sel_list1[ii];
                    for(jj = 0 ; jj < slong2 ; jj++) {
                        j = sel_list2[jj];

                        if(j == i) continue;

/* remove bond i - j from the connectivity matrix */
                        (void)gomp_BreakBond(bond_type, k , i, j);

                    }
                }
            }
            else
                gomp_PrintMessage("?ERROR - no atoms in the selection list 1 or 2");

            break;

        default:
            gomp_PrintMessage("?ERROR - you should not be here in edit_bond \n");
            return(3);
        }
    }

    free(sel_list1);
    free(sel_list2);

    return(0);
}

/************************************************************************/
int gomp_BreakBondM(int    bond_type,
                  const char *text1,
                  const char *text2,
                  const char *text3)
/************************************************************************/
{

    static int i,j,k,si,swop;
    static int ii,kk;
    static char chelp[BUFF_LEN];
    static int ihelp;
    static int *AtomConnectionTable1;
    static int *AtomConnectionTable2;
    static int *sel_list1;
    static int slong1,atom_max;
    static int  kb;

    atom_max  = gomp_GetTotalNumberOfAtoms();

    if(!atom_max) {
        gomp_PrintWARNING("no molecular system is defined");
        return(1);
    }

    sel_list1 = gomp_AllocateIntVector(atom_max);

    for(kb = 0 ; kb < gomp_GetNumMolecStructs() ; kb++) {

        slong1 = gomp_MakeSelectionList(kb , text1,text2,text3,sel_list1);

        if(slong1 < 0) continue;

        if(slong1 > 1) {

            for(ii = 0 ; ii < slong1 ; ii++) {
                i = sel_list1[ii];

                AtomConnectionTable1    = GetConnection[bond_type](kb, i);
                si = AtomConnectionTable1[0];
                for(j = 1 ; j <= si ; j++) {
                    ihelp = AtomConnectionTable1[j];
                    swop = 0;
                    for(k = 0 ; k < slong1 ; k++) {
                        kk = sel_list1[k];
                        if(kk == ihelp) {
                            swop = 1;
                            break;
                        }
                    }
                    if(!swop) {
                        AtomConnectionTable2 = GetConnection[bond_type](kb, ihelp);
                        AtomConnectionTable2[0] = 0;
                    }
                }
                AtomConnectionTable1[0] = 0; 
            }
        }
        else {
            sprintf(chelp,
                    "Structure (%d): no atoms in the selection list 1 or 2 ",(kb+1));
            gomp_PrintWARNING(chelp);
        }
    }
    free(sel_list1);

    return(0);
}

