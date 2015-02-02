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

#include "cell.h"
#include "ldp.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "selection.h"

#include "stdafx.h"

static int  BuildLDPList(void);
static int  DeleteLDPList(void);

static struct LDP_Struct {
    char Atom1[MAX_ATM_NAME_LEN];   /* first atom in list      */
    char Resi1[MAX_RES_NAME_LEN];
    char Segm1[MAX_SEG_NAME_LEN];
    char Atom2[MAX_ATM_NAME_LEN];   /* second atom in list     */
    char Resi2[MAX_RES_NAME_LEN];
    char Segm2[MAX_SEG_NAME_LEN];
    int  NumAtoms1;              /* number of atoms in list 1 */
    int *AtomIndexList1;        /* Atom index list of entries 1 */
    int  NumAtoms2;              /* number of atoms in list 2 */
    int *AtomIndexList2;        /* Atom index list of entries 2 */
    float *LDParray;             /* floating pointer to the data array */
}  LDP_info;

static struct {
    float min1;
    float max1;
    float min2;
    float max2;
    float min3;
    float max3;
    float min4;
    float max4;
    float rest;
} LDPBounds = { 3.5 , 4.0 , 4.0 , 6.0 , 6.0 , 8.0 , 8.0 ,
                               10.0 , 10.0};

#ifdef ENABLE_GRAPHICS
/***********************************************************************/
int gomp_BuildLDParray()
/***********************************************************************/
{
    char TempS[BUFF_LEN];
    static float min1,min2,min3,min4;
    static float max1,max2,max3,max4;
    static float rest;

    min1  = LDPBounds.min1;
    max1 = LDPBounds.max1;
    min2  = LDPBounds.min2;
    max2 = LDPBounds.max2;
    min3  = LDPBounds.min3;
    max3 = LDPBounds.max3;
    min4  = LDPBounds.min4;
    max4 = LDPBounds.max4;
    rest  = LDPBounds.rest;

#if defined(OBSOLATE)
    gomp_PrintMessage("The displayed intervals are:");
    sprintf(TempS,"          dist <  %4.2f ",min1);
    gomp_PrintMessage(TempS);
    sprintf(TempS,"%4.2f   < dist <  %4.2f ",min1,max1);
    gomp_PrintMessage(TempS);
    sprintf(TempS,"%4.2f   < dist <  %4.2f ",min2,max2);
    gomp_PrintMessage(TempS);
    sprintf(TempS,"%4.2f   < dist <  %4.2f ",min3,max3);
    gomp_PrintMessage(TempS);
    sprintf(TempS,"%4.2f   < dist <  %4.2f ",min4,max4);
    gomp_PrintMessage(TempS);
    sprintf(TempS,"%4.2f   < dist          ",rest);
    gomp_PrintMessage(TempS);
    gomp_PrintMessage(" ");
#endif

    if(BuildLDPList()) { /* prepare the lists */
        gomp_PrintMessage("?ERROR - can't build the LDP lists");
        return(1);
    }

    sprintf(TempS,"(1) Number of '%s' atoms is: %d",LDP_info.Atom1,
            LDP_info.NumAtoms1);
    gomp_PrintMessage(TempS);

    sprintf(TempS,"(2) Number of '%s' atoms is: %d",LDP_info.Atom2,
            LDP_info.NumAtoms2);
    gomp_PrintMessage(TempS);


    if(LDP_info.NumAtoms1 < 1) {
        gomp_PrintMessage("?WARNING - no atoms in list (1) to be displayed ");
        (void)DeleteLDPList();
        return(1);
    }

    if(LDP_info.NumAtoms2 < 1) {
        gomp_PrintMessage("?WARNING - no atoms in list (2) to be displayed ");
        (void)DeleteLDPList();
        return(1);
    }

    (void)gomp_CalculateLDParraydraw_ldp(
        LDP_info.NumAtoms1 , 
        LDP_info.NumAtoms2 ,
        min1,min2,min3,min4,max1,max2,max3,max4,rest);

    return(0);

}
#endif /* ENABLE_GRAPHICS */
#if 0
/***********************************************************************/
int gomp_CalculateLDParray(int num1,int num2,
                         float min1,float min2,float min3,float min4,
                         float max1,float max2,float max3,float max4,
                         float rest)
/***********************************************************************/
{
    static float dist;
    static int i,ii,j,jj;
    static float tmp1,tmp2,tmp3;
    static const float *XCoord;
    static const float *YCoord;
    static const float *ZCoord;
    static float xboxl,yboxl,zboxl;
    static float boxl1,boxl2,boxl3;

/* Calculate the distance matrix */ 

    XCoord = gomp_GetAtomXCoordPointer(0);
    YCoord = gomp_GetAtomYCoordPointer(0);
    ZCoord = gomp_GetAtomZCoordPointer(0);

/* take the periodic boundaries into account */

    boxl1 = gomp_GetCellA();
    boxl2 = gomp_GetCellB();
    boxl3 = gomp_GetCellC();
    xboxl = 1. / boxl1;
    yboxl = 1. / boxl2;
    zboxl = 1. / boxl3;


    for(i=0 ; i < num1 ; i++) {
        for(j=i+1 ; j < num2 ; j++)   {

            jj = LDP_info.AtomIndexList2[j];
            ii = LDP_info.AtomIndexList1[i];

            tmp1 = (XCoord[ii]-XCoord[jj]);
            tmp2 = (YCoord[ii]-YCoord[jj]);
            tmp3 = (ZCoord[ii]-ZCoord[jj]);

            tmp1 = tmp1 - boxl1 * rint (xboxl * tmp1);
            tmp2 = tmp2 - boxl2 * rint (yboxl * tmp2);
            tmp3 = tmp3 - boxl3 * rint (zboxl * tmp3);

            dist=sqrt( tmp1 * tmp1 + 
                       tmp2 * tmp2 + 
                       tmp3 * tmp3);

            LDP_info.LDParray[i + j * (j - 1) / 2] = dist;
        }
    }

    return(0);
}
#endif
/***********************************************************************/
int BuildLDPList() 
/***********************************************************************/
{

    register int Iloop;

    int *sel_list;

    if(DeleteLDPList()) {
        gomp_PrintMessage("?ERROR - can't build LDP structure");
        return(1);
    }


    sel_list           = gomp_AllocateIntVector(gomp_GetNumAtomsInMolecStruct(0));
    LDP_info.NumAtoms1 = gomp_MakeSelectionList(0 ,
                                              LDP_info.Segm1 , 
                                              LDP_info.Resi1 ,
                                              LDP_info.Atom1 , sel_list);

    if(!LDP_info.NumAtoms1) {
        gomp_PrintMessage("?ERROR - no atoms in selection list (1)");
        free(sel_list);
        return(1);
    }

    LDP_info.AtomIndexList1 = gomp_AllocateIntVector(LDP_info.NumAtoms1);
    if(LDP_info.AtomIndexList1 == NULL) {
        gomp_PrintMessage("?ERROR - can't get space in building the LDP strcuture");
        return(1);
    }

    for(Iloop = 0 ; Iloop < LDP_info.NumAtoms1 ; Iloop++)
        LDP_info.AtomIndexList1[Iloop] = sel_list[Iloop];

    free(sel_list);

    sel_list           = gomp_AllocateIntVector(gomp_GetNumAtomsInMolecStruct(0));
    LDP_info.NumAtoms2 = gomp_MakeSelectionList(0 ,
                                              LDP_info.Segm2 ,
                                              LDP_info.Resi2 ,
                                              LDP_info.Atom2 , sel_list);

    if(!LDP_info.NumAtoms2) {
        gomp_PrintMessage("?ERROR - no atoms in selection list (2)");
        free(sel_list);
        return(1);
    }

    LDP_info.AtomIndexList2 = gomp_AllocateIntVector(LDP_info.NumAtoms2);
    if(LDP_info.AtomIndexList2 == NULL) {
        gomp_PrintMessage("?ERROR - can't get space in building the LDP strcuture");
        return(1);
    }

    for(Iloop = 0 ; Iloop < LDP_info.NumAtoms2 ; Iloop++)
        LDP_info.AtomIndexList2[Iloop] = sel_list[Iloop];

    free(sel_list);

    LDP_info.LDParray = gomp_AllocateFloatVector(
        LDP_info.NumAtoms1 * (LDP_info.NumAtoms1 - 1) / 2);

    if(LDP_info.LDParray == NULL) {
        gomp_PrintMessage("?ERROR - can't get space in building the LDP strcuture");
        return(1);
    }

    return(0);
}
/***********************************************************************/
int DeleteLDPList()
/***********************************************************************/
{
    if(LDP_info.NumAtoms1) free(LDP_info.AtomIndexList1);
    LDP_info.NumAtoms1 = 0;
    if(LDP_info.NumAtoms2) free(LDP_info.AtomIndexList2);
    LDP_info.NumAtoms2 = 0;
    if(LDP_info.NumAtoms1 && LDP_info.NumAtoms2) free(LDP_info.LDParray);
    return(0);
}

/***********************************************************************/
int gomp_PushAtomToLDP(int Watom , const char *segm , const char *resi , const char *atom)
/***********************************************************************/
{
    switch(Watom) {

    case 1: /* atom 1 */
        if(segm[0] == '\0') {
            strncpy(LDP_info.Atom1,"CA",MAX_ATM_NAME_LEN);
            strncpy(LDP_info.Resi1,"*",MAX_RES_NAME_LEN);
            strncpy(LDP_info.Segm1,"*",MAX_SEG_NAME_LEN);
        }
        else {
            strncpy(LDP_info.Atom1,atom,MAX_ATM_NAME_LEN);
            strncpy(LDP_info.Resi1,resi,MAX_RES_NAME_LEN);
            strncpy(LDP_info.Segm1,segm,MAX_SEG_NAME_LEN);
        }
        break;
    case 2: /* atom 2 */
        if(segm[0] == '\0') {
            strncpy(LDP_info.Atom2,"CA",MAX_ATM_NAME_LEN);
            strncpy(LDP_info.Resi2,"*",MAX_RES_NAME_LEN);
            strncpy(LDP_info.Segm2,"*",MAX_SEG_NAME_LEN);
        }
        else {
            strncpy(LDP_info.Atom2,atom,MAX_ATM_NAME_LEN);
            strncpy(LDP_info.Resi2,resi,MAX_RES_NAME_LEN);
            strncpy(LDP_info.Segm2,segm,MAX_SEG_NAME_LEN);
        }
    }

    return(0);
}

#if 0
/***********************************************************************/
int gomp_ChangeLDPBounds( const char *text1 , const char *text2 , const char *text3 ,
                        const char *text4 , const char *text5 , const char *text6 ,
                        const char *text7 , const char *text8 , const char *text9)
/***********************************************************************/
{

    LDPBounds.min1 = atof(text1);
    LDPBounds.max1 = atof(text2);
    LDPBounds.min2 = atof(text3);
    LDPBounds.max2 = atof(text4);
    LDPBounds.min3 = atof(text5);
    LDPBounds.max3 = atof(text6);
    LDPBounds.min4 = atof(text7);
    LDPBounds.max4 = atof(text8);
    LDPBounds.rest = atof(text9);
      
    return(0);
}

/***********************************************************************/
int gomp_ResetLDPBounds()
/***********************************************************************/
{

    LDPBounds.min1 = 3.5;
    LDPBounds.max1 = 4.0;
    LDPBounds.min2 = 4.0;
    LDPBounds.max2 = 6.0;
    LDPBounds.min3 = 6.0;
    LDPBounds.max3 = 8.0;
    LDPBounds.min4 = 8.0;
    LDPBounds.max4 = 10.0;
    LDPBounds.rest = 10.0;
      
    return(0);
}
#endif
/***********************************************************************/
int gomp_GetNumLdpAtoms1()
/***********************************************************************/
{
    return(LDP_info.NumAtoms1);
}
/***********************************************************************/
int gomp_GetNumLdpAtoms2()
/***********************************************************************/
{
    return(LDP_info.NumAtoms2);
}

/***********************************************************************/
const int *gomp_GetLdpAtomList1()
/***********************************************************************/
{
    return(LDP_info.AtomIndexList1);
}
/***********************************************************************/
const int *gomp_GetLdpAtomList2()
/***********************************************************************/
{
    return(LDP_info.AtomIndexList2);
}
#if 0
/***********************************************************************/
int  gomp_ResetLDPstructure()
/***********************************************************************/
{

/*
  if(LDP_info.NumAtoms1 && LDP_info.NumAtoms2) {
  free(LDParray);
  }
*/

    if(LDP_info.NumAtoms1) {
        free(LDP_info.AtomIndexList1);
        LDP_info.NumAtoms1 = 0;
    }

    if(LDP_info.NumAtoms2) {
        free(LDP_info.AtomIndexList2);
        LDP_info.NumAtoms2 = 0;
    }

    return(0);
}
#endif
