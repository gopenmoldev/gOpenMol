/*  

Copyright (c) 1995 - 2004 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2003 by:
Eero Häkkinen

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>

#include "measure.h"
#include "molecule.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "memalloc.h"
#include "printmsg.h"
#include "selection.h"

#include "stdafx.h"

/************************************************************************/
float gomp_CalculateQuatfit(const char *text1,const char *text2,
                          const char *text3,const char *text4,
                          const char *text5,const char *text6,
                          const char *text7,const char *text8,
                          int   no_extra)   
/************************************************************************/
{

    static int   i;
    static float temp1,temp2,temp3;
    static float mmass;
    static int *sel_list1;
    static int *sel_list2;
    static int   slong1,slong2,atom_max;
    static float QuatFitValue;
    static int   StruL1;
    static int   StruL2;

    if(gomp_GetNumMolecStructs() < 2) {
        gomp_PrintMessage("?ERROR - no molecules defined to be fitted");
        return((float)-1.0);
    }

    atom_max = gomp_GetTotalNumberOfAtoms();      

    sel_list1 = gomp_AllocateIntVector(atom_max);
    sel_list2 = gomp_AllocateIntVector(atom_max);

    temp1 = (float)0.0;
    temp2 = (float)0.0;
    temp3 = (float)0.0;
    mmass = (float)0.0;

    if(text1[0] == '\0') {

        StruL1 = 0;
        StruL2 = 1;

        for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(StruL1) ; i++) {
            sel_list1[i] = i;
            sel_list2[i] = i;
        }

        slong1 = gomp_GetNumAtomsInMolecStruct(StruL1);
        slong2 = gomp_GetNumAtomsInMolecStruct(StruL2);

        if(slong1 != slong2) {
            gomp_PrintMessage("?ERROR - number of atoms in list 1 and 2 must be equal");
            free(sel_list1);
            free(sel_list2);
            return((float)-1.0);
        }
    }
    else {

        StruL1 = atoi(text1);
        if(StruL1 < 1 || StruL1 > gomp_GetNumMolecStructs()) {
            gomp_PrintERROR("index to structure #1 out of allowed range");
            return((float)-1.0);
        }

        StruL2 = atoi(text5);
        if(StruL2 < 1 || StruL2 > gomp_GetNumMolecStructs()) {
            gomp_PrintERROR("index to structure #2 out of allowed range");
            return((float)-1.0);
        }

        StruL1 -= 1;
        StruL2 -= 1;

/* select from the first structure (cheat the selection list) */
        slong1 = gomp_MakeSelectionList(StruL1,text2,text3,text4,sel_list1);
        slong2 = gomp_MakeSelectionList(StruL2,text6,text7,text8,sel_list2);

    }

    if(slong1 < 1) {
        gomp_PrintMessage("?ERROR - no atoms selected from list 1");
        free(sel_list1);
        free(sel_list2);
        return((float)-1.0);
    }

    if(slong2 < 1) {
        gomp_PrintMessage("?ERROR - no atoms selected from list 2");
        free(sel_list1);
        free(sel_list2);
        return((float)-1.0);
    }

    if(slong1 != slong2) {
        gomp_PrintMessage("?ERROR - number of atoms in list 1 and 2 must be equal");
        free(sel_list1);
        free(sel_list2);
        return((float)-1.0);
    }

    QuatFitValue = gomp_QuatFit(gomp_GetAtomXCoordPointer(StruL1) ,
                              gomp_GetAtomYCoordPointer(StruL1) ,
                              gomp_GetAtomZCoordPointer(StruL1) ,
                              gomp_GetNumAtomsInMolecStruct(StruL1) ,
                              gomp_GetModifiableAtomXCoordPointer(StruL2) , 
                              gomp_GetModifiableAtomYCoordPointer(StruL2) , 
                              gomp_GetModifiableAtomZCoordPointer(StruL2) ,
                              gomp_GetNumAtomsInMolecStruct(StruL2) ,
                              sel_list1 , slong1 ,
                              sel_list2 , slong2 ,
                              no_extra);

    free(sel_list1);
    free(sel_list2);

    return(QuatFitValue);
}       

