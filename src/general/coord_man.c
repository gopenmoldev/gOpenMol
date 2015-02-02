/*

Copyright (c) 1994 - 2005 by:
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
#include <string.h>
#include "gommath.h"
#include <ctype.h>

#include "coord_man.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "memalloc.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "selection.h"

#include "stdafx.h"

#define Rabs(a)    ( ( a ) > 0.0 ? (a) : -(a))

/* structure to contain a selection */
static struct {
    gom_Plotter *CallbackHandle;
    int  Structure;       /* structure number                */
    int  slong;           /* length of the list              */
    int *sel_list;       /* index to atoms                  */
    float SelCenter[3];   /* geometrical center of selection */
} ManSelection = { NULL , 0 , 0 , NULL, {0.0 , 0.0 , 0.0} };

#if 0
/************************************************************************/
int gomp_TranslateCoordinates(const char *CXtra,
                              const char *CYtra,
                              const char *CZtra,
                              const char *Segment,
                              const char *Residue,
                              const char *Atom)
/************************************************************************/
{
    int *sel_list; /* selection list */
    int  ent_list; /* entries in the selection list */
    int  i,j;
    int atom_list;
    float Xtemp,Ytemp,Ztemp; 
    float *XCoordP;
    float *YCoordP;
    float *ZCoordP;

    atom_list = gomp_GetNumAtomsInMolecStruct(0);

    if(atom_list < 1) {
        gomp_PrintMessage("?ERROR - no atoms defined");
        return(1);
    }

    XCoordP   = gomp_GetModifiableAtomXCoordPointer(0);
    YCoordP  = gomp_GetModifiableAtomYCoordPointer(0);
    ZCoordP = gomp_GetModifiableAtomZCoordPointer(0);

    sel_list = gomp_AllocateIntVector(atom_list);
    ent_list = gomp_MakeSelectionList(0,Segment,Residue,Atom,sel_list);

    Xtemp = atof(CXtra);
    Ytemp = atof(CYtra);
    Ztemp = atof(CZtra);     

    for(i = 0 ; i < ent_list ; i++) {
        j = sel_list[i];
        XCoordP[j] += Xtemp;

        YCoordP[j] += Ytemp;

        ZCoordP[j] += Ztemp;
    }

    free(sel_list);
    return(0);
}
/*

Rotate always around the coordinate center of the selected atoms
     
Changed 1993-09-03  LUL
*/

/************************************************************************/
int gomp_RotateCoordinates1(const char *Caxis,
                            const char *Crot,
                            const char *Segment,
                            const char *Residue,
                            const char *Atom)
/************************************************************************/
{
    static int *sel_list; /* selection list */
    static int  ent_list; /* entries in the selection list */
    static int  i,j,k;
    static int AxisRot;
    static int atom_list;
    static float Piper180 = (float)(M_PI/180.0);
    static float Xtemp,Ytemp,Ztemp; 
    static float *XCtemp,*YCtemp,*ZCtemp;
    static float cosa,sina;
    static float XCcenter;
    static float YCcenter;
    static float ZCcenter;
    static float *XCoordP;
    static float *YCoordP;
    static float *ZCoordP;

    atom_list = gomp_GetNumAtomsInMolecStruct(0);

    if(atom_list < 1) {
        gomp_PrintMessage("?ERROR - no atoms defined");
        return(1);
    }

    XCoordP   = gomp_GetModifiableAtomXCoordPointer(0);
    YCoordP  = gomp_GetModifiableAtomYCoordPointer(0);
    ZCoordP = gomp_GetModifiableAtomZCoordPointer(0);

    sel_list = gomp_AllocateIntVector(atom_list);
    ent_list = gomp_MakeSelectionList(0,Segment,Residue,Atom,sel_list);

    Xtemp = Piper180 * atof(Crot);

#ifdef sgi
    cosa = fcos(Xtemp);
    sina = fsin(Xtemp);
#else
    cosa = cos((double)Xtemp);
    sina = sin((double)Xtemp);
#endif

    XCcenter = YCcenter = ZCcenter = 0.0;

    for(i = 0 ; i < ent_list ; i++) {
        j = sel_list[i];
        XCcenter += XCoordP[j];
        YCcenter += YCoordP[j];
        ZCcenter += ZCoordP[j];
    }

    XCcenter /= (float)ent_list;
    YCcenter /= (float)ent_list;
    ZCcenter /= (float)ent_list;

    AxisRot = 0;
    if(Caxis[0] == 'x' || Caxis[0] == 'X') AxisRot = 1;
    if(Caxis[0] == 'y' || Caxis[0] == 'Y') AxisRot = 2;
    if(Caxis[0] == 'z' || Caxis[0] == 'Z') AxisRot = 3;

    if(AxisRot < 1 || AxisRot > 3) {
        gomp_PrintMessage("?ERROR - unknown axis defined");
        free(sel_list);
        return(1);
    }

    i = atom_list;
    XCtemp   = gomp_AllocateFloatVector(i);
    YCtemp  = gomp_AllocateFloatVector(i);
    ZCtemp = gomp_AllocateFloatVector(i);

    for(i = 0 ; i < ent_list ; i++) {

        k = sel_list[i];

        switch(AxisRot) {

        case 1:

            Ytemp = YCoordP[k] - YCcenter;
            Ztemp = ZCoordP[k] - ZCcenter;

            YCtemp[i] =  cosa * Ytemp - sina * Ztemp;
            ZCtemp[i] =  sina * Ytemp + cosa * Ztemp;
            break;

        case 2:

            Xtemp = XCoordP[k] - XCcenter;
            Ztemp = ZCoordP[k] - ZCcenter;

            XCtemp[i] =  cosa * Xtemp + sina * Ztemp;
            ZCtemp[i] = -sina * Xtemp + cosa * Ztemp;
            break;

        case 3:

            Xtemp = XCoordP[k] - XCcenter;
            Ytemp = YCoordP[k] - YCcenter;

            XCtemp[i] =  cosa * Xtemp - sina * Ytemp;
            YCtemp[i] =  sina * Xtemp + cosa * Ytemp;
            break;

        }
        j++;
    }

    for(i = 0 ; i < ent_list ; i++) {
        k = sel_list[i];

        switch(AxisRot) {

        case 1:
            YCoordP[k] = YCtemp[i] + YCcenter;
            ZCoordP[k] = ZCtemp[i] + ZCcenter;
            break;
        case 2:
            XCoordP[k] = XCtemp[i] + XCcenter;
            ZCoordP[k] = ZCtemp[i] + ZCcenter;
            break;
        case 3:
            XCoordP[k] = XCtemp[i] + XCcenter;
            YCoordP[k] = YCtemp[i] + YCcenter;
            break;
        }
    }


    free(XCtemp);
    free(YCtemp);
    free(ZCtemp);
    free(sel_list);
    return(0);
}

/************************************************************************/
int gomp_RotateCoordinates3(const char *CXrot,
                            const char *CYrot,
                            const char *CZrot,
                            const char *Segment,
                            const char *Residue,
                            const char *Atom)
/************************************************************************/
{
    static int *sel_list; /* selection list */
    static int  ent_list; /* entries in the selection list */
    static int  i,j;
    static int  atom_list;
    static float Piper180 = (float)(M_PI/180.0);
    static float Xtemp,Ytemp,Ztemp; 
    static float *XCtemp,*YCtemp,*ZCtemp;
    static float Fcos,Fsin;
    static float XCcenter;
    static float YCcenter;
    static float ZCcenter;
    static float Fentries;
    static float Temp1;
    static float Temp2;
    static float *XCoordP;
    static float *YCoordP;
    static float *ZCoordP;


    atom_list = gomp_GetNumAtomsInMolecStruct(0);

    if(atom_list < 1) {
        gomp_PrintMessage("?ERROR - no atoms defined");
        return(1);
    }

    XCoordP   = gomp_GetModifiableAtomXCoordPointer(0);
    YCoordP  = gomp_GetModifiableAtomYCoordPointer(0);
    ZCoordP = gomp_GetModifiableAtomZCoordPointer(0);

    sel_list = gomp_AllocateIntVector(atom_list);
    ent_list = gomp_MakeSelectionList(0,Segment,Residue,Atom,sel_list);

    Xtemp = Piper180 * atof(CXrot);
    Ytemp = Piper180 * atof(CYrot);
    Ztemp = Piper180 * atof(CZrot);     

    XCtemp   = gomp_AllocateFloatVector(ent_list);
    YCtemp  = gomp_AllocateFloatVector(ent_list);
    ZCtemp = gomp_AllocateFloatVector(ent_list);

    XCcenter = YCcenter = ZCcenter = 0.0;

    for(i = 0 ; i < ent_list ; i++) {
        j = sel_list[i];
        XCcenter += XCoordP[j];
        YCcenter += YCoordP[j];
        ZCcenter += ZCoordP[j];
    }

    Fentries  = 1./((float)ent_list);
    XCcenter *= Fentries;
    YCcenter *= Fentries;
    ZCcenter *= Fentries;

/* rotate around x - axis */

    if(Rabs(Xtemp) > 0.001) {
        Fcos = cos(Xtemp);
        Fsin = sin(Xtemp);

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];

            Temp1 = YCoordP[j] - YCcenter;
            Temp2 = ZCoordP[j] - ZCcenter;

            YCtemp[i] =  Fcos * Temp1 - Fsin * Temp2;
            ZCtemp[i] =  Fsin * Temp1 + Fcos * Temp2;
        }

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];
            YCoordP[j] = YCtemp[i]+YCcenter;
            ZCoordP[j] = ZCtemp[i]+ZCcenter;
        }
    }

/* rotate around y - axis */

    if(Rabs(Ytemp) > 0.001) {
        Fcos = cos(Ytemp);
        Fsin = sin(Ytemp);

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];

            Temp1 = XCoordP[j] - XCcenter;
            Temp2 = ZCoordP[j] - ZCcenter;

            XCtemp[i] =  Fcos * Temp1 + Fsin * Temp2;
            ZCtemp[i] = -Fsin * Temp1 + Fcos * Temp2;
        }

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];
            XCoordP[j] = XCtemp[i]+XCcenter;
            ZCoordP[j] = ZCtemp[i]+ZCcenter;
        }
    }

/* rotate around z - axis */

    if(Rabs(Ztemp) > 0.001) {
        Fcos = cos(Ztemp);
        Fsin = sin(Ztemp);

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];

            Temp1 = XCoordP[j] - XCcenter;
            Temp2 = YCoordP[j] - YCcenter;

            XCtemp[i] =  Fcos * Temp1 - Fsin * Temp2;
            YCtemp[i] =  Fsin * Temp1 + Fcos * Temp2;
        }

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];
            XCoordP[j] = XCtemp[i] + XCcenter;
            YCoordP[j] = YCtemp[i] + YCcenter;
        }
    }


    free(sel_list);
    free(XCtemp);
    free(YCtemp);
    free(ZCtemp);
    return(0);
}
#endif
/************************************************************************/
int gomp_TranslateCoordinatesX(float CXtra, float CYtra, float CZtra)
/************************************************************************/
{
    static const int *sel_list; /* selection list */
    static int  ent_list; /* entries in the selection list */
    static int  i,j;
    static float Xtemp,Ytemp,Ztemp; 
    static float *XCoordP;
    static float *YCoordP;
    static float *ZCoordP;
    static int    atom_list;

    atom_list = gomp_GetTotalNumberOfAtoms();

    if(atom_list < 1) {
        gomp_PrintMessage("?ERROR - no atoms defined");
        return(1);
    }

    if(ManSelection.slong < 1) {
        gomp_PrintWARNING("no atom list selected");
        return(1);
    }

    XCoordP   = gomp_GetModifiableAtomXCoordPointer(ManSelection.Structure);
    YCoordP  = gomp_GetModifiableAtomYCoordPointer(ManSelection.Structure);
    ZCoordP = gomp_GetModifiableAtomZCoordPointer(ManSelection.Structure);

    sel_list = ManSelection.sel_list;
    ent_list = ManSelection.slong;

    Xtemp = CXtra;
    Ytemp = CYtra;
    Ztemp = CZtra;     

    for(i = 0 ; i < ent_list ; i++) {
        j = sel_list[i];

        XCoordP[j]  +=  Xtemp;

        YCoordP[j]  +=  Ytemp;

        ZCoordP[j]  +=  Ztemp;
    }

    (void)gomp_MoveSelectionCenter(Xtemp , Ytemp , Ztemp);

    return(0);
}
/*

Rotate always around the coordinate center of the selected atoms
     
Changed 1993-09-03  LUL
*/

/************************************************************************/
int gomp_RotateCoordinates1X(float Crot , char Caxis)
/************************************************************************/
{
    static const int *sel_list; /* selection list */
    static int  ent_list; /* entries in the selection list */
    static int  i,j,k;
    static int atom_list;
    static int AxisRot;
    static const float *XYZCcenter;
    static float Piper180 = (float)(M_PI/180.0);
    static float Xtemp,Ytemp,Ztemp; 
    static float *XCtemp,*YCtemp,*ZCtemp;
    static float cosa,sina;
    static float XCcenter;
    static float YCcenter;
    static float ZCcenter;
    static float *XCoordP;
    static float *YCoordP;
    static float *ZCoordP;


    atom_list = gomp_GetTotalNumberOfAtoms();

    if(atom_list < 1) {
        gomp_PrintMessage("?ERROR - no atoms defined");
        return(1);
    }

    if(ManSelection.slong < 1) {
        gomp_PrintWARNING("no atom list selected");
        return(1);
    }

    AxisRot = 0;
    if(Caxis == 'x' || Caxis == 'X') AxisRot = 1;
    if(Caxis == 'y' || Caxis == 'Y') AxisRot = 2;
    if(Caxis == 'z' || Caxis == 'Z') AxisRot = 3;

    if(AxisRot < 1 || AxisRot > 3) {
        gomp_PrintMessage("?ERROR - unknown axis defined");
        return(1);
    }

    XCoordP   = gomp_GetModifiableAtomXCoordPointer(ManSelection.Structure);
    YCoordP  = gomp_GetModifiableAtomYCoordPointer(ManSelection.Structure);
    ZCoordP = gomp_GetModifiableAtomZCoordPointer(ManSelection.Structure);

    sel_list =  ManSelection.sel_list;
    ent_list =  ManSelection.slong;

    Xtemp = Piper180 * Crot;

#ifdef sgi
    cosa = fcos(Xtemp);
    sina = fsin(Xtemp);
#else
    cosa = cos((double)Xtemp);
    sina = sin((double)Xtemp);
#endif

    XYZCcenter = gomp_GetSelectionCenter();
    XCcenter = XYZCcenter[0];
    YCcenter = XYZCcenter[1];
    ZCcenter = XYZCcenter[2];

    i = atom_list;
    XCtemp   = gomp_AllocateFloatVector(i);
    YCtemp  = gomp_AllocateFloatVector(i);
    ZCtemp = gomp_AllocateFloatVector(i);

    for(i = 0 ; i < ent_list ; i++) {

        k = sel_list[i];

        switch(AxisRot) {

        case 1:

            Ytemp = YCoordP[k] - YCcenter;
            Ztemp = ZCoordP[k] - ZCcenter;
            YCtemp[i] =  cosa * Ytemp - sina * Ztemp;
            ZCtemp[i] =  sina * Ytemp + cosa * Ztemp;
            break;

        case 2:

            Xtemp = XCoordP[k] - XCcenter;
            Ztemp = ZCoordP[k] - ZCcenter;
            XCtemp[i] =  cosa * Xtemp + sina * Ztemp;
            ZCtemp[i] = -sina * Xtemp + cosa * Ztemp;
            break;

        case 3:

            Xtemp = XCoordP[k] - XCcenter;
            Ytemp = YCoordP[k] - YCcenter;
            XCtemp[i] =  cosa * Xtemp - sina * Ytemp;
            YCtemp[i] =  sina * Xtemp + cosa * Ytemp;
            break;
        }
        j++;
    }

    for(i = 0 ; i < ent_list ; i++) {
        k = sel_list[i];

        switch(AxisRot) {

        case 1:
            YCoordP[k] = YCtemp[i] + YCcenter;
            ZCoordP[k] = ZCtemp[i] + ZCcenter;
            break;
        case 2:
            XCoordP[k] = XCtemp[i] + XCcenter;
            ZCoordP[k] = ZCtemp[i] + ZCcenter;
            break;
        case 3:
            XCoordP[k] = XCtemp[i] + XCcenter;
            YCoordP[k] = YCtemp[i] + YCcenter;
            break;
        }
    }


    free(XCtemp);
    free(YCtemp);
    free(ZCtemp);
    return(0);
}
#if 0
/************************************************************************/
int gomp_RotateCoordinates3X(float CXrot,float CYrot,float CZrot)
/************************************************************************/
{
    static const int *sel_list; /* selection list */
    static int  ent_list; /* entries in the selection list */
    static int  i,j;
    static int  atom_list;
    static float Piper180 = (float)(M_PI/180.0);
    static float Xtemp,Ytemp,Ztemp; 
    static float *XCtemp,*YCtemp,*ZCtemp;
    static float Fcos,Fsin;
    static float XCcenter;
    static float YCcenter;
    static float ZCcenter;
    static float Fentries;
    static float Temp1;
    static float Temp2;
    static float *XCoordP;
    static float *YCoordP;
    static float *ZCoordP;


    atom_list = gomp_GetTotalNumberOfAtoms();

    if(atom_list < 1) {
        gomp_PrintMessage("?ERROR - no atoms defined");
        return(1);
    }

    if(ManSelection.slong < 1) {
        gomp_PrintWARNING("no atom list selected");
        return(1);
    }

    XCoordP   = gomp_GetModifiableAtomXCoordPointer(ManSelection.Structure);
    YCoordP  = gomp_GetModifiableAtomYCoordPointer(ManSelection.Structure);
    ZCoordP = gomp_GetModifiableAtomZCoordPointer(ManSelection.Structure);

    ent_list = ManSelection.slong;
    sel_list = ManSelection.sel_list;

    Xtemp = Piper180 * CXrot;
    Ytemp = Piper180 * CYrot;
    Ztemp = Piper180 * CZrot;     

    XCtemp = gomp_AllocateFloatVector(ent_list);
    YCtemp = gomp_AllocateFloatVector(ent_list);
    ZCtemp = gomp_AllocateFloatVector(ent_list);

    XCcenter = YCcenter = ZCcenter = 0.0;

    for(i = 0 ; i < ent_list ; i++) {
        j = sel_list[i];
        XCcenter += XCoordP[j];
        YCcenter += YCoordP[j];
        ZCcenter += ZCoordP[j];
    }

    Fentries  = 1.0 / ((float)ent_list);
    XCcenter *= Fentries;
    YCcenter *= Fentries;
    ZCcenter *= Fentries;

/* rotate around x - axis */

    if(Rabs(Xtemp) > 0.001) {
        Fcos = cos(Xtemp);
        Fsin = sin(Xtemp);

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];

            Temp1     =  YCoordP[j] - YCcenter;
            Temp2     =  ZCoordP[j] - ZCcenter;

            YCtemp[i] =  Fcos * Temp1 - Fsin * Temp2;
            ZCtemp[i] =  Fsin * Temp1 + Fcos * Temp2;
        }

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];
            YCoordP[j] = YCtemp[i] + YCcenter;
            ZCoordP[j] = ZCtemp[i] + ZCcenter;
        }
    }

/* rotate around y - axis */

    if(Rabs(Ytemp) > 0.001) {
        Fcos = cos(Ytemp);
        Fsin = sin(Ytemp);

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];

            Temp1     =  XCoordP[j] - XCcenter;
            Temp2     =  ZCoordP[j] - ZCcenter;

            XCtemp[i] =  Fcos * Temp1 + Fsin * Temp2;
            ZCtemp[i] = -Fsin * Temp1 + Fcos * Temp2;
        }

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];
            XCoordP[j] = XCtemp[i] + XCcenter;
            ZCoordP[j] = ZCtemp[i] + ZCcenter;
        }
    }

/* rotate around z - axis */

    if(Rabs(Ztemp) > 0.001) {
        Fcos = cos(Ztemp);
        Fsin = sin(Ztemp);

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];

            Temp1     =  XCoordP[j] - XCcenter;
            Temp2     =  YCoordP[j] - YCcenter;

            XCtemp[i] =  Fcos * Temp1 - Fsin * Temp2;
            YCtemp[i] =  Fsin * Temp1 + Fcos * Temp2;
        }

        for(i = 0 ; i < ent_list ; i++) {
            j = sel_list[i];
            XCoordP[j] = XCtemp[i] + XCcenter;
            YCoordP[j] = YCtemp[i] + YCcenter;
        }
    }


    free(XCtemp);
    free(YCtemp);
    free(ZCtemp);
    return(0);
}
#endif
/************************************************************************/
int  gomp_SetSelectionCenter(float XCenter, float YCenter, float ZCenter)
/************************************************************************/
{
    ManSelection.SelCenter[0]   = XCenter;
    ManSelection.SelCenter[1]  = YCenter;
    ManSelection.SelCenter[2] = ZCenter;

    return(0);
}
/************************************************************************/
int    gomp_MoveSelectionCenter(float XMove, float YMove, float ZMove) 
/************************************************************************/
{
    ManSelection.SelCenter[0]    += XMove;
    ManSelection.SelCenter[1]   += YMove;
    ManSelection.SelCenter[2]  += ZMove;

    return(0);
}
/************************************************************************/
const float *gomp_GetSelectionCenter()
/************************************************************************/
{
    return(ManSelection.SelCenter);
}
/************************************************************************/
int gomp_CalcSelectionCenter()
/************************************************************************/
{
    static const int *sel_list; /* selection list */
    static int  ent_list; /* entries in the selection list */
    static int  i,j;
    static int atom_list;
    static float XCcenter;
    static float YCcenter;
    static float ZCcenter;
    static float Fentries;
    static const float *XCoordP;
    static const float *YCoordP;
    static const float *ZCoordP;


    atom_list = gomp_GetTotalNumberOfAtoms();

    if(atom_list < 1) {
        gomp_PrintMessage("?ERROR - no atoms defined");
        return(1);
    }

    if(ManSelection.slong < 1) {
        gomp_PrintWARNING("no atom list selected");
        return(1);
    }

    XCoordP   = gomp_GetAtomXCoordPointer(ManSelection.Structure);
    YCoordP  = gomp_GetAtomYCoordPointer(ManSelection.Structure);
    ZCoordP = gomp_GetAtomZCoordPointer(ManSelection.Structure);

    sel_list =  ManSelection.sel_list;
    ent_list =  ManSelection.slong;

    XCcenter = YCcenter = ZCcenter = 0.0;

    for(i = 0 ; i < ent_list ; i++) {
        j = sel_list[i];
        XCcenter += XCoordP[j];
        YCcenter += YCoordP[j];
        ZCcenter += ZCoordP[j];
    }

    Fentries = 1.0 / ((float)ent_list);

    XCcenter *= Fentries;
    YCcenter *= Fentries;
    ZCcenter *= Fentries;

    (void)gomp_SetSelectionCenter(XCcenter , YCcenter , ZCcenter);

    return(0);
}
#if 0
/************************************************************************/
int    gomp_InvalidateAtomSelectionGraphics()
/************************************************************************/
{
    if( ManSelection.CallbackHandle )
        gomp_InvalidatePlotterState(ManSelection.CallbackHandle);

    return(0);
}
#endif
/************************************************************************/
int    gomp_UpdateSelectionList(int Which , int ListEnt , const int *ListPointer)
/************************************************************************/
{
    static int i;

    free(ManSelection.sel_list);

    ManSelection.sel_list =  gomp_AllocateIntVector(ListEnt);

    for(i = 0 ; i < ListEnt ; i++)
        ManSelection.sel_list[i] = ListPointer[i];

    ManSelection.slong     = ListEnt;
    ManSelection.Structure = Which;

    if( !ManSelection.CallbackHandle )
        ManSelection.CallbackHandle = gomp_RegisterPlotter(
            gomp_PlotSelectedAtoms,NULL,
            PLOTTER_NAME_ATOMSELECTION,PLOTTER_ORDER_ATOMSELECTION);

    return(0);
}

/************************************************************************/
int    gomp_DeleteSelectionList()
/************************************************************************/
{
    free(ManSelection.sel_list);
    ManSelection.sel_list  = NULL;
    ManSelection.slong     = 0;
    ManSelection.Structure = 0;
    gomp_SetSelectionCenter( 0.0 , 0.0 , 0.0 );
    
    gomp_UnregisterPlotter(ManSelection.CallbackHandle);
    ManSelection.CallbackHandle = NULL;

    return(0);
}

/************************************************************************/
int    gomp_IsSelectionListEmpty()
/************************************************************************/
{
    return(ManSelection.slong==0);
}

/************************************************************************/
int    gomp_GetSelectionListLength(int Which)
/************************************************************************/
{
    if(ManSelection.Structure != Which)
        return(0);
    return(ManSelection.slong);
}

/************************************************************************/
const int *gomp_GetSelectionList(int Which)
/************************************************************************/
{
    if(ManSelection.Structure != Which)
        return(NULL);
    return(ManSelection.sel_list);
}
