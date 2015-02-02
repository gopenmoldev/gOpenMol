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
#include "gommath.h"
#include <string.h>
#include <ctype.h>
#include <tcl.h>

#include "colouring.h"
#include "correl.h"
#include "gomclipbrd.h"
#include "gommonitor.h"
#include "math_oper.h"
#include "measure.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "plot_prop.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "trajectory.h"
#include "tclutils.h"
#include "selection.h"

#include "stdafx.h"

#define DIST_LIST   0
#define ANGLE_LIST  1

#define MAX_FLOAT  1.e+35    /* This max float is just a guess */
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define IABS(a)    ( ( a ) > 0   ? (a) : -(a))
#define RABS(a)    ( ( a ) > 0.0 ? (a) : -(a))
#define SCALE(a,b,c)     scaleO(a,b,c) /* use own scale function */

/* manipulation list */
/*
#define     DAVERAGE      1
#define     SQUARE        2
#define     COS           3
#define     COS2          4
#define     SQRT          5
#define     DINITIAL      6
#define     COPYNR        7
#define     ADDNR         8
#define     LOG           9
#define     EXP          10
#define     POWERREAL    11
#define     MULTREAL     12
#define     DIVIDEREAL   13
#define     SHIFTREAL    14
#define     DMIN         15
#define     ABS          16
#define     DIVFIRST     17
#define     DIVMAXIMUM   18
#define     PSPECTRUM    19
#define     ZERO         20
*/
#define DIST_TYPE  0
#define ANG_TYPE   1
#define TORS_TYPE  2

static int    MonitorLength = 0;
     
/* distance monitor array */
static int    dist_len= 0;  /* length of distance array (in fact 2 times) */
static gom_Plotter  *dist_callback_handle = NULL;
static int *dist_vec       = NULL;
static int *dist_vec_type  = NULL;
static float *dist_vec_color = NULL;

/* angle monitor array */
static int    ang_len = 0;       /* length of angle array (3 times) */
static gom_Plotter  *ang_callback_handle = NULL;
static int *ang_vec       = NULL;
static int *ang_vec_type  = NULL;
static float *ang_vec_color = NULL;

/* torsion monitor array */
static int    tors_len= 0;       /* length of torsion array (4 times) */
static gom_Plotter  *tors_callback_handle = NULL;
static int *tors_vec       = NULL;
static int *tors_vec_type  = NULL;
static float *tors_vec_color = NULL;

static int    cc_len  = 0;       /* length of coordinate centre list */
static float *cc_vec = NULL;
static int    mc_len  = 0;       /* length of mass centre list       */
static float *mc_vec = NULL;
static int    num_obs = 0;
static float *obs_vec = NULL;
static int    obs_vec_point = 0;

/* distance */
static int         NumberDistanceSeries;
static TimeSeries *DistanceSeries;

/* angle */
static int         NumberAngleSeries;
static TimeSeries *AngleSeries;

/* torsion */
static int         NumberTorsionSeries;
static TimeSeries *TorsionSeries;

/* correlation data structure */
static corr_info_t corr_info;

/*
  DAVErage  ;    Q(t) = Q(t) - <Q(t)>
  SQUAre  ;      Q(t) = Q(t) ** 2
  COS  ;         Q(t) = cos(Q(t))  
  COS2  ;        Q(t) = 3*cos(Q(t))**2 - 1  
  SQRT  ;        Q(t) = sqrt(Q(t)) 
  DINItial  ;    Q(t) = Q(t) - Q(0) 
  COPY nr ;      Q(t) = Q2(t)  
  ADD nr ;       Q(t) = Q(t) + Q2(t)   
  LOG ;          Q(t) = log(Q(t)) 
  EXP ;          Q(t) = exp(Q(t)) 
  POWer real ;   Q(t) = Q(t) ** real 
  MULT real ;    Q(t) = real * Q(t) 
  DIVIde real ;  Q(t) = Q(t) / real 
  SHIFt real ;   Q(t) = Q(t) + real  
  DMIN ;         Q(t) = Q(t) - Q(min)
  ABS ;          Q(t) = abs(Q(t))
  DIVFirst ;     Q(t) = Q(t) / Q(0) 
  DIVMaximum ;   Q(t) = Q(t) /max(Q(t)) 
  PSPEctrum nr ; Power spectrum
  ZERO ;         Q(t) = 0.0

*/

/**************************************************************************/
static int CalcVecAv(const float *vec_list,int numbers, int alt)  
    /* calc average for the list */
/*
  const float *vec_list;     list of values 
  int numbers;         number of elements in the list 
  int alt;             if = 0 normal list if > 0 it's an angle list 
*/
/**************************************************************************/
{
    register int i;
    register float vec_sum,devi;
    static   float vmin,vmax;
    static   float veci,pia2;
    static   int imin,imax;
    static   float vecx,vecy,vecxi,vecyi,vecxii,vecyii;
    static   float d2ij,d2jk,d2ik;
    static   int swap_flag;
    static   char OutText[BUFF_LEN];
    static   float pia = (float)M_PI;

    pia2 = pia / 180.0;
    vmin =   1.e+20f;
    vmax =      0.0;
    imin = -1;
    imax = -1;

    vecx = 0.0;
    vecy = 0.0;

    vec_sum=0.0;
    devi=0.0;
    swap_flag = 0;

    switch(alt) {

    case 0:
        for(i=0;i<numbers;i++) vec_sum += vec_list[i];
        vec_sum = vec_sum / (float)numbers;
        break;

    default:

        for(i=0;i<numbers;i++) {
            veci = vec_list[i];
            devi = veci * pia2;
            vecx += cos((double)devi);
            vecy += sin((double)devi);
        }
        vecx /= (float)numbers;
        vecy /= (float)numbers;

        vec_sum = asin((double)(vecy / sqrt(vecx*vecx + vecy*vecy))) / pia * 180.0;
        swap_flag = 0;
        if(vecx < 0.0) {
            vec_sum = 180.0 - vec_sum;
            if(vec_sum > 180.0) swap_flag = 1;
        } 
        break;
    }

    devi = 0.0;

    switch(alt) {

    case 0:
        for(i=0;i<numbers;i++)  {
            veci = vec_list[i] - vec_sum;
            devi += veci * veci;
            if(RABS(veci) < RABS(vmin)) {
                vmin = veci;
                imin = i + 1;
            }

            if(RABS(veci) > RABS(vmax)) {
                vmax = veci;
                imax = i + 1;
            }
        }
        break;

    default:
        devi = 0.0;
        d2ij = vecx*vecx + vecy*vecy;
        for(i = 0 ; i < numbers ; i++)  {
            veci = vec_list[i] * pia2;
            vecxii=cos((double)veci);
            vecyii=sin((double)veci);
            d2jk = vecxii*vecxii + vecyii*vecyii; 
            vecxi = vecx - vecxii;
            vecyi = vecy - vecyii;
            d2ik = vecxi*vecxi + vecyi*vecyi;
            veci = 0.5 *(d2ij + d2jk - d2ik) / sqrt(d2ij*d2jk);
            veci = (veci >  1.0 ?  1.0 : veci);
            veci = (veci < -1.0 ? -1.0 : veci);
            veci = acos((double)veci) / pia * 180.0;

            devi += veci * veci;
   
            if(RABS(veci) < RABS(vmin)) {
                vmin = veci;
                imin = i + 1;
            }

            if(RABS(veci) > RABS(vmax)) {
                vmax = veci;
                imax = i + 1;
            }
        }
        break;
    }

    if(swap_flag) vec_sum = vec_sum - 360.0;

    if(numbers > 1)
        devi = devi /(float)(numbers-1);

    devi = sqrt(devi);

    sprintf(OutText,"Av. value , max. diff.  (entry nr '%d') , min. diff.  \
(entry nr '%d') , st. dev.",imax,imin);
    gomp_PrintMessage(OutText);
    sprintf(OutText," %f     %f                 %f                 %f ",vec_sum,vmax,vmin,devi);
    gomp_PrintMessage(OutText);

    return(0);
}

/**********************************************************************/
int gomp_SelectDistArray(const char *Seg1, const char *Res1, const char *Atm1,
                       const char *Seg2, const char *Res2, const char *Atm2,
                       const char *Type, const char *Colour)
/**********************************************************************/
{
    static int    i,j,si,sj,loop;
    static int *sel_list1;
    static int    slong1;
    static int *sel_list2;
    static int    slong2,atom_max;
    static float  dist,tmp1,tmp2,tmp3;
    static int    cross;
    static const int *res1;
    char OutText[BUFF_LEN];
    char TAtm[MAX_ATM_NAME_LEN];
    char TRes[MAX_RES_NAME_LEN];
    char TSeg[MAX_SEG_NAME_LEN];
    int  Wstr;
    const float *x;
    const float *y;
    const float *z;
    float  RedC;
    float  GreenC;
    float  BlueC;
    int    LineType;
    
    Wstr      = 0;
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);
    res1      = gomp_GetAtomResNum1Pointer(Wstr);

    atom_max  = gomp_GetNumAtomsInMolecStruct(0);
    sel_list1 = gomp_AllocateIntVector(atom_max);
    sel_list2 = gomp_AllocateIntVector(atom_max);

    slong1 = gomp_MakeSelectionList(0,Seg1,Res1,Atm1,sel_list1);
    slong2 = gomp_MakeSelectionList(0,Seg2,Res2,Atm2,sel_list2);

    if(slong1 > 0 && slong2 > 0) {

        if(*Type != (char)NULL)
            LineType = atoi(Type);
        else
            LineType = 1;

        if(*Colour != (char)NULL) {

            if(gomp_ColourName2RGB(Colour , &RedC , &GreenC , &BlueC)) {
                sprintf(OutText,"can't resolve colour name '%s', will use default colour",Colour);
                gomp_PrintERROR(OutText);
                RedC    = 1.0;
                GreenC  = 0.0;
                BlueC   = 1.0;
            }
        } else {
            RedC    = 1.0;
            GreenC  = 0.0;
            BlueC   = 1.0;
        }

        cross = 0;

/*     *dist_len = 0; */
        loop = 0;
        for(i = 0 ; i < slong1 ; i++ ) {
            si = sel_list1[i];
            for(j = 0 ; j < slong2 ; j++) {
                sj = sel_list2[j];

                if(cross) {
                    if(res1[si] != res1[sj])
                        continue;
                }

                loop++;
                tmp1 = (x[si]-x[sj]);
                tmp2 = (y[si]-y[sj]);
                tmp3 = (z[si]-z[sj]);

                dist = sqrt( tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3);

                strncpy(TAtm,gomp_GetAtomAtmName(Wstr , si),MAX_ATM_NAME_LEN);
                strncpy(TRes,gomp_GetAtomResName(Wstr , si),MAX_RES_NAME_LEN);
                strncpy(TSeg,gomp_GetAtomSegName(Wstr , si),MAX_SEG_NAME_LEN);

                sprintf(OutText,"Distance >%.4s<>%.4s(%d)<>%.4s<  -  >%.4s<>%.4s(%d)<>%.4s<",
                        TSeg,TRes,gomp_GetAtomResNum1(Wstr , si),TAtm,
                        gomp_GetAtomSegName(Wstr , sj),gomp_GetAtomResName(Wstr , sj),
                        gomp_GetAtomResNum1(Wstr , sj),gomp_GetAtomAtmName(Wstr , sj));
                gomp_PrintMessage(OutText);
                sprintf(OutText," %f ",dist);
                gomp_PrintMessage(OutText);

                gomp_ExpandDistVec(1);

                dist_vec[(dist_len - 2)    ] = si;
                dist_vec[(dist_len - 2) + 1] = sj;

                dist_vec_type[((dist_len - 2) / 2)]          = LineType;

                dist_vec_color[3 * ((dist_len - 2) / 2)    ] = RedC;
                dist_vec_color[3 * ((dist_len - 2) / 2) + 1] = GreenC;
                dist_vec_color[3 * ((dist_len - 2) / 2) + 2] = BlueC;
            }
        }
    }
    else
        gomp_PrintMessage("?ERROR - no atoms in the selection list 1 or 2 ");

    sprintf(OutText,"Selected %d distance(s) from the list ",loop);
    gomp_PrintMessage(OutText);

    free(sel_list1);
    free(sel_list2);

    return(0);
}

/**********************************************************************/
int gomp_SelectAngArray(const char *Seg1, const char *Res1, const char *Atm1,
                      const char *Seg2, const char *Res2, const char *Atm2,
                      const char *Seg3, const char *Res3, const char *Atm3,
                      const char *Type, const char *Colour)
/**********************************************************************/
{

    static int i,j,k,si,sj,sk,loop,atom_max;
    static int *sel_list1;
    static int    slong1;
    static int *sel_list2;
    static int    slong2;
    static int *sel_list3;
    static int    slong3;
    static float  angle;
    static int    cross;
    static const int *res1;
    char OutText[BUFF_LEN];
    char TAtm[MAX_ATM_NAME_LEN];
    char TRes[MAX_RES_NAME_LEN];
    char TSeg[MAX_SEG_NAME_LEN];
    char TAtm1[MAX_ATM_NAME_LEN];
    char TRes1[MAX_RES_NAME_LEN];
    char TSeg1[MAX_SEG_NAME_LEN];
    int  Wstr;
    const float *x;
    const float *y;
    const float *z;
    float  RedC;
    float  GreenC;
    float  BlueC;
    int    LineType;

    Wstr      = 0;
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);
    res1      = gomp_GetAtomResNum1Pointer(Wstr);

    atom_max  = gomp_GetNumAtomsInMolecStruct(0);
    sel_list1 = gomp_AllocateIntVector(atom_max);
    sel_list2 = gomp_AllocateIntVector(atom_max);
    sel_list3 = gomp_AllocateIntVector(atom_max);

    slong1 = gomp_MakeSelectionList(0,Seg1,Res1,Atm1,sel_list1);
    slong2 = gomp_MakeSelectionList(0,Seg2,Res2,Atm2,sel_list2);
    slong3 = gomp_MakeSelectionList(0,Seg3,Res3,Atm3,sel_list3);

    if(slong1 > 0 && slong2 > 0 && slong3 > 0) {

        if(*Type != (char)NULL) {
            LineType = atoi(Type);
        }
        else {
            LineType = 2;
        }

        if(*Colour != (char)NULL) {

            if(gomp_ColourName2RGB(Colour , &RedC , &GreenC , &BlueC)) {
                sprintf(OutText,"can't resolve colour name '%s', will use default colour",Colour);
                gomp_PrintERROR(OutText);
                RedC    = 1.0;
                GreenC  = 0.0;
                BlueC   = 1.0;
            }
        } else {
            RedC    = 1.0;
            GreenC  = 0.0;
            BlueC   = 1.0;
        }

        cross = 0;

/*     *dist_len = 0; */
        loop = 0;
        for(i = 0 ; i < slong1 ; i++ ) {
            si = sel_list1[i];
            for(j = 0 ; j < slong2 ; j++) {
                sj = sel_list2[j];
                for(k = 0 ; k < slong3 ; k++) {
                    sk = sel_list3[k];

                    if(cross) {
                        if(res1[si] != res1[sj] ||
                           res1[si] != res1[sk] || 
                           res1[sj] != res1[sk]) continue;
                    }

                    loop++;

                    gomp_BondAngle(x[si], y[si] , z[si] ,
                                 x[sj], y[sj] , z[sj] ,
                                 x[sk], y[sk] , z[sk] , &angle);

                    strncpy(TAtm,gomp_GetAtomAtmName(Wstr , si),MAX_ATM_NAME_LEN);
                    strncpy(TRes,gomp_GetAtomResName(Wstr , si),MAX_RES_NAME_LEN);
                    strncpy(TSeg,gomp_GetAtomSegName(Wstr , si),MAX_SEG_NAME_LEN);

                    strncpy(TAtm1,gomp_GetAtomAtmName(Wstr , sj),MAX_ATM_NAME_LEN);
                    strncpy(TRes1,gomp_GetAtomResName(Wstr , sj),MAX_RES_NAME_LEN);
                    strncpy(TSeg1,gomp_GetAtomSegName(Wstr , sj),MAX_SEG_NAME_LEN);

                    sprintf(OutText,"Angle >%.4s<>%.4s(%d)<>%.4s<->%.4s<>%.4s(%d)<>%.4s<->%.4s<>%.4s(%d)<>%.4s<",
                            TSeg ,TRes ,gomp_GetAtomResNum1(Wstr , si),TAtm,
                            TSeg1,TRes1,gomp_GetAtomResNum1(Wstr , sj),TAtm1,
                            gomp_GetAtomSegName(Wstr , sk),gomp_GetAtomResName(Wstr , sk),
                            gomp_GetAtomResNum1(Wstr , sk),gomp_GetAtomAtmName(Wstr , sk));
                    gomp_PrintMessage(OutText);

                    sprintf(OutText," %f ",180.*angle/M_PI);
                    gomp_PrintMessage(OutText);

                    gomp_ExpandAngVec(1);

                    ang_vec[(ang_len - 3)    ] = si;
                    ang_vec[(ang_len - 3) + 1] = sj;
                    ang_vec[(ang_len - 3) + 2] = sk;

                    ang_vec_type[((ang_len - 3) / 3)]             = LineType;

                    ang_vec_color[3 * ((ang_len - 3) / 3)    ]    = RedC;
                    ang_vec_color[3 * ((ang_len - 3) / 3) + 1]    = GreenC;
                    ang_vec_color[3 * ((ang_len - 3) / 3) + 2]    = BlueC;
                }

            }
        }
    }
    else
        gomp_PrintMessage("?ERROR - no atoms in selection list 1,2 or 3"); 

    sprintf(OutText,"Selected %d angle(s) from the list ",loop);
    gomp_PrintMessage(OutText);

    free(sel_list1);
    free(sel_list2);
    free(sel_list3);

    return(0);
}

/**********************************************************************/
int gomp_SelectTorsArray(const char *Seg1, const char *Res1, const char *Atm1,
                       const char *Seg2, const char *Res2, const char *Atm2,
                       const char *Seg3, const char *Res3, const char *Atm3,
                       const char *Seg4, const char *Res4, const char *Atm4,
                       const char *Type, const char *Colour)
/**********************************************************************/
{

    static int    i,j,k,l,si,sj,sk,sl,loop,atom_max;
    static int *sel_list1;
    static int    slong1;
    static int *sel_list2;
    static int    slong2;
    static int *sel_list3;
    static int    slong3;
    static int *sel_list4;
    static int    slong4;
    static float  angle;
    static int    cross;
    static const int *res1;
    char OutText[BUFF_LEN];
    char TAtm[MAX_ATM_NAME_LEN];
    char TRes[MAX_RES_NAME_LEN];
    char TSeg[MAX_SEG_NAME_LEN];
    char TAtm1[MAX_ATM_NAME_LEN];
    char TRes1[MAX_RES_NAME_LEN];
    char TSeg1[MAX_SEG_NAME_LEN];
    char TAtm2[MAX_ATM_NAME_LEN];
    char TRes2[MAX_RES_NAME_LEN];
    char TSeg2[MAX_SEG_NAME_LEN];
    int  Wstr;
    const float *x;
    const float *y;
    const float *z;
    float  RedC;
    float  GreenC;
    float  BlueC;
    int    LineType;

    Wstr      = 0;
    x         = gomp_GetAtomXCoordPointer(Wstr);
    y         = gomp_GetAtomYCoordPointer(Wstr);
    z         = gomp_GetAtomZCoordPointer(Wstr);
    res1      = gomp_GetAtomResNum1Pointer(Wstr);


    atom_max = gomp_GetNumAtomsInMolecStruct(0);

    sel_list1 = gomp_AllocateIntVector(atom_max);
    sel_list2 = gomp_AllocateIntVector(atom_max);
    sel_list3 = gomp_AllocateIntVector(atom_max);
    sel_list4 = gomp_AllocateIntVector(atom_max);

    slong1 = gomp_MakeSelectionList(0,Seg1,Res1,Atm1,sel_list1);
    slong2 = gomp_MakeSelectionList(0,Seg2,Res2,Atm2,sel_list2);
    slong3 = gomp_MakeSelectionList(0,Seg3,Res3,Atm3,sel_list3);
    slong4 = gomp_MakeSelectionList(0,Seg4,Res4,Atm4,sel_list4);

    if(slong1 > 0 && slong2 > 0 && slong3 > 0 && slong4 > 0) {

        if(*Type != (char)NULL) {
            LineType = atoi(Type);
        }
        else {
            LineType = 3;
        }

        if(*Colour != (char)NULL) {

            if(gomp_ColourName2RGB(Colour , &RedC , &GreenC , &BlueC)) {
                sprintf(OutText,"can't resolve colour name '%s', will use default colour",Colour);
                gomp_PrintERROR(OutText);
                RedC    = 1.0;
                GreenC  = 0.0;
                BlueC   = 1.0;
            }
        } else {
            RedC    = 1.0;
            GreenC  = 0.0;
            BlueC   = 1.0;
        }

        cross = 0;

/*     *dist_len = 0; */
        loop = 0;
        for(i = 0 ; i < slong1 ; i++ ) {
            si = sel_list1[i];
            for(j = 0 ; j < slong2 ; j++) {
                sj = sel_list2[j];
                for(k = 0 ; k < slong3 ; k++) {
                    sk = sel_list3[k];
                    for(l = 0 ; l < slong4 ; l++) {
                        sl = sel_list4[l];

                        if(cross) {
                            if(res1[si] != res1[sj] || 
                               res1[si] != res1[sk] || 
                               res1[si] != res1[sl] || 
                               res1[sj] != res1[sk] ||
                               res1[sj] != res1[sl] || 
                               res1[sk] != res1[sl]) continue;
                        }

                        loop++;

                        gomp_floDihedAngle(x[si], y[si] , z[si] ,
                                      x[sj], y[sj] , z[sj] ,
                                      x[sk], y[sk] , z[sk] , 
                                      x[sl], y[sl] , z[sl] , &angle);

                        strncpy(TAtm,gomp_GetAtomAtmName(Wstr , si),MAX_ATM_NAME_LEN);
                        strncpy(TRes,gomp_GetAtomResName(Wstr , si),MAX_RES_NAME_LEN);
                        strncpy(TSeg,gomp_GetAtomSegName(Wstr , si),MAX_SEG_NAME_LEN);

                        strncpy(TAtm1,gomp_GetAtomAtmName(Wstr , sj),MAX_ATM_NAME_LEN);
                        strncpy(TRes1,gomp_GetAtomResName(Wstr , sj),MAX_RES_NAME_LEN);
                        strncpy(TSeg1,gomp_GetAtomSegName(Wstr , sj),MAX_SEG_NAME_LEN);

                        strncpy(TAtm2,gomp_GetAtomAtmName(Wstr , sk),MAX_ATM_NAME_LEN);
                        strncpy(TRes2,gomp_GetAtomResName(Wstr , sk),MAX_RES_NAME_LEN);
                        strncpy(TSeg2,gomp_GetAtomSegName(Wstr , sk),MAX_SEG_NAME_LEN);

                        sprintf(OutText,"Torsion >%.4s<>%.4s(%d)<>%.4s<->%.4s<>%.4s(%d)<>%.4s<->%.4s<>%.4s(%d)<>%.4s<->%.4s<>%.4s(%d)<>%.4s<",
                                TSeg ,TRes ,gomp_GetAtomResNum1(Wstr , si),TAtm,
                                TSeg1,TRes1,gomp_GetAtomResNum1(Wstr , sj),TAtm1,
                                TSeg2,TRes2,gomp_GetAtomResNum1(Wstr , sk),TAtm2,
                                gomp_GetAtomSegName(Wstr , sl),gomp_GetAtomResName(Wstr , sl),
                                gomp_GetAtomResNum1(Wstr , sl),gomp_GetAtomAtmName(Wstr , sl));
                        gomp_PrintMessage(OutText);

                        sprintf(OutText," %f ",180.*angle/M_PI);
                        gomp_PrintMessage(OutText);

                        gomp_ExpandTorsVec(1);

                        tors_vec[(tors_len - 4)    ] = si;
                        tors_vec[(tors_len - 4) + 1] = sj;
                        tors_vec[(tors_len - 4) + 2] = sk;
                        tors_vec[(tors_len - 4) + 3] = sl;

                        tors_vec_type[((tors_len - 4) / 4)]             = LineType;

                        tors_vec_color[3 * ((tors_len - 4) / 4)    ]    = RedC;
                        tors_vec_color[3 * ((tors_len - 4) / 4) + 1]    = GreenC;
                        tors_vec_color[3 * ((tors_len - 4) / 4) + 2]    = BlueC;
                    }

                }
            }
        }
    }
    else
        gomp_PrintMessage("?ERROR - no atoms in the selection list 1,2,3 or 4 ");

    sprintf(OutText,"Selected %d torsion angle(s) from the list ",loop);
    gomp_PrintMessage(OutText);

    free(sel_list1);
    free(sel_list2);
    free(sel_list3);
    free(sel_list4);

    return(0);
}
/************************************************************************/
int    gomp_ExpandDistVec(int NewLength)
/************************************************************************/
{
    if(!MonitorLength) {
        dist_vec       = gomp_AllocateIntVector(2 * NewLength);
        dist_vec_type  = gomp_AllocateIntVector(    NewLength);
        dist_vec_color = gomp_AllocateFloatVector(3 * NewLength);

        MonitorLength = NewLength;
        dist_len      = 2;
    }
    else {
        if(MonitorLength != NewLength) { /* reset all monitor arrays */
            (void)gomp_ResetMonitorData();
        }
        if(!dist_len) {
            dist_vec       = gomp_AllocateIntVector(2 * NewLength);
            dist_vec_type  = gomp_AllocateIntVector(    NewLength);
            dist_vec_color = gomp_AllocateFloatVector(3 * NewLength);
            MonitorLength = NewLength;
            dist_len      = 2;
        }
        else { 
            dist_len      += 2;
            dist_vec       = gomp_ReallocateIntVector(dist_vec       , NewLength 
                                           *      dist_len);
            dist_vec_type  = gomp_ReallocateIntVector(dist_vec_type  , NewLength
                                           *     (dist_len / 2));
            dist_vec_color = gomp_ReallocateFloatVector(dist_vec_color , NewLength
                                           * 3 * (dist_len / 2));
            MonitorLength  = NewLength;
        }
    }

    return(0);
}

/************************************************************************/
int    gomp_ExpandAngVec(int NewLength)
/************************************************************************/
{
    if(!MonitorLength) {
        ang_vec       = gomp_AllocateIntVector(3 * NewLength);
        ang_vec_type  = gomp_AllocateIntVector(    NewLength);
        ang_vec_color = gomp_AllocateFloatVector(3 * NewLength);

        MonitorLength = NewLength;
        ang_len       = 3;
    }
    else {
        if(MonitorLength != NewLength) { /* reset all monitor arrays */
            (void)gomp_ResetMonitorData();
        }
        if(!ang_len) {
            ang_vec       = gomp_AllocateIntVector(3 * NewLength);
            ang_vec_type  = gomp_AllocateIntVector(    NewLength);
            ang_vec_color = gomp_AllocateFloatVector(3 * NewLength);

            MonitorLength = NewLength;
            ang_len       = 3;
        }
        else { 
            ang_len      += 3;
            ang_vec       = gomp_ReallocateIntVector(ang_vec       , NewLength 
                                          * ang_len);
            ang_vec_type  = gomp_ReallocateIntVector(ang_vec_type  , NewLength 
                                          *     (ang_len / 3));
            ang_vec_color = gomp_ReallocateFloatVector(ang_vec_color , NewLength 
                                          * 3 * (ang_len / 3));
            MonitorLength = NewLength;
        }
    }

    return(0);
}
/************************************************************************/
int    gomp_ExpandTorsVec(int NewLength)
/************************************************************************/
{
    if(!MonitorLength) {
        tors_vec       = gomp_AllocateIntVector(4 * NewLength);
        tors_vec_type  = gomp_AllocateIntVector(    NewLength);
        tors_vec_color = gomp_AllocateFloatVector(3 * NewLength);

        MonitorLength  = NewLength;
        tors_len       = 4;
    }
    else {
        if(MonitorLength != NewLength) { /* reset all monitor arrays */
            (void)gomp_ResetMonitorData();
        }
        if(!tors_len) {
            tors_vec       = gomp_AllocateIntVector(4 * NewLength);
            tors_vec_type  = gomp_AllocateIntVector(    NewLength);
            tors_vec_color = gomp_AllocateFloatVector(3 * NewLength);

            MonitorLength  = NewLength;
            tors_len       = 4;
        }
        else { 
            tors_len      += 4;
            tors_vec       = gomp_ReallocateIntVector(tors_vec ,       NewLength 
                                           * tors_len);
            tors_vec_type  = gomp_ReallocateIntVector(tors_vec_type ,  NewLength 
                                           *     (tors_len / 4));
            tors_vec_color = gomp_ReallocateFloatVector(tors_vec_color , NewLength 
                                           * 3 * (tors_len / 4));
            MonitorLength  = NewLength;
        }
    }

    return(0);
}
#if 0
/************************************************************************/
int    gomp_InvalidateMonitorGraphics()
/************************************************************************/
{
    if( dist_callback_handle )
        gomp_InvalidatePlotterState(dist_callback_handle);

    if( ang_callback_handle )
        gomp_InvalidatePlotterState(ang_callback_handle);

    if( tors_callback_handle )
        gomp_InvalidatePlotterState(tors_callback_handle);

    return(0);
}
#endif
/************************************************************************/
int    gomp_ResetMonitorData()
/************************************************************************/
{

    MonitorLength = 0;
     
    gomp_UnregisterPlotter(dist_callback_handle);
    gomp_UnregisterPlotter(ang_callback_handle);
    gomp_UnregisterPlotter(tors_callback_handle);

    dist_callback_handle = NULL;
    ang_callback_handle  = NULL;
    tors_callback_handle = NULL;

    dist_len= 0;
    if(dist_vec)
        free(dist_vec);
    if(dist_vec_type)
        free(dist_vec_type);
    if(dist_vec_color)
        free(dist_vec_color);

    dist_vec       = NULL;
    dist_vec_type  = NULL;
    dist_vec_color = NULL;

    ang_len = 0;
    if(ang_vec)
        free(ang_vec);
    ang_vec       =  NULL;
    if(ang_vec_type)
        free(ang_vec_type);
    ang_vec_type  =  NULL;
    if(ang_vec_color)
        free(ang_vec_color);
    ang_vec_color =  NULL;

    tors_len= 0;
    if(tors_vec)
        free(tors_vec);
    tors_vec = NULL;
    if(tors_vec_type)
        free(tors_vec_type);
    tors_vec_type = NULL;
    if(tors_vec_color)
        free(tors_vec_color);
    tors_vec_color = NULL;

    cc_len  = 0;
    if(cc_vec)
        free(cc_vec);
    cc_vec =   NULL;

    mc_len  = 0;
    if(mc_vec)
        free(mc_vec);
    mc_vec =   NULL;

    num_obs = 0;
    if(obs_vec)
        free(obs_vec);
    obs_vec =  NULL;

    obs_vec_point = 0;

    return(0);
}

/************************************************************************/
int    gomp_ResetMonitorDistanceData()
/************************************************************************/
{

    dist_len     = 0;

    gomp_UnregisterPlotter(dist_callback_handle);
    dist_callback_handle = NULL;

    if(dist_vec)
        free(dist_vec);
    if(dist_vec_type)
        free(dist_vec_type);
    if(dist_vec_color)
        free(dist_vec_color);

    dist_vec       = NULL;
    dist_vec_type  = NULL;
    dist_vec_color = NULL;


    return(0);
}
/************************************************************************/
int    gomp_ResetMonitorAngleData()
/************************************************************************/
{

    ang_len     = 0;

    gomp_UnregisterPlotter(ang_callback_handle);
    ang_callback_handle = NULL;

    if(ang_vec)
        free(ang_vec);
    ang_vec =  NULL;
    if(ang_vec_type)
        free(ang_vec_type);
    ang_vec_type =  NULL;
    if(ang_vec_color)
        free(ang_vec_color);
    ang_vec_color =  NULL;

    return(0);
}
/************************************************************************/
int    gomp_ResetMonitorTorsionData()
/************************************************************************/
{

    tors_len     = 0;

    gomp_UnregisterPlotter(tors_callback_handle);
    tors_callback_handle = NULL;

    if(tors_vec)
        free(tors_vec);
    tors_vec = NULL;
    if(tors_vec_type)
        free(tors_vec_type);
    tors_vec_type = NULL;
    if(tors_vec_color)
        free(tors_vec_color);
    tors_vec_color = NULL;

    return(0);
}

/**************************************************************************/
int gomp_SetDistMonitor(int State)
/**************************************************************************/
{
    if( State ) {
        if( !dist_callback_handle )
            gomp_RegisterPlotter(
                gomp_PlotMonitorDistance,NULL,
                PLOTTER_NAME_MONITOR_DISTANCE,
                PLOTTER_ORDER_MONITOR_DISTANCE);
    }
    else {
        gomp_UnregisterPlotter(dist_callback_handle);
        dist_callback_handle = NULL;
    }

    return(0);
}
/**************************************************************************/
int gomp_GetDistMonitor()
/**************************************************************************/
{
    return(dist_callback_handle!=NULL);
}
/**************************************************************************/
int  gomp_GetDistMonitorSamples()
/**************************************************************************/
{
    return(dist_len);
}
/**************************************************************************/
const int *gomp_GetDistMonitorSamplesList()
/**************************************************************************/
{
    return(dist_vec);
}
/**************************************************************************/
int gomp_SetAngMonitor(int State)
/**************************************************************************/
{
    if( State ) {
        if( !ang_callback_handle )
            gomp_RegisterPlotter(
                gomp_PlotMonitorAngle,NULL,
                PLOTTER_NAME_MONITOR_ANGLE,
                PLOTTER_ORDER_MONITOR_ANGLE);
    }
    else {
        gomp_UnregisterPlotter(ang_callback_handle);
        ang_callback_handle = NULL;
    }

    return(0);
}
/**************************************************************************/
int gomp_GetAngMonitor()
/**************************************************************************/
{
    return(ang_callback_handle!=NULL);
}
/**************************************************************************/
int  gomp_GetAngMonitorSamples()
/**************************************************************************/
{
    return(ang_len);
}
/**************************************************************************/
const int *gomp_GetAngMonitorSamplesList()
/**************************************************************************/
{
    return(ang_vec);
}
/**************************************************************************/
int gomp_SetTorsMonitor(int State)
/**************************************************************************/
{
    if( State ) {
        if( !tors_callback_handle )
            gomp_RegisterPlotter(
                gomp_PlotMonitorTorsion,NULL,
                PLOTTER_NAME_MONITOR_TORSION,
                PLOTTER_ORDER_MONITOR_TORSION);
    }
    else {
        gomp_UnregisterPlotter(tors_callback_handle);
        tors_callback_handle = NULL;
    }

    return(0);
}
/**************************************************************************/
int gomp_GetTorsMonitor()
/**************************************************************************/
{
    return(tors_callback_handle!=NULL);
}
/**************************************************************************/
int  gomp_GetTorsMonitorSamples()
/**************************************************************************/
{
    return(tors_len);
}
/**************************************************************************/
const int *gomp_GetTorsMonitorSamplesList()
/**************************************************************************/
{
    return(tors_vec);
}

/**************************************************************************/
int         gomp_ExpandeDistanceSeries(int Numbers)
/**************************************************************************/
{
    if(!NumberDistanceSeries) {
        DistanceSeries = gomp_AllocateVoidVector(sizeof(TimeSeries));
        DistanceSeries[0].Value    = gomp_AllocateFloatVector(Numbers);
        DistanceSeries[0].Length   = Numbers;
        NumberDistanceSeries       = 1;
    }
    else {
        printf("%d\n",DistanceSeries[0].Length);
        NumberDistanceSeries++;
        DistanceSeries = gomp_ReallocateVoidVector(
            DistanceSeries , NumberDistanceSeries *
            sizeof(*DistanceSeries));
        DistanceSeries[NumberDistanceSeries - 1].Value    
            = gomp_AllocateFloatVector(Numbers);
        DistanceSeries[NumberDistanceSeries - 1].Length
            =           Numbers;
    }

    return(NumberDistanceSeries);
}
/**************************************************************************/
int         gomp_DeleteDistanceSeries()
/**************************************************************************/
{
    int i;

    if(!NumberDistanceSeries) return(0);          /* no series defined */

    for(i = 0 ; i < NumberDistanceSeries ; i++)
        free(DistanceSeries[i].Value);

    free(DistanceSeries);
    NumberDistanceSeries = 0;

    return(0);
}
/**************************************************************************/
static int         ExpandeAngleSeries(int Numbers)
/**************************************************************************/
{
    if(!NumberAngleSeries) {
        AngleSeries = gomp_AllocateVoidVector(sizeof(TimeSeries));
        AngleSeries[0].Value    = gomp_AllocateFloatVector(Numbers);
        AngleSeries[0].Length   = Numbers;
        NumberAngleSeries       = 1;
    }
    else {
        NumberAngleSeries++;
        AngleSeries     = gomp_ReallocateVoidVector(
            AngleSeries , NumberAngleSeries *
            sizeof(*AngleSeries));
        AngleSeries[NumberAngleSeries - 1].Value    
            = gomp_AllocateFloatVector(Numbers);
        AngleSeries[NumberAngleSeries - 1].Length
            =           Numbers;
    }

    return(NumberAngleSeries);
}
/**************************************************************************/
int         gomp_DeleteAngleSeries()
/**************************************************************************/
{
    int i;

    if(!NumberAngleSeries) return(0);          /* no series defined */

    for(i = 0 ; i < NumberAngleSeries ; i++)
        free(AngleSeries[i].Value);

    free(AngleSeries);
    NumberAngleSeries = 0;

    return(0);
}
/**************************************************************************/
static int ExpandTorsionSeries(int Numbers)
/**************************************************************************/
{
    if(!NumberTorsionSeries) {
        TorsionSeries = gomp_AllocateVoidVector(sizeof(TimeSeries));
        TorsionSeries[0].Value    = gomp_AllocateFloatVector(Numbers);
        TorsionSeries[0].Length   = Numbers;
        NumberTorsionSeries       = 1;
    }
    else {
        NumberTorsionSeries++;
        TorsionSeries = gomp_ReallocateVoidVector(
            TorsionSeries , NumberTorsionSeries *
            sizeof(*TorsionSeries));
        TorsionSeries[NumberTorsionSeries - 1].Value    
            = gomp_AllocateFloatVector(Numbers);
        TorsionSeries[NumberTorsionSeries - 1].Length
            =           Numbers;

    }

    return(NumberTorsionSeries);
}
/**************************************************************************/
int         gomp_DeleteTorsionSeries()
/**************************************************************************/
{
    int i;

    if(!NumberTorsionSeries) return(0);          /* no series defined */

    for(i = 0 ; i < NumberTorsionSeries ; i++)
        free(TorsionSeries[i].Value);

    free(TorsionSeries);
    NumberTorsionSeries = 0;

    return(0);
}
/**************************************************************************/
int   gomp_StoreValues2DistanceSeries(int Which, int Numbers, const float *Values)
/**************************************************************************/
{
    return(0);
}
/*
  corr_info.corr_val = vector(2 * gomp_check_pow_2(FramesInSet()));
  pre_correl(obs_vec_point,obs_vec,FramesInSet(),corr_info.corr_va\
*/
/***************************************************************************/
int    gomp_CalculateCorrelation(int Type, int Index1, int Index2)
/***************************************************************************/
{
    if(gomp_CheckTimeSeriesIndex(Type , Index1))
        return(1);

    if(gomp_CheckTimeSeriesIndex(Type , Index2))
        return(1);

    switch(Type) {

    case DIST_TYPE:
        corr_info.corr_vec1 = DistanceSeries[Index1 - 1].Value;
        corr_info.corr_vec2 = DistanceSeries[Index2 - 1].Value;
        corr_info.corr_obs  = DistanceSeries[Index1 - 1].Length;
        return(gomp_PreCorrel());
        break;
    case ANG_TYPE:
        corr_info.corr_vec1 = AngleSeries[Index1 - 1].Value;
        corr_info.corr_vec2 = AngleSeries[Index2 - 1].Value;
        corr_info.corr_obs  = AngleSeries[Index1 - 1].Length;
        return(gomp_PreCorrel());
        break;
    case TORS_TYPE:
        corr_info.corr_vec1 = TorsionSeries[Index1 - 1].Value;
        corr_info.corr_vec2 = TorsionSeries[Index2 - 1].Value;
        corr_info.corr_obs  = TorsionSeries[Index1 - 1].Length;
        return(gomp_PreCorrel());
        break;
    }

    return(0);
}
/***************************************************************************/
int    gomp_PutCorrelationVectorAddress1(float *Address)
/***************************************************************************/
{
    corr_info.corr_vec1 = Address;

    return(0);
}
/***************************************************************************/
int    gomp_PutCorrelationVectorAddress2(float *Address)
/***************************************************************************/
{
    corr_info.corr_vec2 = Address;

    return(0);
}
/***************************************************************************/
const float *gomp_GetCorrelationResultVector()
/***************************************************************************/
{
    return(corr_info.corr_val);
}
/***************************************************************************/
int    gomp_GetCorrelationResultVectorLength()
/***************************************************************************/
{
    return(corr_info.corr_obs);
}
/***************************************************************************/
int gomp_PreCorrel()
/***************************************************************************/
{

    static float *fft;
    static float *data1;
    static float *data2;
    static float fhelp1,fhelp2;
    static int n,n2;
    static int i;
         
    n     = corr_info.corr_obs;
    n2    = gomp_check_pow_2(n);
/* temprary storage */
    fft   = gomp_AllocateFloatVector(2 * n2);
    data1 = gomp_AllocateFloatVector(n2);
    data2 = gomp_AllocateFloatVector(n2);

    if(corr_info.corr_obs) 
        free(corr_info.corr_val);

    corr_info.corr_val = 
        gomp_AllocateFloatVector(2 * n2);
    corr_info.corr_length = 2 * n2;
      
/* fill data arrays with zero padding if necessary */
    /* first array .... */
    if(corr_info.corr_vec1 == corr_info.corr_vec2) { /* autocorr */

        for(i = 0 ; i < n2 ; i++) {
            if(i < n)   data1[i] = corr_info.corr_vec1[i];
            else
                data1[i] = 0.0;
        }
        /* now handle the other one */
        for(i = 0 ; i < n2 ; i++) {
            if(i < n)  data2[i] = corr_info.corr_vec2[i];
            else
                data2[i] = 0.0;
        }

        gomp_Correl((data1-1),(data2-1),(fft-1),n2,(corr_info.corr_val - 1));

        fhelp1 = corr_info.corr_val[0];
        for(i=0 ; i<n2 ; i++) { /* normalize the result */
            corr_info.corr_val[i] = corr_info.corr_val[i] / fhelp1;
        }
    }
    else { /* cross correlation */

        for(i = 0 ; i < n2 ; i++) {
            if(i < n)   data1[i] = corr_info.corr_vec1[i];
            else
                data1[i] = 0.0;
        }

        gomp_Correl((data1-1),(data1-1),(fft-1),n2,(corr_info.corr_val - 1));

        fhelp1 = corr_info.corr_val[0];

        /* now handle the other one */
        for(i = 0 ; i < n2 ; i++) {
            if(i < n)  data2[i] = corr_info.corr_vec2[i];
            else
                data2[i] = 0.0;
        }

        gomp_Correl((data2-1),(data2-1),(fft-1),n2,(corr_info.corr_val - 1));

        fhelp2 = corr_info.corr_val[0];

        /* do now the cross correlation */

        gomp_Correl((data1-1),(data2-1),(fft-1),n2,(corr_info.corr_val - 1));

        fhelp1 = sqrt(fhelp1*fhelp2);
        for(i=0 ; i<n2 ; i++) { /* normalize the result */
            corr_info.corr_val[i] = corr_info.corr_val[i] / fhelp1;
        }

    }
    free(fft);
    free(data1);
    free(data2);

    return(0);
}
/***************************************************************************/
int   gomp_ManipulateTimeSeries(int Type , int Action , 
                              const char *Text1, const char *Text2 , const char *Text3)
/***************************************************************************/
{
    static int    Mindex;
    static int    Ihelp;
    static int    ndim;
    static float  Fhelp;
    static float *Array1;
    static float *Array2;

    if(Text1[0] == (char)NULL) {
        gomp_PrintERROR("index to list is missing");
        return(1);
    }
    Mindex = atoi(Text1);

    Ihelp  = atoi(Text2);
    Fhelp  = atof(Text2);

    switch(Type) {

    case DIST_TYPE: /* distance series */
        if(gomp_CheckTimeSeriesIndex(Type , Mindex)) {
            return(1);
        }
        ndim   = DistanceSeries[Mindex - 1].Length;
        Array1 = DistanceSeries[Mindex - 1].Value;
        if(Action == COPYNR ||
           Action == ADDNR) {
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("index to second list is missing");
                return(1);
            }
            if(gomp_CheckTimeSeriesIndex(Type , Ihelp)) {
                return(1);
            }
            Array2 = DistanceSeries[Ihelp - 1].Value;
        }
        else
            Array2 = Array1;
        break;
    case ANG_TYPE:
        if(gomp_CheckTimeSeriesIndex(Type , Mindex)) {
            return(1);
        }
        ndim   = AngleSeries[Mindex - 1].Length;
        Array1 = AngleSeries[Mindex - 1].Value;
        if(Action == COPYNR ||
           Action == ADDNR) {
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("index to second list is missing");
                return(1);
            }
            if(gomp_CheckTimeSeriesIndex(Type , Ihelp)) {
                return(1);
            }
            Array2 = AngleSeries[Ihelp - 1].Value;
        }
        else
            Array2 = Array1;
        break;
    case TORS_TYPE:
        if(gomp_CheckTimeSeriesIndex(Type , Mindex)) {
            return(1);
        }
        ndim   = TorsionSeries[Mindex - 1].Length;
        Array1 = TorsionSeries[Mindex - 1].Value;
        if(Action == COPYNR ||
           Action == ADDNR) {
            if(Text2[0] == (char)NULL) {
                gomp_PrintERROR("index to second list is missing");
                return(1);
            }
            if(gomp_CheckTimeSeriesIndex(Type , Ihelp)) {
                return(1);
            }
            Array2 = TorsionSeries[Ihelp - 1].Value;
        }
        else
            Array2 = Array1;
        break;
    }

    return(gomp_Mantime(Array1, Array2, ndim, Action, Fhelp, Ihelp));
}
/***************************************************************************/
int gomp_Mantime(float *obs_vec1, float *obs_vec2, int ndim, 
               int alt, float fnum, int inum)    
    /* manipulate time series */
/*
  float *obs_vec1;      vector with the data1  
  float *obs_vec2;      vector with the data2 
  int    ndim;          number of observations in one "set"
  int    alt;           switch for the different alternatives
  float fnum;           real number needed for some series 
  int   inum;           integer number needed for some series */
/***************************************************************************/
{

    static int i;
    static float fhelp;

    switch(alt) {

    case 1:   /* DAVErage    Q(t) = Q(t) - <Q(t)> */

        /* calculate first average */
        fhelp = gomp_VecSum(obs_vec1,ndim) / ((float) ndim);

        /* do the rest             */
        for(i = 0 ; i < ndim ; i++) obs_vec1[i] -= fhelp;

        break;

    case 2:  /* SQUAre   Q(t) = Q(t) ** 2   */

        for(i = 0 ; i < ndim ; i++) {
            fhelp = obs_vec1[i];
            obs_vec1[i] = fhelp * fhelp;
        }

        break;

    case 3:  /* COS  Q(t) = cos(Q(t))  */

        for(i = 0 ; i < ndim ; i++) obs_vec1[i] = cos(obs_vec1[i]);
        
        break;

    case 4: /* COS2  Q(t) = 3*cos(Q(t))**2 - 1   */

        for(i = 0 ; i < ndim ; i++) {
            fhelp = obs_vec1[i];
            obs_vec1[i] = 3.0 * fhelp*fhelp - 1.;
        }

        break;
      
    case 5:  /* SQRT  Q(t) = sqrt(Q(t))  */

        for(i = 0 ; i < ndim ; i++) obs_vec1[i] = sqrt(obs_vec1[i]);

        break;

    case 6: /* DINItial  Q(t) = Q(t) - Q(0)  */

        fhelp = obs_vec1[0];
        for(i = 0 ; i < ndim ; i++) obs_vec1[i] -= fhelp;

        break;

    case 7: /*  COPY nr Q(t) = Q2(t)   */

        for(i = 0 ; i < ndim ; i++) obs_vec1[i] = obs_vec2[i];

        break;

    case 8: /* ADD nr Q(t) = Q(t) + Q2(t)   */

        for(i = 0 ; i < ndim ; i++) obs_vec1[i] += obs_vec2[i];

        break;

    case 9: /* LOG Q(t) = log(Q(t))  */

        for(i = 0 ; i < ndim ; i++) obs_vec1[i] = log(obs_vec1[i]);

        break;

    case 10: /* EXP Q(t) = exp(Q(t))  */

        for(i = 0 ; i < ndim ; i++) obs_vec1[i] = exp(obs_vec1[i]);

        break;

    case 11: /* POWer real Q(t) = Q(t) ** real  */

        for(i = 0 ; i < ndim ; i++) 
            obs_vec1[i] = pow(obs_vec1[i],fnum);

        break;

    case 12: /* MULT real Q(t) = real * Q(t)  */

        for(i = 0 ; i < ndim ; i++)
            obs_vec1[i] *= fnum;

        break;

    case 13: /* DIVIde real Q(t) = Q(t) / real  */

        for(i = 0 ; i < ndim ; i++) 
            obs_vec1[i] /= fnum;

        break;

    case 14: /* SHIFt real Q(t) = Q(t) + real  */

        for(i = 0 ; i < ndim ; i++) 
            obs_vec1[i] += fnum;

        break;

    case 15: /* DMIN Q(t) = Q(t) - Q(min) */

        /* calc min first */
        fhelp = gomp_FMini(obs_vec1,ndim);
        for(i = 0 ; i < ndim ; i++)
            obs_vec1[i] -= fhelp;

        break;

    case 16: /* ABS Q(t) = abs(Q(t)) */

        for(i = 0 ; i < ndim ; i++) 
            obs_vec1[i] = RABS(obs_vec1[i]);

        break;

    case 17: /* DIVFirst Q(t) = Q(t) / Q(0)  */

        fhelp = obs_vec1[0];
        for(i = 0 ; i < ndim ; i++) obs_vec1[i] /= fhelp;

        break;

    case 18: /* DIVMaximum Q(t) = Q(t) /max(Q(t))  */

        fhelp = gomp_FMaxi(obs_vec1,ndim);
        for(i = 0 ; i < ndim ; i++)
            obs_vec1[i] /= fhelp;

        break;

    case 19: /* ZERO Q(t) = 0.0  */

        for(i = 0 ; i < ndim ; i++) obs_vec1[i] = 0.0;

        break;

    }

    return(0);
}

/***************************************************************************/
int  gomp_FillDistanceSeries()
/***************************************************************************/
{
    static int    i,j;
    static int    iLoop;
    static FILE  *File_p;
    static int    FirstFrame;
    static int    LastFrame;
    static int    StepFrame;
    static int    Frames;
    static int    Swop;
    static int    Wstr;
    static const float *x;
    static const float *y;
    static const float *z;
    static float  value,a,b,c;
    static char   OutText[BUFF_LEN];

    if(!dist_len) {
        gomp_PrintERROR("no distance variables set");
        return(1);
    }

    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintERROR("Number of frames is not defined ");
        return(1);
    }

    if(gomp_GetTrajectoryFileType() == XMOL_TRAJ      ||
       gomp_GetTrajectoryFileType() == GROMOS96A_TRAJ ||
       gomp_GetTrajectoryFileType() == TINKER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FAMBER_TRAJ    ||
       gomp_GetTrajectoryFileType() == FDL_POLY_TRAJ) {
        File_p        = fopen(gomp_GetTrajectoryFileName(),"r");
    } else {
        File_p        = fopen(gomp_GetTrajectoryFileName(),"rb");
    }

    if(File_p    == NULL) {
        sprintf(OutText,"?ERROR - can't open trajectory file : %s\n",
                gomp_GetTrajectoryFileName());
        gomp_PrintERROR(OutText);
        return(1);
    }

    (void)gomp_GetTrajectoryDisplayParams(&FirstFrame ,
                                        &LastFrame  ,
                                        &StepFrame);

    Frames = gomp_GetTrajectoryDisplayFrames();

    Wstr = 0;
    x    = gomp_GetAtomXCoordPointer(Wstr);
    y    = gomp_GetAtomYCoordPointer(Wstr);
    z    = gomp_GetAtomZCoordPointer(Wstr);

/*    SaveCoordinates  */
    (void)gomp_SaveAtomCoords(Wstr);
/*    done .........   */

    (void)gomp_DeleteDistanceSeries();

    for(i = 0 ; i < dist_len ; i = i + 2)
        (void)gomp_ExpandeDistanceSeries(Frames);

    iLoop  =  0;
    for(i  =  FirstFrame;
        i <=  LastFrame ;
        i += StepFrame) {

        if(gomp_GetOneFrame(i , File_p , TRAJ_OLD)) {
            gomp_PrintMessage("?ERROR - can't read frame in fill time series");
            fclose(File_p);
            (void)gomp_DeleteDistanceSeries();
            (void)gomp_GetSavedAtomCoords(Wstr);
            return(1);
        }

        Swop = 0;
        for(j = 0 ; j < dist_len ; j = j + 2) {
            a     = x[dist_vec[j]] - x[dist_vec[j + 1]];
            b     = y[dist_vec[j]] - y[dist_vec[j + 1]];
            c     = z[dist_vec[j]] - z[dist_vec[j + 1]];
            value = sqrt(a * a + b * b + c * c);
            DistanceSeries[Swop].Value[iLoop] = value;
            Swop++;
        }

        iLoop++;
        rewind(File_p);
    }

    NumberDistanceSeries = dist_len / 2;

    for(i = 0 ; i < NumberDistanceSeries ; i++) {
        (void)CalcVecAv(DistanceSeries[i].Value,Frames,DIST_LIST);
    }

    fclose(File_p);

    (void)gomp_GetSavedAtomCoords(Wstr);

    return(0);
}
/***************************************************************************/
int  gomp_FillAngleSeries()
/***************************************************************************/
{
    static int    i,j;
    static int    iLoop;
    static FILE  *File_p;
    static int    FirstFrame;
    static int    LastFrame;
    static int    StepFrame;
    static int    Frames;
    static int    Swop;
    static int    Wstr;
    static const float *x;
    static const float *y;
    static const float *z;
    static float  value;
    static char   OutText[BUFF_LEN];

    if(!ang_len) {
        gomp_PrintERROR("no angle variables set");
        return(1);
    }

    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintERROR("Number of frames is not defined ");
        return(1);
    }

    File_p        = fopen(gomp_GetTrajectoryFileName(),"rb");
    if(File_p    == NULL) {
        sprintf(OutText,"?ERROR - can't open trajectory file : %s\n",
                gomp_GetTrajectoryFileName());
        gomp_PrintERROR(OutText);
        return(1);
    }

    (void)gomp_GetTrajectoryDisplayParams(&FirstFrame ,
                                        &LastFrame  ,
                                        &StepFrame);

    Frames = gomp_GetTrajectoryDisplayFrames();

    Wstr = 0;
    x    = gomp_GetAtomXCoordPointer(Wstr);
    y    = gomp_GetAtomYCoordPointer(Wstr);
    z    = gomp_GetAtomZCoordPointer(Wstr);

/*    SaveCoordinates  */
    (void)gomp_SaveAtomCoords(Wstr);
/*    done .........   */

    (void)gomp_DeleteAngleSeries();
    for(i = 0 ; i < ang_len ; i = i + 3)
        (void)ExpandeAngleSeries(Frames);

    iLoop  =  0;
    for(i  =  FirstFrame;
        i <=  LastFrame ;
        i += StepFrame) {

        if(gomp_GetOneFrame(i , File_p , TRAJ_OLD)) {
            gomp_PrintMessage("?ERROR - can't read frame in fill time series");
            fclose(File_p);
            (void)gomp_DeleteAngleSeries();
            (void)gomp_GetSavedAtomCoords(Wstr);
            return(1);
        }

        Swop = 0;
        for(j = 0 ; j < ang_len ; j = j + 3) {
            gomp_BondAngle(
                x[ang_vec[j    ]] , y[ang_vec[j    ]] , z[ang_vec[j    ]],
                x[ang_vec[j + 1]] , y[ang_vec[j + 1]] , z[ang_vec[j + 1]],
                x[ang_vec[j + 2]] , y[ang_vec[j + 2]] , z[ang_vec[j + 2]],
                &value);
            AngleSeries[Swop].Value[iLoop] = 180.0 * value / M_PI;
            Swop++;
        }
        iLoop++;
        rewind(File_p);
    }

    NumberAngleSeries = ang_len / 3;

    for(i = 0 ; i < NumberAngleSeries ; i++) {
        (void)CalcVecAv(AngleSeries[i].Value,Frames,ANGLE_LIST);
    }

    fclose(File_p);

    (void)gomp_GetSavedAtomCoords(Wstr);

    return(0);
}
/***************************************************************************/
int  gomp_FillTorsionSeries()
/***************************************************************************/
{
    static int    i,j;
    static int    iLoop;
    static FILE  *File_p;
    static int    FirstFrame;
    static int    LastFrame;
    static int    StepFrame;
    static int    Frames;
    static int    Swop;
    static int    Wstr;
    static const float *x;
    static const float *y;
    static const float *z;
    static float  value;
    static char   OutText[BUFF_LEN];

    if(!tors_len) {
        gomp_PrintERROR("no torsion variables set");
        return(1);
    }

    if(gomp_GetNumberOfFrames() < 1) {
        gomp_PrintERROR("Number of frames is not defined ");
        return(1);
    }

    File_p        = fopen(gomp_GetTrajectoryFileName(),"rb");
    if(File_p    == NULL) {
        sprintf(OutText,"?ERROR - can't open trajectory file : %s\n",
                gomp_GetTrajectoryFileName());
        gomp_PrintERROR(OutText);
        return(1);
    }

    (void)gomp_GetTrajectoryDisplayParams(&FirstFrame ,
                                        &LastFrame  ,
                                        &StepFrame);

    Frames = gomp_GetTrajectoryDisplayFrames();

    Wstr = 0;
    x    = gomp_GetAtomXCoordPointer(Wstr);
    y    = gomp_GetAtomYCoordPointer(Wstr);
    z    = gomp_GetAtomZCoordPointer(Wstr);

/*    SaveCoordinates  */
    (void)gomp_SaveAtomCoords(Wstr);
/*    done .........   */

    (void)gomp_DeleteTorsionSeries();
    for(i = 0 ; i < tors_len ; i = i + 4)
        (void)ExpandTorsionSeries(Frames);

    iLoop  =  0;
    for(i  =  FirstFrame;
        i <=  LastFrame ;
        i += StepFrame) {

        if(gomp_GetOneFrame(i , File_p , TRAJ_OLD)) {
            gomp_PrintMessage("?ERROR - can't read frame in fill time series");
            fclose(File_p);
            (void)gomp_DeleteTorsionSeries();
            (void)gomp_GetSavedAtomCoords(Wstr);
            return(1);
        }

        Swop = 0;
        for(j = 0 ; j < tors_len ; j = j + 4) {
            gomp_floDihedAngle(
                x[tors_vec[j    ]] , y[tors_vec[j    ]] , z[tors_vec[j    ]] ,
                x[tors_vec[j + 1]] , y[tors_vec[j + 1]] , z[tors_vec[j + 1]] ,
                x[tors_vec[j + 2]] , y[tors_vec[j + 2]] , z[tors_vec[j + 2]] ,
                x[tors_vec[j + 3]] , y[tors_vec[j + 3]] , z[tors_vec[j + 3]] ,
                &value);
            TorsionSeries[Swop].Value[iLoop] = 180.0 * value / M_PI;
            Swop++;
        }
        iLoop++;
        rewind(File_p);
    }

    NumberTorsionSeries = tors_len / 4;

    for(i = 0 ; i < NumberTorsionSeries ; i++) {
        (void)CalcVecAv(TorsionSeries[i].Value,Frames,ANGLE_LIST);
    }

    fclose(File_p);

    (void)gomp_GetSavedAtomCoords(Wstr);

    return(0);
}
/************************************************************************/
int         gomp_WriteDistanceVector(int Which, const char *FileName)
/************************************************************************/
{
    FILE *File_p;
    int i;
    char OutText[BUFF_LEN];

    if(!NumberDistanceSeries) {
        gomp_PrintERROR("no distance series defined");
        return(1);
    }

    if(Which < 1 || Which > NumberDistanceSeries) {
        gomp_PrintERROR("wrong index to distance array");
        return(1);
    }

    File_p = fopen(FileName,"w");
    if(File_p == NULL) {
        sprintf(OutText,"can't open file '%s' for writing",FileName);
        gomp_PrintERROR(OutText);
        return(1);
    }

    Which--;

    for(i = 0 ; i < DistanceSeries[Which].Length ; i++) 
        fprintf(File_p,"%d  %f\n",(i+1),DistanceSeries[Which].Value[i]);

    fclose(File_p);

    return(0);
}

/************************************************************************/
int         gomp_WriteAngleVector(int Which, const char *FileName)
/************************************************************************/
{
    FILE *File_p;
    int i;
    char OutText[BUFF_LEN];

    if(!NumberAngleSeries) {
        gomp_PrintERROR("no angle series defined");
        return(1);
    }

    if(Which < 1 || Which > NumberAngleSeries) {
        gomp_PrintERROR("wrong index to angle array");
        return(1);
    }

    File_p = fopen(FileName,"w");
    if(File_p == NULL) {
        sprintf(OutText,"can't open file '%s' for writing",FileName);
        gomp_PrintERROR(OutText);
        return(1);
    }

    Which--;

    for(i = 0 ; i < AngleSeries[Which].Length ; i++) 
        fprintf(File_p,"%d  %f\n",(i+1),AngleSeries[Which].Value[i]);

    fclose(File_p);

    return(0);
}

/************************************************************************/
int         gomp_WriteTorsionVector(int Which, const char *FileName)
/************************************************************************/
{
    FILE *File_p;
    int i;
    char OutText[BUFF_LEN];

    if(!NumberTorsionSeries) {
        gomp_PrintERROR("no torsion series defined");
        return(1);
    }

    if(Which < 1 || Which > NumberTorsionSeries) {
        gomp_PrintERROR("wrong index to torsion array");
        return(1);
    }

    File_p = fopen(FileName,"w");
    if(File_p == NULL) {
        sprintf(OutText,"can't open file '%s' for writing",FileName);
        gomp_PrintERROR(OutText);
        return(1);
    }

    Which--;

    for(i = 0 ; i < TorsionSeries[Which].Length ; i++) 
        fprintf(File_p,"%d  %f\n",(i+1),TorsionSeries[Which].Value[i]);

    fclose(File_p);

    return(0);
}
/************************************************************************/
int         gomp_GetNumberDistanceSeries()
/************************************************************************/
{
    return(NumberDistanceSeries);
}
/************************************************************************/
int         gomp_GetNumberAngleSeries()
/************************************************************************/
{
    return(NumberAngleSeries);
}
/************************************************************************/
int         gomp_GetNumberTorsionSeries()
/************************************************************************/
{
    return(NumberTorsionSeries);
}

/************************************************************************/
int gomp_CheckTimeSeriesIndex(int Type , int Mindex)
/************************************************************************/
{
    switch(Type) {

    case DIST_TYPE:
        if(!NumberDistanceSeries) {
            gomp_PrintERROR("no distance series defined");
            return(1);       
        }
        if(Mindex < 1 || Mindex > NumberDistanceSeries) {
            gomp_PrintERROR("list index is out of range");
            return(1);
        }
        return(0);
        break;
    case ANG_TYPE:
        if(!NumberAngleSeries) {
            gomp_PrintERROR("no angle series defined");
            return(1);       
        }
        if(Mindex < 1 || Mindex > NumberAngleSeries) {
            gomp_PrintERROR("list index is out of range");
            return(1);
        }
        return(0);
        break;
    case TORS_TYPE:
        if(!NumberTorsionSeries) {
            gomp_PrintERROR("no torsion series defined");
            return(1);       
        }
        if(Mindex < 1 || Mindex > NumberTorsionSeries) {
            gomp_PrintERROR("list index is out of range");
            return(1);
        }
        return(0);
        break;
    }

    return(0);
}

/************************************************************************/
int    gomp_WriteCorrelationArray(const char *FileName)
/************************************************************************/
{
    FILE  *File_p;
    int    i;
    char   OutText[BUFF_LEN];
    const float *Value;

    if(!gomp_GetCorrelationResultVectorLength()) {
        gomp_PrintERROR("no correlation data defined");
        return(1);
    }

    File_p = fopen(FileName,"w");
    if(File_p == NULL) {
        sprintf(OutText,"can't open file '%s' for writing",FileName);
        gomp_PrintERROR(OutText);
        return(1);
    }

    Value = gomp_GetCorrelationResultVector();

    for(i = 0 ; i < gomp_GetCorrelationResultVectorLength() ; i++) 
        fprintf(File_p,"%d  %f\n",(i+1),Value[i]);

    fclose(File_p);

    return(0);
}

/************************************************************************/
int    gomp_GetDistMonitorType(int Which)
/************************************************************************/
{

    if(!dist_len) {
        gomp_PrintERROR("no distance monitor information is available (GET)");
        return(0);
    }
    if(Which < 0 || Which > ((dist_len / 2) - 1)) {
        gomp_PrintERROR("distance monitor index to type list is out of range (GET)");
        return(0);
    }

    return(dist_vec_type[Which]);
}
/************************************************************************/
int    gomp_GetDistMonitorColor(int Which ,float *RedC, float *GreenC, float *BlueC)
/************************************************************************/
{

    *RedC   = 0.0;
    *GreenC = 0.0;
    *BlueC  = 0.0;

    if(!dist_len) {
        gomp_PrintERROR("no distance monitor information is available (GET)");
        return(1);
    }
    if(Which < 0 || Which > ((dist_len / 2) - 1)) {
        gomp_PrintERROR("distance monitor index to colour list is out of range (GET)");
        return(1);
    }

    *RedC   = dist_vec_color[3 * Which    ];
    *GreenC = dist_vec_color[3 * Which + 1];
    *BlueC  = dist_vec_color[3 * Which + 2];

    return(0);
}
/************************************************************************/
int    gomp_SetDistMonitorType(int Which, int Type)
/************************************************************************/
{
    if(!dist_len) {
        gomp_PrintERROR("no distance monitor information is available (SET)");
        return(1);
    }
    if(Which < 0 || Which > ((dist_len / 2) - 1)) {
        gomp_PrintERROR("distance monitor index type to list is out of range (SET)");
        return(1);
    }
    dist_vec_type[Which] = Type;

    return(0);
}
/************************************************************************/
int    gomp_SetDistMonitorColor(int Which, float RedC, float GreenC, float BlueC)
/************************************************************************/
{

    if(!dist_len) {
        gomp_PrintERROR("no distance monitor information is available (SET)");
        return(1);
    }
    if(Which < 0 || Which > ((dist_len / 2) - 1)) {
        gomp_PrintERROR("distance monitor index to colour list is out of range (SET)");
        return(1);
    }

    dist_vec_color[3 * Which    ]  = RedC;
    dist_vec_color[3 * Which + 1]  = GreenC;
    dist_vec_color[3 * Which + 2]  = BlueC;

    return(0);
}


/* ............. */


/************************************************************************/
int    gomp_GetAngMonitorType(int Which)
/************************************************************************/
{
    if(!ang_len) {
        gomp_PrintERROR("no angle monitor information is available");
        return(0);
    }
    if(Which < 0 || Which > ((ang_len / 3) - 1)) {
        gomp_PrintERROR("angle monitor index to list is out of range");
        return(0);
    }

    return(ang_vec_type[Which]);
}
/************************************************************************/
int    gomp_GetAngMonitorColor(int Which ,float *RedC, float *GreenC, float *BlueC)
/************************************************************************/
{
    *RedC   = 0.0;
    *GreenC = 0.0;
    *BlueC  = 0.0;

    if(!ang_len) {
        gomp_PrintERROR("no angle monitor information is available");
        return(1);
    }
    if(Which < 0 || Which > ((ang_len / 3) - 1)) {
        gomp_PrintERROR("angle monitor index to list is out of range");
        return(1);
    }

    *RedC   = ang_vec_color[3 * Which    ];
    *GreenC = ang_vec_color[3 * Which + 1];
    *BlueC  = ang_vec_color[3 * Which + 2];

    return(0);
}
/************************************************************************/
int    gomp_SetAngMonitorType(int Which, int Type)
/************************************************************************/
{
    if(!ang_len) {
        gomp_PrintERROR("no angle monitor information is available");
        return(1);
    }
    if(Which < 0 || Which > ((ang_len / 3) - 1)) {
        gomp_PrintERROR("angle monitor index to list is out of range");
        return(1);
    }

    ang_vec_type[Which] = Type;

    return(0);
}
/************************************************************************/
int    gomp_SetAngMonitorColor(int Which, float RedC, float GreenC, float BlueC)
/************************************************************************/
{
    if(!ang_len) {
        gomp_PrintERROR("no angle monitor information is available");
        return(1);
    }
    if(Which < 0 || Which > ((ang_len / 3) - 1)) {
        gomp_PrintERROR("angle monitor index to list is out of range");
        return(1);
    }


    ang_vec_color[3 * Which    ]  = RedC;
    ang_vec_color[3 * Which + 1]  = GreenC;
    ang_vec_color[3 * Which + 2]  = BlueC;

    return(0);
}

/* ........... */


/************************************************************************/
int    gomp_GetTorsMonitorType(int Which)
/************************************************************************/
{
    if(!tors_len) {
        gomp_PrintERROR("no torsion monitor information is available");
        return(0);
    }
    if(Which < 0 || Which > ((tors_len / 4) - 1)) {
        gomp_PrintERROR("torsion monitor index to list is out of range");
        return(0);
    }

    return(tors_vec_type[Which]);
}
/************************************************************************/
int    gomp_GetTorsMonitorColor(int Which ,float *RedC, float *GreenC, float *BlueC)
/************************************************************************/
{

    *RedC   = 0.0;
    *GreenC = 0.0;
    *BlueC  = 0.0;

    if(!tors_len) {
        gomp_PrintERROR("no torsion monitor information is available");
        return(1);
    }
    if(Which < 0 || Which > ((tors_len / 4) - 1)) {
        gomp_PrintERROR("torsion monitor index to list is out of range");
        return(1);
    }

    *RedC   = tors_vec_color[3 * Which    ];
    *GreenC = tors_vec_color[3 * Which + 1];
    *BlueC  = tors_vec_color[3 * Which + 2];

    return(0);
}
/************************************************************************/
int    gomp_SetTorsMonitorType(int Which, int Type)
/************************************************************************/
{
    if(!tors_len) {
        gomp_PrintERROR("no torsion monitor information is available");
        return(1);
    }
    if(Which < 0 || Which > ((tors_len / 4) - 1)) {
        gomp_PrintERROR("torsion monitor index to list is out of range");
        return(1);
    }

    tors_vec_type[Which] = Type;

    return(0);
}
/************************************************************************/
int    gomp_SetTorsMonitorColor(int Which, float RedC, float GreenC, float BlueC)
/************************************************************************/
{
    if(!tors_len) {
        gomp_PrintERROR("no torsion monitor information is available");
        return(1);
    }
    if(Which < 0 || Which > ((tors_len / 4) - 1)) {
        gomp_PrintERROR("torsion monitor index to list is out of range");
        return(1);
    }

    tors_vec_color[3 * Which    ]  = RedC;
    tors_vec_color[3 * Which + 1]  = GreenC;
    tors_vec_color[3 * Which + 2]  = BlueC;

    return(0);
}
/************************************************************************/
int   gomp_CopyTimeseries2Clipboard(int Type , int Which)
/************************************************************************/
{
    static int  i;
    static char Text[BUFF_LEN];

    if(gomp_CheckTimeSeriesIndex(Type , Which))
        return(1);

    Which--;

    switch(Type) {
    case DIST_TYPE:
        {
            char *String = NULL;


            for (i = 0 ; i < DistanceSeries[Which].Length ; i++) {
#if defined(WIN32)
                sprintf(Text,"%d  %f\r\n",(i+1),DistanceSeries[Which].Value[i]);
#else
                sprintf(Text,"%d  %f\n",(i+1),DistanceSeries[Which].Value[i]);
#endif

                if(!i) {
                    String = gomp_AllocateCharVector(strlen(Text) + 1);
                    strncpy(String,Text,strlen(Text));
                    String[strlen(Text)] = (char)NULL;
                }
                else {
                    String = gomp_ReallocateCharVector(String , 
                                                    strlen(Text) + strlen(String) + 1);
                    strncat(String,Text,strlen(Text));
                    String[strlen(String)] = (char)NULL;
                }
            }
            (void)gomp_CopyText2Clipboard(String);
            free(String);
        }
        break;
    case ANG_TYPE:
        {
            char *String = NULL;


            for (i = 0 ; i < DistanceSeries[Which].Length ; i++) {
#if defined(WIN32)
                sprintf(Text,"%d  %f\r\n",(i+1),AngleSeries[Which].Value[i]);
#else
                sprintf(Text,"%d  %f\n",(i+1),AngleSeries[Which].Value[i]);
#endif

                if(!i) {
                    String = gomp_AllocateCharVector(strlen(Text) + 1);
                    strncpy(String,Text,strlen(Text));
                    String[strlen(Text)] = (char)NULL;
                }
                else {
                    String = gomp_ReallocateCharVector(String , 
                                                    strlen(Text) + strlen(String) + 1);
                    strncat(String,Text,strlen(Text));
                    String[strlen(String)] = (char)NULL;
                }
            }
            (void)gomp_CopyText2Clipboard(String);
            free(String);
        }
        break;
    case TORS_TYPE:
        {
            char *String = NULL;


            for (i = 0 ; i < DistanceSeries[Which].Length ; i++) {
#if defined(WIN32)
                sprintf(Text,"%d  %f\r\n",(i+1),TorsionSeries[Which].Value[i]);
#else
                sprintf(Text,"%d  %f\n",(i+1),TorsionSeries[Which].Value[i]);
#endif

                if(!i) {
                    String = gomp_AllocateCharVector(strlen(Text) + 1);
                    strncpy(String,Text,strlen(Text));
                    String[strlen(Text)] = (char)NULL;
                }
                else {
                    String = gomp_ReallocateCharVector(String , 
                                                    strlen(Text) + strlen(String) + 1);
                    strncat(String,Text,strlen(Text));
                    String[strlen(String)] = (char)NULL;
                }
            }
            (void)gomp_CopyText2Clipboard(String);
            free(String);
        }
        break;
    }

    return(0);
}

/* ........................... */

/************************************************************************/
int   gomp_CopyCorrelationArray2Clipboard()
/************************************************************************/
{
    static int  i;
    static char Text[BUFF_LEN];

    if(!gomp_GetCorrelationResultVectorLength()) {
        gomp_PrintERROR("no correlation array is available");
        return(1);
    }

    {
        char *String = NULL;
        const float *Value;

        Value  = gomp_GetCorrelationResultVector();

        for (i = 0 ; i < gomp_GetCorrelationResultVectorLength() ; i++) {
#if defined(WIN32)
            sprintf(Text,"%d  %f\r\n",(i+1),Value[i]);
#else
            sprintf(Text,"%d  %f\n",(i+1),Value[i]);
#endif

            if(!i) {
                String = gomp_AllocateCharVector(strlen(Text) + 1);
                strncpy(String,Text,strlen(Text));
                String[strlen(Text)] = (char)NULL;
            }
            else {
                String = gomp_ReallocateCharVector(String , 
                                                strlen(Text) + strlen(String) + 1);
                strncat(String,Text,strlen(Text));
                String[strlen(String)] = (char)NULL;
            }
        }
        (void)gomp_CopyText2Clipboard(String);
        free(String);
    }

    return(0);
}
/************************************************************************/
int    gomp_EditDistVector( int Which, int Newi, int Newj, int NewType, 
                          float NewRed, float NewGreen, float NewBlue)
/************************************************************************/
{
    int ITemp;

/* no dist list available to edit */
    if(!dist_len) return(0);
/* index outside allowed range */
    if(Which < 0 || Which >= dist_len/2) {
        gomp_PrintERROR("distance array index out of allowed range");
        return(1);
    }
    
    if(Newi >= 0)       dist_vec[2*Which]             = Newi;
    if(Newj >= 0)       dist_vec[2*Which + 1]         = Newj;
    if(NewType > 0)     dist_vec_type[Which]          = NewType;
    if(NewRed >= 0.0)   dist_vec_color[3*Which]       = NewRed;
    if(NewGreen >= 0.0) dist_vec_color[3*Which + 1]   = NewGreen;
    if(NewBlue >= 0.0)  dist_vec_color[3*Which + 2]   = NewBlue;

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),"lulCorrectDistanceListWidget");
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulCorrectDistanceListWidget'");
            return(1);
        }
    }

    return(0);
}
/************************************************************************/
int    gomp_EditAngVector( int Which, int Newi, int Newj, int Newk , int NewType, 
                         float NewRed, float NewGreen, float NewBlue)
/************************************************************************/
{
    int ITemp;

/* no dist list available to edit */
    if(!ang_len) return(0);
/* index outside allowed range */
    if(Which < 0 || Which >= ang_len/3) {
        gomp_PrintERROR("angle array index out of allowed range");
        return(1);
    }
    
    if(Newi >= 0)       ang_vec[3*Which]             = Newi;
    if(Newj >= 0)       ang_vec[3*Which + 1]         = Newj;
    if(Newk >= 0)       ang_vec[3*Which + 1]         = Newk;
    if(NewType > 0)     ang_vec_type[Which]          = NewType;
    if(NewRed >= 0.0)   ang_vec_color[3*Which]       = NewRed;
    if(NewGreen >= 0.0) ang_vec_color[3*Which + 1]   = NewGreen;
    if(NewBlue >= 0.0)  ang_vec_color[3*Which + 2]   = NewBlue;

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),"lulCorrectAngleListWidget");
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulCorrectAngleListWidget'");
            return(1);
        }
    }

    return(0);
}
/************************************************************************/
int    gomp_EditTorsVector( int Which, int Newi, int Newj, int Newk , int Newl, int NewType, 
                          float NewRed, float NewGreen, float NewBlue)
/************************************************************************/
{
    int ITemp;
/* no dist list available to edit */
    if(!tors_len) return(0);
/* index outside allowed range */
    if(Which < 0 || Which >= tors_len/4) {
        gomp_PrintERROR("torsion array index out of allowed range");
        return(1);
    }
    
    if(Newi >= 0)       tors_vec[4*Which]             = Newi;
    if(Newj >= 0)       tors_vec[4*Which + 1]         = Newj;
    if(Newk >= 0)       tors_vec[4*Which + 1]         = Newk;
    if(Newl >= 0)       tors_vec[4*Which + 1]         = Newl;
    if(NewType > 0)     tors_vec_type[Which]          = NewType;
    if(NewRed >= 0.0)   tors_vec_color[3*Which]       = NewRed;
    if(NewGreen >= 0.0) tors_vec_color[3*Which + 1]   = NewGreen;
    if(NewBlue >= 0.0)  tors_vec_color[3*Which + 2]   = NewBlue;

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),"lulCorrectTorsionListWidget");
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulCorrectTorsionListWidget'");
            return(1);
        }
    }

    return(0);
}
/************************************************************************/
int    gomp_DeleteDistanceIndex(int Which)
/************************************************************************/
{
    int i;
    int ITemp;

/* no dist list available to edit */
    if(!dist_len) return(0);
/* index outside allowed range */
    if(Which < 1 || Which > dist_len/2) {
        gomp_PrintERROR("distance array index out of allowed range");
        return(1);
    }

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),"lulCleanMonitorDistanceWidget");
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulCleanMonitorDistanceWidget'");
            return(1);
        }
    }

    if(dist_len/2 == 1) {
        gomp_ResetMonitorDistanceData();
    } else {
/* last one */
        if(Which == dist_len/2 - 1) {
        } else {
            for(i = Which ; i < dist_len/2 - 1; i++) {
                dist_vec[2 * i]           = dist_vec[2 * (i + 1)];
                dist_vec[2 * i + 1]       = dist_vec[2 * (i + 1) + 1];
                dist_vec_type[i]          = dist_vec_type[i + 1];
                dist_vec_color[3 * i]     = dist_vec_color[3 * (i + 1)];
                dist_vec_color[3 * i + 1] = dist_vec_color[3 * (i + 1) + 1];
                dist_vec_color[3 * i + 2] = dist_vec_color[3 * (i + 1) + 2];
            }
        }
        dist_len      -= 2;
        dist_vec       = gomp_ReallocateIntVector(dist_vec       , dist_len);
        dist_vec_type  = gomp_ReallocateIntVector(dist_vec_type  , (dist_len / 2));
        dist_vec_color = gomp_ReallocateFloatVector(dist_vec_color , 3 * (dist_len / 2));
    }

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),"lulCorrectDistanceListWidget");
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulCorrectDistanceListWidget'");
            return(1);
        }
    }

    return(0);
}
/************************************************************************/
int    gomp_DeleteAngleIndex(int Which)
/************************************************************************/
{
    int i;
    int ITemp;

/* no dist list available to edit */
    if(!ang_len) return(0);
/* index outside allowed range */
    if(Which < 1 || Which > ang_len/3) {
        gomp_PrintERROR("angle array index out of allowed range");
        return(1);
    }

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),"lulCleanMonitorAngleWidget");
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulCleanMonitorAngleWidget'");
            return(1);
        }
    }

    if(ang_len/3 == 1) {
        gomp_ResetMonitorAngleData();
    } else {
/* last one */
        if(Which == ang_len/3 - 1) {
        } else {
            for(i = Which ; i < ang_len/3 - 1; i++) {
                ang_vec[3 * i]           = ang_vec[3 * (i + 1)];
                ang_vec[3 * i + 1]       = ang_vec[3 * (i + 1) + 1];
                ang_vec[3 * i + 2]       = ang_vec[3 * (i + 1) + 2];
                ang_vec_type[i]          = ang_vec_type[i + 1];
                ang_vec_color[3 * i]     = ang_vec_color[3 * (i + 1)];
                ang_vec_color[3 * i + 1] = ang_vec_color[3 * (i + 1) + 1];
                ang_vec_color[3 * i + 2] = ang_vec_color[3 * (i + 1) + 2];
            }
        }
        ang_len      -= 3;
        ang_vec       = gomp_ReallocateIntVector(ang_vec       , ang_len);
        ang_vec_type  = gomp_ReallocateIntVector(ang_vec_type  , (ang_len / 2));
        ang_vec_color = gomp_ReallocateFloatVector(ang_vec_color , 3 * (ang_len / 2));
    }

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),"lulCorrectAngleListWidget");
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulCorrectAngleListWidget'");
            return(1);
        }
    }

    return(0);
}
/************************************************************************/
int    gomp_DeleteTorsionIndex(int Which)
/************************************************************************/
{
    int i;
    int ITemp;

/* no dist list available to edit */
    if(!tors_len) return(0);
/* index outside allowed range */
    if(Which < 1 || Which > tors_len/4) {
        gomp_PrintERROR("torsion array index out of allowed range");
        return(1);
    }

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),"lulCleanMonitorTorsionWidget");
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulCleanMonitorTorsionWidget'");
            return(1);
        }
    }

    if(tors_len/4 == 1) {
        gomp_ResetMonitorTorsionData();
    } else {
/* last one */
        if(Which == tors_len/4 - 1) {
        } else {
            for(i = Which ; i < tors_len/4 - 1; i++) {
                tors_vec[4 * i]           = tors_vec[4 * (i + 1)];
                tors_vec[4 * i + 1]       = tors_vec[4 * (i + 1) + 1];
                tors_vec[4 * i + 2]       = tors_vec[4 * (i + 1) + 2];
                tors_vec[4 * i + 3]       = tors_vec[4 * (i + 1) + 3];
                tors_vec_type[i]          = tors_vec_type[i + 1];
                tors_vec_color[3 * i]     = tors_vec_color[3 * (i + 1)];
                tors_vec_color[3 * i + 1] = tors_vec_color[3 * (i + 1) + 1];
                tors_vec_color[3 * i + 2] = tors_vec_color[3 * (i + 1) + 2];
            }
        }
        tors_len      -= 4;
        tors_vec       = gomp_ReallocateIntVector(tors_vec       , tors_len);
        tors_vec_type  = gomp_ReallocateIntVector(tors_vec_type  , (tors_len / 4));
        tors_vec_color = gomp_ReallocateFloatVector(tors_vec_color , 3 * (tors_len / 4));
    }

    if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
        ITemp = Tcl_GlobalEval(gomp_GetTclInterp(),"lulCorrectTorsionListWidget");
        if(ITemp != TCL_OK) {
            gomp_PrintERROR("can't execute the script 'lulCorrectTorsionListWidget'");
            return(1);
        }
    }

    return(0);
}
