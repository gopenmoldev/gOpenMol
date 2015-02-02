/*

Copyright (c) 1991 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing , ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved

Enhancements 2003, 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "gomstdio.h"

#include "gommain.h"
#include "math_oper.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"

#include "stdafx.h"

typedef struct {
    float Translate[3]; /* translation in the x-, y- and z-directions */
    float ProjectionM[16];
    float ModelViewM[16];
} ObjectHandling_t;

static ObjectHandling_t ObjectHandling;

typedef struct {
    ObjectHandling_t Curr;
    ObjectHandling_t Save;
} ObjectHandlingStructure;

static const ObjectHandlingStructure DefaultObjectHandlingStructure = {
    {
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 1.0 },
        { 1.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 1.0 }
    },
    {
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 1.0 },
        { 1.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 1.0 }
    }
};

static DataVectorHandler ObjectHandlingMTHandler = {
    gomp_DataVectorCopyInitFunc, /* creator function        */
    NULL,                        /* destructor function     */
    &DefaultObjectHandlingStructure
};

static ObjectHandlingStructure *ObjectHandlingMT;

static int ObjectCenterType = 0;      /* 0 == Global, !0 == local */
static int DoViewingTransformation(int , float * , float *, float *);

/**************************************************************************/
int    gomp_SaveTranslateArray(float Xt, float Yt, float Zt)
/**************************************************************************/
{
    ObjectHandling.Translate[0] = Xt;
    ObjectHandling.Translate[1] = Yt;
    ObjectHandling.Translate[2] = Zt;

    return(0);
}

/**************************************************************************/
int    gomp_SaveProjectionMatrix(const float *Matrix)
/**************************************************************************/
{
    int i;
    for ( i = 0 ; i < 16 ; i++ )
        ObjectHandling.ProjectionM[i] = Matrix[i];

    return(0);
}
/**************************************************************************/
int    gomp_SaveModelViewMatrix(const float *Matrix)
/**************************************************************************/
{
    int i;
    for ( i = 0 ; i < 16 ; i++ )
        ObjectHandling.ModelViewM[i] = Matrix[i];

    return(0);
}

/**************************************************************************/
const float *gomp_GetTranslateArray(void)
/**************************************************************************/
{
    return(ObjectHandling.Translate);
}
/**************************************************************************/
const float *gomp_GetSavedProjectionMatrix(void)
/**************************************************************************/
{
    return(ObjectHandling.ProjectionM);
}
/**************************************************************************/
const float *gomp_GetSavedModelViewMatrix()
/**************************************************************************/
{
    return(ObjectHandling.ModelViewM);
}
/**************************************************************************/
void gomp_ScaleModelView(float scale_x , float scale_y, float scale_z)
/**************************************************************************/
{
    ObjectHandling.ModelViewM[0 ] *= scale_x;
    ObjectHandling.ModelViewM[5 ] *= scale_y;
    ObjectHandling.ModelViewM[10] *= scale_z;
}
/**************************************************************************/
float  gomp_VecSum(const float *vec , int dim)                 /* vector addition */
/**************************************************************************/
{
    register int i;
    static float fhelp;

    fhelp = 0.0;

    for(i = 0 ; i < dim ; i++) fhelp += vec[i];

    return(fhelp);
}

/***********************************************************************/
float gomp_FMaxi(const float *vec , int n) 
    /* Find maximum floating point value in array vec */
/***********************************************************************/
{

    register int i;
    static float flmax;

    flmax = -FLT_MAX;

    for(i = 0 ; i < n ; i++) {

        if(vec[i] > flmax) flmax = vec[i];
    }

    return(flmax);
}

/***********************************************************************/
float gomp_FMini(const float *vec , int n) 
    /* Find minimum floating point value in array vec */
/***********************************************************************/
{

    register int i;
    static float flmin;

    flmin = FLT_MAX;

    for(i = 0 ; i < n ; i++) {

        if(vec[i] < flmin) flmin = vec[i];
    }

    return(flmin);
}
/*
  Very simple routine for integration:

  (1) Check number of points.
  (2) Integrate using the 3/8 Simpson's rule
  (3) a) Use  Simpson's rule for the leftover points or
  b) -"-  Trapezoidal rule for the leftover points.

  Leif Laaksonen 1992

*/
/***********************************************************************/
float gomp_Trapez(const float *Fvalues ,int NumF , float DeltaF)
/*
  const float *Fvalues;    Function values 
  int    NumF;       number of function values
  float  DeltaF;     function step
*/
/***********************************************************************/
{
    int i;
    float sum;

    sum = Fvalues[0] - Fvalues[NumF - 1];

    for(i = 1 ; i < NumF - 1 ; i++) 
        sum += 2.0 * Fvalues[i];

    return(sum = sum * DeltaF / 2.0);
}

/***********************************************************************/
float gomp_Simpson(const float *Fvalues , int NumF , float DeltaF)
/*
  const float *Fvalues;    Function values
  int    NumF;       number of function values
  float  DeltaF;     function step  
*/
/***********************************************************************/
{
    int i;
    float sum;

    sum = Fvalues[0] - Fvalues[NumF - 1];

    for(i = 1 ; i < NumF - 1 ; i += 2) 
        sum += 4.0 * Fvalues[i] + 2.0 * Fvalues[i + 1];

    return(sum = sum * DeltaF / 3.0);
}
/***********************************************************************/
float gomp_Simpson38(const float *Fvalues , int NumF , float DeltaF)
/*
  const float *Fvalues;    Function values 
  int    NumF;       number of function values 
  float  DeltaF;     function step  
*/
/***********************************************************************/
{
    int i;
    int RealNumF;
    int LeftOver;
    register float sum;

    RealNumF = 3 * ((NumF - 1) / 3) + 1;

    LeftOver = NumF - RealNumF;

    sum = 3.0 * (Fvalues[0] - Fvalues[RealNumF - 1]);

    for(i = 1 ; i < RealNumF - 1 ; i += 3) {
        sum += 9.0 * Fvalues[i] + 
            9.0 * Fvalues[i + 1] + 
            6.0 * Fvalues[i + 2];
    }

    sum *=  DeltaF / 8.0;

    switch(LeftOver) {

    case 0:
        break;
    case 1:
        sum += gomp_Trapez(&Fvalues[RealNumF - 1] , 2 , DeltaF);
        break;
    case 2:
        sum += gomp_Simpson(&Fvalues[RealNumF - 1] , 3 , DeltaF);
        break;
    default:
        printf("LeftOver: %d\n",LeftOver);
        printf("***** ERROR *****\n");
        gomp_Exit(1);
    }

    return(sum);
}

/**************************************************************************/
int    gomp_SaveTranslateArrayMT(int Which , float Xt, float Yt, float Zt)
/**************************************************************************/
{
    ObjectHandlingMT[Which].Curr.Translate[0] = Xt;
    ObjectHandlingMT[Which].Curr.Translate[1] = Yt;
    ObjectHandlingMT[Which].Curr.Translate[2] = Zt;

    return(0);
}

/**************************************************************************/
int    gomp_SaveProjectionMatrixMT(int Which , const float *Matrix)
/**************************************************************************/
{
    int i;
    for ( i = 0 ; i < 16 ; i++ )
        ObjectHandlingMT[Which].Curr.ProjectionM[i] = Matrix[i];

    return(0);
}
/**************************************************************************/
int    gomp_SaveModelViewMatrixMT(int Which , const float *Matrix)
/**************************************************************************/
{
    int i;
    for ( i = 0 ; i < 16 ; i++ )
        ObjectHandlingMT[Which].Curr.ModelViewM[i] = Matrix[i];

    return(0);
}

/**************************************************************************/
const float *gomp_GetTranslateArrayMT(int Which)
/**************************************************************************/
{
    return(ObjectHandlingMT[Which].Curr.Translate);
}
/**************************************************************************/
const float *gomp_GetSavedProjectionMatrixMT(int Which)
/**************************************************************************/
{
    return(ObjectHandlingMT[Which].Curr.ProjectionM);
}
/**************************************************************************/
const float *gomp_GetSavedModelViewMatrixMT(int Which)
/**************************************************************************/
{
    return(ObjectHandlingMT[Which].Curr.ModelViewM);
}
/**************************************************************************/
void gomp_ScaleModelViewMT(int Which , float scale_x , float scale_y, float scale_z)
/**************************************************************************/
{
    ObjectHandlingMT[Which].Curr.ModelViewM[0 ] *= scale_x;
    ObjectHandlingMT[Which].Curr.ModelViewM[5 ] *= scale_y;
    ObjectHandlingMT[Which].Curr.ModelViewM[10] *= scale_z;
}
/**************************************************************************/
int    gomp_AppendSpaceModelViewMatrixMT(size_t Count)
/**************************************************************************/
{
    /* There is usually need for just one structure so there is no
     * point to allocate memory in advance. Let's use Exact
     * version.
     */
    if ( ! gomp_DataVectorCreateOrAppendExact(
             &ObjectHandlingMT,&ObjectHandlingMTHandler,Count) )
        /* Allocatation failed. Error is printed. */
        return(1);
    return(0);
}
/**************************************************************************/
int    gomp_RemoveSpaceModelViewMatrixMT(int Which, size_t Count)
/**************************************************************************/
{
    gomp_DataVectorRemove(&ObjectHandlingMT,Which,Count);

    return(0);
}
/**************************************************************************/
int    gomp_SetObjectCenterType(int Value)
/**************************************************************************/
{
    ObjectCenterType = Value;

    return(0);
}
/**************************************************************************/
int    gomp_GetObjectCenterType(void)
/**************************************************************************/
{
    return(ObjectCenterType);
}

/************************************************************************/
int DoViewingTransformation(int Which, float *x , float *y, float *z)
/************************************************************************/
{
    static float         xt,yt,zt;


    if((Which < 0) || (Which > gomp_GetNumMolecStructs() - 1)) {
        gomp_PrintERROR("structure number outside allowed range");
        return(1);
    }

    xt = ObjectHandlingMT[Which].Curr.ModelViewM[0 ] * (*x) + 
         ObjectHandlingMT[Which].Curr.ModelViewM[4 ] * (*y) + 
         ObjectHandlingMT[Which].Curr.ModelViewM[8 ] * (*z) + 
         ObjectHandlingMT[Which].Curr.ModelViewM[12];
    yt = ObjectHandlingMT[Which].Curr.ModelViewM[1 ] * (*x) + 
         ObjectHandlingMT[Which].Curr.ModelViewM[5 ] * (*y) + 
         ObjectHandlingMT[Which].Curr.ModelViewM[9 ] * (*z) + 
         ObjectHandlingMT[Which].Curr.ModelViewM[13];
    zt = ObjectHandlingMT[Which].Curr.ModelViewM[2 ] * (*x) + 
         ObjectHandlingMT[Which].Curr.ModelViewM[6 ] * (*y) + 
         ObjectHandlingMT[Which].Curr.ModelViewM[10] * (*z) + 
         ObjectHandlingMT[Which].Curr.ModelViewM[14];
    
    *x = xt;
    *y = yt;
    *z = zt;

    return(0);
}

/************************************************************************/
int gomp_DoViewingTransformationOverStructures(int What)
/************************************************************************/
{
    register float *XC,*YC,*ZC;
    static int NumStruct,i,j;
    static  float Identity[16] = {
        1.0 , 0.0 , 0.0 , 0.0,
        0.0 , 1.0 , 0.0 , 0.0,
        0.0 , 0.0 , 1.0 , 0.0,
        0.0 , 0.0 , 0.0 , 1.0};

    NumStruct = gomp_GetNumMolecStructs();   

    if(What < 0) {
/* Structure loop */
        for(i = 0 ; i < NumStruct ; i++) {

            XC = gomp_GetModifiableAtomXCoordPointer(i);
            YC = gomp_GetModifiableAtomYCoordPointer(i);
            ZC = gomp_GetModifiableAtomZCoordPointer(i);
            for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++ ) { /* main loop over atoms */

                (void)DoViewingTransformation( i, &XC[j] , &YC[j] , &ZC[j] );

            }

            (void)gomp_SaveModelViewMatrixMT(   i, Identity);

        }
    } else {

        if((What < 0) || (What >  gomp_GetNumMolecStructs() - 1)) {
            gomp_PrintERROR("structure index outside allowed range");
            return(1);
        }

        i        = What;
        XC       = gomp_GetModifiableAtomXCoordPointer(i);
        YC       = gomp_GetModifiableAtomYCoordPointer(i);
        ZC       = gomp_GetModifiableAtomZCoordPointer(i);

        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++ ) { /* main loop over atoms */

            (void)DoViewingTransformation( i, &XC[j] , &YC[j] , &ZC[j] );

        }

        (void)gomp_SaveModelViewMatrixMT(   i, Identity);
    }


    return(0);

}
/************************************************************************/
int gomp_PushModelViewingData()
/************************************************************************/
{
    static int i;

    for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {
        /*
        To copy or not to copy ProjectionM?
        ObjectHandlingMT[i].Save = ObjectHandlingMT[i].Curr;
        */
        memcpy(
            ObjectHandlingMT[i].Save.Translate,
            ObjectHandlingMT[i].Curr.Translate,
            sizeof( ObjectHandlingMT[i].Save.Translate ) );
        memcpy(
            ObjectHandlingMT[i].Save.ModelViewM,
            ObjectHandlingMT[i].Curr.ModelViewM,
            sizeof( ObjectHandlingMT[i].Save.ModelViewM ) );
    }

    return(0);
}
/************************************************************************/
int gomp_PopModelViewingData()
/************************************************************************/
{
    static int i;

    for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {
        /*
        To copy or not to copy ProjectionM?
        ObjectHandlingMT[i].Save = ObjectHandlingMT[i].Curr;
        */
        memcpy(
            ObjectHandlingMT[i].Curr.Translate,
            ObjectHandlingMT[i].Save.Translate,
            sizeof( ObjectHandlingMT[i].Curr.Translate ) );
        memcpy(
            ObjectHandlingMT[i].Curr.ModelViewM,
            ObjectHandlingMT[i].Save.ModelViewM,
            sizeof( ObjectHandlingMT[i].Curr.ModelViewM ) );
    }

    return(0);
}

/************************************************************************/
int gomp_CalcPlaneFrom3Points(const float *p1 , const float *p2 , const float *p3 , 
                            float *A  , float *B  , float *C , float *D)
/************************************************************************/
{

    float AB[3];
    float AC[3];

    AB[0] = p1[0] - p2[0];
    AB[1] = p1[1] - p2[1];
    AB[2] = p1[2] - p2[2];

    AC[0] = p1[0] - p3[0];
    AC[1] = p1[1] - p3[1];
    AC[2] = p1[2] - p3[2];

    *A =  (AB[1] * AC[2] - AC[1] * AB[2]);
    *B = -(AB[0] * AC[2] - AC[0] * AB[2]);
    *C =  (AB[0] * AC[1] - AC[0] * AB[1]);
    *D = -((*A) * p1[0] + (*B) * p1[1] + (*C) * p1[2]);

    return(0);

}
