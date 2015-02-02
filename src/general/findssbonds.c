/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2003 by:
Eero Häkkinen

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "bond.h"
#include "memalloc.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "selection.h"

#include "stdafx.h"

#define MAX_FLOAT  1.e+35    /* This max float is just a guess */
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define Rabs(a)    ( ( a ) > 0   ? (a) : -(a))
#define Fabs(a)    ( ( a ) > 0.0 ? (a) : -(a))
#define SCALE(a,b,c)     scaleO(a,b,c) /* use own scale function */


#define FIX  1        /* 1 flexible , 0 is stiff */

#define DEFSSDIST2   4.5    /* default S-S search distance to power 2 */

/***********************************************************************/
int gomp_FindSSbondsAll(float SDistance)
    /* SDistance is search distance */
/***********************************************************************/
{

    static int i,j,k,l,m,n;                 /* loop index     */
    static int *SSlist;                     /* pointer list   */
    static int *InStruct;                   /* structure pointer */
    static int Incr=10;                     /* list increment */
    static int RIndx;
    static int TempStruct;
    static float SSdist2 = DEFSSDIST2;
    static float   dist;

    static char TRes[MAX_RES_NAME_LEN];
    static char TAtm[MAX_ATM_NAME_LEN];

    static char OutText[BUFF_LEN];
    static char Tmp1[BUFF_LEN];
    static char Tmp2[BUFF_LEN];
    static char Tmp3[BUFF_LEN];
    static char Tmp4[BUFF_LEN];

    static const float *x;
    static const float *y;
    static const float *z;

    if(gomp_GetNumMolecStructs() < 1) {
        gomp_PrintMessage("?ERROR - No atoms defined ");
        return(1);
    }

    TempStruct = 0;

    if(SDistance > 0.01) SSdist2 = SDistance*SDistance;

    SSlist = gomp_AllocateIntVector(Incr);
    if(SSlist == NULL) {
        gomp_PrintMessage("?ERROR - Can't allocate space in 'FindSSbonds'");
        return(1);
    }

    InStruct = gomp_AllocateIntVector(Incr);
    if(InStruct == NULL) {
        gomp_PrintMessage("?ERROR - Can't allocate space in 'FindSSbonds'");
        return(1);
    }

/* Look for the CYS residues and SG atoms (capital letters!!!!) ... */

    RIndx = Incr;
    for(i = 0 ; i < gomp_GetNumMolecStructs() ; i++) {

        k = 0;
        if(!gomp_GetSelectedStructure(i)) continue;

        for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) {

            strncpy(TRes,gomp_GetAtomResName(i , j),MAX_RES_NAME_LEN);
/*        gomp_toller(TRes); */
            strncpy(TAtm,gomp_GetAtomAtmName(i , j),MAX_ATM_NAME_LEN);
/*        gomp_toller(TAtm); */

            if((!strncmp(TRes,"CYS",3)) && (!strncmp(TAtm,"SG",2))) {

                if(k < RIndx) { 
                    SSlist[k] = j;
                    InStruct[k] = i;
                    k++;
                }
                else {
                    RIndx = RIndx + Incr;
                    SSlist = realloc(SSlist , RIndx * sizeof(int));

                    if(SSlist == NULL) { 
                        gomp_PrintMessage("?ERROR - Can't allocate space in 'FindSSbonds'");
                        return(1);
                    }

                    InStruct = realloc(InStruct , RIndx * sizeof(int));

                    if(InStruct == NULL) { 
                        gomp_PrintMessage("?ERROR - Can't allocate space in 'FindSSbonds'");
                        return(1);
                    }

                    SSlist[k] = j;
                    InStruct[k] = i;
                    k++;
                }
            }

        }
        if(k > 0) { /* yes I did find some */
            sprintf(OutText,"SS - bonds for molecule nr: '%d'",(i+1));
            gomp_PrintMessage(OutText);

            for(j = 0 ; j < k ; j++) {

                l = SSlist[j];
                sprintf(OutText,"(%d) - '%.4s' CYS(%d) SG(%d)",
                        (j+1),gomp_GetAtomSegName(i , l),gomp_GetAtomResNum1(i , l),(l+1));
                gomp_PrintMessage(OutText);
            }

            if(k < 2) { /* only one CYS no bonds possible */
                return(0);
            }

            strncpy(Tmp3,"SG",2);
            Tmp3[2] = '\0';
            strncpy(Tmp4,"SG",2);
            Tmp4[2] = '\0';

            x = gomp_GetAtomXCoordPointer(i);
            y = gomp_GetAtomYCoordPointer(i);
            z = gomp_GetAtomZCoordPointer(i);

            for(j = 0 ; j < (k - 1) ; j++) {

                l = SSlist[j];
                strncpy(TRes,gomp_GetAtomSegName(i , l),MAX_SEG_NAME_LEN);
/*               gomp_toller(TRes); */

                for(m = j+1 ; m < k ; m++)  {

                    n = SSlist[m];

                    strncpy(TAtm,gomp_GetAtomSegName(i , n),MAX_SEG_NAME_LEN);
/*               gomp_toller(TAtm); */
               
                    if(strncmp(TRes,TAtm,MAX_SEG_NAME_LEN)) continue;

                    dist = (x[l] - x[n])*(x[l] - x[n]) +
                        (y[l] - y[n])*(y[l] - y[n]) +
                        (z[l] - z[n])*(z[l] - z[n]);

                    if(dist < SSdist2) {
                        sprintf(OutText,"Distance < %f ==> Bonding between",sqrt(SSdist2));
                        gomp_PrintMessage(OutText);
                        sprintf(OutText,"(%d) - '%.4s' CYS(%d) SG(%d)",
                                (j+1),gomp_GetAtomSegName(i , l),gomp_GetAtomResNum1(i , l),(l+1));
                        gomp_PrintMessage(OutText);
                        sprintf(OutText,"(%d) - '%.4s' CYS(%d) SG(%d)",
                                (m+1),gomp_GetAtomSegName(i , n),gomp_GetAtomResNum1(i , n),(n+1));
                        gomp_PrintMessage(OutText);
                        sprintf(Tmp1,"%d",gomp_GetAtomResNum1(i , l));
                        sprintf(Tmp2,"%d",gomp_GetAtomResNum1(i , n));
                        if(InStruct[j] == InStruct[m]) {
                            (void)gomp_EditBondI(0,1,i,TRes,Tmp1,Tmp3,TRes,Tmp2,Tmp4);
                        }
                    }
                }
            }            
        }
        else {
            gomp_PrintMessage("No SS-bonds found (is this a protein?)");
        }

    }


    free(InStruct);
    free(SSlist);

    return(0);
}

/***********************************************************************/
int gomp_FindSSbondsI(int Which , float SDistance)
    /* SDistance is search distance */
/***********************************************************************/
{

    static int i,j,k,l,m,n;                 /* loop index     */
    static int *SSlist;                     /* pointer list   */
    static int *InStruct;                   /* structure pointer */
    static int Incr=10;                     /* list increment */
    static int RIndx;
    static int TempStruct;
    static float SSdist2 = DEFSSDIST2;
    static float   dist;

    static char TRes[MAX_RES_NAME_LEN];
    static char TAtm[MAX_ATM_NAME_LEN];

    static char OutText[BUFF_LEN];
    static char Tmp1[BUFF_LEN];
    static char Tmp2[BUFF_LEN];
    static char Tmp3[BUFF_LEN];
    static char Tmp4[BUFF_LEN];

    static const float *x;
    static const float *y;
    static const float *z;

    if(gomp_GetNumMolecStructs() < 1) {
        gomp_PrintMessage("?ERROR - No atoms defined ");
        return(1);
    }

    TempStruct = 0;

    if(SDistance > 0.01) SSdist2 = SDistance*SDistance;

    SSlist = gomp_AllocateIntVector(Incr);
    if(SSlist == NULL) {
        gomp_PrintMessage("?ERROR - Can't allocate space in 'FindSSbonds'");
        return(1);
    }

    InStruct = gomp_AllocateIntVector(Incr);
    if(InStruct == NULL) {
        gomp_PrintMessage("?ERROR - Can't allocate space in 'FindSSbonds'");
        return(1);
    }

/* Look for the CYS residues and SG atoms (capital letters!!!!) ... */

    RIndx = Incr;

    if(Which < 0 || Which > (gomp_GetNumMolecStructs() - 1)) {
        gomp_PrintERROR("molecule index out of range");
        return(1);
    }

    i = Which;
    k = 0;

    for(j = 0 ; j < gomp_GetNumAtomsInMolecStruct(i) ; j++) {

        strncpy(TRes,gomp_GetAtomResName(i , j),MAX_RES_NAME_LEN);
/*        gomp_toller(TRes); */
        strncpy(TAtm,gomp_GetAtomAtmName(i , j),MAX_ATM_NAME_LEN);
/*        gomp_toller(TAtm); */

        if((!strncmp(TRes,"CYS",3)) && (!strncmp(TAtm,"SG",2))) {

            if(k < RIndx) { 
                SSlist[k] = j;
                InStruct[k] = i;
                k++;
            }
            else {
                RIndx = RIndx + Incr;
                SSlist = realloc(SSlist , RIndx * sizeof(int));

                if(SSlist == NULL) { 
                    gomp_PrintMessage("?ERROR - Can't allocate space in 'FindSSbonds'");
                    return(1);
                }

                InStruct = realloc(InStruct , RIndx * sizeof(int));

                if(InStruct == NULL) { 
                    gomp_PrintMessage("?ERROR - Can't allocate space in 'FindSSbonds'");
                    return(1);
                }

                SSlist[k] = j;
                InStruct[k] = i;
                k++;
            }
        }

    }
    if(k > 0) { /* yes I did find some */
        sprintf(OutText,"SS - bonds for molecule nr: '%d'",(i+1));
        gomp_PrintMessage(OutText);

        for(j = 0 ; j < k ; j++) {

            l = SSlist[j];
            sprintf(OutText,"(%d) - '%.4s' CYS(%d) SG(%d)",
                    (j+1),gomp_GetAtomSegName(i , l),gomp_GetAtomResNum1(i , l),(l+1));
            gomp_PrintMessage(OutText);
        }

        if(k < 2) { /* only one CYS no bonds possible */
            return(0);
        }

        strncpy(Tmp3,"SG",2);
        Tmp3[2] = '\0';
        strncpy(Tmp4,"SG",2);
        Tmp4[2] = '\0';

        x = gomp_GetAtomXCoordPointer(i);
        y = gomp_GetAtomYCoordPointer(i);
        z = gomp_GetAtomZCoordPointer(i);

        for(j = 0 ; j < (k - 1) ; j++) {

            l = SSlist[j];
            strncpy(TRes,gomp_GetAtomSegName(i , l),MAX_SEG_NAME_LEN);
/*               gomp_toller(TRes); */

            for(m = j+1 ; m < k ; m++)  {

                n = SSlist[m];

                strncpy(TAtm,gomp_GetAtomSegName(i , n),MAX_SEG_NAME_LEN);
/*               gomp_toller(TAtm); */
               
                if(strncmp(TRes,TAtm,MAX_SEG_NAME_LEN)) continue;

                dist = (x[l] - x[n])*(x[l] - x[n]) +
                    (y[l] - y[n])*(y[l] - y[n]) +
                    (z[l] - z[n])*(z[l] - z[n]);

                if(dist < SSdist2) {
                    sprintf(OutText,"Distance < %f ==> Bonding between",sqrt(SSdist2));
                    gomp_PrintMessage(OutText);
                    sprintf(OutText,"(%d) - '%.4s' CYS(%d) SG(%d)",
                            (j+1),gomp_GetAtomSegName(i , l),gomp_GetAtomResNum1(i , l),(l+1));
                    gomp_PrintMessage(OutText);
                    sprintf(OutText,"(%d) - '%.4s' CYS(%d) SG(%d)",
                            (m+1),gomp_GetAtomSegName(i , n),gomp_GetAtomResNum1(i , n),(n+1));
                    gomp_PrintMessage(OutText);
                    sprintf(Tmp1,"%d",gomp_GetAtomResNum1(i , l));
                    sprintf(Tmp2,"%d",gomp_GetAtomResNum1(i , n));
                    if(InStruct[j] == InStruct[m]) {
                        (void)gomp_EditBondI(0,1,i,TRes,Tmp1,Tmp3,TRes,Tmp2,Tmp4);
                    }
                }
            }
        }            
    }
    else {
        gomp_PrintMessage("No SS-bonds found (is this a protein?)");
    }

    free(InStruct);
    free(SSlist);

    return(0);
}
