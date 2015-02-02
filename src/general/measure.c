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

#include <ctype.h>
#include <float.h>
#include "gommath.h"
#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>

#include "gomstring.h"
#include "memalloc.h"
#include "measure.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "selection.h"

#include "stdafx.h"

#define SIGN(a,b)  ((b >= 0.0) ? (a) : -(a))

static void FormatAtomInfo(char *target,int Wstr,int IAtom) {
    sprintf(target,"%s:%s(%d):%s(%d)",
            gomp_GetAtomSegName(Wstr,IAtom),
            gomp_GetAtomResName(Wstr,IAtom),gomp_GetAtomResNum1(Wstr,IAtom),
            gomp_GetAtomAtmName(Wstr,IAtom),IAtom+1);
}

static void FormatCoordinateInfo(char *target,float xc,float yc,float zc) {
    sprintf(target,"x=%f,y=%f,z=%f",xc,yc,zc);
}

/************************************************************************/
float gomp_CalculateDistance(
    const char *Text1,    /* Atom             "atom 1" */
    const char *Text2,    /* Residue                   */
    const char *Text3,    /* Segement                  */
    const char *Text4,    /* Atom             "atom 2" */
    const char *Text5,    /* Residue                   */
    const char *Text6)    /* Segement                  */
/************************************************************************/
{
    static char Info1[BUFF_LEN],Info2[BUFF_LEN];
    static char OutText[BUFF_LEN];
    int *sel_list1;
    int slong1;
    int *sel_list2;
    int slong2;
    int i,j,si,sj;
    float dist;
    int atom_max;
    int FirstFloat  = 0;
    int SecondFloat = 0;
    float TempX1,TempY1,TempZ1;
    float TempX2,TempY2,TempZ2;

    const float *Xc, *Yc, *Zc;
    float  sumx,sumy,sumz;
    const float *TranslateXYZ;

    dist = FLT_MAX;

    sumx = sumy = sumz = 0.0;
    TranslateXYZ = gomp_GetTranslateArray();
    sumx   = TranslateXYZ[0];
    sumy  = TranslateXYZ[1];
    sumz = TranslateXYZ[2];

    atom_max = gomp_GetNumAtomsInMolecStruct(0);

    Xc       = gomp_GetAtomXCoordPointer(0);
    Yc       = gomp_GetAtomYCoordPointer(0);
    Zc       = gomp_GetAtomZCoordPointer(0);

    FirstFloat  = gomp_IsStringAFloat(Text1);
    SecondFloat = gomp_IsStringAFloat(Text4);

/* neither is const float */
    if(!FirstFloat && !SecondFloat) {
        atom_max = gomp_GetNumAtomsInMolecStruct(0);

        sel_list1 = gomp_AllocateIntVector(atom_max);
        slong1    = gomp_MakeSelectionList(0,Text1,Text2,Text3,sel_list1);     
        sel_list2 = gomp_AllocateIntVector(atom_max);
        slong2    = gomp_MakeSelectionList(0,Text4,Text5,Text6,sel_list2);     

        if(slong1 > 0 && slong2 > 0) {
            for(i = 0 ; i < slong1 ; i++) {
                si = sel_list1[i];

                FormatAtomInfo(Info1,0,si);

                for(j = 0 ; j < slong2 ; j++) { 
                    sj = sel_list2[j];

                    dist = sqrt(pow((Xc[si]-Xc[sj]),2.)+
                                pow((Yc[si]-Yc[sj]),2.)+
                                pow((Zc[si]-Zc[sj]),2.));

                    FormatAtomInfo(Info2,0,sj);
                    sprintf(OutText,"Distance %s - %s:\n %f",
                            Info1,Info2,dist);
                    gomp_PrintMessage(OutText);
                }
            }
        }
        else
            gomp_PrintMessage("?ERROR - no atoms in selection list 1 or 2");

        free(sel_list1);
        free(sel_list2);
        return(dist);
    }

/* first is const float */
    if(FirstFloat && !SecondFloat) {

        TempX1 = atof(Text1);
        TempY1 = atof(Text2);
        TempZ1 = atof(Text3);

        sel_list2 = gomp_AllocateIntVector(atom_max);
        slong2    = gomp_MakeSelectionList(0,Text4,Text5,Text6,sel_list2);     

        if(slong2 > 0) {
            FormatCoordinateInfo(Info1,TempX1,TempY1,TempZ1);

            for(j = 0 ; j < slong2 ; j++) { 
                sj = sel_list2[j];

                dist = sqrt(pow((TempX1-Xc[sj] - sumx),2.)+
                            pow((TempY1-Yc[sj] - sumy),2.)+
                            pow((TempZ1-Zc[sj] - sumz),2.));

                FormatAtomInfo(Info2,0,sj);
                sprintf(OutText,"Distance %s - %s:\n %f",
                        Info1,Info2,dist);
                gomp_PrintMessage(OutText);
            }
        }
        else
            gomp_PrintMessage("?ERROR - no atoms in selection list 1 or 2");

        free(sel_list2);
        return(dist);
    }

/* second is const float */
    if(!FirstFloat && SecondFloat) {

        TempX2 = atof(Text4);
        TempY2 = atof(Text5);
        TempZ2 = atof(Text6);

        sel_list1 = gomp_AllocateIntVector(atom_max);
        slong1    = gomp_MakeSelectionList(0,Text1,Text2,Text3,sel_list1);     

        if(slong1 > 0) {
            FormatCoordinateInfo(Info2,TempX2,TempY2,TempZ2);
            
            for(j = 0 ; j < slong1 ; j++) { 
                si = sel_list1[j];

                dist = sqrt(pow((TempX2-Xc[si] - sumx),2.)+
                            pow((TempY2-Yc[si] - sumy),2.)+
                            pow((TempZ2-Zc[si] - sumz),2.));

                FormatAtomInfo(Info1,0,si);
                sprintf(OutText,"Distance %s - %s:\n %f",
                        Info1,Info2,dist);
                gomp_PrintMessage(OutText);
            }
        }
        else
            gomp_PrintMessage("?ERROR - no atoms in selection list 1 or 2");

        free(sel_list1);
        return(dist);
    }

/* both are floats      */
    if(FirstFloat && SecondFloat) {

        TempX1     = atof(Text1);
        TempX2     = atof(Text4);
        TempY1   = atof(Text2);
        TempY2   = atof(Text5);
        TempZ1 = atof(Text3);
        TempZ2 = atof(Text6);

        dist = sqrt(pow((TempX1 - TempX2), 2.) + 
                    pow((TempY1 - TempY2), 2.) + 
                    pow((TempZ1 - TempZ2), 2.));

        FormatCoordinateInfo(Info1,TempX1,TempY1,TempZ1);
        FormatCoordinateInfo(Info2,TempX2,TempY2,TempZ2);
        sprintf(OutText,"Distance %s - %s:\n %f",
                Info1,Info2,dist);
        gomp_PrintMessage(OutText);
    }
    else
        gomp_PrintMessage("?ERROR - no atoms in selection list 1 or 2");
    
    return(dist);
}



/************************************************************************/
float gomp_CalculateAngle(const char *Text1,
                          const char *Text2,
                          const char *Text3,
                          const char *Text4,
                          const char *Text5,
                          const char *Text6,
                          const char *Text7,
                          const char *Text8,
                          const char *Text9)
/************************************************************************/
{
    const char *Text[3][3];
    static char Info[3][BUFF_LEN];
    static char OutText[BUFF_LEN];
    int *sel_list[3];
    int  slong[3];
    int  loop[2];
    int  level,i,k,sk;
    double x[3],y[3],z[3];
    float angle;
    int atom_max;
    const float *XCoordP;
    const float *YCoordP;
    const float *ZCoordP;     

    Text[0][0]=Text1; Text[0][1]=Text2; Text[0][2]=Text3;
    Text[1][0]=Text4; Text[1][1]=Text5; Text[1][2]=Text6;
    Text[2][0]=Text7; Text[2][1]=Text8; Text[2][2]=Text9;

    angle = FLT_MAX;
    atom_max = gomp_GetNumAtomsInMolecStruct(0);

    XCoordP   = gomp_GetAtomXCoordPointer(0);
    YCoordP   = gomp_GetAtomYCoordPointer(0);
    ZCoordP   = gomp_GetAtomZCoordPointer(0);
     
    for( i=0; i<3; i++ ) {
        if(gomp_IsStringAFloat(Text[i][0])) {
            x[i] = atof(Text[i][0]);
            y[i] = atof(Text[i][1]);
            z[i] = atof(Text[i][2]);

            sel_list[i] = NULL;
            slong[i]    = 1;
        
            FormatCoordinateInfo(Info[i],x[i],y[i],z[i]);
        }
        else {
            sel_list[i] = gomp_AllocateIntVector(atom_max);
            slong[i]    = gomp_MakeSelectionList(0,
                                               Text[i][0],Text[i][1],Text[i][2],sel_list[i]);

            if( slong[i] <= 0 ) {
                sprintf(OutText,"?ERROR - no atoms in the selection list %d",i+1);
                gomp_PrintMessage(OutText);

                for( i--; i>=0; i-- )
                    free(sel_list[i]);
                return angle;
            }
        }
    }

    /* Two loops:                                   */
    /*     loop[0] from 0 to (slong[0]-1)           */
    /*         loop[1] from 0 to (slong[1]-1)       */
    level = 0;
    loop[level] = 0;
    while( level >= 0 ) {

        if( level < 2 ) {
            i = loop[level];

            if( sel_list[level] ) {
                x[level] = XCoordP[sel_list[level][i]];
                y[level] = YCoordP[sel_list[level][i]];
                z[level] = ZCoordP[sel_list[level][i]];

                FormatAtomInfo(Info[level],0,sel_list[level][i]);
            }

            /* Continue to inner loop. */
            level++;
            if( level < 2 ) {
                loop[level] = 0;
                continue;
            }
        }

        /* Inner most loop */
        if( sel_list[2] ) {
            for( k=0; k<slong[2]; k++ ) {
                sk = sel_list[2][k];

                (void)gomp_BondAngle(x[0],y[0],z[0],
                                   x[1],y[1],z[1],
                                   XCoordP[sk],YCoordP[sk],ZCoordP[sk],
                                   &angle);

                FormatAtomInfo(Info[2],0,sk);
                sprintf(OutText,"Angle %s - %s - %s:\n %f",
                        Info[0],Info[1],Info[2],(float)(180.*angle/M_PI));
                gomp_PrintMessage(OutText);
            }
        }
        else {
            (void)gomp_BondAngle(x[0],y[0],z[0],
                               x[1],y[1],z[1],
                               x[2],y[2],z[2],
                               &angle);

            FormatCoordinateInfo(Info[2],x[2],y[2],z[2]);
            sprintf(OutText,"Angle %s - %s - %s:\n %f",
                    Info[0],Info[1],Info[2],angle);
            gomp_PrintMessage(OutText);
        }

        /* Exit loops. */
        level--;
        while( level >= 0 ) {
            loop[level]++;
            if( loop[level] < slong[level] )
                /* We can continue on this loop level. */
                break;
            /* Still exiting */
            level--;
        }
    }

    for( i=0; i<3; i++ )
        free(sel_list[i]);

    return(180.*angle/M_PI);
}

/************************************************************************/
float gomp_CalculateTorsionAngle(const char *Text1,
                                 const char *Text2,
                                 const char *Text3,
                                 const char *Text4,
                                 const char *Text5,
                                 const char *Text6,
                                 const char *Text7,
                                 const char *Text8,
                                 const char *Text9,
                                 const char *Text10,
                                 const char *Text11,
                                 const char *Text12)
/************************************************************************/
{
    const char * Text[4][3];
    static char Info[4][BUFF_LEN];
    static char OutText[BUFF_LEN];
    int *sel_list[4];
    int  slong[4];
    int  loop[3];
    int  level,i,l,sl;
    double x[4],y[4],z[4];
    float angle;
    int atom_max;
    const float *XCoordP;
    const float *YCoordP;
    const float *ZCoordP;     
    
    Text[0][0]=Text1; Text[0][1]=Text2; Text[0][2]=Text3;
    Text[1][0]=Text4; Text[1][1]=Text5; Text[1][2]=Text6;
    Text[2][0]=Text7; Text[2][1]=Text8; Text[2][2]=Text9;
    Text[3][0]=Text10;Text[3][1]=Text11;Text[3][2]=Text12;

    angle = FLT_MAX;

    atom_max = gomp_GetNumAtomsInMolecStruct(0);

    XCoordP  = gomp_GetAtomXCoordPointer(0);
    YCoordP  = gomp_GetAtomYCoordPointer(0);
    ZCoordP  = gomp_GetAtomZCoordPointer(0);

    for( i=0; i<4; i++ ) {
        if(gomp_IsStringAFloat(Text[i][0])) {
            x[i] = atof(Text[i][0]);
            y[i] = atof(Text[i][1]);
            z[i] = atof(Text[i][2]);

            sel_list[i] = NULL;
            slong[i]    = 1;
        
            FormatCoordinateInfo(Info[i],x[i],y[i],z[i]);
        }
        else {
            sel_list[i] = gomp_AllocateIntVector(atom_max);
            slong[i]    = gomp_MakeSelectionList(0,
                                               Text[i][0],Text[i][1],Text[i][2],sel_list[i]);

            if( slong[i] <= 0 ) {
                sprintf(OutText,"?ERROR - no atoms in the selection list %d",i+1);
                gomp_PrintMessage(OutText);

                for( i--; i>=0; i-- )
                    free(sel_list[i]);
                return angle;
            }
        }
    }

    /* Three loops:                                 */
    /*     loop[0] from 0 to (slong[0]-1)           */
    /*         loop[1] from 0 to (slong[1]-1)       */
    /*             loop[2] from 0 to (slong[2]-1)   */
    level = 0;
    loop[level] = 0;
    while( level >= 0 ) {

        if( level < 3 ) {
            i = loop[level];

            if( i == 0 ) {
                /* Loop initialization */
                if( sel_list[level] ) {
                    x[level] = XCoordP[sel_list[level][i]];
                    y[level] = YCoordP[sel_list[level][i]];
                    z[level] = ZCoordP[sel_list[level][i]];

                    FormatAtomInfo(Info[level],0,sel_list[level][i]);
                }
            }

            /* Continue to inner loop. */
            level++;
            if( level < 3 ) {
                loop[level] = 0;
                continue;
            }
        }

        /* Inner most loop */
        if( sel_list[3] ) {
            for( l=0; l<slong[3]; l++ ) {
                sl = sel_list[3][l];

                (void)gomp_floDihedAngle(x[0],y[0],z[0],
                                    x[1],y[1],z[1],
                                    x[2],y[2],z[2],
                                    XCoordP[sl],YCoordP[sl],ZCoordP[sl],
                                    &angle);

                FormatAtomInfo(Info[3],0,sl);
                sprintf(OutText,"Torsion %s - %s - %s - %s:\n %f",
                        Info[0],Info[1],Info[2],Info[3],(float)(180.*angle/M_PI));
                gomp_PrintMessage(OutText);
            }
        }
        else {
            (void)gomp_floDihedAngle(x[0],y[0],z[0],
                                x[1],y[1],z[1],
                                x[2],y[2],z[2],
                                x[3],y[3],z[3],
                                &angle);

            FormatCoordinateInfo(Info[3],x[3],y[3],z[3]);
            sprintf(OutText,"Torsion %s - %s - %s - %s:\n %f",
                    Info[0],Info[1],Info[2],Info[3],angle);
            gomp_PrintMessage(OutText);
        }

        /* Exit loops. */
        level--;
        while( level >= 0 ) {
            loop[level]++;
            if( loop[level] < slong[level] )
                /* We can continue on this loop level. */
                break;
            /* Still exiting */
            level--;
        }
    }

    for( i=0; i<4; i++ )
        free(sel_list[i]);

    return(180.*angle/M_PI);
}


/*
  This routine calculates the torsion angle for the atoms ia , ib , ic , id.

  This routine is kindly provided by Dr. Florian Muller-Plathe.

  1990
*/

/***********************************************************************/
void gomp_floDihedAngle(float iax , float iay , float iaz ,
                   float ibx , float iby , float ibz ,
                   float icx , float icy , float icz ,
                   float idx , float idy , float idz ,
                   float *angle)
/***********************************************************************/
{
    static float eabx,eaby,eabz;
    static float ebcx,ebcy,ebcz;
    static float ecdx,ecdy,ecdz;
    static float abbcx,abbcy,abbcz;
    static float dccbx,dccby,dccbz;
    static float abcdx,abcdy,abcdz,signum;
    static float cosdel,xrcd;
    static float rab,rbc,rcd,xrbc,xrab,cosb,phib,xsinb,cosc,phic,xsinc;
         
/*        bond lengths and unit vectors    */

    eabx = ibx - iax;
    eaby = iby - iay;
    eabz = ibz - iaz;

    rab = sqrt (eabx * eabx + eaby * eaby + eabz * eabz);
    xrab = 1.0 / rab;
    eabx = eabx * xrab;
    eaby = eaby * xrab;
    eabz = eabz * xrab;
    ebcx = icx - ibx;
    ebcy = icy - iby;
    ebcz = icz - ibz;

    rbc = sqrt (ebcx * ebcx + ebcy * ebcy + ebcz * ebcz);
    xrbc = 1.0 / rbc;
    ebcx = ebcx * xrbc;
    ebcy = ebcy * xrbc;
    ebcz = ebcz * xrbc;
    ecdx = idx - icx;
    ecdy = idy - icy;
    ecdz = idz - icz;

    rcd = sqrt (ecdx * ecdx + ecdy * ecdy + ecdz * ecdz);
    xrcd = 1.0 / rcd;
    ecdx = ecdx * xrcd;
    ecdy = ecdy * xrcd;
    ecdz = ecdz * xrcd;

/*
  cross and dot products between unit vectors, and bond (!)
  angles
*/
    abbcx = eaby * ebcz - eabz * ebcy;
    abbcy = eabz * ebcx - eabx * ebcz;
    abbcz = eabx * ebcy - eaby * ebcx;
    cosb = - (eabx * ebcx + eaby * ebcy + eabz * ebcz);
    phib = acos(cosb);
    xsinb = 1.0 / sin(phib);
    dccbx = ecdy * ebcz - ecdz * ebcy;
    dccby = ecdz * ebcx - ecdx * ebcz;
    dccbz = ecdx * ebcy - ecdy * ebcx;
    cosc = - (ecdx * ebcx + ecdy * ebcy + ecdz * ebcz);
    phic = acos(cosc);
    xsinc = 1.0 / sin(phic);

/*        torsional angle    */

    abcdx = - ( abbcy * dccbz - abbcz * dccby );
    abcdy = - ( abbcz * dccbx - abbcx * dccbz );
    abcdz = - ( abbcx * dccby - abbcy * dccbx );
    signum = SIGN(1.0,
                  (abcdx * ebcx + abcdy * ebcy + abcdz * ebcz) );
    cosdel = - (abbcx * dccbx + abbcy * dccby + abbcz * dccbz)
        *  xsinb * xsinc;

    if(cosdel < -1.0) cosdel = -1.0;
    if(cosdel >  1.0) cosdel =  1.0;

    *angle = signum * acos (cosdel);

}


/*   Subroutine bangle          */
/*   Calculate the angle between atoms i,j and k
     using c2=a2+b2-2ab*cos(angle)
     and calculate the bond lengths i-j and j-k.

     leif laaksonen   1989 and 1995

*/

/***********************************************************************/
void gomp_BondAngle( float xi, float yi, float zi,
                   float xj, float yj, float zj,
                   float xk, float yk, float zk ,
                   float *angijk )  /* calculate bond angle */
/***********************************************************************/
{
    static double d2ij,d2jk,d2ik,tmp1,tmp2,tmp3,temp;

    tmp1 = (xi - xj);
    tmp2 = (yi - yj);
    tmp3 = (zi - zj);

    d2ij= tmp1 * tmp1 +
        tmp2 * tmp2 +
        tmp3 * tmp3;

    tmp1 = (xj - xk);
    tmp2 = (yj - yk);
    tmp3 = (zj - zk);

    d2jk= tmp1 * tmp1 +
        tmp2 * tmp2 +
        tmp3 * tmp3;

    tmp1 = (xi - xk);
    tmp2 = (yi - yk);
    tmp3 = (zi - zk);

    d2ik= tmp1 * tmp1 +
        tmp2 * tmp2 +
        tmp3 * tmp3;

    temp =  0.5*(d2ij+d2jk-d2ik)/sqrt(d2ij*d2jk);

    *angijk= ( temp >  1. ?  1. : temp);
    *angijk= ( temp < -1. ? -1. : temp);

    *angijk = acos( *angijk ) ;

}




