// altest.cpp : Defines the entry point for the console application.
//

#include "altest.h"

namespace Plugin {
namespace Symmetry {

void FillMats(double **M1,double **M2,atomrec* ar)
{
    int i,j;
    const float *Str1x,*Str1y,*Str1z;
    const float *Str2x,*Str2y,*Str2z;

    Str1x=gom_GetAtomXCoordPointer(ar->Str1);
    Str1y=gom_GetAtomYCoordPointer(ar->Str1);
    Str1z=gom_GetAtomZCoordPointer(ar->Str1);
    Str2x=gom_GetAtomXCoordPointer(ar->Str2);
    Str2y=gom_GetAtomYCoordPointer(ar->Str2);
    Str2z=gom_GetAtomZCoordPointer(ar->Str2);
    for(j=0;j<ar->NAtom;j++)
    {
        M1[0][j]=(double)Str1x[ar->AtomList1[j]];
        M1[1][j]=(double)Str1y[ar->AtomList1[j]];
        M1[2][j]=(double)Str1z[ar->AtomList1[j]];
        M2[0][j]=(double)Str2x[ar->AtomList2[j]];
        M2[1][j]=(double)Str2y[ar->AtomList2[j]];
        M2[2][j]=(double)Str2z[ar->AtomList2[j]];
    }
    if((ar->FixedAtom1>=0)&&(ar->FixedAtom2>=0))
    {
        for(j=0;j<ar->NAtom;j++)
        {
            if(j!=ar->FixedAtom1)
            {
                for(i=0;i<3;i++)
                {
                    M1[i][j]-=M1[i][ar->FixedAtom1];
                }
            }
            if(j!=ar->FixedAtom2)
            {
                for(i=0;i<3;i++)
                {
                    M2[i][j]-=M2[i][ar->FixedAtom2];
                }
            }
        }
        for(i=0;i<3;i++)
        {
            M1[i][ar->FixedAtom1]=0.0;
            M2[i][ar->FixedAtom2]=0.0;
        }
    }
}

void SetAlignedStructs(double **R,double *t,atomrec* ar)
{
    int i,j;
    int NA;
    double x,y,z;
    float fx,fy,fz;
    float *Str2x,*Str2y,*Str2z;

    NA=gom_GetNumAtomsInMolecStruct(ar->Str1);
    Str2x=gom_GetModifiableAtomXCoordPointer(ar->Str1);
    Str2y=gom_GetModifiableAtomYCoordPointer(ar->Str1);
    Str2z=gom_GetModifiableAtomZCoordPointer(ar->Str1);
    if((ar->FixedAtom1>=0)&&(ar->FixedAtom2>=0))
    {
        fx=-Str2x[ar->FixedAtom2];
        fy=-Str2y[ar->FixedAtom2];
        fz=-Str2z[ar->FixedAtom2];
        for(i=0;i<NA;i++)
        {
            Str2x[i]+=fx;
            Str2y[i]+=fy;
            Str2z[i]+=fz;
        }
    }
    else
    {
        fx=(float)t[0];
        fy=(float)t[1];
        fz=(float)t[2];
    }
    for(j=0;j<NA;j++)
    {
        x=(R[0][0]*Str2x[j]+R[0][1]*Str2y[j]+R[0][2]*Str2z[j]);
        y=(R[1][0]*Str2x[j]+R[1][1]*Str2y[j]+R[1][2]*Str2z[j]);
        z=(R[2][0]*Str2x[j]+R[2][1]*Str2y[j]+R[2][2]*Str2z[j]);
        if((ar->FixedAtom1<0)||(ar->FixedAtom2<0))
        {
            Str2x[j]=(float)x+fx;
            Str2y[j]=(float)y+fy;
            Str2z[j]=(float)z+fz;
        }
        else
        {
            Str2x[j]=(float)x;
            Str2y[j]=(float)y;
            Str2z[j]=(float)z;
        }

    }
    if((ar->FixedAtom1>=0)&&(ar->FixedAtom2>=0))
    {
        NA=gom_GetNumAtomsInMolecStruct(ar->Str2);
        Str2x=gom_GetModifiableAtomXCoordPointer(ar->Str2);
        Str2y=gom_GetModifiableAtomYCoordPointer(ar->Str2);
        Str2z=gom_GetModifiableAtomZCoordPointer(ar->Str2);
        fx=Str2x[ar->FixedAtom1];
        fy=Str2y[ar->FixedAtom1];
        fz=Str2z[ar->FixedAtom1];
        for(j=0;j<NA;j++)
        {
            Str2x[j]-=fx;
            Str2y[j]-=fy;
            Str2z[j]-=fz;
        }
    }
}


void AlignStructs(atomrec* ar)
{
    double **A,**B;
    double **U,**V;
    double *S;
    double *muA;
    double *muB;
    double **SigmaXY;
    double **SigmaXYT;
    double sigmaA,sigmaB;
    double d[3],e[3];
    double det;
    int i,j,k;
    int NPts=ar->NAtom;

    A=(double**)malloc(3*sizeof(double*));
    B=(double**)malloc(3*sizeof(double*));
    U=(double**)malloc(3*sizeof(double*));
    V=(double**)malloc(3*sizeof(double*));
    SigmaXY=(double**)malloc(3*sizeof(double*));
    SigmaXYT=(double**)malloc(3*sizeof(double*));
    S=(double*)malloc(3*sizeof(double));
    muA=(double*)malloc(3*sizeof(double));
    muB=(double*)malloc(3*sizeof(double));
    for(i=0;i<3;i++)
    {
        A[i]=(double*)malloc(NPts*sizeof(double));
        B[i]=(double*)malloc(NPts*sizeof(double));
        U[i]=(double*)malloc(3*sizeof(double));
        V[i]=(double*)malloc(3*sizeof(double));
        SigmaXY[i]=(double*)malloc(3*sizeof(double));
        SigmaXYT[i]=(double*)malloc(3*sizeof(double));
    }
//  printf("Loading matrices\n");
    FillMats(A,B,ar);
    for(i=0;i<3;i++)
    {
        muA[i]=muB[i]=0.0;
    }
    for(i=0;i<NPts;i++)
    {
        for(j=0;j<3;j++)
        {
            muA[j]+=A[j][i];
            muB[j]+=B[j][i];
        }
    }
    for(i=0;i<3;i++)
    {
        muA[i]/=(double)NPts;
        muB[i]/=(double)NPts;
        for(j=0;j<3;j++)
            SigmaXY[i][j]=SigmaXYT[i][j]=0.0;
    }
    sigmaA=sigmaB=0.0;
    for(i=0;i<NPts;i++)
    {
        for(j=0;j<3;j++)
        {
            sigmaA+=(A[j][i]-muA[j])*(A[j][i]-muA[j]);
            sigmaB+=(B[j][i]-muB[j])*(B[j][i]-muB[j]);
            for(k=0;k<3;k++)
            {
                SigmaXY[k][j]+=(A[j][i]-muA[j])*(B[k][i]-muB[k]);
            }
        }
    }
    sigmaA/=(double)NPts;
    sigmaB/=(double)NPts;
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            SigmaXY[i][j]/=(double)NPts;
//          printf("%5.3f ",SigmaXY[i][j]);
        }
//      printf("\n");
    }
//  PauseForAction();
    BiDiag(SigmaXY,3,3,U,V);
/*  for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            printf("%6.4f ",U[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            printf("%6.4f ",V[i][j]);
        }
        printf("\n");
    }
    printf("\n");*/
    GolubKahan(SigmaXY,3,3,U,V);
/*  printf("Back from GK.\n");
    PauseForAction();*/
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            SigmaXYT[i][j]=0.0;
            for(k=0;k<3;k++)
                SigmaXYT[i][j]+=U[i][k]*V[j][k];
            printf("%6.4f ",U[i][j]);
        }
        printf("\t");
        for(j=0;j<3;j++)
        {
            printf("%6.4f ",V[i][j]);
        }
        printf("\n");
    }
    det=SigmaXYT[0][0]*(SigmaXYT[1][1]*SigmaXYT[2][2]-SigmaXYT[1][2]*SigmaXYT[2][1])
        -SigmaXYT[0][1]*(SigmaXYT[1][0]*SigmaXYT[2][2]-SigmaXYT[1][2]*SigmaXYT[2][0])
        +SigmaXYT[0][2]*(SigmaXYT[1][0]*SigmaXYT[2][1]-SigmaXYT[1][1]*SigmaXYT[2][0]);
//  PauseForAction();
    if(det<0.0)
    {
        for(i=0;i<3;i++)
        {
            SigmaXYT[i][2]=-SigmaXYT[i][2];
        }
    }
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            printf("%6.4f ",SigmaXYT[i][j]);
        }
        printf("\n");
    }
//  PauseForAction();
    if((ar->FixedAtom1>=0)&&(ar->FixedAtom2>=0))
    {
        d[0]=d[1]=d[2]=0.0;
    }
    else
    {
        for(i=0;i<3;i++)
        {
            S[i]=0.0;
            for(j=0;j<3;j++)
                S[i]+=SigmaXYT[i][j]*muA[j];
            d[i]=muB[i]-S[i];
        }
    }
//  printf("%g %g %g\n",d[0],d[1],d[2]);
    for(i=0;i<NPts;i++)
    {
        for(j=0;j<3;j++)
        {
            e[j]=0.0;
            for(k=0;k<3;k++)
                e[j]+=SigmaXYT[j][k]*A[k][i];
            e[j]+=d[j];
        }
        for(j=0;j<3;j++)
            A[j][i]=e[j];
    }
    for(i=0;i<NPts;i++)
    {
        for(j=0;j<3;j++)
        {
            printf("%6.4f\t%6.4f\n",A[j][i],B[j][i]);
        }
//      PauseForAction();
    }
    SetAlignedStructs(SigmaXYT,d,ar);
    printf("Freeing allocated space.\n");
    for(i=0;i<3;i++)
    {
        free((void*)A[i]);
        free((void*)B[i]);
        free((void*)U[i]);
        free((void*)V[i]);
        free((void*)SigmaXY[i]);
        free((void*)SigmaXYT[i]);
    }
    printf("Freeing matrix space.\n");
    free((void*)A);
    free((void*)B);
    free((void*)U);
    free((void*)V);
    free((void*)SigmaXY);
    free((void*)SigmaXYT);
    printf("Freeing vector space.\n");
    free((void*)S);
    free((void*)muA);
    free((void*)muB);
}

} // namespace Symmetry
} // namespace Plugin
