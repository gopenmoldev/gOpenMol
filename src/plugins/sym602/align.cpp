#include <math.h>
//typedef _complex complex;
#include "align.h"

namespace Plugin {
namespace Symmetry {

#define pi 3.14159265359

void PrintMatrix(complex** C,int Order,char* format)
{
    char frmtstr[32];
    int i,j;

    sprintf(frmtstr,"(%s,%s) ",format,format);
    for(i=0;i<Order;i++)
    {
        for(j=0;j<Order;j++)
        {
            printf(frmtstr,C[i][j].Re(),C[i][j].Im());
        }
        printf("\n");
    }
    printf("\n");
}

void MatMul(complex** c1,complex **c2,int Order,complex **cp)
{
    int i,j,k;

    for(i=0;i<Order;i++)
    {
        for(j=0;j<Order;j++)
        {
            cp[i][j]=Complex(0.0,0.0);
            for(k=0;k<Order;k++)
                cp[i][j]+=c1[i][k]*c2[k][j];
        }
    }
}

void MatHMul(complex** c1h,complex **c2,int Order,complex **cp)
{
    int i,j,k;

    for(i=0;i<Order;i++)
    {
        for(j=0;j<Order;j++)
        {
            cp[i][j]=Complex(0.0,0.0);
            for(k=0;k<Order;k++)
                cp[i][j]+=c1h[k][i]^c2[k][j];
        }
    }
}

void MatMulH(complex** c1,complex **c2h,int Order,complex **cp)
{
    int i,j,k;

    for(i=0;i<Order;i++)
    {
        for(j=0;j<Order;j++)
        {
            cp[i][j]=Complex(0.0,0.0);
            for(k=0;k<Order;k++)
                cp[i][j]+=c2h[j][k]^c1[i][k];
        }
    }
}


void Householder(double* x,double *v,int n,double* beta)
{
    double sigma,mu;
    int i;

    sigma=0.0;

    for(i=1;i<n;i++)
    {
        sigma+=x[i]*x[i];
        v[i]=x[i];
    }
    if(::fabs(sigma)<1.e-9)
    {
            if(n==1)
                *beta=1.0;
            else
                *beta=0.0;
    }
    else
    {
        mu=sqrt(x[0]*x[0]+sigma);
        if(x[0]<=0.0)
            v[0]=x[0]-mu;
        else
            v[0]=-sigma/(x[0]+mu);
        *beta=(2.0*v[0]*v[0])/(sigma+v[0]*v[0]);
        sigma=v[0];
        for(i=0;i<n;i++)
            v[i]/=sigma;
    }
}

complex Householder(complex* x,complex *v,int n,complex* beta)
{
    double mag;
    complex eit;
    complex eit2;
    double normx;
    complex normvp,normvm;
    complex rv;
    int i;

    mag=x[0].magnitude();
    eit=x[0]/mag;
/*  if(eit.imaginary>0.0)
    {
        eit2.imaginary=sqrt(0.5*(1.0-eit.real));
        eit2.real=sqrt(0.5*(1.0+eit.real));
    }
    else
    {
        eit2.imaginary=-sqrt(0.5*(1.0-eit.real));
        eit2.real=sqrt(0.5*(1.0+eit.real));
    }
    normvp=eit2*eit2;*/
    normx=0.0;
    for(i=0;i<n;i++)
        normx+=(x[i]^x[i]).Re();
    normx=sqrt(normx);
    if(x[0].Re()<0.0)
        normx=-normx;
    for(i=1;i<n;i++)
    {
        v[i]=x[i];
    }
    v[0]=x[0]+normx;
    for(i=0;i<n;i++)
        v[i]/=(x[0]+normx);
    *beta=(x[0]+normx)/normx;
//  *beta=beta->Cj();
    return Complex(-1.0,0.0);
}



int BiDiag(double **A,int m,int n,double **U,double **V)
{
    double *tvec,*thvec;
    double beta,tem;
    double tem2;
    int i,j,k,l;
    int nmax;

    if(m<n)
    {
        printf("Matrix must have m>=n!.  Try transposing matrix.\n");
        return -1;
    }
    if(U!=NULL)
    {
        for(i=0;i<m;i++)
        {
            for(j=i+1;j<m;j++)
                U[i][j]=U[j][i]=0.0;
            U[i][i]=1.0;
        }
    }
    if(V!=NULL)
    {
        for(i=0;i<n;i++)
        {
            for(j=i+1;j<n;j++)
                V[i][j]=V[j][i]=0.0;
            V[i][i]=1.0;
        }
    }
    tvec=(double*)malloc(m*sizeof(double));
    thvec=(double*)malloc(m*sizeof(double));
    if(m==n)
        nmax=n-1;
    else
        nmax=n;
    for(j=0;j<nmax;j++)
    {
            for(i=j;i<m;i++)
            {
                tvec[i-j]=A[i][j];
            }
            Householder(tvec,thvec,m-j,&beta);
            for(k=0;k<n;k++)
            {
                tem=0.0;
                for(l=j;l<m;l++)
                {
                    tem+=A[l][k]*thvec[l-j];
                }
                for(i=j;i<m;i++)
                {
                    A[i][k]-=beta*thvec[i-j]*tem;
                }
            }
            if(U!=NULL)
            {
                for(k=0;k<m;k++)
                {
                    tem=0.0;
                    for(l=j;l<m;l++)
                    {
                        tem+=U[k][l]*thvec[l-j];
                    }
                    for(i=j;i<m;i++)
                    {
                        U[k][i]-=beta*thvec[i-j]*tem;
                    }
                }
            }
/*          for(i=0;i<m;i++)
            {
                for(k=0;k<m;k++)
                {
                    printf("%5.3lf ",U[i][k]);
                }
                printf("\t");
                for(k=0;k<n;k++)
                {
                    printf("%5.3lf ",A[i][k]);
                }
                printf("\n");
            }
            printf("\n");
            PauseForAction();*/
            if(j<(n-2))
            {
                if(::fabs(A[j][j+1])>1.e-12)
                {
                    for(i=j+1;i<n;i++)
                        tvec[i-j-1]=A[j][i];
                    Householder(tvec,thvec,n-j-1,&beta);
                    for(i=0;i<m;i++)
                    {
                        tem=0.0;
                        tem2=0.0;
                        for(l=j+1;l<n;l++)
                        {
                            tem+=A[i][l]*thvec[l-j-1];
                            if((V!=NULL)&&(i<n))
                                tem2+=V[i][l]*thvec[l-j-1];
                        }
                        for(k=j+1;k<n;k++)
                        {
                            A[i][k]-=beta*tem*thvec[k-j-1];
                            if((V!=NULL)&&(i<n))
                                V[i][k]-=beta*tem2*thvec[k-j-1];
                        }
                    }
/*                  for(i=0;i<m;i++)
                    {
                        for(k=0;k<n;k++)
                        {
                            printf("%5.3lf\t",A[i][k]);
                        }
                        printf("\n");
                    }
                    printf("\n");*/
                }
            }
        //PauseForAction();
    }
    return 0;
}

void Givens(double a,double b,double trigs[])
{
    double tau;

    if(b==0)
    {
        trigs[0]=1.0;
        trigs[1]=0.0;
    }
    else
    {
        if(::fabs(b)>::fabs(a))
        {
            tau=-a/b;
            trigs[1]=1.0/sqrt(1.0+tau*tau);
            trigs[0]=trigs[1]*tau;
        }
        else
        {
            tau=-b/a;
            trigs[0]=1.0/sqrt(1.0+tau*tau);
            trigs[1]=trigs[0]*tau;
        }
    }
}

void Givens(complex a,complex b,complex trigs[])
{
    complex tau;

    if(a==Complex(0.0,0.0))
    {
        trigs[0]=Complex(0.0,0.0);
        trigs[1]=Complex(1.0,0.0);
    }
    else
    {
        tau=b/a;
        trigs[0]=1.0/sqrt(1.0+(tau^tau).Re());
        trigs[1]=-tau*trigs[0];
    }
}

void GivC(double **M,int m,int n,int i,int j,double **U=NULL,double **V=NULL)
{
    double y,z,g[2];
    double tau1,tau2;
    int h;

    y=M[i][n-1];
    z=M[j][n-1];
    Givens(y,z,g);
    for(h=0;h<n;h++)
    {
        tau1=M[i][h];
        tau2=M[j][h];
        M[i][h]=g[0]*tau1-g[1]*tau2;
        M[j][h]=g[1]*tau1+g[0]*tau2;
    }
    if(U!=NULL)
    {
        for(h=0;h<m;h++)
        {
            tau1=U[h][i];
            tau2=U[h][j];
            U[h][i]=g[0]*tau1-g[1]*tau2;
            U[h][j]=g[1]*tau1+g[0]*tau2;
        }
    }
}

void GivR(double **M,int m,int n,int i,int j,double **U=NULL,double **V=NULL)
{
    double y,z,g[2];
    double tau1,tau2;
    int h;

    y=M[i][i];
    z=M[i][j];
    Givens(y,z,g);
    for(h=0;h<m;h++)
    {
        tau1=M[h][i];
        tau2=M[h][j];
        M[h][i]=g[0]*tau1-g[1]*tau2;
        M[h][j]=g[1]*tau1+g[0]*tau2;
    }
    if(V!=NULL)
    {
        for(h=0;h<n;h++)
        {
            tau1=V[h][i];
            tau2=V[h][j];
            V[h][i]=g[0]*tau1-g[1]*tau2;
            V[h][j]=g[1]*tau1+g[0]*tau2;
        }
    }
}

void GolubKahan(double **B,int m,int n,double **U,double **V)
{
/* B should be upper bidiagonal.  */
    int i,j,k;
    int p,q;
    int np,nq;
    double lam1,lam2,mu;
    double tau1,tau2;
    double givels[2];
    double eps=1.e-17;
    double tnn,a,b,c;
    double y,z;

    q=0;
    while(q!=n)
    {
/*      for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                printf("%7.5le ",B[i][j]);
            }
            printf("\n");
        }
        printf("\n");*/
        for(i=0;i<=n-2;i++)
        {
            if(::fabs(B[i][i+1])<eps*(::fabs(B[i][i])+::fabs(B[i+1][i+1])))
                B[i][i+1]=0.0;
        }
        nq=n-1;
        while((B[nq-1][nq]==0.0)&&(nq>1))
            nq--;
        if(nq==1)
            if(B[0][1]==0.0)
                nq--;
        q=n-nq;
        if(q>0)
            a=1.0;
        if(q!=n)
        {
            np=nq-1;
            nq+=1;
            while((::fabs(B[np][np+1])>0.0)&&(np>0))
                np--;
            if(np==0)
                if(::fabs(B[0][1])==0.0)
                    np++;
            p=n-np;
        }
        else
        {
            p=0;
            np=0;
            nq=0;
        }
//      printf("%i %i\n",p,q);
//      PauseForAction();
        if(q<n)
        {
            k=-1;
            for(i=np;i<nq;i++)
            {
                if(::fabs(B[i][i])<eps)
                {
                    B[i][i]=0.0;
                    k=i;
                }
            }
            if(k!=-1)
            {
                if(k==n-1)
                {
                    for(j=n-2;j>=0;j--)
                        GivC(B,m,n,k-1,j);
                }
                else
                {
                    for(j=k+1;j<n;j++)
                        GivR(B,m,n,k,j);
                }
            }
            else
            {
                if(nq-np>2)
                {
                    tnn=B[nq-1][nq-1]*B[nq-1][nq-1]+B[nq-2][nq-1]*B[nq-2][nq-1];
                    a=B[nq-2][nq-2]*B[nq-2][nq-2]+B[nq-3][nq-2]+B[nq-3][nq-2];
                    b=B[nq-2][nq-2]*B[nq-2][nq-1];
                    c=B[nq-2][nq-1]*B[nq-2][nq-1]+B[nq-1][nq-1]*B[nq-1][nq-1];
                }
                else
                {
                    tnn=B[nq-1][nq-1]*B[nq-1][nq-1]+B[nq-2][nq-1]*B[nq-2][nq-1];
                    a=B[nq-2][nq-2]*B[nq-2][nq-2];
                    b=B[nq-2][nq-2]*B[nq-2][nq-1];
                    c=B[nq-2][nq-1]*B[nq-2][nq-1]+B[nq-1][nq-1]*B[nq-1][nq-1];
                }
                lam1=0.5*(a+c+sqrt((a+c)*(a+c)-4.0*(a*c-b*b)));
                lam2=0.5*(a+c-sqrt((a+c)*(a+c)-4.0*(a*c-b*b)));
                if(::fabs(tnn-lam1)>::fabs(tnn-lam2))
                    mu=lam2;
                else
                    mu=lam1;
                y=B[np][np]*B[np][np]-mu;
                z=B[np][np+1]*B[np][np];
                for(k=np;k<nq-1;k++)
                {
                    Givens(y,z,givels);
                    for(j=0;j<m;j++)
                    {
                        tau1=B[j][k];
                        tau2=B[j][k+1];
                        B[j][k]=givels[0]*tau1-givels[1]*tau2;
                        B[j][k+1]=givels[1]*tau1+givels[0]*tau2;
                    }
                    if(V!=NULL)
                    {
                        for(j=0;j<n;j++)
                        {
                            tau1=V[j][k];
                            tau2=V[j][k+1];
                            V[j][k]=givels[0]*tau1-givels[1]*tau2;
                            V[j][k+1]=givels[1]*tau1+givels[0]*tau2;
                        }
                    }
                    y=B[k][k];
                    z=B[k+1][k];
                    Givens(y,z,givels);
                    for(j=0;j<n;j++)
                    {
                        tau1=B[k][j];
                        tau2=B[k+1][j];
                        B[k][j]=givels[0]*tau1-givels[1]*tau2;
                        B[k+1][j]=givels[1]*tau1+givels[0]*tau2;
                    }
                    if(U!=NULL)
                    {
                        for(j=0;j<m;j++)
                        {
                            tau1=U[j][k];
                            tau2=U[j][k+1];
                            U[j][k]=givels[0]*tau1-givels[1]*tau2;
                            U[j][k+1]=givels[1]*tau1+givels[0]*tau2;
                        }
                    }
                    if(k<nq-2)
                    {
                        y=B[k][k+1];
                        z=B[k][k+2];
                    }
                }
            }
        }
    }
/* There is one little thing to check here.  We want to force all the elements of the diagonal
matrix (now stored in B) to be >=0.  We thus loop over the diagonal elements of B and invert the
corresponding column of U (so that UFDVT=B; U<-UF) */
    for(i=0;i<n;i++)
    {
        if(B[i][i]<0.0)
        {
            B[i][i]=-B[i][i];
            for(k=0;k<m;k++)
            {
                U[k][i]=-U[k][i];
            }
        }
    }
}

int SVD(double **A,int m,int n,double **U,double **V)
{
    double **Af;
    double *d;
    double *e;
    int i;
    int j;

    Af=(double**)malloc((m+n)*sizeof(double*));
    d=(double*)malloc((m+n)*sizeof(double));
    e=(double*)malloc((m+n)*sizeof(double));
    for(i=0;i<m+n;i++)
        Af[i]=(double*)malloc((m+n)*sizeof(double));
    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
            Af[i][j]=0.0;
        for(j=0;j<n;j++)
            Af[i][j+m]=A[i][j];
        for(j=0;j<n;j++)
            Af[i+m][j]=A[j][i];
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            Af[i+m][j+m]=0.0;
        }
    }
    for(i=0;i<n+m;i++)
    {
        for(j=0;j<n+m;j++)
        {
            printf("%6.4f ",Af[i][j]);
        }
        printf("\n");
    }
    tred2(Af,m+n,d,e);
    tqli(d,e,m+n,Af);
    for(i=0;i<n;i++)
    {
        for(j=0;j<n+m;j++)
        {
            if(d[j]>0.0)
                printf("%6.4f ",Af[i][j]*0.7071);
            else
                printf("%6.4f ",-Af[i][j]*0.7071);
        }
        printf("\n");
    }
    for(i=0;i<m;i++)
    {
        for(j=0;j<n+m;j++)
        {
                printf("%6.4f ",Af[i+n][j]*0.7071);
        }
        printf("\n");
    }
    for(i=0;i<n+m;i++)
        free((void*)Af[i]);
    free((void*)d);
    free((void*)e);
    free((void*)Af);
    return 0;
}

void UpperTriang(complex** A,int Order,complex** Q)
{
    int i,j,k,l;
    complex* tvec;
    complex* thvec;
    complex tem;
    complex beta;

    tvec=(complex*)malloc(Order*sizeof(complex));
    thvec=(complex*)malloc(Order*sizeof(complex));
    if(Q!=NULL)
    {
        for(i=0;i<Order;i++)
        {
            Q[i][i]=Complex(1.0,0.0);
            for(j=i+1;j<Order;j++)
            {
                Q[i][j]=Complex(0.0,0.0);
                Q[j][i]=Complex(0.0,0.0);
            }
        }
    }
    for(j=0;j<Order;j++)
    {
        for(i=j;i<Order;i++)
            tvec[i-j]=A[i][j];
        Householder(tvec,thvec,Order-j,&beta);
        for(k=0;k<Order;k++)
        {
            tem=0.0;
            for(l=j;l<Order;l++)
            {
                tem+=A[l][k]^thvec[l-j];
            }
            for(i=j;i<Order;i++)
            {
                A[i][k]-=tem^thvec[i-j]*beta;
            }
        }
        if(Q!=NULL)
        {
            for(k=0;k<Order;k++)
            {
                tem=0.0;
                for(l=j;l<Order;l++)
                {
                    tem+=Q[k][l]*thvec[l-j];
                }
                for(i=j;i<Order;i++)
                {
                    Q[k][i]-=thvec[i-j]^tem*beta;
                }
            }
        }
    }
}

void TriDiag(complex** A,int Order,complex** Q)
{
    int i,j,k,l;
    complex* tvec;
    complex* thvec;
    complex* p;
    complex* w;
    complex q;
    complex tem;
    complex sgn;
    complex beta;

    tvec=(complex*)malloc(Order*sizeof(complex));
    thvec=(complex*)malloc(Order*sizeof(complex));
    p=(complex*)malloc(Order*sizeof(complex));
    w=(complex*)malloc(Order*sizeof(complex));
    if(Q!=NULL)
    {
        for(i=0;i<Order;i++)
        {
            Q[i][i]=Complex(1.0,0.0);
            for(j=i+1;j<Order;j++)
            {
                Q[i][j]=Complex(0.0,0.0);
                Q[j][i]=Complex(0.0,0.0);
            }
        }
    }
    for(k=0;k<Order-1;k++)
    {
        for(i=k+1;i<Order;i++)
            tvec[i-k-1]=A[i][k];
        sgn=Householder(tvec,thvec,Order-k-1,&beta);
        for(i=k+1;i<Order;i++)
        {
//          p[i-k-1]=Complex(0.0,0.0);
            w[i-k-1]=Complex(0.0,0.0);
            for(j=k+1;j<Order;j++)
            {
//              p[i-k-1]+=A[j][i]^thvec[j-k-1];
                w[i-k-1]+=A[i][j]*thvec[j-k-1];
            }
//          p[i-k-1]*=beta;
            w[i-k-1]*=beta;
        }
        tem=Complex(0.0,0.0);
        q=Complex(0.0,0.0);
        for(i=k+1;i<Order;i++)
        {
            q+=(w[i-k-1]^thvec[i-k-1]);
            tem+=A[i][k]^A[i][k];
        }
        q*=beta;
        for(i=k+1;i<Order;i++)
        {
            if(i==k+1)
            {
                A[i][k]=sgn*Complex(sqrt(tem.Re()),0.0);
            }
            else
            {
                A[i][k]=Complex(0.0,0.0);
            }
            A[k][i]=A[i][k];
            for(j=i;j<Order;j++)
            {
                A[i][j]-=(w[j-k-1]^thvec[i-k-1])+(thvec[j-k-1]^w[i-k-1])-(q*(thvec[j-k-1]^thvec[i-k-1]));
                if(i!=j)
                    A[j][i]=A[i][j].Cj();
            }
        }
        if(Q!=NULL)
        {
            for(j=0;j<Order;j++)
            {
                tem=0.0;
                for(l=k+1;l<Order;l++)
                {
                    tem+=Q[l][j]^thvec[l-k-1];
                }
                for(i=k+1;i<Order;i++)
                {
                    Q[i][j]-=beta*(tem^thvec[i-k-1]);
                }
            }
        }
    }
}
    
void TQRImp(complex** C,int Order,complex** Q)
{
    int i,j,k;
    int p,q;
    int np,nq;
    complex lam1,lam2,mu;
    complex tau1,tau2;
    complex cs[2];
    complex d,x,z;
    double eps=1.e-15;

    q=0;
    while(q!=Order)
    {
        for(i=0;i<=Order-2;i++)
        {
            if(fabs(C[i][i+1])<eps*(fabs(C[i][i])+fabs(C[i+1][i+1])))
            {
                C[i][i+1]=Complex(0.0,0.0);
                C[i+1][i]=Complex(0.0,0.0);
            }
        }
        nq=Order-1;
        while((C[nq-1][nq]==Complex(0.0,0.0))&&(nq>1))
            nq--;
        if(nq==1)
            if(C[0][1]==Complex(0.0,0.0))
                nq--;
        q=Order-nq;
        if(q!=Order)
        {
            np=nq-1;
            nq+=1;
            while((fabs(C[np][np+1])>0.0)&&(np>0))
                np--;
            if(np==0)
                if(fabs(C[0][1])==0.0)
                    np++;
            p=Order-np;
        }
        else
        {
            p=0;
            np=0;
            nq=0;
        }
        if(q!=Order)
        {
            d=(C[nq-2][nq-2]-C[nq-1][nq-1])*0.5;
            mu=C[nq-1][nq-1];
            x=C[np][np]-mu;
            z=C[np+1][np];
            for(k=np;k<nq-1;k++)
            {
                Givens(x,z,cs);
                for(i=0;i<Order;i++)
                {
                    tau1=C[k][i];
                    tau2=C[k+1][i];
                    C[k][i]=(cs[0]*tau1)-(cs[1]^tau2);
                    C[k+1][i]=(cs[1]*tau1)+(cs[0]*tau2);
                }
                for(i=0;i<Order;i++)
                {
                    tau1=C[i][k];
                    tau2=C[i][k+1];
                    C[i][k]=(cs[0]*tau1)-(cs[1]*tau2);
                    C[i][k+1]=(cs[1]^tau1)+(cs[0]*tau2);
                }
            }
            if(k<nq-2)
            {
                x=C[k+1][k];
                z=C[k+2][k];
            }
        }
    }
}

} // namespace Symmetry
} // namespace Plugin
