/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
/* Contains matrix routines for operation on matrices stored as one-dimensional arrays.
Everything that didn't come from Numerical Recipes is a brute force implementation of
Golub and Van Loan, and may not work for all sets of input.*/
#include "matrices.h"

namespace Plugin {
namespace Fchk {
    
static double sqrarg;
#define SQR(a)((sqrarg=(a))==0.0?0.0:sqrarg*sqrarg)
#define SIGN(a,b) ((b)>=0.0?fabs(a):-fabs(a))
//These are all pretty self-explanatory

void MatMul(double* Ma,double *Mb,int ma,int na,int mb,int nb,double* Mp)
{
    int i,j,k;

    if(na!=mb)
        Mp[0]=-9.99e99;
    for(i=0;i<ma;i++)
    {
        for(j=0;j<nb;j++)
        {
            Mp[i*nb+j]=0.0;
            for(k=0;k<na;k++)
                Mp[i*nb+j]+=Ma[i*na+k]*Mb[k*nb+j];
        }
    }
}

void MatPr(double *Ma,int ma,int na)
{
    int i,j;

    for(i=0;i<ma;i++)
    {
        for(j=0;j<na;j++)
        {
            printf("%g ",Ma[i*na+j]);
        }
        printf("\n");
    }
}

void MatTr(double *Ma,double *Mat,int ma)
{
    int i,j;

    for(i=0;i<ma;i++)
    {
        for(j=0;j<ma;j++)
        {
            Mat[j*ma+i]=Ma[i*ma+j];
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
    if(fabs(sigma)<1.e-9)
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
        if(fabs(b)>fabs(a))
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

int Cholesky(double* C,int Order,double *G)
{
/* Calculates the Cholesky decomposition GG(H)=C of the Hermitian matrix C. Returns -1 if C is not positive 
    definite; -2 if C is not Hermitian; otherwise returns 0. If G is NULL on input, C is replaced with its 
    Cholesky decomposition, otherwise the decomposition is returned in G, leaving C unharmed.*/
    int i,j,k;

    if(G!=NULL)
    {
        for(i=0;i<Order;i++)
        {
            for(j=0;j<Order;j++)
            {
                G[i*Order+j]=C[i*Order+j];
            }
        }
    }
    else
    {
        G=C;
    }
    for(k=0;k<Order;k++)
    {
        if(G[k*Order+k]<=0.0)
            return -1;
        G[k*Order+k]=sqrt(G[k*Order+k]);
        for(i=k+1;i<Order;i++)
            G[i*Order+k]/=G[k*Order+k];
        for(j=k+1;j<Order;j++)
        {
            for(i=j;i<Order;i++)
            {
                G[i*Order+j]-=G[j*Order+k]*G[i*Order+k];
            }
        }
        for(j=0;j<k;j++)
            G[j*Order+k]=0.0;
    }
    return 0;
}

int BiDiag(double *A,int m,int n,double *U,double *V)
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
                U[i*m+j]=U[j*m+i]=0.0;
            U[i*m+i]=1.0;
        }
    }
    if(V!=NULL)
    {
        for(i=0;i<n;i++)
        {
            for(j=i+1;j<n;j++)
                V[i*n+j]=V[j*n+i]=0.0;
            V[i*n+i]=1.0;
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
                tvec[i-j]=A[i*n+j];
            }
            Householder(tvec,thvec,m-j,&beta);
            for(k=0;k<n;k++)
            {
                tem=0.0;
                for(l=j;l<m;l++)
                {
                    tem+=A[l*n+k]*thvec[l-j];
                }
                for(i=j;i<m;i++)
                {
                    A[i*n+k]-=beta*thvec[i-j]*tem;
                }
            }
            if(U!=NULL)
            {
                for(k=0;k<m;k++)
                {
                    tem=0.0;
                    for(l=j;l<m;l++)
                    {
                        tem+=U[k*m+l]*thvec[l-j];
                    }
                    for(i=j;i<m;i++)
                    {
                        U[k*m+i]-=beta*thvec[i-j]*tem;
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
                if(fabs(A[j*n+j+1])>1.e-12)
                {
                    for(i=j+1;i<n;i++)
                        tvec[i-j-1]=A[j*n+i];
                    Householder(tvec,thvec,n-j-1,&beta);
                    for(i=0;i<m;i++)
                    {
                        tem=0.0;
                        tem2=0.0;
                        for(l=j+1;l<n;l++)
                        {
                            tem+=A[i*n+l]*thvec[l-j-1];
                            if((V!=NULL)&&(i<n))
                                tem2+=V[i*n+l]*thvec[l-j-1];
                        }
                        for(k=j+1;k<n;k++)
                        {
                            A[i*n+k]-=beta*tem*thvec[k-j-1];
                            if((V!=NULL)&&(i<n))
                                V[i*n+k]-=beta*tem2*thvec[k-j-1];
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

void GivC(double *M,int m,int n,int i,int j,double *U=NULL,double *V=NULL)
{
    double y,z,g[2];
    double tau1,tau2;
    int h;

    y=M[i*n+n-1];
    z=M[j*n+n-1];
    Givens(y,z,g);
    for(h=0;h<n;h++)
    {
        tau1=M[i*n+h];
        tau2=M[j*n+h];
        M[i*n+h]=g[0]*tau1-g[1]*tau2;
        M[j*n+h]=g[1]*tau1+g[0]*tau2;
    }
    if(U!=NULL)
    {
        for(h=0;h<m;h++)
        {
            tau1=U[h*m+i];
            tau2=U[h*m+j];
            U[h*m+i]=g[0]*tau1-g[1]*tau2;
            U[h*m+j]=g[1]*tau1+g[0]*tau2;
        }
    }
}

void GivR(double *M,int m,int n,int i,int j,double *U=NULL,double *V=NULL)
{
    double y,z,g[2];
    double tau1,tau2;
    int h;

    y=M[i*n+i];
    z=M[i*n+j];
    Givens(y,z,g);
    for(h=0;h<m;h++)
    {
        tau1=M[h*n+i];
        tau2=M[h*n+j];
        M[h*n+i]=g[0]*tau1-g[1]*tau2;
        M[h*n+j]=g[1]*tau1+g[0]*tau2;
    }
    if(V!=NULL)
    {
        for(h=0;h<n;h++)
        {
            tau1=V[h*n+i];
            tau2=V[h*n+j];
            V[h*n+i]=g[0]*tau1-g[1]*tau2;
            V[h*n+j]=g[1]*tau1+g[0]*tau2;
        }
    }
}

void GolubKahan(double *B,int m,int n,double *U,double *V)
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
/*      for(i=0;i<n-1;i++)
        {
                printf("%7.5le ",B[i*n+i+1]);
        }
        printf("\n");
        PauseForAction();*/
        for(i=0;i<=n-2;i++)
        {
            if(fabs(B[i*n+i+1])<eps*(fabs(B[i*n+i])+fabs(B[(i+1)*n+i+1])))
                B[i*n+i+1]=0.0;
        }
        nq=n-1;
        while((B[(nq-1)*n+nq]==0.0)&&(nq>1))
            nq--;
        if(nq==1)
            if(B[1]==0.0)
                nq--;
        q=n-nq;
        if(q>0)
            a=1.0;
        if(q!=n)
        {
            p=0;
            nq+=1;
            while(fabs(B[p*n+p+1])==0.0)
                p++;
/*          if(p==0)
                if(fabs(B[1])==0.0)
                    p++;*/
            np=p;
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
                if(fabs(B[i*n+i])<eps)
                {
                    B[i*n+i]=0.0;
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
                    tnn=B[(nq-1)*n+nq-1]*B[(nq-1)*n+nq-1]+B[(nq-2)*n+nq-1]*B[(nq-2)*n+nq-1];
                    a=B[(nq-2)*n+nq-2]*B[(nq-2)*n+nq-2]+B[(nq-3)*n+nq-2]+B[(nq-3)*n+nq-2];
                    b=B[(nq-2)*n+nq-2]*B[(nq-2)*n+nq-1];
                    c=B[(nq-2)*n+nq-1]*B[(nq-2)*n+nq-1]+B[(nq-1)*n+nq-1]*B[(nq-1)*n+nq-1];
                }
                else
                {
                    tnn=B[(nq-1)*n+nq-1]*B[(nq-1)*n+nq-1]+B[(nq-2)*n+nq-1]*B[(nq-2)*n+nq-1];
                    a=B[(nq-2)*n+nq-2]*B[(nq-2)*n+nq-2];
                    b=B[(nq-2)*n+nq-2]*B[(nq-2)*n+nq-1];
                    c=B[(nq-2)*n+nq-1]*B[(nq-2)*n+nq-1]+B[(nq-1)*n+nq-1]*B[(nq-1)*n+nq-1];
                }
                lam1=0.5*(a+c+sqrt((a+c)*(a+c)-4.0*(a*c-b*b)));
                lam2=0.5*(a+c-sqrt((a+c)*(a+c)-4.0*(a*c-b*b)));
                if(fabs(tnn-lam1)>fabs(tnn-lam2))
                    mu=lam2;
                else
                    mu=lam1;
                y=B[np*n+np]*B[np*n+np]-mu;
                z=B[np*n+np+1]*B[np*n+np];
                for(k=np;k<nq-1;k++)
                {
                    Givens(y,z,givels);
                    for(j=0;j<m;j++)
                    {
                        tau1=B[j*n+k];
                        tau2=B[j*n+k+1];
                        B[j*n+k]=givels[0]*tau1-givels[1]*tau2;
                        B[j*n+k+1]=givels[1]*tau1+givels[0]*tau2;
                    }
                    if(V!=NULL)
                    {
                        for(j=0;j<n;j++)
                        {
                            tau1=V[j*n+k];
                            tau2=V[j*n+k+1];
                            V[j*n+k]=givels[0]*tau1-givels[1]*tau2;
                            V[j*n+k+1]=givels[1]*tau1+givels[0]*tau2;
                        }
                    }
                    y=B[k*n+k];
                    z=B[(k+1)*n+k];
                    Givens(y,z,givels);
                    for(j=0;j<n;j++)
                    {
                        tau1=B[k*n+j];
                        tau2=B[(k+1)*n+j];
                        B[k*n+j]=givels[0]*tau1-givels[1]*tau2;
                        B[(k+1)*n+j]=givels[1]*tau1+givels[0]*tau2;
                    }
                    if(U!=NULL)
                    {
                        for(j=0;j<m;j++)
                        {
                            tau1=U[j*m+k];
                            tau2=U[j*m+k+1];
                            U[j*m+k]=givels[0]*tau1-givels[1]*tau2;
                            U[j*m+k+1]=givels[1]*tau1+givels[0]*tau2;
                        }
                    }
                    if(k<nq-2)
                    {
                        y=B[k*n+k+1];
                        z=B[k*n+k+2];
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
        if(B[i*n+i]<0.0)
        {
            B[i*n+i]=-B[i*n+i];
            for(k=0;k<m;k++)
            {
                U[k*m+i]=-U[k*m+i];
            }
        }
    }
}

void MatSqrtSVD(double* Mat,double* Mat12,int n)
{
    double *U;
    double *G;
    double *Ut;
    double *St;

    G=(double*)malloc(n*n*sizeof(double));
    U=(double*)malloc(n*n*sizeof(double));
    Ut=(double*)malloc(n*n*sizeof(double));
    St=(double*)malloc(n*n*sizeof(double));
    Cholesky(Mat,n,G);
    BiDiag(G,n,n,U);
    MatPr(G,n,n);
    GolubKahan(G,n,n,U);
    MatMul(U,G,n,n,n,n,St);
    MatTr(U,Ut,n);
    MatMul(St,Ut,n,n,n,n,Mat12);
    free((void*)U);
    free((void*)G);
    free((void*)Ut);
    free((void*)St);
}

double pythag(double a,double b)
{
    double aba,abb;

    aba=fabs(a);
    abb=fabs(b);
    if(aba>abb)
        return aba*sqrt(1.0+SQR(abb/aba));
    else
        return (abb==0.0?0.0:abb*sqrt(1.0+SQR(aba/abb)));
}

void tred2(double *a,int n,double *d,double *e)
{
        int l,k,j,i;
        double scale,hh,h,g,f;

        for(i=n-1;i>=1;i--)
        {
                l=i-1;
                h=scale=0.0;
                if(l>0)
                {
                        for(k=0;k<=l;k++)
                                scale+=fabs(a[i*n+k]);
                        if(scale==0.0)
                                e[i]=a[i*n+l];
                        else
                        {
                                for(k=0;k<=l;k++)
                                {
                                        a[i*n+k]/=scale;
                                        h+=a[i*n+k]*a[i*n+k];
                                }
                                f=a[i*n+l];
                                g=(f>=0.0?-sqrt(h):sqrt(h));
                                e[i]=scale*g;
                                h-=f*g;
                                a[i*n+l]=f-g;
                                f=0.0;
                                for(j=0;j<=l;j++)
                                {
                                        a[j*n+i]=a[i*n+j]/h;
                                        g=0.0;
                                        for(k=0;k<=j;k++)
                                                g+=a[j*n+k]*a[i*n+k];
                                        for(k=j+1;k<=l;k++)
                                                g+=a[k*n+j]*a[i*n+k];
                                        e[j]=g/h;
                                        f+=e[j]*a[i*n+j];
                                }
                                hh=f/(h+h);
                                for(j=0;j<=l;j++)
                                {
                                        f=a[i*n+j];
                                        e[j]=g=e[j]-hh*f;
                                        for(k=0;k<=j;k++)
                                                a[j*n+k]-=(f*e[k]+g*a[i*n+k]);
                                }
                        }
                }
                else
                        e[i]=a[i*n+l];
                d[i]=h;
        }
        d[0]=0.0;
        e[0]=0.0;
        for(i=0;i<n;i++)
        {
                l=i-1;
                if(d[i])
                {
                        for(j=0;j<=l;j++)
                        {
                                g=0.0;
                                for(k=0;k<=l;k++)
                                        g+=a[i*n+k]*a[k*n+j];
                                for(k=0;k<=l;k++)
                                        a[k*n+j]-=g*a[k*n+i];
                        }
                }
                d[i]=a[i*n+i];
                a[i*n+i]=1.0;
                for(j=0;j<=l;j++)
                        a[j*n+i]=a[i*n+j]=0.0;
        }
}

int tqli(double *d,double *e,int n,double *z)
{
        int m,l,iter,i,k;
        double s,r,p,g,f,dd,c,b;

        for(i=1;i<n;i++)
                e[i-1]=e[i];
        e[n-1]=0.0;
        for(l=0;l<n;l++)
        {
                iter=0;
                do
                {
                        for(m=l;m<n-1;m++)
                        {
                                dd=fabs(d[m])+fabs(d[m+1]);
                                if((fabs(e[m])+dd)==dd)
                                        break;
                        }
                        if(m!=l)
                        {
                                if(iter++ == 30)
                                {
                                        printf("Error--too many iterations in tqli\n");
                                        return -1;
                                }
                                g=(d[l+1]-d[l])/(2.0*e[l]);
                                r=pythag(g,1.0);
                                g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
                                s=c=1.0;
                                p=0.0;
                                for(i=m-1;i>=l;i--)
                                {
                                        f=s*e[i];
                                        b=c*e[i];
                                        e[i+1]=(r=pythag(f,g));
                                        if(r==0.0)
                                        {
                                                d[i+1]-=p;
                                                e[m]=0.0;
                                                break;
                                        }
                                        s=f/r;
                                        c=g/r;
                                        g=d[i+1]-p;
                                        r=(d[i]-g)*s+2.0*c*b;
                                        d[i+1]=g+(p=s*r);
                                        g=c*r-b;
                                        for(k=0;k<n;k++)
                                        {
                                                f=z[k*n+i+1];
                                                z[k*n+i+1]=s*z[k*n+i]+c*f;
                                                z[k*n+i]=c*z[k*n+i]-s*f;
                                        }
                                }
                                if((r==0.0)&&(i>=l))
                                        continue;
                                d[l]-=p;
                                e[l]=g;
                                e[m]=0.0;
                        }
                }while(m!=l);
        }
/*      for(i=0;i<n;i++)
        {
            for(k=0;k<n;k++)
            {
                printf("%g ",z[i][k]);
            }
            printf("\n");
        }*/
        return 0;
}

void MatSqrtS(double* M,double* M12,int n)
{
    double *Mm;
    double *d;
    double *e;
    int i,j,k;
    double tsqr;

    Mm=(double*)malloc(n*n*sizeof(double));
    d=(double*)malloc(n*sizeof(double));
    e=(double*)malloc(n*sizeof(double));
    for(i=0;i<n*n;i++)
    {
        Mm[i]=M[i];
        M12[i]=0.0;
    }
    tred2(Mm,n,d,e);
    tqli(d,e,n,Mm);
    for(k=0;k<n;k++)
    {
        tsqr=sqrt(d[k]);
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                M12[i*n+j]+=tsqr*Mm[i*n+k]*Mm[j*n+k];
            }
        }
    }
    free((void*)Mm);
    free((void*)d);
    free((void*)e);
}

} // namespace Fchk
} // namespace Plugin
