// gullick.cpp : Defines the entry point for the console application.
//
/* Copyright 2002-2003 Kevin J. Boyd and the University of New Orleans.  
Permission granted to distribute and modify this code for personal, educational, and
research use. */

/*Implements Gullickson's method for constrained or weighted nonlinear least squares.
    We'll see if this sucker works. */
#include "gullick.h"
#define matacc(q1,q2,q3,q4) (q1)[(q2)*(q4)+(q3)]
#define mtol 1.e30           


/*maintenance routines. Mostly for debugging.*/
void MatMul(double* Ma,double *Mb,int ma,int na,int mb,int nb,double* Mp)
{
    int i,j,k;

    if(na!=mb)
        Mp[0]=-9.99e99;
    for(i=0;i<ma;i++)
    {
        for(j=0;j<nb;j++)
        {
            matacc(Mp,i,j,nb)=0.0;
            for(k=0;k<na;k++)
                matacc(Mp,i,j,nb)+=matacc(Ma,i,k,na)*matacc(Mb,k,j,nb);
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
            printf("%g ",matacc(Ma,i,j,na));
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
            matacc(Mat,j,i,ma)=matacc(Ma,i,j,ma);
        }
    }
}


void MatVec(double* Ma,double *V,int N,int M,double *RVec)
{
    int i,j;

    for(i=0;i<N;i++)
    {
        RVec[i]=0.0;
        for(j=0;j<M;j++)
            RVec[i]+=matacc(Ma,i,j,N)*V[j];
    }
}

/*Modified QR decomposition*/
int RPR(double *R,double *mu,double *qvec,int* pcol,int m,int n,double *Q)
{
    int i,j,h;
    int k;
    int kprev;
    int tempi;
    int maxid;
    double tempd;
    double *n12;
    double *gamma;
    double quote;
    double gammamax;
    double sbeta,sbetainv;
    double s;
    double *work;
    int info;

    n12=(double*)malloc(m*sizeof(double));
    gamma=(double*)malloc(n*sizeof(double));
    work=(double*)malloc(m*sizeof(double));
    kprev=0;
    for(i=0;i<n;i++)
        pcol[i]=i;
    for(i=0;i<m;i++)
        qvec[i]=0.0;
    if(Q!=NULL)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<m;j++)
            {
                if(i==j)
                    matacc(Q,i,j,m)=1.0;
                else
                    matacc(Q,i,j,m)=0.0;
            }
        }
    }
    for(k=0;k<n;k++)
    {
        for(i=0;i<m;i++)
            n12[i]=1.0;
        for(i=k;i<m;i++)
        {
            if((mu[i]!=0)&&(fabs(mu[k])<9.e99))
                n12[i]=mu[k]/mu[i];
        }
        if((k==0)||((mu[kprev]==0.0)&&(mu[k]>0.0)))
        {
            quote=9.e99;
        }
        else
        {
            if(mu[k]==mu[kprev])
                quote=1.0;
            else
                quote=mu[k]/mu[kprev];
        }
        if(quote>mtol)
        {
            for(j=k;j<n;j++)
            {
                gamma[j]=0.0;
                for(h=k;h<m;h++)
                    gamma[j]+=(n12[h]*matacc(R,h,j,n))*(n12[h]*matacc(R,h,j,n));
            }
            kprev=k;
        }
        else
        {
            if(mu[k]==mu[k-1])
                quote=1.0;
            else
                quote=mu[k]/mu[k-1];
            for(j=k;j<n;j++)
                gamma[j]=quote*quote*(gamma[j]-matacc(R,k-1,j,n)*matacc(R,k-1,j,n));
        }
        gammamax=0.0;
        maxid=k;
        for(i=k;i<n;i++)
        {
            if(gamma[i]>gammamax)
            {
                maxid=i;
                gammamax=gamma[i];
            }
        }
        tempd=gamma[k];
        gamma[k]=gamma[maxid];
        gamma[maxid]=tempd;
        tempi=pcol[k];
        pcol[k]=pcol[maxid];
        pcol[maxid]=tempi;
        for(i=0;i<m;i++)
        {
            tempd=matacc(R,i,k,n);
            matacc(R,i,k,n)=matacc(R,i,maxid,n);
            matacc(R,i,maxid,n)=tempd;
        }
        sbeta=0.0;
        for(i=k;i<m;i++)
        {
            qvec[i]=n12[i]*matacc(R,i,k,n);
            sbeta+=qvec[i]*qvec[i];
        }
        sbeta=sqrt(sbeta);
        for(i=k;i<m;i++)
            qvec[i]*=n12[i];
        if(sbeta==0.0)
        {
            info=k-1;
            free((void*)n12);
            free((void*)gamma);
            return info;
        }
        else
        {
            if(matacc(R,k,k,n)!=0.0)
            {
                if(matacc(R,k,k,n)<0.0)
                    sbeta*=-1.0;
            }
            sbetainv=1.0/sbeta;
            for(i=k;i<m;i++)
            {
                matacc(R,i,k,n)*=sbetainv;
                qvec[i]*=sbetainv;
            }
            matacc(R,k,k,n)+=1.0;
            qvec[k]=matacc(R,k,k,n);
            for(j=k+1;j<n;j++)
            {
                s=0.0;
                for(i=k;i<m;i++)
                    s-=qvec[i]*matacc(R,i,j,n);
                s/=qvec[k];
                for(i=k;i<m;i++)
                    matacc(R,i,j,n)+=s*matacc(R,i,k,n);
            }
            if(Q!=NULL)
            {
                for(i=0;i<m;i++)
                {
                    work[i]=0.0;
                    for(j=k;j<m;j++)
                        work[i]+=matacc(Q,i,j,m)*matacc(R,j,k,n);
                }
                for(i=0;i<m;i++)
                {
                    for(j=k;j<m;j++)
                    {
                        matacc(Q,i,j,m)-=work[i]*qvec[j]/qvec[k];
                    }
                }
            }
            matacc(R,k,k,n)=-sbeta;
        }
    }
    info=n;
    free((void*)work);
    free((void*)n12);
    free((void*)gamma);
    return info;
}
            

/*Modified QR solver.  Takes decomposed matrix and all */
int RPRSolve(double *A,double *x,double *b,double *mu,int* pcol,double *qvec,double *lambda,int m,int n,int p)
{
    double *v;
    double *Qinvb;
    int i,j;
    double ni;
    double temp;
    int h;
    double *alfa;
    double *pt;
    int allnonzero;


    v=(double*)malloc(m*sizeof(double));
    alfa=(double*)malloc(m*sizeof(double));
    Qinvb=(double*)malloc(m*p*sizeof(double));
    pt=(double*)malloc(m*p*sizeof(double));
    for(i=0;i<m;i++)
    {
        v[i]=0.0;
        for(h=0;h<p;h++)
        {
            matacc(Qinvb,i,h,p)=matacc(b,i,h,p);
            matacc(lambda,i,h,p)=0.0;
        }
    }
    for(j=0;j<n;j++)
    {
        for(i=j+1;i<m;i++)
        {
            if((mu[i]==0)||(mu[j]>=1.0e99))
                ni=1.0;
            else
                ni=mu[j]/mu[i];
            v[i]=ni*ni*matacc(A,i,j,n);
        }
        temp=matacc(A,j,j,n);
        v[j]=qvec[j];
        matacc(A,j,j,n)=qvec[j];
        for(h=0;h<p;h++)
        {
            alfa[h]=0.0;
            for(i=j;i<m;i++)
            {
                alfa[h]-=matacc(Qinvb,i,h,p)*v[i];
            }
            alfa[h]/=qvec[j];
        }
        for(i=j;i<m;i++)
        {
            for(h=0;h<p;h++)
            {
                matacc(Qinvb,i,h,p)+=matacc(A,i,j,n)*alfa[h];
            }
        }
        matacc(A,j,j,n)=temp;
    }
    for(i=0;i<n;i++)
    {
        for(h=0;h<p;h++)
        {
            matacc(x,i,h,p)=matacc(Qinvb,i,h,p);
        }
    }
    for(i=n-1;i>=0;i--)
    {
        if(matacc(A,i,i,n)!=0.0)
        {
            for(h=0;h<p;h++)
                matacc(x,i,h,p)/=matacc(A,i,i,n);
            if(i!=0)
            {
                for(j=0;j<i;j++)
                {
                    for(h=0;h<p;h++)
                        matacc(x,j,h,p)-=matacc(A,j,i,n)*matacc(x,i,h,p);
                }
            }
        }
        else
        {
            return i-1;
        }
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<p;j++)
        {
            matacc(pt,pcol[i],j,p)=matacc(x,i,j,p);
        }
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<p;j++)
        {
            matacc(x,i,j,p)=matacc(pt,i,j,p);
        }
    }
    allnonzero=1;
    for(j=0;j<n;j++)
        for(h=0;h<p;h++)
            matacc(lambda,j,h,p)=0.0;   
    for(j=n;j<m;j++)
    {
        if(mu[j]==0.0)
            allnonzero=0;
    }
    if(allnonzero)
    {
        for(j=n;j<m;j++)
        {
            for(h=0;h<p;h++)
                matacc(lambda,j,h,p)=matacc(Qinvb,j,h,p)/(mu[j]*mu[j]);
        }
        for(j=n-1;j>=0;j--)
        {
            for(i=j+1;i<m;i++)
            {
                if((mu[i]==0.0)||(mu[j]>=9.e99))
                {
                    ni=1.0;
                }
                else
                {
                    ni=mu[j]/mu[i];
                }
                v[i]=ni*ni*matacc(A,i,j,n);
            }
            temp=matacc(A,j,j,n);
            v[j]=qvec[j];
            matacc(A,j,j,n)=qvec[j];
            for(h=0;h<p;h++)
            {
                alfa[h]=0.0;
                for(i=j;i<m;i++)
                {
                    alfa[h]-=matacc(lambda,i,h,p)*matacc(A,i,j,n);
                }
                alfa[h]/=qvec[j];
            }
            for(i=j;i<m;i++)
            {
                for(h=0;h<p;h++)
                {
                    matacc(lambda,i,h,p)+=v[i]*alfa[h];
                }
            }
            matacc(A,j,j,n)=temp;
        }
    }
    else
    {
        printf("Lambda is undefined.\n");
    }
    free((void*)pt);
    free((void*)Qinvb);
    free((void*)alfa);
    free((void*)v);
    return n;
}

/*LU decomposition done in place */
void LUDecomp(double *M,int N,int* pc)
{
    int i,k,m;
    double mmax;
    double temp;

    for(k=0;k<N-1;k++)
    {
        m=k;
        mmax=fabs(M[m*N+k]);
        for(i=k;i<N;i++)
        {
            if(fabs(M[i*N+k])>mmax)
            {
                mmax=fabs(M[m*N+k]);
                m=i;
            }
        }
        pc[k]=m;
        for(i=k;i<N;i++)
        {
            temp=M[k*N+i];
            M[k*N+i]=M[m*N+i];
            M[m*N+i]=temp;
        }
        if(M[k*N+k]!=0.0)
        {
            for(i=k+1;i<N;i++)
            {
                M[i*N+k]/=M[k*N+k];
                for(m=k+1;m<N;m++)
                    M[i*N+m]-=M[i*N+k]*M[k*N+m];
            }
        }
    }
}

/*Back substitution into LU decomposed matrix.  b is destroyed! */
void LUSolve(double* LM,double* b,double *x,int N,int* pc)
{
    double temp;
    int i,j,k;

    for(k=0;k<N-1;k++)
    {
        temp=b[k];
        b[k]=b[pc[k]];
        b[pc[k]]=temp;
        for(i=k+1;i<N;i++)
            b[i]-=b[k]*LM[i*N+k];
    }
    for(k=N-1;k>=0;k--)
    {
        x[k]=b[k];
        for(j=k+1;j<N;j++)
        {
            x[k]-=LM[k*N+j]*x[j];
        }
        x[k]/=LM[k*N+k];
    }
}
            
/*If the Jacobian is ill-conditioned the modified QR algorithm can't find a unique answer.
Regularization augments the problem with additional equations and finds a solution.  I think
this is Tikhonov, but it may be either a bastardization thereof or some brute-force kludge which
I made work at some point.*/
void Regularize(double *J,double* b,double *x,int M,int N,double h)
{
    int i,j,k;
    double* Jaug;
    double *xx;
    double *yy;
    double *yt;
    double *yaug;
    double xadj;
    int* pc;

    Jaug=(double*)malloc(N*N*sizeof(double));
    xx=(double*)malloc(N*sizeof(double));
    yy=(double*)malloc(M*sizeof(double));
    yaug=(double*)malloc(N*sizeof(double));
    yt=(double*)malloc(M*sizeof(double));
    pc=(int*)malloc(N*sizeof(int));
    for(i=0;i<N;i++)
    {
        yy[i]=b[i];
        xx[i]=0.0;
        for(j=0;j<N;j++)
        {
            Jaug[i*N+j]=0.0;
            for(k=0;k<M;k++)
            {
                Jaug[i*N+j]+=J[k*N+i]*J[k*N+j];
            }
            if(i==j)
            {
                Jaug[i*N+j]+=h;
            }
        }
    }
    for(i=N;i<M;i++)
        yy[i]=b[i];
    LUDecomp(Jaug,N,pc);
    do
    {
        for(i=0;i<N;i++)
        {
            yaug[i]=0.0;
            for(k=0;k<M;k++)
                yaug[i]+=J[k*N+i]*yy[k];
        }
        LUSolve(Jaug,yaug,yt,N,pc);
        xadj=0.0;
        for(i=0;i<N;i++)
        {
            xx[i]+=yt[i];
            xadj+=yt[i]*yt[i];
        }
        for(i=0;i<M;i++)
        {
            yt[i]=0.0;
            for(k=0;k<N;k++)
                yt[i]+=J[i*N+k]*xx[k];
            yy[i]=b[i]-yt[i];
        }
        printf("%g %g\n",h,xadj);
    }while(xadj>1.0e-6);
    for(i=0;i<N;i++)
        x[i]=xx[i];
    free((void*)Jaug);
    free((void*)xx);
    free((void*)yy);
    free((void*)yaug);
    free((void*)yt);
    free((void*)pc);
}
