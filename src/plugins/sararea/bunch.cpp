#include "bunch.h"

namespace Plugin {
namespace Surfarea {
    
double Cutoff;

void SISolver(int* P,double *d,double *e,double *L,double *x,double *b,int Order)
{
/*solves the system.  d is the diagonal, e the superdiagonal (if applicable) 
  P is the permutation matrix, and L is lower triangular.  For an initial 
  matrix A, the Bunch-Kaufman algorithm calculated P,L,d, and e such that if
  D is the matrix with diagonal d and sub- and superdiagonals e, then
  PAP*=LDL*.  P is stored in permutation format; L is stored in diagonal 
  format.*/
    double *bp;
    double *z;
    double *w;
    double *y;
    double frac;
    int i,j;

    bp=(double*)malloc(Order*sizeof(double));
    z=(double*)malloc(Order*sizeof(double));
    w=(double*)malloc(Order*sizeof(double));
    y=(double*)malloc(Order*sizeof(double));
    for(i=0;i<Order;i++)
        bp[i]=b[P[i]];
//Now we backsolve Lz=bp.
    z[0]=bp[0]/L[0];
    for(i=1;i<Order;i++)
    {
        for(j=0;j<i;j++)
            bp[i]-=L[j+Order*(i-j)-(i-j)*(i-j-1)/2]*z[j];
        z[i]=bp[i]/L[i];
    }
    i=0;
    while(i<Order)
    {
        if(e[i]==0.0)
        {
            if(fabs(d[i])<1.e-9)
                d[i]=1.e-9;
            w[i]=z[i]/d[i];
            i+=1;
        }
        else
        {
            if(fabs(d[i])>fabs(e[i]))
            {
                frac=e[i]/d[i];
                z[i+1]-=z[i]*frac;
                w[i+1]=z[i+1]/(d[i+1]-e[i]*frac);
                w[i]=(z[i]-e[i]*w[i+1])/d[i];
            }
            else
            {
                frac=d[i]/e[i];
                z[i]-=z[i+1]*frac;
                w[i+1]=z[i]/(e[i]-d[i+1]*frac);
                w[i]=(z[i+1]-d[i+1]*w[i+1])/e[i];
            }
            i+=2;
        }
    }
    y[Order-1]=w[Order-1]/L[Order-1];
    for(i=Order-2;i>=0;i--)
    {
        for(j=Order-1;j>i;j--)
            w[i]-=L[i+Order*(j-i)-(j-i)*(j-i-1)/2]*y[j];
        y[i]=w[i]/L[i];
    }
    for(i=0;i<Order;i++)
        x[P[i]]=y[i];
    free((void*)bp);
    free((void*)w);
    free((void*)z);
    free((void*)y);
}


void IncrInert(double e,int Inertia[])
{
    if(fabs(e)>Cutoff)
    {
        if(e>0.0)
            Inertia[0]+=1;
        else
            Inertia[1]+=1;
    }
    else
    {
        Inertia[2]+=1;
    }
}

void GJSolve(double **E,double **x,double **b,int NCol)
{
    int i,j;
    double frac;
    double **b2;

    b2=(double**)malloc(2*sizeof(double*));
    for(i=0;i<2;i++)
    {
        b2[i]=(double*)malloc(NCol*sizeof(double));
        for(j=0;j<NCol;j++)
            b2[i][j]=b[i][j];
    }
    if(fabs(E[0][0])>=fabs(E[1][0]))
    {
        frac=E[1][0]/E[0][0];
        for(i=0;i<2;i++)
        {
            E[1][i]-=frac*E[0][i];
        }
        for(i=0;i<NCol;i++)
        {
            b2[1][i]-=frac*b2[0][i];
        }
        for(i=0;i<NCol;i++)
        {
            x[1][i]=b2[1][i]/E[1][1];
            x[0][i]=(b2[0][i]-E[0][1]*x[1][i])/E[0][0];
        }
    }
    else
    {
        frac=E[0][0]/E[1][0];
        for(i=0;i<2;i++)
        {
            E[0][i]-=frac*E[1][i];
        }
        for(i=0;i<NCol;i++)
        {
            b2[0][i]-=frac*b2[1][i];
        }
        for(i=0;i<NCol;i++)
        {
            x[1][i]=b2[0][i]/E[0][1];
            x[0][i]=(b2[1][i]-E[1][1]*x[1][i])/E[1][0];
        }
    }
    for(i=0;i<2;i++)
        free((void*)b2[i]);
    free((void*)b2);
}

void GJSolve(double *E,double *x,double *b,int NCol)
{
    int i;
    double frac;
    double *b2;

    b2=(double*)malloc(2*NCol*sizeof(double*));
    for(i=0;i<2*NCol;i++)
    {
        b2[i]=b[i];
    }
    if(fabs(E[0])>=fabs(E[2]))
    {
        frac=E[2]/E[0];
        for(i=0;i<2;i++)
        {
            E[2+i]-=frac*E[i];
        }
        for(i=0;i<NCol;i++)
        {
            b2[2+i]-=frac*b2[i];
        }
        for(i=0;i<NCol;i++)
        {
            x[2+i]=b2[2+i]/E[3];
            x[i]=(b2[i]-E[1]*x[2+i])/E[0];
        }
    }
    else
    {
        frac=E[0]/E[2];
        for(i=0;i<2;i++)
        {
            E[i]-=frac*E[2+i];
        }
        for(i=0;i<NCol;i++)
        {
            b2[i]-=frac*b2[2+i];
        }
        for(i=0;i<NCol;i++)
        {
            x[2+i]=b2[i]/E[1];
            x[i]=(b2[2+i]-E[3]*x[2+i])/E[2];
        }
    }
    free((void*)b2);
}



void BunchKaufman(double **Mat,int Order,int Inertia[],double *d,double *e,int *P,double *L,int Mag)
{
    double **E;
    double **C;
    double **Omega;
    double temp;
    int i,j,k,l;
    int r;
    int tint;
    double lambda;
    double sigma;
    double alpha=0.64038820320220756872767623199676; //(1+sqrt(17))/8
    
    Cutoff=1.0;
    for(i=0;i<Mag;i++)
        Cutoff*=0.1;
    E=(double**)malloc(2*sizeof(double*));
    C=(double**)malloc(Order*sizeof(double*));
    Omega=(double**)malloc(Order*sizeof(double*));
    for(i=0;i<Order;i++)
    {
        if(i<2)
            E[i]=(double*)malloc(2*sizeof(double));
        C[i]=(double*)malloc((Order-2)*sizeof(double));
        Omega[i]=(double*)malloc((Order-2)*sizeof(double));
        if(L!=NULL)
            L[i]=1.0;
        if(P!=NULL)
            P[i]=i;
        if(e!=NULL)
            e[i]=0.0;
    }
    if(L!=NULL)
    {
        for(i=Order;i<Order*(Order+1)/2;i++)
            L[i]=0.0;
    }
    for(i=0;i<3;i++)
        Inertia[i]=0;
    k=0;
    while(k<Order)
    {
        lambda=-1.0;
        for(j=k+1;j<Order;j++)
        {
            if(fabs(Mat[j][k])>lambda)
            {
                lambda=fabs(Mat[j][k]);
                r=j;
            }
        }
        if((lambda==0.0)||(lambda==-1.0))
        {
            E[0][0]=Mat[k][k];
            IncrInert(E[0][0],Inertia);
            if(d!=NULL)
                d[k]=E[0][0];
            k+=1;
        }
        else
        {
            if(fabs(Mat[k][k])>=alpha*lambda)
            {
                E[0][0]=Mat[k][k];
                if(d!=NULL)
                    d[k]=E[0][0];
                if(L!=NULL)
                {
                    for(i=0;i<Order;i++)
                    {
                        if(i>k)
                            L[k+Order*(i-k)-(i-k)*(i-k-1)/2]=Mat[i][k]/E[0][0];
                    }
                }
                for(i=k+1;i<Order;i++)
                {
                    for(j=k+1;j<Order;j++)
                    {
                        Mat[i][j]-=Mat[i][k]*Mat[k][j]/E[0][0];
                    }
                }
                IncrInert(E[0][0],Inertia);
                k+=1;
            }
            else
            {
                sigma=-1.0;
                for(j=k;j<Order;j++)
                {
                    if(j!=r)
                    {
                        if(fabs(Mat[j][r])>sigma)
                            sigma=fabs(Mat[j][r]);
                    }
                }
                if(alpha*lambda*lambda<=sigma*fabs(Mat[k][k]))
                {
                    E[0][0]=Mat[k][k];
                    if(L!=NULL)
                    {
                        for(i=0;i<Order;i++)
                        {
                            if(i>k)
                                L[k+Order*(i-k)-(i-k)*(i-k-1)/2]=Mat[i][k]/E[0][0];
                        }
                    }
                    for(i=k+1;i<Order;i++)
                    {
                        {
                            for(j=k+1;j<Order;j++)
                            {
                                Mat[i][j]-=Mat[i][k]*Mat[k][j]/E[0][0];
                            }
                        }
                    }
                    IncrInert(E[0][0],Inertia);
                    k+=1;
                }
                else
                {
                    if(fabs(Mat[r][r])>=alpha*sigma)
                    {
                        for(j=k;j<Order;j++)
                        {
                            temp=Mat[k][j];
                            Mat[k][j]=Mat[r][j];
                            Mat[r][j]=temp;
                        }
                        for(j=k;j<Order;j++)
                        {
                            temp=Mat[j][k];
                            Mat[j][k]=Mat[j][r];
                            Mat[j][r]=temp;
                        }
                        if(P!=NULL)
                        {
                            tint=P[k];
                            P[k]=P[r];
                            P[r]=tint;
                        }
                        E[0][0]=Mat[k][k];
                        if(d!=NULL)
                            d[k]=E[0][0];
                        if(L!=NULL)
                        {
                            for(i=0;i<k;i++)
                            {
                                temp=L[i+(k-i)*Order-(k-i)*(k-i-1)/2];
                                L[i+(k-i)*Order-(k-i)*(k-i-1)/2]=L[i+(r-i)*Order-(r-i)*(r-i-1)/2];
                                L[i+(r-i)*Order-(r-i)*(r-i-1)/2]=temp;
                            }
                            for(i=0;i<Order;i++)
                            {
                                if(i>k)
                                    L[k+Order*(i-k)-(i-k)*(i-k-1)/2]=Mat[i][k]/E[0][0];
                            }
                        }
                        for(i=k+1;i<Order;i++)
                        {
                            for(j=k+1;j<Order;j++)
                            {
                                Mat[i][j]-=Mat[i][k]*Mat[k][j]/E[0][0];
                            }   
                        }
                        IncrInert(E[0][0],Inertia);
                        k+=1;
                    }
                    else  //We need to do a 2x2 pivot here
                    {
                        for(j=k;j<Order;j++)
                        {
                            temp=Mat[k+1][j];
                            Mat[k+1][j]=Mat[r][j];
                            Mat[r][j]=temp;
                        }
                        for(j=k;j<Order;j++)
                        {
                            temp=Mat[j][k+1];
                            Mat[j][k+1]=Mat[j][r];
                            Mat[j][r]=temp;
                        }
                        if(P!=NULL)
                        {
                            tint=P[k+1];
                            P[k+1]=P[r];
                            P[r]=tint;
                        }
                        if(L!=NULL)
                        {
                            for(i=0;i<k+1;i++)
                            {
                                temp=L[i+(k+1-i)*Order-(k+1-i)*(k-i)/2];
                                L[i+(k+1-i)*Order-(k+1-i)*(k-i)/2]=L[i+(r-i)*Order-(r-i)*(r-i-1)/2];
                                L[i+(r-i)*Order-(r-i)*(r-i-1)/2]=temp;
                            }
                        }
                        for(i=0;i<2;i++)
                        {
                            for(j=0;j<2;j++)
                                E[i][j]=Mat[i+k][j+k];
                            for(j=0;j<Order-k-2;j++)
                            {
                                C[i][j]=Mat[i+k][k+j+2];
                            }
                        }
                        if(d!=NULL)
                        {
                            d[k]=E[0][0];
                            d[k+1]=E[1][1];
                        }
                        if(e!=NULL)
                            e[k]=E[0][1];
//                      printf("Trying 2x2. %i\n",k);
                        GJSolve(E,Omega,C,Order-k-2);
                        for(i=k+2;i<Order;i++)
                        {
                            for(j=k+2;j<Order;j++)
                            {
                                temp=0.0;
                                for(l=0;l<2;l++)
                                {
                                    temp+=C[l][i-k-2]*Omega[l][j-k-2];
                                }
                                Mat[i][j]-=temp;
                            }
                        }
                        if(L!=NULL)
                        {
                            for(i=0;i<Order;i++)
                            {
                                if(i>k+1)
                                {
                                    L[k+(i-k)*Order-(i-k)*(i-k-1)/2]=Omega[0][i-k-2];
                                    L[k+1+(i-k-1)*Order-(i-k-1)*(i-k-2)/2]=Omega[1][i-k-2];
//                                  Ll[k+2+(i-k-2)*Order-(i-k-2)*(i-k-3)/2]=Omega[1][i];
                                }
                            }
                        }
                        Inertia[0]+=1;
                        Inertia[1]+=1;
                        k+=2;
                    }
                }
            }
        }
#ifdef dbgbk
        if(L!=NULL)
        {
            for(i=0;i<Order;i++)
            {
                for(j=0;j<Order;j++)
                {
                    if(j>i)
                        temp=0;
                    else
                        temp=L[j+(i-j)*Order-(i-j)*(i-j-1)/2];
                    printf("%6.4f ",temp);
                }
                printf("\n");
            }
            printf("\n");
            for(i=0;i<Order;i++)
            {
                for(j=0;j<Order;j++)
                {
                    printf("%6.4f ",Mat[i][j]);
                }
                printf("\n");
            }
            printf("\n");
//          PauseForAction();
        }
#endif
    }
    for(i=0;i<Order;i++)
    {
        if(i<2)
            free((void*)E[i]);
        free((void*)C[i]);
        free((void*)Omega[i]);
    }
    free((void*)E);
    free((void*)C);
    free((void*)Omega);
}

void BunchKaufman(double *Mat,int Order,int Inertia[],double *d,double *e,int *P,double *L,int Mag)
{
    double *E;
    double *C;
    double *Omega;
    double temp;
    int i,j,k,l;
    int r;
    int tint;
    double lambda;
    double sigma;
    double alpha=0.64038820320220756872767623199676; //(1+sqrt(17))/8
    
    Cutoff=1.0;
    for(i=0;i<Mag;i++)
        Cutoff*=0.1;
    E=(double*)malloc(4*sizeof(double*));
    C=(double*)malloc(Order*(Order-2)*sizeof(double*));
    Omega=(double*)malloc(Order*(Order-2)*sizeof(double*));
    for(i=0;i<Order;i++)
    {
        if(L!=NULL)
            L[i]=1.0;
        if(P!=NULL)
            P[i]=i;
        if(e!=NULL)
            e[i]=0.0;
    }
    if(L!=NULL)
    {
        for(i=Order;i<Order*(Order+1)/2;i++)
            L[i]=0.0;
    }
    for(i=0;i<3;i++)
        Inertia[i]=0;
    k=0;
    while(k<Order)
    {
        lambda=-1.0;
        for(j=k+1;j<Order;j++)
        {
            if(fabs(Mat[Order*j+k])>lambda)
            {
                lambda=fabs(Mat[Order*j+k]);
                r=j;
            }
        }
        if((lambda==0.0)||(lambda==-1.0))
        {
            E[0]=Mat[k*Order+k];
            IncrInert(E[0],Inertia);
            if(d!=NULL)
                d[k]=E[0];
            k+=1;
        }
        else
        {
            if(fabs(Mat[k*Order+k])>=alpha*lambda)
            {
                E[0]=Mat[k*Order+k];
                if(d!=NULL)
                    d[k]=E[0];
                if(L!=NULL)
                {
                    for(i=0;i<Order;i++)
                    {
                        if(i>k)
                            L[k+Order*(i-k)-(i-k)*(i-k-1)/2]=Mat[i*Order+k]/E[0];
                    }
                }
                for(i=k+1;i<Order;i++)
                {
                    for(j=k+1;j<Order;j++)
                    {
                        Mat[i*Order+j]-=Mat[i*Order+k]*Mat[k*Order+j]/E[0];
                    }
                }
                IncrInert(E[0],Inertia);
                k+=1;
            }
            else
            {
                sigma=-1.0;
                for(j=k;j<Order;j++)
                {
                    if(j!=r)
                    {
                        if(fabs(Mat[j*Order+r])>sigma)
                            sigma=fabs(Mat[j*Order+r]);
                    }
                }
                if(alpha*lambda*lambda<=sigma*fabs(Mat[k*Order+k]))
                {
                    E[0]=Mat[k*Order+k];
                    if(L!=NULL)
                    {
                        for(i=0;i<Order;i++)
                        {
                            if(i>k)
                                L[k+Order*(i-k)-(i-k)*(i-k-1)/2]=Mat[i*Order+k]/E[0];
                        }
                    }
                    for(i=k+1;i<Order;i++)
                    {
                        {
                            for(j=k+1;j<Order;j++)
                            {
                                Mat[i*Order+j]-=Mat[i*Order+k]*Mat[k*Order+j]/E[0];
                            }
                        }
                    }
                    IncrInert(E[0],Inertia);
                    k+=1;
                }
                else
                {
                    if(fabs(Mat[r*Order+r])>=alpha*sigma)
                    {
                        for(j=k;j<Order;j++)
                        {
                            temp=Mat[k*Order+j];
                            Mat[k*Order+j]=Mat[r*Order+j];
                            Mat[r*Order+j]=temp;
                        }
                        for(j=k;j<Order;j++)
                        {
                            temp=Mat[j*Order+k];
                            Mat[j*Order+k]=Mat[j*Order+r];
                            Mat[j*Order+r]=temp;
                        }
                        if(P!=NULL)
                        {
                            tint=P[k];
                            P[k]=P[r];
                            P[r]=tint;
                        }
                        E[0]=Mat[k*Order+k];
                        if(d!=NULL)
                            d[k]=E[0];
                        if(L!=NULL)
                        {
                            for(i=0;i<k;i++)
                            {
                                temp=L[i+(k-i)*Order-(k-i)*(k-i-1)/2];
                                L[i+(k-i)*Order-(k-i)*(k-i-1)/2]=L[i+(r-i)*Order-(r-i)*(r-i-1)/2];
                                L[i+(r-i)*Order-(r-i)*(r-i-1)/2]=temp;
                            }
                            for(i=0;i<Order;i++)
                            {
                                if(i>k)
                                    L[k+Order*(i-k)-(i-k)*(i-k-1)/2]=Mat[i*Order+k]/E[0];
                            }
                        }
                        for(i=k+1;i<Order;i++)
                        {
                            for(j=k+1;j<Order;j++)
                            {
                                Mat[i*Order+j]-=Mat[i*Order+k]*Mat[k*Order+j]/E[0];
                            }   
                        }
                        IncrInert(E[0],Inertia);
                        k+=1;
                    }
                    else  //We need to do a 2x2 pivot here
                    {
                        for(j=k;j<Order;j++)
                        {
                            temp=Mat[(k+1)*Order+j];
                            Mat[(k+1)*Order+j]=Mat[r*Order+j];
                            Mat[r*Order+j]=temp;
                        }
                        for(j=k;j<Order;j++)
                        {
                            temp=Mat[j*Order+k+1];
                            Mat[j*Order+k+1]=Mat[j*Order+r];
                            Mat[j*Order+r]=temp;
                        }
                        if(P!=NULL)
                        {
                            tint=P[k+1];
                            P[k+1]=P[r];
                            P[r]=tint;
                        }
                        if(L!=NULL)
                        {
                            for(i=0;i<k+1;i++)
                            {
                                temp=L[i+(k+1-i)*Order-(k+1-i)*(k-i)/2];
                                L[i+(k+1-i)*Order-(k+1-i)*(k-i)/2]=L[i+(r-i)*Order-(r-i)*(r-i-1)/2];
                                L[i+(r-i)*Order-(r-i)*(r-i-1)/2]=temp;
                            }
                        }
                        for(i=0;i<2;i++)
                        {
                            for(j=0;j<2;j++)
                                E[i*2+j]=Mat[(i+k)*Order+j+k];
                            for(j=0;j<Order-k-2;j++)
                            {
                                C[i*(Order-2)+j]=Mat[(i+k)*Order+k+j+2];
                            }
                        }
                        if(d!=NULL)
                        {
                            d[k]=E[0];
                            d[k+1]=E[3];
                        }
                        if(e!=NULL)
                            e[k]=E[1];
//                      printf("Trying 2x2. %i\n",k);
                        GJSolve(E,Omega,C,Order-k-2);
                        for(i=k+2;i<Order;i++)
                        {
                            for(j=k+2;j<Order;j++)
                            {
                                temp=0.0;
                                for(l=0;l<2;l++)
                                {
                                    temp+=C[l*(Order-2)+i-k-2]*Omega[l*(Order-2)+j-k-2];
                                }
                                Mat[i*Order+j]-=temp;
                            }
                        }
                        if(L!=NULL)
                        {
                            for(i=0;i<Order;i++)
                            {
                                if(i>k+1)
                                {
                                    L[k+(i-k)*Order-(i-k)*(i-k-1)/2]=Omega[i-k-2];
                                    L[k+1+(i-k-1)*Order-(i-k-1)*(i-k-2)/2]=Omega[Order-2+i-k-2];
//                                  Ll[k+2+(i-k-2)*Order-(i-k-2)*(i-k-3)/2]=Omega[1][i];
                                }
                            }
                        }
                        Inertia[0]+=1;
                        Inertia[1]+=1;
                        k+=2;
                    }
                }
            }
        }
#ifdef dbgbk
        if(L!=NULL)
        {
            for(i=0;i<Order;i++)
            {
                for(j=0;j<Order;j++)
                {
                    if(j>i)
                        temp=0;
                    else
                        temp=L[j+(i-j)*Order-(i-j)*(i-j-1)/2];
                    printf("%6.4f ",temp);
                }
                printf("\n");
            }
            printf("\n");
            for(i=0;i<Order;i++)
            {
                for(j=0;j<Order;j++)
                {
                    printf("%6.4f ",Mat[i*Order+j]);
                }
                printf("\n");
            }
            printf("\n");
//          PauseForAction();
        }
#endif
    }
    free((void*)E);
    free((void*)C);
    free((void*)Omega);
}

} // namespace Surfarea
} // namespace Plugin
