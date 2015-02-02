#include "matrix.h"

namespace Plugin {
namespace Symmetry {

static double sqrarg;
#define SQR(a)((sqrarg=(a))==0.0?0.0:sqrarg*sqrarg)
#define SIGN(a,b) ((b)>=0.0?::fabs(a):-::fabs(a))
//These are all pretty self-explanatory

double dotp(double* v1,double* v2,int len)
{
    double rv=0.0;
    int i;

    for(i=0;i<len;i++)
    {
        rv+=v1[i]*v2[i];
    }
    return rv;
}

complex dotp(complex* v1,complex* v2,int len)
{
    complex rv=Complex(0.0,0.0);
    int i;

    for(i=0;i<len;i++)
    {
        rv+=v1[i]^v2[i];
    }
    return rv;
}



double pythag(complex a,complex b)
{
        complex aba,abb;

        aba=a^a;
        abb=b^b;
        if(aba.Re()>abb.Re())
                return sqrt(aba.Re())*sqrt(1.0+abb.Re()/aba.Re());
        else
                return (abb==0.0?0.0:sqrt(abb.Re())*sqrt(1.0+aba.Re()/abb.Re()));
}

double pythag(double a,double b)
{
    double aba,abb;

    aba=::fabs(a);
    abb=::fabs(b);
    if(aba>abb)
        return aba*sqrt(1.0+SQR(abb/aba));
    else
        return (abb==0.0?0.0:abb*sqrt(1.0+SQR(aba/abb)));
}

void tred2(double **a,int n,double *d,double *e)
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
                                scale+=::fabs(a[i][k]);
                        if(scale==0.0)
                                e[i]=a[i][l];
                        else
                        {
                                for(k=0;k<=l;k++)
                                {
                                        a[i][k]/=scale;
                                        h+=a[i][k]*a[i][k];
                                }
                                f=a[i][l];
                                g=(f>=0.0?-sqrt(h):sqrt(h));
                                e[i]=scale*g;
                                h-=f*g;
                                a[i][l]=f-g;
                                f=0.0;
                                for(j=0;j<=l;j++)
                                {
                                        a[j][i]=a[i][j]/h;
                                        g=0.0;
                                        for(k=0;k<=j;k++)
                                                g+=a[j][k]*a[i][k];
                                        for(k=j+1;k<=l;k++)
                                                g+=a[k][j]*a[i][k];
                                        e[j]=g/h;
                                        f+=e[j]*a[i][j];
                                }
                                hh=f/(h+h);
                                for(j=0;j<=l;j++)
                                {
                                        f=a[i][j];
                                        e[j]=g=e[j]-hh*f;
                                        for(k=0;k<=j;k++)
                                                a[j][k]-=(f*e[k]+g*a[i][k]);
                                }
                        }
                }
                else
                        e[i]=a[i][l];
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
                                        g+=a[i][k]*a[k][j];
                                for(k=0;k<=l;k++)
                                        a[k][j]-=g*a[k][i];
                        }
                }
                d[i]=a[i][i];
                a[i][i]=1.0;
                for(j=0;j<=l;j++)
                        a[j][i]=a[i][j]=0.0;
        }
}

void tred2ne(double **a,int n,double *d,double *e)
{
        int l,k,j,i;
        double scale,hh,h,g,f;

        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                printf("%g ",a[i][j]);
            }
            printf("\n");
        }
        for(i=n-1;i>=1;i--)
        {
                l=i-1;
                h=scale=0.0;
                if(l>0)
                {
                        for(k=0;k<=l;k++)
                                scale+=::fabs(a[i][k]);
                        if(scale==0.0)
                                e[i]=a[i][l];
                        else
                        {
                                for(k=0;k<=l;k++)
                                {
                                        a[i][k]/=scale;
                                        h+=a[i][k]*a[i][k];
                                }
                                f=a[i][l];
                                g=(f>=0.0?-sqrt(h):sqrt(h));
                                e[i]=scale*g;
                                h-=f*g;
                                a[i][l]=f-g;
                                f=0.0;
                                for(j=0;j<=l;j++)
                                {
//                                      a[j][i]=a[i][j]/h;
                                        g=0.0;
                                        for(k=0;k<=j;k++)
                                                g+=a[j][k]*a[i][k];
                                        for(k=j+1;k<=l;k++)
                                                g+=a[k][j]*a[i][k];
                                        e[j]=g/h;
                                        f+=e[j]*a[i][j];
                                }
                                hh=f/(h+h);
                                for(j=0;j<=l;j++)
                                {
                                        f=a[i][j];
                                        e[j]=g=e[j]-hh*f;
                                        for(k=0;k<=j;k++)
                                                a[j][k]-=(f*e[k]+g*a[i][k]);
                                }
                        }
                }
                else
                        e[i]=a[i][l];
                d[i]=h;
        }
//      d[0]=0.0;
        e[0]=0.0;
        for(i=0;i<n;i++)
        {
/*              l=i-1;
                if(d[i])
                {
                        for(j=0;j<=l;j++)
                        {
                                g=0.0;
                                for(k=0;k<=l;k++)
                                        g+=a[i][k]*a[k][j];
                                for(k=0;k<=l;k++)
                                        a[k][j]-=g*a[k][i];
                        }
                }*/
                d[i]=a[i][i];
/*              a[i][i]=1.0;
                for(j=0;j<=l;j++)
                        a[j][i]=a[i][j]=0.0;*/
        }
}

int tqli(double *d,double *e,int n,double **z)
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
                                dd=::fabs(d[m])+::fabs(d[m+1]);
                                if((::fabs(e[m])+dd)==dd)
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
                                                f=z[k][i+1];
                                                z[k][i+1]=s*z[k][i]+c*f;
                                                z[k][i]=c*z[k][i]-s*f;
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

int tqline(double *d,double *e,int n,double **z)
{
        int m,l,iter,i;
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
                                dd=::fabs(d[m])+::fabs(d[m+1]);
                                if((::fabs(e[m])+dd)==dd)
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
/*                                      for(k=0;k<n;k++)
                                        {
                                                f=z[k][i+1];
                                                z[k][i+1]=s*z[k][i]+c*f;
                                                z[k][i]=c*z[k][i]-s*f;
                                        }*/
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

int choldc(double **a,int n,double p[])
{
/* Cholesky decomposition of a real matrix */
    int i,j,k;
    double sum;

    for(i=0;i<n;i++)
    {
        for(j=i;j<n;j++)
        {
            for(sum=a[i][j],k=i-1;k>=0;k--)
                sum-=a[i][k]*a[j][k];
            if(i==j)
            {
                if(sum<=0.0)
                {
                    printf("Nonphysical overlap matrix\n");
                    p[i]=1.e-99;
                    return -1;
                }
                else
                    p[i]=sqrt(sum);
            }
            else
                a[j][i]=sum/p[i];
        }
    }
    return 0;
}

int ccholdc(complex **a,int n,double p[])
{
/* Cholesky decomposition of a complex matrix */
    int i,j,k;
    complex sum;

    for(i=0;i<n;i++)
    {
        for(j=i;j<n;j++)
        {
            for(sum=a[i][j],k=i-1;k>=0;k--)
                sum-=a[i][k]*a[j][k];
            if(i==j)
            {
                if((sum.Re()<=0.0)||(sum.Im()!=0.0))
                {
                    printf("Nonphysical overlap matrix\n");
                    p[i]=1.e-99;
                    return -1;
                }
                else
                    p[i]=sqrt(sum.Re());
            }
            else
                a[j][i]=sum/p[i];
        }
    }
    return 0;
}

int ccholsl(complex **a,int n,double p[],complex b[],complex x[])
{
/*Cholesky solution of a complex matrix */
    int i,k;
    complex sum;

    for(i=0;i<n;i++)
    {
        for(sum=b[i],k=i-1;k>=0;k--)
            sum-=a[i][k]*x[k];
        x[i]=sum/p[i];
    }
    for(i=n-1;i>=0;i--)
    {
        for(sum=x[i],k=i+1;k<n;k++)
            sum-=a[k][i]*x[k];
        x[i]=sum/p[i];
    }
    return 0;
}

void geneigset(double **a,double **b,double **aeff,double **lt,int n)
{
/* sets up the cholesky decomposition */
    double **atem;
    double **chk;
    double *p;
    double sum;
    int i;
    int j;
    int k;

    atem=(double**)malloc(n*sizeof(double*));
    chk=(double**)malloc(n*sizeof(double*));
    p=(double*)malloc(n*sizeof(double));
    for(i=0;i<n;i++)
    {
        atem[i]=(double*)malloc(n*sizeof(double));
        chk[i]=(double*)malloc(n*sizeof(double));
    }
    for(i=0;i<n;i++)
        for(j=i;j<n;j++)
            lt[j][i]=lt[i][j]=b[i][j];
    if(choldc(lt,n,p)<0)
    {
        tred2(b,n,chk[1],chk[2]);
        tqli(chk[1],chk[2],n,b);
        for(i=0;i<n;i++)
            if(chk[1][i]<=0.0)
                printf("Offending eigenvalue %i %g \n",i,chk[1][i]);
    }
    for(i=0;i<n;i++)
    {
        lt[i][i]=p[i];
        for(j=i+1;j<n;j++)
            lt[i][j]=lt[j][i];
    }
/*  printf("ltLtt\n");
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            sum=0.0;
            for(k=0;(k<=i)&&(k<=j);k++)
            {
                sum+=lt[i][k]*lt[k][j];
            }
            printf("%g ",sum);
        }
        printf("\n");
    }*/
                
    for(i=0;i<n;i++)
    {   //loop over vectors
        for(j=0;j<n;j++)
        {
            sum=a[i][j];
            for(k=0;k<j;k++)
                sum-=atem[i][k]*lt[k][j];
            atem[i][j]=sum/lt[j][j];
        }
    }
/*  printf("Atem\n");
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%g ",atem[i][j]);
        }
        printf("\n");
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            sum=0.0;
            for(k=0;(k<=j)&&(k<=i);k++)
            {
                sum+=atem[i][k]*lt[k][j];
            }
            printf("%g ",sum);
        }
        printf("\n");
    }*/
    for(i=0;i<n;i++)
    {
        for(j=i;j<n;j++)
        {
            sum=atem[i][j];
            for(k=0;k<i;k++)
                sum-=aeff[k][j]*lt[i][k];
            aeff[i][j]=sum/lt[i][i];
            if(i!=j)
                aeff[j][i]=aeff[i][j];
        }
    }
    for(i=0;i<n;i++)
        free((void*)atem[i]);
    free((void*)atem);
    free((void*)p);
}

void reseteigs(double **evc,double **lt,int n)
{
/* Back transform of the Cholesky thing */
    double sum;
    int i;
    int j;
    int k;

    for(i=0;i<n;i++)
    {   //loop over eigenvectors
        for(j=n-1;j>=0;j--)
        {   //backsubstitution loop
            sum=evc[j][i];
            for(k=j+1;k<n;k++)
            {
                sum-=evc[k][i]*lt[j][k];
            }
            evc[j][i]=sum/lt[j][j];
        }
    }
//Since lt is not orthogonal, the transformed eigenvectors are not of unit length.  We renormalize
// them here.
    for(i=0;i<n;i++)
    {
        sum=0.0;
        for(j=0;j<n;j++)
            sum+=evc[j][i]*evc[j][i];
        sum=sqrt(sum);
        for(j=0;j<n;j++)
            evc[j][i]/=sum;
    }
}



void eigsrt(double d[],double **v,int n)
{
    int i,j,k;
    double p;

    for(i=0;i<n-1;i++)
    {
        p=d[k=i];
        for(j=i+1;j<n;j++)
            if(d[j]>=p)
                p=d[k=j];
        if(k!=i)
        {
            d[k]=d[i];
            d[i]=p;
            for(j=0;j<n;j++)
            {
                p=v[j][i];
                v[j][i]=v[j][k];
                v[j][k]=p;
            }
        }
    }
}

void aeigsrt(double d[],double **v,int n)
{
/*Sorts by magnitude (absolute value)*/
    int i,j,k;
    double p;

    for(i=0;i<n-1;i++)
    {
        p=::fabs(d[k=i]);
        for(j=i+1;j<n;j++)
            if(::fabs(d[j])>=p)
                p=::fabs(d[k=j]);
        if(k!=i)
        {
            p=d[k];
            d[k]=d[i];
            d[i]=p;
            for(j=0;j<n;j++)
            {
                p=v[j][i];
                v[j][i]=v[j][k];
                v[j][k]=p;
            }
        }
    }
}

void ceigsrt(double d[],complex **v,int n)
{
/* sorts complex eigenvectors */
    int i,j,k;
    double p;
    complex pp;

    for(i=0;i<n-1;i++)
    {
        p=d[k=i];
        for(j=i+1;j<n;j++)
            if(d[j]>=p)
                p=d[k=j];
        if(k!=i)
        {
            p=d[k];
            d[k]=d[i];
            d[i]=p;
            for(j=0;j<n;j++)
            {
                pp=v[j][i];
                v[j][i]=v[j][k];
                v[j][k]=pp;
            }
        }
    }
}

void veigsrt(double d[],double **v,int n,complex **w,complex* r,int o)
{
/* Sorts by relationship to chosen vector */
    int i,j,k;
    double p;
    complex* tvec;
    double* dvec;

    tvec=(complex*)malloc(o*sizeof(complex));
    dvec=(double*)malloc(n*sizeof(double));
    for(i=0;i<n;i++)
    {
        for(j=0;j<o;j++)
        {
            tvec[j]=0.0;
            for(k=0;k<n;k++)
            {
                tvec[j]=w[k][j]*v[k][i];
            }
            dvec[i]=dotp(tvec,r,o).magnitude();
        }
    }
    for(i=0;i<n-1;i++)
    {
        p=dvec[k=i];
        for(j=i+1;j<n;j++)
            if(d[j]>=p)
                p=dvec[k=j];
        if(k!=i)
        {
            p=dvec[k];
            dvec[k]=dvec[i];
            dvec[i]=p;
            p=d[k];
            d[k]=d[i];
            d[i]=p;
            for(j=0;j<n;j++)
            {
                p=v[i][j];
                v[i][j]=v[k][j];
                v[k][j]=p;
            }
        }
    }
    free((void*)tvec);
    free((void*)dvec);
}

int aeq(double x,double y)
{
/* Close enough? */
    if(::fabs(x-y)<1.e-9)
        return 1;
    return 0;
}

void SortCRots(double** evecs,double* evals,int n)
{
/* Arranges complex eigenvalues to be paired (i,i+1) columns are paired eigenvectors*/
    int i,j,k;
    complex c1,c2,c3;
    double** cc1;
    double p;
    double ang;

    cc1=(double**)malloc(n*sizeof(double*));
    for(i=0;i<n;i++)
    {
        cc1[i]=(double*)malloc(2*sizeof(double));
    }
    for(i=0;i<2*n-1;i++)
    {
        for(j=i+1;j<2*n;j++)
        {
            if(aeq(evals[i],evals[j]))
            {
                for(k=0;k<n;k++)
                {
                    c1=Complex(evecs[k][i],evecs[k+n][i]);
                    c2=Complex(evecs[k][j],evecs[k+n][j]);
                    if(!aeq(c1.magnitude(),c2.magnitude()))
                    {
                        k=n;
                        cc1[0][0]=-9.0e99;
                    }
                    else
                    {
                        if(aeq(c1.magnitude(),0.0))
                        {
                            cc1[k][0]=0.0;
                            cc1[k][1]=0.0;
                        }
                        else
                        {
                            c3=c2/c1;
                            c3.Argand(cc1[k]);
                        }
                    }
                }
                if(cc1[0][0]!=-9.0e99)
                {
                    ang=9.0e99;
                    for(k=0;k<n;k++)
                    {
                        if((!aeq(cc1[k][0],1.0))&&(!aeq(cc1[k][0],0.0)))
                        {
                            cc1[0][0]=-9.0e99;
                            k=n;
                        }
                        else
                        {
                            if(aeq(cc1[k][0],1.0))
                            {
                                if(ang==9.0e99)
                                {
                                    printf("Attempting rotation angle of %g\n",cc1[k][1]);
                                    ang=cc1[k][1];
                                }
                                else
                                {
                                    if(!aeq(cc1[k][1],ang))
                                    {
                                        cc1[0][0]=-9.0e99;
                                        k=n;
                                    }
                                }
                            }
                        }
                    }
                }
                if(cc1[0][0]!=-9.0e99)
                {
                    printf("%i rotation of %i; angle is %g\n",i,j,ang);
//The two vectors are complex rotations of each other.  Put them next to each other
                    for(k=0;k<2*n;k++)
                    {
                        p=evecs[k][i+1];
                        evecs[k][i+1]=evecs[k][j];
                        evecs[k][j]=p;
                    }
                    p=evals[i+1];
                    evals[i+1]=evals[i];
                    evals[i]=p;
                    j=2*n; //don't risk another swap
                }
            }
        }
    }
}

} // namespace Symmetry
} // namespace Plugin
