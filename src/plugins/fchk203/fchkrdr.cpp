// fchkrdr.cpp : Defines the entry point for the console application.
//

#include "fchkrdr.h"

namespace Plugin {
namespace Fchk {
    
#define pi 3.1415926535897932384626433832795
#define sq2 1.4142135623730950488016887242097
#define symprec 1.e-4
#define NumSobPts 50
#define dbg 1

extern double** FullCoeffs;
extern GaussianOrbital* Basis;
extern GaussianOrbital* SBasisS;
extern GaussianOrbital* SBasisM;
extern GaussAtom* AtomList;
extern int NumGBasis;
extern int fun1,fun2;

extern int NumContShells;
extern int NumUncontShells;
extern int NumIndepFun;
extern int NumAtoms;

extern double* AlphaEnergies;
extern double* BetaEnergies;
extern double* AlphaCoeffs;
extern double* BetaCoeffs;
extern ContractedShell* ContShells;
extern char** AlfOrbID;
extern char** BetaOrbID;

extern int norb;

double absc[30];
double weights[30];

int Order=3;
double Yml[NumSobPts*3];

double SHConstrain(double *Parms)
{
    double rt;
    int i;

    rt=0.0;
    for(i=1;i<Order+1;i++)
        rt+=Parms[i]*Parms[i];
    rt=1.0*fabs(sqrt(rt)-1.0);
    return rt;
}

void dSHConstrain(double *x,double *Parms,double *dFdP)
{
    double rt=0.0;
    int i;

    for(i=1;i<Order+1;i++)
        rt+=Parms[i]*Parms[i];
    rt=sqrt(rt);
    if(rt==1.0)
        return;
    if(rt<1.0)
    {
        for(i=1;i<Order+1;i++)
            dFdP[i]-=1.0*Parms[i]/rt;
    }
    else
    {
        for(i=1;i<Order+1;i++)
            dFdP[i]+=1.0*Parms[i]/rt;
    }
}

double SHFit(double *x,double *Parms)
{
    double coeffs[6];
    double pf=Parms[0];
    double rval=0.0;
    int xi=(int)(*x);
    int i;

    for(i=0;i<Order;i++)
    {
        coeffs[i]=Parms[i+1];
    }
    for(i=0;i<Order;i++)
    {
        rval+=coeffs[i]*Yml[3*xi+i];
    }
    return pf*rval;
}

void dSHFit(double *x,double *Parms,double *dFdP)
{
    double coeffs[6];
    double pf=Parms[0];
    double rval=0.0;
    int xi=(int)(*x);
    int i;

    for(i=0;i<Order;i++)
    {
        coeffs[i]=Parms[i+1];
    }
    for(i=0;i<Order;i++)
    {
        rval+=coeffs[i]*Yml[3*xi+i];
        dFdP[i+1]=pf*Yml[3*xi+i];
    }
    dFdP[0]=rval;
}

int GaussElim(double **M,int *p,int n)
{
    int j,k,r;
    double temp;
    int d=1;

    for(k=0;k<n-1;k++)
    {
        temp=0.0;
        for(j=k;j<n;j++)
        {
            if(fabs(M[j][k])>temp)
            {
                p[k]=j;
                temp=fabs(M[j][k]);
            }
        }
        if(p[k]!=k)
        {
            d*=-1;
            for(j=k;j<n;j++)
            {
                temp=M[k][j];
                M[k][j]=M[p[k]][j];
                M[p[k]][j]=temp;
            }
        }
        if(M[k][k]!=0.0)
        {
            for(r=k+1;r<n;r++)
            {
                M[r][k]=M[r][k]/M[k][k];
                for(j=k+1;j<n;j++)
                    M[r][j]-=M[r][k]*M[k][j];
            }
        }
    }
    return d;
}

void GaussBackSubst(double **M,double *x,double *b,int *p,int N)
{
    int k,j;
    double temp;

    for(k=0;k<N-1;k++)
    {
        temp=b[p[k]];
        b[p[k]]=b[k];
        b[k]=temp;
        for(j=k+1;j<N;j++)
        {
            b[j]-=b[k]*M[j][k];
        }
    }
    for(k=N-1;k>=0;k--)
    {
        for(j=N-1;j>k;j--)
        {
            b[k]-=M[k][j]*x[j];
        }
        x[k]=b[k]/M[k][k];
    }
}

void MSolve(double **M,double *x,double *b,int Order)
{
    int *p;

    p=(int*)malloc(Order*sizeof(int));
    GaussElim(M,p,Order);
    GaussBackSubst(M,x,b,p,Order);
    free((void*)p);
}

double levmar(double **X,double *y,double *weight,int NPts,
            double *Parms,int *vParms,int NParms,
            double (*Func)(double *x,double *Parms),
            void (*dFunc)(double *x,double *Parms,double *dFdP),
            double (*Constraints)(double* Parms),
            void (*dConstraints)(double*,double*,double*))
{
/* Performs a nonlinear least-squares (Levenberg-Marquardt) fit.  The
fitting function is passed as Func, and maps R^(N)->R.  The to-be-fit
x-vectors (X in R^N) are passed as the vectors X[i], where i goes from
0 to NPts-1.  The corresponding to-be-fit data are in y[0..NPts-1].  
The weights are passed explicitly.  vParms flags which Parameters are to
be adjusted.  NParms is the total number of parameters. dFunc returns
the gradient of Func with respect to the parameter vector in dFdP.*/
    double lambda;
    double *beta;
    double **alpha;
    double *dFdP;
    double *yhat;
    double *yhat2;
    double *tParms;
    double chio;
    double chi2;
    int i,j,k;
    int iter=0;
    int count=0;
    int converged=0;

    beta=(double*)malloc(NParms*sizeof(double));
    alpha=(double**)malloc(NParms*sizeof(double*));
    dFdP=(double*)malloc(NParms*sizeof(double));
    tParms=(double*)malloc(NParms*sizeof(double));
    yhat=(double*)malloc(NPts*sizeof(double));
    yhat2=(double*)malloc(NPts*sizeof(double));
    for(i=0;i<NParms;i++)
        alpha[i]=(double*)malloc(NParms*sizeof(double));
    chio=0.0;
    for(i=0;i<NPts;i++)
    {
        yhat[i]=Func(X[i],Parms);
        chio+=weight[i]*(yhat[i]-y[i])*(yhat[i]-y[i]);
    }
    if(Constraints!=NULL)
        chio+=Constraints(Parms);
    lambda=0.001;
    while(!converged)
    {
#ifdef dbg
        printf("%i %g\n",++iter,chio);
        printf("\nCurrent Parameters\n");
        for(i=0;i<NParms;i++)
        {
            printf("%g ",Parms[i]);
        }
        printf("\n"); 
#endif
        for(i=0;i<NParms;i++)
        {
            beta[i]=0.0;
            for(j=0;j<NParms;j++)
                alpha[i][j]=0.0;
        }
        for(k=0;k<NPts;k++)
        {
            dFunc(X[k],Parms,dFdP);
            if(dConstraints!=NULL)
                dConstraints(X[k],Parms,dFdP);
            for(i=0;i<NParms;i++)
            {
                beta[i]+=(y[k]-yhat[k])*(weight[k])*dFdP[i];
                for(j=0;j<NParms;j++)
                {
                    alpha[i][j]+=weight[k]*dFdP[i]*dFdP[j];
                }
            }
        }
        for(i=0;i<NParms;i++)
        {
#ifdef dbg
            printf("%g ",beta[i]);
#endif
            alpha[i][i]*=(1.0+lambda);
        }
#ifdef dbg
        printf("\nLambda=%g\n",lambda);
#endif
        MSolve(alpha,tParms,beta,NParms);
        for(i=0;i<NParms;i++)
        {
#ifdef dbg
            printf("%g ",tParms[i]);
#endif
            if(vParms[i])
                tParms[i]+=Parms[i];
            else
                tParms[i]=Parms[i];
#ifdef dbg
            printf("%g\n",tParms[i]);
#endif
        }
        chi2=0.0;
        for(i=0;i<NPts;i++)
        {
            yhat2[i]=Func(X[i],tParms);
            chi2+=weight[i]*(yhat2[i]-y[i])*(yhat2[i]-y[i]);
        }
        if(Constraints!=NULL)
            chi2+=Constraints(tParms);
#ifdef dbg
        printf("Current: %g Best: %g\n",chi2,chio);
        PauseForAction();
#endif
        if(chi2>chio)
        {
            lambda*=10.0;
        }
        else
        {
            lambda*=0.1;
            for(i=0;i<NParms;i++)
            {
                yhat[i]=yhat2[i];
                Parms[i]=tParms[i];
//              printf("%g ",Parms[i]);
            }
//          printf("\n");
//          PauseForAction();
            if(((chio/chi2-1.0)<1.e-3)||((chio-chi2)<1.e-2))
            {
                if(++count>2)
                    converged=1;
            }
            chio=chi2;
        }
    }
    for(i=0;i<NParms;i++)
        free((void*)alpha[i]);
    free((void*)alpha);
    free((void*)beta);
    free((void*)yhat);
    free((void*)yhat2);
    free((void*)tParms);
    free((void*)dFdP);
    return chi2;
}

double gammaln(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.1800917947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp-=(x+0.5)*log(tmp);
    ser=1.000000000190015;
    for(j=0;j<=5;j++)
        ser+=cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

double aslg(int l,int m,double x)
{
    double fact,pll,pmm,pmmp1,somx2;
    int i,ll;

    pmm=1.0;
    if(m>0)
    {
        somx2=sqrt((1.0-x)*(1.0-x));
        fact=1.0;
        for(i=1;i<=m;i++)
        {
            pmm*=-fact*somx2;
            fact+=2.0;
        }
    }
    if(l==m)
        return pmm;
    else
    {
        pmmp1=x*((double)(2*m+1))*pmm;
        if(l==(m+1))
            return pmmp1;
        else
        {
            for(ll=m+2;ll<=l;ll++)
            {
                pll=(x*(double)(2*ll+1)*pmmp1-(double)(ll+m)*pmm)/(double)(ll-m+1);
                pmm=pmmp1;
                pmmp1=pll;
            }
            return pll;
        }
    }
}

double fact(int N)
{
    int i;
    double v;

    if((N==0)||(N==1))
        return 1.0;
    v=1.0;
    for(i=2;i<=N;i++)
        v*=(double)i;
    return v;
}

double SphHarm(int l,int m,double cth,double phi)
{
    double iphi=phi;
    double Plm;
    double Norm;
    double prefact=1.0;
    double eiphi;

    if(m>0)
    {
        eiphi=sin((double)m*phi);
    }
    else 
    {
        if(m==0)
            eiphi=1.0;
        else
            eiphi=cos((double)(-m)*phi);
    }
    if(m<0)
    {
        m=-m;
        if(m%2!=0)
            prefact=-1.0;
    }
    Norm=(double)(2*l+1)*fact(l-m)/(4.0*pi*fact(l+m));
    Norm=sqrt(Norm);
    Plm=aslg(l,m,cth);
    return eiphi*prefact*Norm*Plm;
}

double EvalLaguerre(int order,double alpha,double xval,double* dlag)
{
    double rval;
    double rvalm1;
    double rvalp1;
    int j;

    rval=1.0;
    rvalm1=0.0;
    j=0;
    while(j<order)
    {
        rvalp1=((-xval+2.0*(double)j+alpha+1.0)*rval-((double)j+alpha)*rvalm1)/((double)j+1.0);
        rvalm1=rval;
        rval=rvalp1;
        j+=1;
    }
    if(dlag!=NULL)
    {
        *dlag=((double)order*rval-((double)order+alpha)*rvalm1)/xval;
    }
    return rval;
}

double FitLaguerre::Eval(double x)
{
    int i;
    double rv=0.0;

    for(i=0;i<NPoly;i++)
    {
        rv+=coeffs[i]*EvalLaguerre(i,0.0,x);
    }
    rv*=exp(-rscale*x);
    return rv;
}

void SetGauLagPar(double *a,double *w,int order,double alpha)
{
    int i,its;
    double v,dv;
    double z,z1,ai;

    for(i=0;i<order;i++)
    {
        if(i==0)
        {
            z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*(double)order+1.8*alpha);
        }
        else
        {
            if(i==1)
            {
                z+=(15.0+6.25*alpha)/(1.0+0.9*alpha+2.5*(double)order);
            }
            else
            {
                ai=(double)i-1.0;
                z+=((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alpha/
                    (1.0+3.5*ai))*(z-a[i-2])/(1.0+0.3*alpha);
            }
        }
        for(its=0;its<30;its++)
        {
            v=EvalLaguerre(order,alpha,z,&dv);
            z1=z;
            z-=v/dv;
            if(fabs(z-z1)<1.e-14)
                break;
        }
        a[i]=z;
        z1=EvalLaguerre(order,alpha,z,&v);
        dv=EvalLaguerre(order-1,alpha,z);
        w[i]=-exp(gammaln(alpha+(double)order)-gammaln((double)order))/(v*dv*(double)order);
    }
}

double func(double cl)
{
    int i,j;
    double FVal;

    for(i=0;i<NumGBasis;i++)
        SBasisM[i].Center[0]=cl;
    FVal=0.0;
    for(i=0;i<NumGBasis;i++)
    {
        for(j=0;j<NumGBasis;j++)
        {
            FVal+=SBasisS[i].NormCoeff*SBasisM[j].NormCoeff*FullCoeffs[i][fun1]*FullCoeffs[j][fun2]*oldrecurse::OverlapIntegral(SBasisS[i],SBasisM[j]);
        }
    }
//  printf("%g %g\n",cl,FVal*exp(cl));
    return FVal;
}


double LagFit(int NMin,int NMax,double aver,double maxer,double* xx,double* yy,int numpts,double xscal,FitLaguerre* FL)
{
    int i,j;
    double fx;
    double coeff[25];
    int NLP;
    double errtrm;
    double avgerr=1.e3,maxerr=0.0;

    NLP=NMin;
    for(j=0;j<=NLP;j++)
    {
        coeff[j]=0.0;
        for(i=0;i<30;i++)
        {
            fx=func(absc[i])*exp(xscal*absc[i])*EvalLaguerre(j,0.0,absc[i]);
            coeff[j]+=fx*weights[i];
        }
//      printf("%i %g\n",j,coeff[j]);
    }
    while((NLP<=NMax)&&((avgerr>aver)||(maxerr>maxer)))
    {
        coeff[NLP]=0.0;
        for(i=0;i<30;i++)
        {
            fx=func(absc[i])*exp(xscal*absc[i])*EvalLaguerre(NLP,0.0,absc[i]);
            coeff[NLP]+=fx*weights[i];
        }
//      printf("%i %g\n",NLP,coeff[NLP]);
        avgerr=0.0;
        maxerr=0.0;
        for(i=0;i<numpts;i++)
        {
            fx=0.0;
            for(j=0;j<NLP;j++)
            {
                fx+=coeff[j]*EvalLaguerre(j,0.0,xx[i]);
            }
            fx*=exp(-xscal*xx[i]);
//          printf("%g %g %g %g\n",xx[i],fx,func(xx[i]),yy[i]);
            errtrm=func(xx[i])-fx;
            avgerr+=errtrm*errtrm;
            if(fabs(errtrm)>maxerr)
                maxerr=fabs(errtrm);
        }
        avgerr=sqrt(avgerr)/numpts;
//      printf("%i %g\t%g\n",NLP,avgerr,maxerr);
//      PauseForAction();
        if((avgerr>aver)||(maxerr>maxer))
            NLP+=1;
    }
//  PauseForAction();
    if(FL!=NULL)
    {
        FL->NPoly=NLP;
        for(i=0;i<NLP;i++)
            FL->coeffs[i]=coeff[i];
        FL->averror=avgerr;
        FL->maxerror=maxerr;
        FL->rscale=xscal;
    }
    return (0.7*maxerr+0.3*avgerr)*(double)NLP;
}

FitLaguerre DetermineFit(int NP,double* r,double* iv,double tol1,double tol2)
{
/*This fit is complicated by the need to scale the independent values.  We are fitting
    a set of polynomials which are orthogonal with the weight function exp(-x); this means
    we need to scale r appropriately.*/
    FitLaguerre FL;
    double* xt;
    double* yt;
    double xsum,ysum,xysum,xxsum;
    double expg,dexpg;
    double res,resh1,resh2,resh3;
    double expgh,expgl,expgt;
    int i;
    int ntpts;
    
//First we need to make an initial guess for our scale factor.  We take the last NP/2 points
//and fit an exponential decay.
    SetGauLagPar(absc,weights,30,0.0);
/*  res=0.0;
    for(i=0;i<30;i++)
    {
        res+=weights[i];
        printf("%g\t%g\n",absc[i],weights[i]);
    }
    printf("%g\n",res);
    PauseForAction();*/
    ntpts=NP/2;
    if(ntpts<3)
        ntpts=3;
    xt=(double*)malloc(ntpts*sizeof(double));
    yt=(double*)malloc(ntpts*sizeof(double));
    xsum=xxsum=xysum=ysum=0.0;
    for(i=0;i<ntpts;i++)
    {
        xt[i]=r[NP-i-1];
        yt[i]=log(fabs(iv[NP-i-1]));
        xsum+=xt[i];
        xxsum+=xt[i]*xt[i];
        xysum+=xt[i]*yt[i];
        ysum+=yt[i];
    }
    expg=(xsum*ysum-(double)ntpts*xysum)/((double)ntpts*xxsum-(xsum*xsum));
    dexpg=(expg>1.e-1?1.0e-4:1.0e-3*expg);
    printf("%g %g\n",expg,dexpg);
    expgl=expg-3.0*dexpg;
    expgh=expg+3.0*dexpg;
    resh1=LagFit(3,25,tol1,tol2,r,iv,NP,expgl);
    res=LagFit(3,25,tol1,tol2,r,iv,NP,expg);
    resh2=LagFit(3,25,tol1,tol2,r,iv,NP,expgh);
    while(resh1<res)
    {
        res=resh1;
        expg=expgl;
        expgl-=3.0*dexpg;
        resh1=LagFit(3,25,tol1,tol2,r,iv,NP,expgl);
        printf("Retrying...%g %g %g %g\n",resh1,res,expgl,expg);
    }
    while(resh2<res)
    {
        resh1=res;
        res=resh2;
        expgl=expg;
        expg=expgh;
        expgh+=(dexpg*=5.0);
        resh2=LagFit(3,25,tol1,tol2,r,iv,NP,expgh);
        printf("Looking to pass minimum..%g %g %g %g\n",resh2,res,expgh,expg);
    }
    while(resh2>1.0)
    {
        expgh=0.5*(expg+expgh);
        resh2=LagFit(3,25,tol1,tol2,r,iv,NP,expgh);
        printf("Closing the gap. %g %g\n",expgh,resh2);
    }
    printf("%g %g %g\n",expgl,expg,expgh);
    printf("%g %g %g\n",resh1,res,resh2);
    while(expgh-expgl>1.e-6)
    {
        if(resh1<resh2)
        {
            expgt=expgl+0.618*(expgh-expgl);
            resh3=LagFit(3,25,tol1,tol2,r,iv,NP,expgh);
        }
        else
        {
            expgt=expgl+0.382*(expgh-expgl);
            resh3=LagFit(3,25,tol1,tol2,r,iv,NP,expgt);
        }
        if(expgt>expg)
        {
            if(resh3<res)
            {
                resh1=res;
                expgl=expg;
                res=resh3;
                expg=expgt;
            }
            else
            {
                resh2=resh3;
                expgh=expgt;
            }
        }
        else
        {
            if(resh3<res)
            {
                resh2=res;
                expgh=expg;
                res=resh3;
                expg=expgt;
            }
            else
            {
                resh1=resh3;
                expgl=expgt;
            }
        }
        printf("%g %g %g\n",expgl,expg,expgh);
        printf("%g %g %g\n",resh1,res,resh2);
    }
//  PauseForAction();
    res=LagFit(3,25,tol1,tol2,r,iv,NP,expg,&FL);
    printf("Best exponent %g; using %i polynomials\n",FL.rscale,FL.NPoly);
    printf("Avg. error %g Maximum error %g\n",FL.averror,FL.maxerror);
    free((void*)xt);
    free((void*)yt);
    return FL;
}





char ExtractLetter(char* orb)
{
    int i;
    char tval;

    i=0;
    while(1)
    {
        tval=orb[i++];
//      printf("%c",tval);
        if(tval>'9')
            return tval;
    }
}

void Rotate(int axis,double* vec)
{
    double tmp;
//rotates by 90 degrees about axis
    switch(axis)
    {
    case 0: tmp=vec[1];
            vec[1]=-vec[2];
            vec[2]=tmp;
            break;
    case 1: tmp=vec[0];
            vec[0]=vec[2];
            vec[2]=-tmp;
            break;
    case 2: tmp=vec[0];
            vec[0]=-vec[1];
            vec[1]=tmp;
            break;
    default:tmp=0.0;
    }
}

void Reflect(int ax1,int ax2,double* vec)
{
    int i;

    for(i=0;i<3;i++)
    {
        if((i!=ax1)&&(i!=ax2))
            vec[i]=-vec[i];
    }
}

void Invert(double* vec)
{
    vec[0]=-vec[0];
    vec[1]=-vec[1];
    vec[2]=-vec[2];
}

void Reflect45(int ax1,int ax2,double* vec)
{
    double tmp;

    if(ax2>ax1)
    {
        tmp=vec[ax1];
        vec[ax1]=vec[ax2];
        vec[ax2]=tmp;
    }
    else
    {
        tmp=-vec[ax1];
        vec[ax1]=-vec[ax2];
        vec[ax2]=tmp;
    }
}

int GetSymmetryDesignation(int orbid)
{
    double val[3];
    int i;
    double TestPoint[3]={1.0,1.0,1.0};

    val[0]=0.0;
    val[1]=0.0;
    for(i=0;i<NumGBasis;i++)
    {
        val[0]+=FullCoeffs[i][orbid]*Basis[i].Evaluate(TestPoint);
    }
    Invert(TestPoint);
    for(i=0;i<NumGBasis;i++)
    {
        val[1]+=FullCoeffs[i][orbid]*Basis[i].Evaluate(TestPoint);
    }
    if(fabs(val[1]-val[0])>symprec)  //for us, lack of inversion symmetry means p
        return 1;
    val[1]=0.0;
    Reflect(0,1,TestPoint);
    for(i=0;i<NumGBasis;i++)
    {
        val[1]+=FullCoeffs[i][orbid]*Basis[i].Evaluate(TestPoint);
    }
    if(fabs(val[1]-val[0])>symprec)  //for us, lack of mirror symmetry means d
        return 2;
    val[1]=0.0;
    Reflect(0,2,TestPoint);
    for(i=0;i<NumGBasis;i++)
    {
        val[1]+=FullCoeffs[i][orbid]*Basis[i].Evaluate(TestPoint);
    }
    if(fabs(val[1]-val[0])>symprec)  //for us, lack of mirror symmetry means d
        return 2;
//if we're here its either dz2 or s  two rotations will tell
    val[1]=0.0;
    Rotate(0,TestPoint);
    for(i=0;i<NumGBasis;i++)
    {
        val[1]+=FullCoeffs[i][orbid]*Basis[i].Evaluate(TestPoint);
    }
    if(fabs(val[1]-val[0])>symprec)  //for us, lack of rotational symmetry means d
        return 2;
    val[1]=0.0;
    Rotate(1,TestPoint);
    for(i=0;i<NumGBasis;i++)
    {
        val[1]+=FullCoeffs[i][orbid]*Basis[i].Evaluate(TestPoint);
    }
    if(fabs(val[1]-val[0])>symprec)  //for us, lack of rotational symmetry means d
        return 2;
    return 0;                   //totally symmetric is an s
}

void InsertID(char* orb,char* msg)
{
    int i;
    int j=0;

    i=0;
    while(orb[i++]!='(');
    while(msg[j]!='\000')
    {
        orb[i++]=msg[j++];
    }
    orb[i++]=')';
    orb[i]='\000';
}

void GetProperpDesignation(char* orb,int orbid)
{
    int SobDim=-2;
    double ThPh[2];
    double costheta,sintheta,cosphi,sinphi,phi;
    double Psi[NumSobPts];
    double psiav=0.0;
    double r[3];
    double prm[4];
    double weights[NumSobPts];
    int i,j;
    double err;
    double **Count;
    int iprm[4]={1,1,1,1};

    Count=(double**)malloc(NumSobPts*sizeof(double*));
    for(i=0;i<NumSobPts;i++)
    {
        Count[i]=(double*)malloc(sizeof(double));
        Count[i][0]=(double)i;
        weights[i]=1.0;
    }
    sobseq(&SobDim,ThPh);
    SobDim*=-1;
    psiav=-1.0;
    for(i=0;i<NumSobPts;i++)
    {
        sobseq(&SobDim,ThPh);
        costheta=2*ThPh[0]-1.0;
        sintheta=sqrt(1.0-costheta*costheta);
        phi=2.0*pi*ThPh[1];
        cosphi=cos(phi);
        sinphi=sin(phi);
        r[0]=1.5*sintheta*cosphi;
        r[1]=1.5*sintheta*sinphi;
        r[2]=1.5*costheta;
        for(j=-1;j<2;j++)
            Yml[3*i+j+1]=SphHarm(1,j,cosphi,phi);
        Psi[i]=0.0;
        for(j=0;j<NumGBasis;j++)
            Psi[i]+=FullCoeffs[j][orbid]*Basis[j].Evaluate(r);
        if(fabs(Psi[i])>psiav)
            psiav=fabs(Psi[i]);
    }
    prm[0]=psiav;
    prm[1]=prm[2]=0.6;
    prm[3]=0.529150262212918118100323150727852;
    Order=3;
    printf("%g\n",psiav);
    PauseForAction();
/*  err=levmar(Count,Psi,weights,NumSobPts,prm,iprm,4,SHFit,dSHFit);
    psiav=0.0;
    for(j=0;j<3;j++)
        psiav+=prm[i+1]*prm[i+1];
    psiav=1.0/sqrt(psiav);
    for(j=0;j<3;j++)
        psiav*=psiav;
    printf("%g %g %g %g %g\n",prm[0],prm[1],prm[2],prm[3],err);*/
    err=levmar(Count,Psi,weights,NumSobPts,prm,iprm,4,SHFit,dSHFit,SHConstrain,dSHConstrain);
    printf("%g %g %g %g %g\n",prm[0],prm[1],prm[2],prm[3],err);
    PauseForAction();
    for(i=0;i<NumSobPts;i++)
        free((void*)Count[i]);
    free((void*)Count);
}

void GetProperdDesignation(char* orb,int orbid)
{
//  printf("Checking d symmetries\n");
}

void WritePlot(char* fn,int orbid)
{
    double r[3];
    double val;
    float fval;
    FILE* of;
    int i,j,k,n;
    int ival;

    of=fopen(fn,"wb");
    ival=3;
    fwrite(&ival,1,sizeof(int),of);
    ival=25;
    fwrite(&ival,1,sizeof(int),of);
    ival=41;
    fwrite(&ival,1,sizeof(int),of);
    ival=41;
    fwrite(&ival,1,sizeof(int),of);
    ival=41;
    fwrite(&ival,1,sizeof(int),of);
    fval=-8.0;
    fwrite(&fval,1,sizeof(float),of);
    fval=8.0;
    fwrite(&fval,1,sizeof(float),of);
    fval=-8.0;
    fwrite(&fval,1,sizeof(float),of);
    fval=8.0;
    fwrite(&fval,1,sizeof(float),of);
    fval=-8.0;
    fwrite(&fval,1,sizeof(float),of);
    fval=8.0;
    fwrite(&fval,1,sizeof(float),of);
    for(k=0;k<41;k++)
    {
        printf("%i\n",k);
        r[2]=-8.0+8.0*(double)k/20.0;
        for(j=0;j<41;j++)
        {
            r[1]=-8.0+8.0*(double)j/20.0;
            for(i=0;i<41;i++)
            {
                r[0]=-8.0+8.0*(double)i/20.0;
                val=0.0;
                for(n=0;n<NumGBasis;n++)
                {
                    val+=FullCoeffs[n][orbid]*Basis[n].Evaluate(r);
                }
                fval=(float)val;
                fwrite(&fval,1,sizeof(float),of);
            }
        }
    }
    fclose(of);
}
    
 
} // namespace Fchk
} // namespace Plugin
