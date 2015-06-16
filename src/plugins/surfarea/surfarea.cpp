// surfarea.cpp : Defines the entry point for the DLL application.
//
/* Copyright 2002-2003 Kevin J. Boyd and the University of New Orleans.  
Permission granted to distribute and modify this code for personal, educational, and
research use. */
#include <sys/types.h>
#include <inttypes.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tcl.h>
#include <gopenmolext.h>
#include "gullick.h"

#define pi 3.1415926535897932384626433832795
#define TOL 1.e-8
#define ATOL 1.e-5
#define MAXMSTEP 1000
#define MAXASTEP 10000
#define Nn 3
#define Mm 4
#define SQ3O3 0.577350269189625764509148780501957
#define R_CUT_RATIO 0.35
#define MAX_NEST 35

//prototype for function which Tcl will see
static int CalcTriangleAreas(ClientData ,Tcl_Interp *,int ,const char **);
static int SetContourData(ClientData ,Tcl_Interp *,int ,const char **);

DYNEXPORT_C int Surfarea_Init(Tcl_Interp *Interp)
//declare exported command
{
    printf("Creating Tcl Extensions.\n");
Tcl_CreateCommand(Interp,"triangles",CalcTriangleAreas,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"sacontour",SetContourData,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_PkgProvide(Interp,"Sarea","1.0");
    printf("Finished.\n");
return(TCL_OK);
}

/*Make these have file scope so we don't have to pass eleventeen variables all the time*/
double bdotc,bdotd,cdotd;
double bdotn,cdotn,ddotn;
char gh;

/*Contour data */
double cmin[3];
double dx[3];
float* ContData;
int nx,ny,nz;

/*Second contour data */
double cmin2[3];
double dx2[3];
float* ContData2;
int nx2,ny2,nz2;

/*Used to limit depth of recursion in triangel area calculation*/
double P_AREA_CUTOFF=2.e-5;


int SetContourData(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    if(nv<12)
        return TCL_ERROR;
    if(atoi(argv[1])==1)
    {
        cmin[0]=atof(argv[2]);
        dx[0]=(atof(argv[3])-cmin[0]);
        cmin[1]=atof(argv[4]);
        dx[1]=(atof(argv[5])-cmin[1]);
        cmin[2]=atof(argv[6]);
        dx[2]=(atof(argv[7])-cmin[2]);
        nx=atoi(argv[8]);
        ny=atoi(argv[9]);
        nz=atoi(argv[10]);
        ContData=(float*)strtol(argv[11], NULL, 10);
    }
    else
    {
        cmin[0]=atof(argv[2]);
        dx[0]=(atof(argv[3])-cmin[0]);
        cmin[1]=atof(argv[4]);
        dx[1]=(atof(argv[5])-cmin[1]);
        cmin[2]=atof(argv[6]);
        dx[2]=(atof(argv[7])-cmin[2]);
        nx=atoi(argv[8]);
        ny=atoi(argv[9]);
        nz=atoi(argv[10]);
        ContData=(float*)strtol(argv[11], NULL, 10);
    }
    return TCL_OK;
}

void FillDots(double* b,double* c,double* d,double* n)
{
/*These don't change during the fitting, so olnly calculate them once. */

    int i;

    bdotc=bdotd=cdotd=0.0;
    bdotn=cdotn=ddotn=0.0;
    for(i=0;i<3;i++)
    {
        bdotc+=b[i]*c[i];
        bdotd+=b[i]*d[i];
        cdotd+=c[i]*d[i];
        bdotn+=b[i]*n[i];
        cdotn+=c[i]*n[i];
        ddotn+=d[i]*n[i];
    }
}

double CalcOmega(double* lmn)
{
    return lmn[0]*lmn[0]+lmn[1]*lmn[1]+lmn[2]*lmn[2]+2.0*lmn[0]*lmn[1]*bdotc+2.0*lmn[0]*lmn[2]*bdotd
            +2.0*lmn[1]*lmn[2]*cdotd;  //length of the vector l*b+m*c+m*d
}

void CalcTheta(double* lmn,double* r)
{
    double t1,t2,t3;

    r[0]=CalcOmega(lmn)-1.0;  //length of the normal should be one
    r[1]=lmn[0]*bdotn+lmn[1]*cdotn+lmn[2]*ddotn-1.0;  //N.n==1 (calculated normal is the triangle plane normal
    t1=lmn[0]+bdotc*lmn[1]+bdotd*lmn[2];  //N.b
    t2=lmn[0]*bdotc+lmn[1]+lmn[2]*cdotd;  //N.c
    t3=lmn[0]*bdotd+lmn[1]*cdotd+lmn[2];  //N.d
    r[2]=t2/t1-1.0;     
    r[3]=t3/t1-1.0;
}

void CalcDTheta(double* lmn,double* J)
{
/* Calculates the Jacobian for the equations above. */
    double t1,t2,t3;

    J[0]=2.0*(lmn[0]+lmn[1]*bdotc+lmn[2]*bdotd);
    J[1]=2.0*(lmn[1]+lmn[0]*bdotc+lmn[2]*cdotd);
    J[2]=2.0*(lmn[2]+lmn[0]*bdotd+lmn[1]*cdotd);
    J[3]=bdotn;
    J[4]=cdotn;
    J[5]=ddotn;
    t1=lmn[0]+bdotc*lmn[1]+bdotd*lmn[2];
    t2=lmn[0]*bdotc+lmn[1]+lmn[2]*cdotd;
    t3=lmn[0]*bdotd+lmn[1]*cdotd+lmn[2];
    J[6]=bdotc/t1-t2/(t1*t1);
    J[7]=1.0/t1-t2*bdotc/(t1*t1);
    J[8]=cdotd/t1-t2*bdotd/(t1*t1);
    J[9]=bdotd/t1-t3/(t1*t1);
    J[10]=cdotd/t1-t3*bdotc/(t1*t1);
    J[11]=1.0/t1-t3*bdotd/(t1*t1);
}

double steplength(double *lmn,double *dir,double *wei)
{
//Determine the steplength in the Gauss-Newton evaluation.  This is very primitive.
    int i;
    double alpha=1.0;
    double rv;
    double rv2;
    double xx[(Nn)];
    double fv[(Mm)];

    CalcTheta(lmn,fv);
    rv=0.0;
    for(i=0;i<(Mm);i++)
    {
        rv+=wei[i]*fv[i]*fv[i];
    }
    while(alpha>=0.3)
    {
        for(i=0;i<(Nn);i++)
        {
            xx[i]=lmn[i]+alpha*dir[i];
        }
        CalcTheta(xx,fv);
        rv2=0.0;
        for(i=0;i<(Mm);i++)
        {
            rv2+=wei[i]*fv[i]*fv[i];
        }
        if(((rv-rv2)>0.4)||(fabs((rv-rv2)/rv)>0.1))
            return alpha;
        else
            alpha*=0.5;
    }
    return alpha;
}

void GenGuess(double* b,double* c,double *d,double* n,double* xg)
{
/*find the values of l,m, and n which make the normal exactly the normal to the triangle plane.  This is the starting
    point for the least squares solution. */
    double M[9];
    double sv[3];
    int pc[3];
    int i;
    double cond;

    for(i=0;i<3;i++)
    {
        M[3*i]=b[i];
        M[3*i+1]=c[i];
        M[3*i+2]=d[i];
        sv[i]=n[i];
    }
    LUDecomp(M,3,pc);
    cond=M[0]*M[4]*M[8];
    if(fabs(cond)<1.e-6)
    {
        printf("We have a conditioning problem! %g\n",cond);
        xg[0]=1.0;
        xg[1]=0.0;
        xg[2]=0.0;
    }
    else
    {
        LUSolve(M,sv,xg,3,pc);
    }
}

void CalcLMN(double* b,double* c,double* d,double* n,double* nret)
{
/* This calculates the normal at the triangle midpoint.  It uses Gullikson's method for performing constrained least-squares
fits, which is based on a modified QR decomposition.  Somewhere I have a reference for this; one of the other authors is
Per-Ake Wedin.  This isn't their best method, since there were aspects of that which I didn't figure out.  It was published
in one of the SIAM journals.  The fitting is actually performed on the variables (l,m,n), where
    N=l*b+m*c+n*d, 
with b, c, and d being the vertex normals.
*/
    double xp[Nn];  //the fitting variables
    double wei[Mm];
    double mwei[Mm];
    double Jc[(Mm)*(Nn)];
    double Jj[(Mm)*(Nn)];
    double fv[Mm];
    int i,j,k;
    int pc[(Nn)];
    double solvec[(Nn)];
    double qvec[Mm];
    double lambda[Mm];
    double temp;
    double dvec[(Mm)];
    double alpha;
    int p;
    double ht=0.1;

    FillDots(b,c,d,n);
    k=0;
    temp=10.0;
    mwei[0]=0.0;  //This is a constraint
    mwei[1]=1.0;  //equally weight the other three
    mwei[2]=1.0;
    mwei[3]=1.0;
    GenGuess(b,c,d,n,xp);//Initial guess for (l,m,n)
    while((k<3501)&&(temp>1.e-20))
    {
        CalcTheta(xp,fv);
        for(i=0;i<Mm;i++)
            fv[i]*=-1.0;
        CalcDTheta(xp,Jc);
        for(i=0;i<(Nn)*(Mm);i++)
            Jj[i]=Jc[i];
        p=0;    
        p=RPR(Jc,mwei,qvec,pc,Mm,(Nn));
//      printf("step %i  Jp=%g rank=%i\n",k,temp,p);
        if(p<(Nn))
        {
            printf("Regularizing....%i\n",k);
            printf("%g %g %g\n",n[0],n[1],n[2]);
            printf("%g %g %g\n",b[0],b[1],b[2]);
            printf("%g %g %g\n",c[0],c[1],c[2]);
            printf("%g %g %g\n",d[0],d[1],d[2]);
/*          gh=getch();
            for(i=0;i<(Nn)*(Mm);i++)
                Jc[i]=Jj[i];
            Regularize(Jc,fv,solvec,Mm,(Nn));
            temp=0.0;
            for(i=0;i<(Nn);i++)
                temp+=solvec[i]*solvec[i];
            while(temp==0)
            {
                ht*=2.0;
                Regularize(Jc,fv,solvec,Mm,(Nn),ht);
                temp=0.0;
                for(i=0;i<(Nn);i++)
                    temp+=solvec[i]*solvec[i];
            }
            ht=0.1;*/
/*If we're here, we probably have degenerate normals.  We'll cheat here:
adopt an alternate tack and pretend that the "sphere" is really a plane.  Find the normals which are
degenerate and return one of them.  We need to look for a problem later, though*/
            temp=0.0;
            temp+=(b[1]*c[2]-b[2]*c[1])*(b[1]*c[2]-b[2]*c[1]);
            temp+=(b[2]*c[0]-b[0]*c[2])*(b[2]*c[0]-b[0]*c[2]);
            temp+=(b[0]*c[1]-b[1]*c[0])*(b[0]*c[1]-b[1]*c[0]);
            if(temp<1.e-18)
            {
                xp[0]=1.0;
                xp[1]=0.0;
                xp[2]=0.0;
            }
            else
            {
                xp[2]=1.0;
                xp[0]=0.0;
                xp[1]=0.0;
            }
            solvec[0]=0.0;
            solvec[1]=0.0;
            solvec[2]=0.0;
        }
        else
        {
            RPRSolve(Jc,solvec,fv,mwei,pc,qvec,lambda,Mm,(Nn),1);
        }
        temp=0.0;
        for(i=0;i<Mm;i++)
        {
            dvec[i]=0.0;
            for(j=0;j<(Nn);j++)
            {
                dvec[i]+=Jj[(Nn)*i+j]*solvec[j];
            }
            temp+=dvec[i]*dvec[i];
        }
        alpha=steplength(xp,solvec,wei);
        for(i=0;i<(Nn);i++)
        {
            xp[i]+=alpha*solvec[i];
        }
        k+=1;
    }
    CalcTheta(xp,fv);
    printf("%g %g %g %g\n",fv[0],fv[1],fv[2],fv[3]);
    temp=0.0;
    for(i=0;i<3;i++)
    {
        nret[i]=xp[0]*b[i]+xp[1]*c[i]+xp[2]*d[i];
        temp+=nret[i]*nret[i];
    }
    temp=sqrt(temp);
    if(temp>1.1)
    {
        printf("Problem in Normal Calculation.");
    }
    for(i=0;i<3;i++)
    {
        nret[i]/=temp;
    }
}


double MAX3(double a,double b,double c)
{
    double r1,r2,r3;

    r1=fabs(a/b-1.0);
    r2=fabs(b/c-1.0);
    r3=fabs(c/a-1.0);
    if(r1>r2)
    {
        if(r1>r3)
            return r1;
        else
            return r3;
    }
    else
    {
        if(r2>r3)
            return r2;
        else
            return r3;
    }
}

void FindCircumcenter(double** vertex,double* midpt)
{   
    int cnvrg=0;
    int i,j;
    int stepcnt;
    double r[3];
    double dotp,nlen;
    double ww[3];
    double mx3;

    stepcnt=0;
    while(!cnvrg)  //find the midpoint of the triangle
    {
        for(i=0;i<3;i++)
        {
            r[i]=0.0;
            for(j=0;j<3;j++)
            {
                r[i]+=(midpt[j]-vertex[i][j])*(midpt[j]-vertex[i][j]);
            }
            r[i]=sqrt(r[i]);
            if(r[i]>3.0)
                cnvrg=0;
        }
        cnvrg=((fabs(r[0]/r[1]-1.0)<TOL)&&(fabs(r[1]/r[2]-1.0)<TOL)&&(fabs(r[2]/r[0]-1.0)<TOL));
        if(!cnvrg)
        {
            dotp=(r[0]+r[1]+r[2])/3.0;
            nlen=0.0;
            for(i=0;i<3;i++)
            {
                ww[i]=dotp-r[i];
                nlen+=ww[i]*ww[i];
            }
            mx3=MAX3(r[0],r[1],r[2]);
            if(mx3>1.0)
                mx3=1.0;
            nlen=0.5*mx3/sqrt(nlen);
            for(j=0;j<3;j++)
            {
                for(i=0;i<3;i++)
                {
                    midpt[j]+=ww[i]*(midpt[j]-vertex[i][j])*nlen;
                }
            }
        }
        if(stepcnt++>MAXMSTEP)
            cnvrg=1;
    }
}

double InterpolateContourValue(double* x)
{
    int xi[3];
    double ddx[3];
    double dv[8];
    int i;
    int offset;

    for(i=0;i<3;i++)
    {
        xi[i]=(int)((x[i]-cmin[i])/dx[i]);
        ddx[i]=x[i]-(xi[i]*dx[i])-cmin[i];
    }
    offset=xi[0]+nx*(xi[1]+ny*xi[2]);
    dv[0]=(double)ContData[offset];
    ++xi[2];
    offset=xi[0]+nx*(xi[1]+ny*xi[2]);
    dv[1]=(double)ContData[offset];
    --xi[2];++xi[1];
    offset=xi[0]+nx*(xi[1]+ny*xi[2]);
    dv[2]=(double)ContData[offset];
    ++xi[2];
    offset=xi[0]+nx*(xi[1]+ny*xi[2]);
    dv[3]=(double)ContData[offset];
    --xi[2];--xi[1];++xi[0];
    offset=xi[0]+nx*(xi[1]+ny*xi[2]);
    dv[4]=(double)ContData[offset];
    ++xi[2];
    offset=xi[0]+nx*(xi[1]+ny*xi[2]);
    dv[5]=(double)ContData[offset];
    --xi[2];++xi[1];
    offset=xi[0]+nx*(xi[1]+ny*xi[2]);
    dv[6]=(double)ContData[offset];
    ++xi[2];
    offset=xi[0]+nx*(xi[1]+ny*xi[2]);
    dv[7]=(double)ContData[offset];
    return dv[0]+dv[1]*ddx[2]*ddx[2]+dv[2]*ddx[1]*ddx[1]+dv[3]*ddx[1]*ddx[2]
        +dv[4]*ddx[0]*ddx[0]+dv[5]*ddx[0]*ddx[2]+dv[6]*ddx[0]*ddx[1]+dv[7]*ddx[0]*ddx[1]*ddx[2];
}

double InterpolateAveragedValue(double* x)
{
    int xi[3];
    double ddx[3];
    double dv[8];
    int i;
    int offset;

    for(i=0;i<3;i++)
    {
        xi[i]=(int)((x[i]-cmin2[i])/dx2[i]);
        ddx[i]=x[i]-(xi[i]*dx2[i])-cmin2[i];
    }
    offset=xi[0]+nx2*(xi[1]+ny2*xi[2]);
    dv[0]=(double)ContData2[offset];
    ++xi[2];
    offset=xi[0]+nx2*(xi[1]+ny2*xi[2]);
    dv[1]=(double)ContData2[offset];
    --xi[2];++xi[1];
    offset=xi[0]+nx2*(xi[1]+ny2*xi[2]);
    dv[2]=(double)ContData2[offset];
    ++xi[2];
    offset=xi[0]+nx2*(xi[1]+ny2*xi[2]);
    dv[3]=(double)ContData2[offset];
    --xi[2];--xi[1];++xi[0];
    offset=xi[0]+nx2*(xi[1]+ny2*xi[2]);
    dv[4]=(double)ContData2[offset];
    ++xi[2];
    offset=xi[0]+nx2*(xi[1]+ny2*xi[2]);
    dv[5]=(double)ContData2[offset];
    --xi[2];++xi[1];
    offset=xi[0]+nx2*(xi[1]+ny2*xi[2]);
    dv[6]=(double)ContData2[offset];
    ++xi[2];
    offset=xi[0]+nx2*(xi[1]+ny2*xi[2]);
    dv[7]=(double)ContData2[offset];
    return dv[0]+dv[1]*ddx[2]*ddx[2]+dv[2]*ddx[1]*ddx[1]+dv[3]*ddx[1]*ddx[2]
        +dv[4]*ddx[0]*ddx[0]+dv[5]*ddx[0]*ddx[2]+dv[6]*ddx[0]*ddx[1]+dv[7]*ddx[0]*ddx[1]*ddx[2];
}

void FindApex(double* start,double* dir,double* result,double threshold)
{
    double alpha=0.0;
    double dalpha;
    int cnvrg=0;
    double cv,cv2;
    double dr;
    int i;


    dalpha=0.001;
    while(!cnvrg)
    {
        for(i=0;i<3;i++)
            result[i]=start[i]+alpha*dir[i];
        cv=InterpolateContourValue(result);
        dalpha=0.001*dx[0];
        cv2=cv;
        while(fabs(cv2-cv)<1.e-9)
        {
            for(i=0;i<3;i++)
                result[i]=start[i]+(alpha+dalpha)*dir[i];
            cv2=InterpolateContourValue(result);
            if(fabs(cv2-cv)<1.e-9)
                dalpha*=2.0;
        }
        dr=(cv2-cv)/dalpha;
        dalpha=(cv-threshold)/dr;
//      printf("%g %g %g %g %g\n",dr,dalpha,alpha,cv,threshold);
        if((fabs(dalpha/alpha)<1.e-14)||(fabs(cv-threshold)<1.e-14))
            cnvrg=1;
        else
            alpha-=dalpha;
    }
    for(i=0;i<3;i++)
        result[i]=start[i]+alpha*dir[i];
//  printf("%g %g %g\n",result[0],result[1],result[2]);
}

int nestlevel=0;

void TriangleStuff(double **vertex,double **normals,char* retstr,int Avg)
{
    double midpt[3];
    double edge1[3],edge2[3];
    double midptnorm[3];
    double rmnorm[3];
    double r[3];
    double nlen;
    double dotp;
    int i,j,k;
    int cnvrg=0;
    double sa,sb,sg,ca,cb,cg,tb;
    double sA,sB,sC,cA,cB,cC,tB;
    double A,B,C;
    double ww[3];
    double ravg;
    double parea,sarea;
    double h[3];
    double havg;
    double qbar,delq,delr;
    double dA;
    double minr,maxr;
    double shv[3][3],shn[3][3];
    double **nn;
    double **nv;
    double nmidpt[3];
    double apex[3];
    double threshold;
    double qval,qv;
    double pa,sfa;
    double dda[3];
    double nt[3];
    FILE* of;

    if(!nestlevel)
    {
        of=fopen("trigfile.txt","w");
        for(i=0;i<3;i++)
            fprintf(of,"%g %g %g\n",vertex[i][0],vertex[i][1],vertex[i][2]);
        fclose(of);
    }
    ++nestlevel;
    for(i=0;i<3;i++)
    {
        edge1[i]=vertex[1][i]-vertex[0][i];
        edge2[i]=vertex[2][i]-vertex[0][i];
        nlen=0.0;
        midpt[i]=0.0;
        for(j=0;j<3;j++)
        {
            nlen+=normals[i][j]*normals[i][j];
            midpt[i]+=vertex[j][i];
        }
        midpt[i]/=3.0;
        nlen=1.0/sqrt(nlen);
        for(j=0;j<3;j++)
            normals[i][j]*=nlen;
    }
    rmnorm[0]=ww[0]=edge1[1]*edge2[2]-edge1[2]*edge2[1];  //find the normal to the triangle
    rmnorm[1]=ww[1]=edge1[2]*edge2[0]-edge1[0]*edge2[2];
    rmnorm[2]=ww[2]=edge1[0]*edge2[1]-edge1[1]*edge2[0];
    parea=0.5*sqrt(ww[0]*ww[0]+ww[1]*ww[1]+ww[2]*ww[2]);
    if(nestlevel==1)
        P_AREA_CUTOFF=parea*0.1;
    for(i=0;i<3;i++)
        rmnorm[i]*=0.5/parea;  //normalize the normal
    k=0;
    for(i=0;i<3;i++)
    {
        dotp=0.0;
        for(j=0;j<3;j++)
            dotp+=rmnorm[j]*normals[i][j];
//      printf("%g ",dotp);
        if(dotp<0.0)
            ++k;
    }
    if(k>1)
    {
        for(i=0;i<3;i++)
        {
            rmnorm[i]*=-1.0;
        }
    }
//  printf("\n");
    FindCircumcenter(vertex,midpt);
    CalcLMN(normals[0],normals[1],normals[2],rmnorm,midptnorm);  //find the LS normal at the midpoint
/*  qbar=(r[0]+r[1]+r[2])/3.0;  //this is now nonsense.  Leave it here to be improved later.
    for(i=0;i<3;i++)
        ww[i]=r[i]-qbar;
    qbar=sqrt(1.0-qbar);
    delq=sqrt((ww[0]*ww[0]+ww[1]*ww[1]+ww[2]*ww[2])*0.5);
    delq*=0.5/qbar;*/
    delq=0.0;
//  printf("%g %g %g\n",midptnorm[0],midptnorm[1],midptnorm[2]);
    for(i=0;i<3;i++)
    {
        dotp=0.0;
        nlen=0.0;
        for(j=0;j<3;j++)
        {
            nlen+=(midpt[j]-vertex[i][j])*(midpt[j]-vertex[i][j]);
            dotp+=normals[i][j]*midptnorm[j];
        }
        nlen=sqrt(nlen);
        r[i]=nlen/sqrt(1-dotp*dotp);
        h[i]=sqrt(nlen*nlen+r[i]*dotp*r[i]*dotp)-r[i]*dotp;
//      printf("%g %g %g %g\n",nlen,dotp,r[i],h[i]);
    }
    maxr=minr=r[0];
    for(i=1;i<3;i++)
    {
        if(r[i]>maxr)maxr=r[i];
        if(r[i]<minr)minr=r[i];
    }
    ravg=(r[0]+r[1]+r[2])/3.0;  //use ravg.  We could also use the maximum r.
    delr=0.0;
    for(i=0;i<3;i++)
        delr+=(r[i]-ravg)*(r[i]-ravg);
    delr=sqrt(0.5*delr);    //this is an error estimate for the radius.  We could also use max-min.
    if(!Avg)
    {
        havg=(h[0]+h[1]+h[2])/3.0;
//      printf("%g %g %g %g\n",h[0],h[1],h[2],havg);
        for(i=0;i<3;i++)
            nmidpt[i]=midpt[i]-(ravg-havg)*midptnorm[i];   //translate to put the center of the sphere at the origin
    }
    else
    {
        threshold=InterpolateContourValue(vertex[0]);
        FindApex(midpt,midptnorm,h,threshold);
        qval=InterpolateAveragedValue(apex);
        for(i=0;i<3;i++)
            nmidpt[i]=h[i]-ravg*midptnorm[i];
    }
    for(i=0;i<3;i++)
    {
        nlen=0.0;
        for(j=0;j<3;j++)
        {
            shv[i][j]=vertex[i][j]-nmidpt[j];
            shn[i][j]=shv[i][j];        //if the center of the sphere is at the origin, the normal is parallel to the position
            nlen+=shn[i][j]*shn[i][j];
        }
//      printf("%g\n",sqrt(nlen));
        nlen=1.0/sqrt(nlen);
        for(j=0;j<3;j++)
            shn[i][j]*=nlen;
    }
    nlen=shn[0][0]*(shn[1][1]*shn[2][2]-shn[2][1]*shn[1][2])
        -shn[0][1]*(shn[1][0]*shn[2][2]-shn[2][0]*shn[1][2])
        +shn[0][2]*(shn[1][0]*shn[2][1]-shn[1][1]*shn[2][0]);  //this may be legacy
    cg=0.0;ca=0.0;cb=0.0;
    for(j=0;j<3;j++)
    {
        ca+=shn[0][j]*shn[1][j];
        cb+=shn[1][j]*shn[2][j];
        cg+=shn[2][j]*shn[0][j];
    }
    if((((maxr-minr)/minr)<R_CUT_RATIO)&&((1.0-ca*ca)>0.0)&&((1.0-cb*cb)>0.0)&&((1.0-cg*cg)>0.0))  //make sure that none of the angles is =0.  
    {
        sa=sqrt(1.0-ca*ca);
        sb=sqrt(1.0-cb*cb);
        sg=sqrt(1.0-cg*cg);
        tb=cb/sb;
        cA=(ca-cb*cg)/(sb*sg);
        sA=sqrt(1.0-cA*cA);
        cB=(cb*sg-sb*cg*cA)/sa;
        sB=sqrt(1.0-cB*cB);
        tB=cB/sB;
        sC=sA/sa*sg;
        if(ca!=0.0)
            cC=(sa*tb-sC*tB)/ca;
        else
            cC=sqrt(1.0-sC*sC);
        A=atan2(sA,cA);
        B=atan2(sB,cB);
        C=atan2(sC,cC);
//  printf("%g %g %g %g\n",A,B,C,A+B+C-pi);
//  printf("%g %g %g\n",r[0]*r[0]*(A+B+C-pi),r[1]*r[1]*(A+B+C-pi),r[2]*r[2]*(A+B+C-pi));
        sarea=(A+B+C-pi)*ravg*ravg; 
//  printf("%g %g\n",delq,delr);
        dA=sqrt(ravg*ravg*(A+B+C-pi)*(A+B+C-pi)*delr*delr*4.0+ravg*ravg*ravg*ravg*delq*delq);
        qval*=sarea;
/*      printf("%g %g %g\n",vertex[0][0],vertex[0][1],vertex[0][2]);
        printf("%g %g %g\n",vertex[1][0],vertex[1][1],vertex[1][2]);
        printf("%g %g %g\n",vertex[2][0],vertex[2][1],vertex[2][2]);
        printf("%g %g %g\n",midpt[0],midpt[1],midpt[2]);
        printf("%g %g %g\n",normals[0][0],normals[0][1],normals[0][2]);
        printf("%g %g %g\n",normals[1][0],normals[1][1],normals[1][2]);
        printf("%g %g %g\n",normals[2][0],normals[2][1],normals[2][2]);
        printf("%g %g %g\n",midptnorm[0],midptnorm[1],midptnorm[2]);
        printf("%g %g %g %g\n",r[0],r[1],r[2],ravg);
        printf("%g %g %g %g\n",A,B,C,(A+B+C-pi));*/
        printf("%i Areas: %g %g %g %g\n",nestlevel,parea,sarea,dA,(sarea-parea)/sarea*100);
//      gh=getch();
    }
    else
    {
//if one of the angles is 0, then the "sphere" is really a plane.
        if(parea<P_AREA_CUTOFF)
        {
            sarea=parea;
            dA=0.0;
            printf("Small triangle threshold reached.  Area=%g\n",parea);
        }
        else if(nestlevel>=MAX_NEST)
        {
            sarea=parea;
            dA=0.0;
            printf("Maximum depth of mesh refinement reached.  Area=%g\n",parea);
        }
        else if((maxr-minr)/minr<=R_CUT_RATIO)
        {
            printf("Degenerate spherical triangle! ");
            if((1.0-ca*ca)>0.0)
                printf("Angle a=0\n");
            else if ((1.0-ca*ca)>0.0)
                printf("Angle b=0\n");
            else
                printf("Angle c=0\n");
            sarea=parea;
            dA=0.0;
        }
        else
        {
            printf("Increasing mesh Area=%g\n",parea);
//          printf("Radii inconsistent by %g (%g %g)\n",(maxr-minr)/minr,maxr,minr);
            nv=(double**)malloc(3*sizeof(double*));
            nn=(double**)malloc(3*sizeof(double*));
            for(i=0;i<3;i++)
            {
                nv[i]=(double*)malloc(3*sizeof(double));
                nn[i]=(double*)malloc(3*sizeof(double));
            }
            threshold=InterpolateContourValue(vertex[0]);
            FindApex(midpt,midptnorm,apex,threshold);
            of=fopen("trigfile.txt","a");
            fprintf(of,"%g %g %g\n",apex[0],apex[1],apex[2]);
            fclose(of);
            for(i=0;i<3;i++)
            {
                edge1[i]=apex[i]-vertex[0][i];
                edge2[i]=apex[i]-vertex[1][i];
            }
            nt[0]=edge1[1]*edge2[2]-edge1[2]*edge2[1];
            nt[1]=edge1[2]*edge2[0]-edge1[0]*edge2[2];
            nt[2]=edge1[0]*edge2[1]-edge1[1]*edge2[1];
            dotp=0.0;
            nlen=0.0;
            for(i=0;i<3;i++)
            {
                dotp+=nt[i]*midptnorm[i];
                nlen+=nt[i]*nt[i];
            }
            nlen=1.0/sqrt(nlen);
            for(i=0;i<3;i++)
            {
                if(dotp<0.0)
                    nn[0][i]=-nt[i]*nlen;
                else
                    nn[0][i]=nt[i]*nlen;
            }
            for(i=0;i<3;i++)
            {
                edge1[i]=apex[i]-vertex[0][i];
                edge2[i]=apex[i]-vertex[2][i];
            }
            nt[0]=edge1[1]*edge2[2]-edge1[2]*edge2[1];
            nt[1]=edge1[2]*edge2[0]-edge1[0]*edge2[2];
            nt[2]=edge1[0]*edge2[1]-edge1[1]*edge2[1];
            dotp=0.0;
            nlen=0.0;
            for(i=0;i<3;i++)
            {
                dotp+=nt[i]*midptnorm[i];
                nlen+=nt[i]*nt[i];
            }
            nlen=1.0/sqrt(nlen);
            for(i=0;i<3;i++)
            {
                if(dotp<0.0)
                    nn[0][i]-=nt[i]*nlen;
                else
                    nn[0][i]+=nt[i]*nlen;
            }
            for(i=0;i<3;i++)
            {
                edge1[i]=apex[i]-vertex[1][i];
                edge2[i]=apex[i]-vertex[2][i];
            }
            nt[0]=edge1[1]*edge2[2]-edge1[2]*edge2[1];
            nt[1]=edge1[2]*edge2[0]-edge1[0]*edge2[2];
            nt[2]=edge1[0]*edge2[1]-edge1[1]*edge2[1];
            dotp=0.0;
            nlen=0.0;
            for(i=0;i<3;i++)
            {
                dotp+=nt[i]*midptnorm[i];
                nlen+=nt[i]*nt[i];
            }
            nlen=1.0/sqrt(nlen);
            for(i=0;i<3;i++)
            {
                if(dotp<0.0)
                    nn[0][i]-=nt[i]*nlen;
                else
                    nn[0][i]+=nt[i]*nlen;
            }
            for(i=0;i<3;i++)
            {
                nv[0][i]=apex[i];
                nv[1][i]=vertex[0][i];
                nv[2][i]=vertex[1][i];
                nn[1][i]=normals[0][i];
                nn[2][i]=normals[1][i];
            }
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf %lf",&parea,&sarea,dda,&qv);
            qval+=qv*sarea;
            for(i=0;i<3;i++)
            {
                nv[1][i]=vertex[0][i];
                nv[2][i]=vertex[2][i];
                nn[1][i]=normals[0][i];
                nn[2][i]=normals[2][i];
            }
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf",&pa,&sfa,dda+1,&qv);
            qval+=qv*sfa;
            parea+=pa;sarea+=sfa;
            for(i=0;i<3;i++)
            {
                nv[1][i]=vertex[1][i];
                nv[2][i]=vertex[2][i];
                nn[1][i]=normals[1][i];
                nn[2][i]=normals[2][i];
            }
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf",&pa,&sfa,dda+1,qv);
            qval+=qv*sfa;
            parea+=pa;sarea+=sfa;
            dA=sqrt(dda[0]*dda[0]+dda[1]*dda[1]+dda[2]*dda[2]);
            for(i=0;i<3;i++)
            {
                free((void*)nv[i]);
                free((void*)nn[i]);
            }
            free((void*)nv);
            free((void*)nn);
        }
    }
    --nestlevel;
    printf("%i\n",nestlevel);
    sprintf(retstr,"%g %g %g %g %g\n",parea,sarea,(sarea-parea)/sarea*100,dA,qval);//return our result
}

int CalcTriangleAreas(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    char tstr[80];
    double **vertices;
    double **normals;
    int i;
    int Avg;


    if(nv<19)
        return TCL_ERROR;  //not enough stuff came in so bail
    Avg=atoi(argv[19]);
    vertices=(double**)malloc(3*sizeof(double*));//set up the arrays
    normals=(double**)malloc(3*sizeof(double*));
    for(i=0;i<3;i++)
    {
        vertices[i]=(double*)malloc(3*sizeof(double));
        normals[i]=(double*)malloc(3*sizeof(double));
    }
    vertices[0][0]=atof(argv[1]);  //fill the arrays
    vertices[0][1]=atof(argv[2]);
    vertices[0][2]=atof(argv[3]);
    vertices[1][0]=atof(argv[4]);
    vertices[1][1]=atof(argv[5]);
    vertices[1][2]=atof(argv[6]);
    vertices[2][0]=atof(argv[7]);
    vertices[2][1]=atof(argv[8]);
    vertices[2][2]=atof(argv[9]);
    normals[0][0]=atof(argv[10]);
    normals[0][1]=atof(argv[11]);
    normals[0][2]=atof(argv[12]);
    normals[1][0]=atof(argv[13]);
    normals[1][1]=atof(argv[14]);
    normals[1][2]=atof(argv[15]);
    normals[2][0]=atof(argv[16]);
    normals[2][1]=atof(argv[17]);
    normals[2][2]=atof(argv[18]);
    TriangleStuff(vertices,normals,tstr,Avg);//calculate the results; the result is listed in tstr
    for(i=0;i<3;i++)
    {
        free((void*)vertices[i]);  //give back the memory
        free((void*)normals[i]);
    }
    free((void*)vertices);
    free((void*)normals);
    Tcl_SetResult(ti,tstr,TCL_VOLATILE);  //return the result to Tcl
    return TCL_OK;
}
