// surfarea.cpp : Defines the entry point for the DLL application.
//
/* Copyright 2002-2003 Kevin J. Boyd and the University of New Orleans.  
Permission granted to distribute and modify this code for personal, educational, and
research use. */

#include <tcl.h>
#include <gopenmolext.h>
#include "bunch.h"
#include "gullick.h"

namespace Plugin {
namespace Surfarea {
    
/*Numerical constants*/
#define pi 3.1415926535897932384626433832795
#define SQ3O3 0.577350269189625764509148780501957
#define OneThird 0.3333333333333333333333333333333333333333
/*Convergence criteria*/
#define TOL 1.e-8
#define ATOL 1.e-5
#define MAXMSTEP 1000
#define MAXASTEP 10000
/*System size for NLLS*/
#define Nn 3
#define Mm 5
/*Constants for Integration*/
#define A0 0.797426985353087322398025276170
#define A1 0.10286507323456338800987361915
#define A2 0.3333333333333333333333333333333333333333
#define A3 0.470142064105115089770441209513
#define A4 0.059715871789769520459117580973
#define CRUDE 1
/*Output flags for debugging*/
#define WriteSomeTriangleDetails 1

//prototype for function which Tcl will see
static int CalcTriangleAreas(ClientData ,Tcl_Interp *,int ,const char **);
static int MakeEllipsoidalMesh(ClientData ,Tcl_Interp *,int ,const char **);
static int SetTriangleParms(ClientData ,Tcl_Interp *,int ,const char **);
static int SetContourData(ClientData ,Tcl_Interp *,int ,const char **);
static int ResetTriangleCount(ClientData, Tcl_Interp *,int,const char **);
static int InterpolateMesh(ClientData, Tcl_Interp *,int,const char **);
static int SumMesh(ClientData, Tcl_Interp *,int,const char **);
static int DiffMesh(ClientData, Tcl_Interp *,int,const char **);
static int ProdMesh(ClientData, Tcl_Interp *,int,const char **);
static int QuotMesh(ClientData, Tcl_Interp *,int,const char **);

DYNEXPORT_C int Surfarea_Init(Tcl_Interp *Interp)
//declare exported command
{
    printf("Creating Tcl Extensions.\n");
Tcl_CreateCommand(Interp,"triangles",CalcTriangleAreas,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"ellipse",MakeEllipsoidalMesh,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"trigparm",SetTriangleParms,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"sacontour",SetContourData,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"resettrig",ResetTriangleCount,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"interpmesh",InterpolateMesh,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"addmesh",SumMesh,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"submesh",DiffMesh,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"mulmesh",ProdMesh,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"divmesh",QuotMesh,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_PkgProvide(Interp,"Sarea","1.0");
    printf("Finished.\n");
return(TCL_OK);
}

/*Make these have file scope so we don't have to pass eleventeen variables all the time*/
double bdotc,bdotd,cdotd;
double bdotn,cdotn,ddotn;
char gh;

/*Contour data */
double contval;
double cmin[3];
double dx[3];
float* ContData;
int nx[3];

/*Second contour data */
double cmin2[3];
double dx2[3];
float* ContData2;
int nx2[3];

/*Used to limit depth of recursion in triangle area calculation*/
double P_AREA_CUTOFF=2.e-7; //minimum triangle size
double P_CUT_SCALE=0.1;     //scales triangle area to initial triangle size
double R_CUT_RATIO=0.35;    //how close the calculated radii must be to accept the sphere
int    MAX_NEST=3;          //maximum recursion depth

int SetContourData(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
/*Load the contour values into the places we want them for interpolation.  If the first value
is one, its the surface contour data, if not its the contour to be averaged.*/
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
        nx[0]=atoi(argv[8]);
        nx[1]=atoi(argv[9]);
        nx[2]=atoi(argv[10]);
        dx[0]/=(double)(nx[0]-1);
        dx[1]/=(double)(nx[1]-1);
        dx[2]/=(double)(nx[2]-1);
        ContData=(float*)atoi(argv[11]);
        if(nv<13)
            contval=1.e-3;
        else
            contval=atof(argv[12]);
    }
    else
    {
        cmin2[0]=atof(argv[2]);
        dx2[0]=(atof(argv[3])-cmin2[0]);
        cmin2[1]=atof(argv[4]);
        dx2[1]=(atof(argv[5])-cmin2[1]);
        cmin2[2]=atof(argv[6]);
        dx2[2]=(atof(argv[7])-cmin2[2]);
        nx2[0]=atoi(argv[8]);
        nx2[1]=atoi(argv[9]);
        nx2[2]=atoi(argv[10]);
        dx2[0]/=(double)(nx2[0]-1);
        dx2[1]/=(double)(nx2[1]-1);
        dx2[2]/=(double)(nx2[2]-1);
        ContData2=(float*)atoi(argv[11]);
    }
    return TCL_OK;
}

void FillDots(double* b,double* c,double* d,double* n)
{
/*These don't change during the fitting, so only calculate them once. */

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
/*The constraint that the length of the normal be unity.*/
    return lmn[0]*lmn[0]+lmn[1]*lmn[1]+lmn[2]*lmn[2]+2.0*lmn[0]*lmn[1]*bdotc+2.0*lmn[0]*lmn[2]*bdotd
            +2.0*lmn[1]*lmn[2]*cdotd;  //length of the vector l*b+m*c+m*d
}

void CalcTheta(double* lmn,double* r)
{
/*Calculates the constraint and fitting equations.*/
    double t1,t2,t3;

    r[0]=CalcOmega(lmn)-1.0;  //length of the normal should be one
    r[1]=(lmn[0]*bdotn+lmn[1]*cdotn+lmn[2]*ddotn)/CalcOmega(lmn)-1.0;  //N.n==1 (calculated normal is the triangle plane normal
    t1=lmn[0]+bdotc*lmn[1]+bdotd*lmn[2];  //N.b
    t2=lmn[0]*bdotc+lmn[1]+lmn[2]*cdotd;  //N.c
    t3=lmn[0]*bdotd+lmn[1]*cdotd+lmn[2];  //N.d
/*  r[2]=t2/t1-1.0;     
    r[3]=t3/t1-1.0;*/
    r[2]=t1-t2;
    r[3]=t1-t3;
    r[4]=t3-t2;
}

void CalcDTheta(double* lmn,double* J)
{
/* Calculates the Jacobian for the equations above. */
    double t1,t2,t3;
    double omega;

    omega=CalcOmega(lmn);
    J[0]=2.0*(lmn[0]+lmn[1]*bdotc+lmn[2]*bdotd);
    J[1]=2.0*(lmn[1]+lmn[0]*bdotc+lmn[2]*cdotd);
    J[2]=2.0*(lmn[2]+lmn[0]*bdotd+lmn[1]*cdotd);
    J[3]=bdotn/omega-(lmn[0]*bdotn+lmn[1]*cdotn+lmn[2]*ddotn)/(omega*omega)*J[0];
    J[4]=cdotn/omega-(lmn[0]*bdotn+lmn[1]*cdotn+lmn[2]*ddotn)/(omega*omega)*J[1];
    J[5]=ddotn/omega-(lmn[0]*bdotn+lmn[1]*cdotn+lmn[2]*ddotn)/(omega*omega)*J[2];
    t1=lmn[0]+bdotc*lmn[1]+bdotd*lmn[2];
    t2=lmn[0]*bdotc+lmn[1]+lmn[2]*cdotd;
    t3=lmn[0]*bdotd+lmn[1]*cdotd+lmn[2];
    J[6]=1.0-bdotc;
    J[7]=bdotc-1.0;
    J[8]=bdotd-cdotd;
    J[9]=1.0-bdotd;
    J[10]=bdotc-cdotd;
    J[11]=bdotd-1.0;
    J[12]=bdotd-bdotc;
    J[13]=cdotd-1.0;
    J[14]=1.0-cdotd;
/*  J[6]=bdotc/t1-t2/(t1*t1);
    J[7]=1.0/t1-t2*bdotc/(t1*t1);
    J[8]=cdotd/t1-t2*bdotd/(t1*t1);
    J[9]=bdotd/t1-t3/(t1*t1);
    J[10]=cdotd/t1-t3*bdotc/(t1*t1);
    J[11]=1.0/t1-t3*bdotd/(t1*t1);*/
}

void CalcDdTheta(double* lmn,double* DDmat)
{
/* Second derivatives are needed only for constraints, not fitting equations.*/
    DDmat[0]=DDmat[4]=DDmat[8]=1.0;
    DDmat[1]=DDmat[3]=bdotc;
    DDmat[2]=DDmat[6]=bdotd;
    DDmat[5]=DDmat[7]=cdotd;
}

double steplength(double *lmn,double *dir,double *J,double *wei)
{
//Determine the steplength in the Gauss-Newton evaluation.  This is very primitive.
    double omega;
    double alpha;
    double s[(Mm)];
    double fv[(Mm)];
    double tfs,sds,fds;
    int i,j;

    CalcTheta(lmn,fv);
    sds=fds=0.0;
    for(i=0;i<(Mm);i++)
    {
        s[i]=0.0;
        for(j=0;j<(Nn);j++)
        {
            s[i]+=J[i*(Nn)+j]*dir[j];
        }
        if(i>0)
        {
            fds+=fv[i]*s[i];
            sds+=fv[i]*s[i];
        }
    }
    if(fds<0.0)
    {
        if((tfs=fv[0]*s[0])>0.0)
            omega=-fds/tfs;
        else
            omega=0.0;
    }
    else
    {
        if((tfs=fv[0]*s[0])<0.0)
            omega=-fds/tfs;
        else
            omega=0.0;
    }
    fds+=omega*s[0]*fv[0];
    sds+=omega*s[0]*s[0];
    alpha=fds/sds;
    if(alpha<0.05)
        return 0.05;
    if(alpha>1.0)
        return 1.0;
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
/*      cond=sqrt(1.0/(3.0+2.0*(bdotc+cdotd+bdotd)));
        xg[0]=xg[1]=xg[2]=cond;*/
    }
}

void FindAverage(double* b,double* c,double *d,double* n,double* nr)
{
/*When all else fails, evenly weight the three vertex normals and use that.*/
    double midpt;
    int i;

    midpt=sqrt(1.0/(3.0+2.0*(bdotc+cdotd+bdotd)));
    for(i=0;i<3;i++)
        nr[i]=midpt*(b[i]+c[i]+d[i]);
}

    

void GetG(double** Gmat,double* vec,double* lambda)
{
/*Calculates the second derivative matrix for the gNR method.*/
    double TGMat[9];
    int i,j;

    CalcDdTheta(vec,TGMat);
    for(i=0;i<3;i++)
    {
        for(j=0;j<3;j++)
        {
            Gmat[(5+i)][5+j]=2.0*lambda[0]*TGMat[3*i+j];
        }
    }
}

void CalcLMN(double* b,double* c,double* d,double* n,double* nret)
{
/* This calculates the normal at the triangle midpoint.  It uses Gullikson's method for performing 
constrained least-squares fits, which is based on a modified QR decomposition.  The fitting is done by
either a Gauss-Newton method (when far from the solution) or a generalized Newton-Raphson method.
Cf. M. Gullikson, I. Soderkvist, P.-A. Wedin, Algorithms for Constrained and Weighted Nonlinear Least 
Squares, SIAM J. Optim. 7 (1997) 208-224.

  The fitting is actually performed on the variables (l,m,n), where
    N=l*b+m*c+n*d, 
with b, c, and d being the vertex normals.

*/
    double xp[Nn];  //the fitting variables
    double wei[Mm]={1000.0,100.0,1.0,1.0,1.0};
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
    double **Ourmat;
    double Ourvec[(Mm+Nn)];
    int Inertia[3];
    double e[(Mm+Nn)];
    int P[(Mm+Nn)];
    double L[(Mm+Nn)*(Mm+Nn)];
    double rho[(Mm+Nn)];
    double dd[(Mm+Nn)];
    int gnr;


    Ourmat=(double**)malloc((Mm+Nn)*sizeof(double*));
    for(i=0;i<(Mm+Nn);i++)
        Ourmat[i]=(double*)malloc((Mm+Nn)*sizeof(double));
    FillDots(b,c,d,n);
    k=0;
    temp=10.0;
    mwei[0]=0.0;  //This is a constraint
    mwei[1]=1.0;  //equally weight the other four
    mwei[2]=1.0;
    mwei[3]=1.0;
    mwei[4]=1.0;
    GenGuess(b,c,d,n,xp);//Initial guess for (l,m,n)
    while((k<51)&&(temp>1.e-8))
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
            CalcTheta(xp,fv);
//          printf("%g %g %g %g %g\n",fv[0],fv[1],fv[2],fv[3],fv[4]);
            gnr=0;
            if(temp<fabs(lambda[0]))
            {
/*If we can try a gNR step, do it here.*/
                for(i=0;i<5;i++)
                {
                    for(j=0;j<5;j++)
                    {
                        if(i==j)
                            Ourmat[i][i]=mwei[i];
                        else
                            Ourmat[i][j]=0.0;
                    }
                    for(j=0;j<3;j++)
                    {
                        Ourmat[i][j+5]=Jj[5*i+j];
                        Ourmat[j+5][i]=Jj[5*i+j];
                        GetG(Ourmat,xp,lambda);
                    }
                    Ourvec[i]=-fv[i];
                }
                for(i=5;i<8;i++)
                    Ourvec[i]=0.0;
/*              for(i=0;i<8;i++)
                {
                    for(j=0;j<8;j++)
                    {
                        printf("%g ",Ourmat[i][j]);
                    }
                    printf("\t%g\n",Ourvec[i]);
                }
                gh=getch();*/
/*The matrix is indefinite.  We use Bunch-Kaufman decomposition, which also gives the inertia of the
matrix.  If there are zero eigenvalues, then we can't use gNR, and must fall through to the GN step.*/
                BunchKaufman(Ourmat,8,Inertia,dd,e,P,L,(Mm+Nn));
                if(Inertia[2]==0)
                {
                    gnr=1;
                    SISolver(P,dd,e,L,rho,Ourvec,(Mm+Nn));
                    for(i=0;i<3;i++)
                        solvec[i]=rho[i+5];
                    alpha=1.0;
//                  printf("Using gNR %g %g %g\n",temp,lambda[0],alpha);
                }
            }
            if(!gnr)
            {
//              printf("Using GN %i %i %i %g %g\n",Inertia[0],Inertia[1],Inertia[2],temp,lambda[0]);
                alpha=steplength(xp,Jj,solvec,wei);
//              printf("Steplength is %g\n",alpha);
            }
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
        for(i=0;i<(Nn);i++)
        {
            xp[i]+=alpha*solvec[i];
        }
        k+=1;
    }
    CalcTheta(xp,fv);
//  printf("%g %g %g %g %g FV\n",fv[0],fv[1],fv[2],fv[3],fv[4]);
    temp=0.0;
    for(i=0;i<3;i++)
    {
        nret[i]=xp[0]*b[i]+xp[1]*c[i]+xp[2]*d[i];
        temp+=nret[i]*nret[i];
    }
    temp=sqrt(temp);
    if(temp>1.1)
    {
//      printf("Taking average....\n");
        FindAverage(b,c,d,n,nret);
//      gh=getch();
    }
    else
    {
        for(i=0;i<3;i++)
        {
            nret[i]/=temp;
        }
    }
    for(i=0;i<(Mm+Nn);i++)
        free((void*)Ourmat[i]);
    free((void*)Ourmat);
}


double MAX3(double a,double b,double c)
{
/*Returns the maximum ratio of sides.*/
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
/*Finds the circumcenter of the triangle.*/
    int cnvrg=0;
    int i,j;
    int stepcnt;
    double r[3];
    double dotp,nlen;
    double ww[3];
    double mx3;
    double edge[2][3];
    double norm[3];
    double bisect[2][3];
    double lmidpt[2][3];

/*From 9th grade geometry, construct the circumcenter as the intersection of the perpendicular
    bisectors of two (or more) sides of the triangle*/
    for(i=0;i<3;i++)
    {
        edge[0][i]=vertex[1][i]-vertex[0][i];
        edge[1][i]=vertex[2][i]-vertex[1][i];
        lmidpt[0][i]=0.5*(vertex[1][i]+vertex[0][i]);
        lmidpt[1][i]=0.5*(vertex[2][i]+vertex[1][i]);
    }
    norm[0]=edge[0][1]*edge[1][2]-edge[0][2]*edge[1][1];
    norm[1]=edge[0][2]*edge[1][0]-edge[0][0]*edge[1][2];
    norm[2]=edge[0][0]*edge[1][1]-edge[0][1]*edge[1][0];
    for(j=0;j<2;j++)
    {
        bisect[j][0]=edge[j][1]*norm[2]-edge[j][2]*norm[1];
        bisect[j][1]=edge[j][2]*norm[0]-edge[j][0]*norm[2];
        bisect[j][2]=edge[j][0]*norm[1]-edge[j][1]*norm[0];
    }
    mx3=(lmidpt[1][0]-lmidpt[0][0])*bisect[0][1]+(lmidpt[0][1]-lmidpt[1][1])*bisect[0][0];
    mx3/=(bisect[1][1]*bisect[0][0]-bisect[1][0]*bisect[0][1]);
    for(i=0;i<3;i++)
    {
        midpt[i]=lmidpt[1][i]+mx3*bisect[1][i];
    }
/*Now test and refine if necessary.*/
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
    if(stepcnt>MAXMSTEP)
        printf("Maximum number of iterations in FindCircumcenter!\n");
}

inline int CDOffset(int xx,int yy,int zz)
{
/*Calculates the offset into the data array*/
    return xx+(nx[0]*(yy+nx[1]*zz));
}

double InterpolateContourValue(double* x,int* aa=NULL,double* dd=NULL)
{
/*Interpolates the value of the contour function at the point x */
    int xi[3];
    double ddx[3];
    double v[8];
    double a[8];
    static double M[64]={1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,
                         1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0,
                         1.0,-1.0,1.0,-1.0,-1.0,1.0,-1.0,1.0,
                         1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,
                         1.0,1.0,1.0,-1.0,1.0,-1.0,-1.0,-1.0,
                         1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0,
                         1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,
                         1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    static int decomp=0;
    static int pc[8];
    int i;

    if(!decomp)
    {
        LUDecomp(M,8,pc);
        decomp=1;
    }
    for(i=0;i<3;i++)
    {
        xi[i]=(int)((x[i]-cmin[i])/dx[i]);
        if((xi[i]<0)||(xi[i]>=nx[i]))
            printf("Bad thing %i %i %g %g %g.\n",xi[i],i,x[i],cmin[i],dx[i]);
        if(aa)
            aa[i]=xi[i];
        ddx[i]=2.0*(((x[i]-(xi[i]*dx[i])-cmin[i])/dx[i])-0.5);
        if(dd)
            dd[i]=ddx[i];
    }
    v[0]=ContData[CDOffset(xi[0],xi[1],xi[2])];
    v[1]=ContData[CDOffset(xi[0]+1,xi[1],xi[2])];
    v[2]=ContData[CDOffset(xi[0],xi[1]+1,xi[2])];
    v[3]=ContData[CDOffset(xi[0],xi[1],xi[2]+1)];
    v[4]=ContData[CDOffset(xi[0]+1,xi[1]+1,xi[2])];
    v[5]=ContData[CDOffset(xi[0]+1,xi[1],xi[2]+1)];
    v[6]=ContData[CDOffset(xi[0],xi[1]+1,xi[2]+1)];
    v[7]=ContData[CDOffset(xi[0]+1,xi[1]+1,xi[2]+1)];
    LUSolve(M,v,a,8,pc);
    return a[0]+ddx[0]*(a[1]+ddx[1]*(a[4]+ddx[2]*a[7])+ddx[2]*a[5])
           +ddx[1]*(a[2]+ddx[2]*a[6])+ddx[2]*a[3];
//  printf("%g %g %g\n",ddx[0],ddx[1],ddx[2]);
}

inline int CDOffset2(int xx,int yy,int zz)
{
/*Calculates the index in the one dimensional data array*/
    return xx+(nx2[0]*(yy+nx2[1]*zz));
}

double InterpolateAveragedValue(double* x)
{
/*Interpolates the value which is being averaged over the surface.*/
    int xi[3];
    double ddx[3];
    double v[8];
    double a[8];
    static double M[64]={1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,
                         1.0,1.0,-1.0,-1.0,-1.0,-1.0,1.0,1.0,
                         1.0,-1.0,1.0,-1.0,-1.0,1.0,-1.0,1.0,
                         1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,
                         1.0,1.0,1.0,-1.0,1.0,-1.0,-1.0,-1.0,
                         1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0,
                         1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,
                         1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    static int decomp=0;
    static int pc[8];
    int i;

    if(!decomp)
        LUDecomp(M,8,pc);
    for(i=0;i<3;i++)
    {
        xi[i]=(int)((x[i]-cmin2[i])/dx2[i]);
        ddx[i]=2.0*(((x[i]-(xi[i]*dx2[i])-cmin2[i])/dx2[i])-0.5);
    }
    v[0]=ContData2[CDOffset2(xi[0],xi[1],xi[2])];
    v[1]=ContData2[CDOffset2(xi[0]+1,xi[1],xi[2])];
    v[2]=ContData2[CDOffset2(xi[0],xi[1]+1,xi[2])];
    v[3]=ContData2[CDOffset2(xi[0],xi[1],xi[2]+1)];
    v[4]=ContData2[CDOffset2(xi[0]+1,xi[1]+1,xi[2])];
    v[5]=ContData2[CDOffset2(xi[0]+1,xi[1],xi[2]+1)];
    v[6]=ContData2[CDOffset2(xi[0],xi[1]+1,xi[2]+1)];
    v[7]=ContData2[CDOffset2(xi[0]+1,xi[1]+1,xi[2]+1)];
    LUSolve(M,v,a,8,pc);
    return a[0]+ddx[0]*(a[1]+ddx[1]*(a[4]+ddx[2]*a[7])+ddx[2]*a[5])
           +ddx[1]*(a[2]+ddx[2]*a[6])+ddx[2]*a[3];
//  printf("%g %g %g\n",ddx[0],ddx[1],ddx[2]);
}


//  triangle counting parameters
int nestlevel=0;
int ntriangles=0;
int exptype;

void RotateDir(double* dir1,double* dir2,int k)
{
    const double F2=0.994987437106619954734479821001206;
    const double rotate[4][9]={{1.0,0.0,0.0,0.0,0.1,F2,0.0,-F2,0.1},{1.0,0.0,0.0,0.0,0.1,-F2,0.0,F2,0.1},
                                {1.0,0.0,0.0,0.0,F2,0.1,0.0,-0.1,F2},{1.0,0.0,0.0,0.0,F2,-0.1,0.0,0.1,F2}};
    int i,j;

    for(i=0;i<3;i++)
    {
        dir2[i]=0.0;
        for(j=0;j<3;j++)
            dir2[i]+=rotate[k][3*i+j]*dir1[j];
    }
}

double FindOpt(double* start,double* dir,double* result,double alphamin,double alphamax)
{
    double alpha=0.0,alphat=0.0;
    double v1,v2,v3,vt;
    double testpt[3];
    const double frac=0.618;
    int i;
    int cnvrg=0;
    double dir2[3];


    v2=fabs(InterpolateContourValue(start)-contval);
    for(i=0;i<3;i++)
        testpt[i]=start[i]+alphamin*dir[i];
    v1=fabs(InterpolateContourValue(testpt)-contval);
    for(i=0;i<3;i++)
        testpt[i]=start[i]+alphamax*dir[i];
    v3=fabs(InterpolateContourValue(testpt)-contval);
    vt=v2;
    while(!cnvrg)
    {
        if(v3<v1)
        {
            alphat=frac*alphamin+(1.0-frac)*alpha;
            for(i=0;i<3;i++)
                testpt[i]=start[i]+alphat*dir[i];
            vt=fabs(InterpolateContourValue(testpt)-contval);
            if(vt<v2)
            {
                alphamax=alpha;
                v3=v2;
                alpha=alphat;
                v2=vt;
            }
            else
            {
                alphamin=alphat;
                v1=vt;
            }
        }
        else
        {
            alphat=frac*alphamax+(1.0-frac)*alpha;
            for(i=0;i<3;i++)
                testpt[i]=start[i]+alphat*dir[i];
            vt=fabs(InterpolateContourValue(testpt)-contval);
            if(vt<v2)
            {
                alphamin=alpha;
                v1=v2;
                alpha=alphat;
                v2=vt;
            }
            else
            {
                alphamax=alphat;
                v3=vt;
            }
        }
//      printf("%g %g %g %g\n",alphamin,alpha-alphamin,alphat-alphamin,alphamax-alphamin);
        cnvrg=(alphamax-alphamin<1.e-14);
        if(!cnvrg)
        {
            if((alphamax-alpha<1.e-14)||(alpha-alphamin<1.e-14))
                alpha=0.5*(alphamax+alphamin);
        }
    }
    for(i=0;i<3;i++)
        result[i]=start[i]+alpha*dir[i];
    return alpha;
}

double WriteApexPath(double* start,double* dir,double alphamin,double alphamax)
{
    static int FirstWritten=0;
    double alpha=alphamin;
    double dalpha=(alphamax-alphamin)*0.001;
    double point[3];
    double result;
    double rdir[3];
    double rresult;
/*  double rsc[3];
    int a[3];*/
    int i;
    int k=0;
    int dirfound=0;
#ifdef WriteUnsuccessfulPaths
    FILE* of;


    if(FirstWritten)
    {
        of=fopen("apexfile.txt","a");
    }
    else
    {
        of=fopen("apexfile.txt","w");
        FirstWritten=1;
    }
    fprintf(of,"%g %g\n",alphamin,alphamax);
#endif
    while(alpha<=(alphamax+0.2*dalpha))
    {
        for(i=0;i<3;i++)
            point[i]=start[i]+alpha*dir[i];
        result=InterpolateContourValue(point);
#ifdef WriteUnsuccessfulPaths
        fprintf(of,"%g %g %i\n",alpha,result,exptype);
#endif
        alpha+=dalpha;
    }
    alpha=FindOpt(start,dir,point,alphamin,alphamax);
    result=InterpolateContourValue(point)-contval;
#ifdef WriteUnsuccessfulPaths
    fprintf(of,"\n%g %g\n\n",alpha,result);
#endif
    while((k<3)&&(!dirfound))
    {
        RotateDir(dir,rdir,k);
        for(i=0;i<3;i++)
            point[i]=start[i]+alpha*rdir[i];
        rresult=InterpolateContourValue(point)-contval;
        if(rresult*result<0.0)
        {
            dirfound=1;
            for(i=0;i<3;i++)
                dir[i]=rdir[i];
#ifdef WriteUnsuccessfulPaths
            fprintf(of,"%g %g\n\n",alpha,rresult);
#endif
        }
        else
        {
            k++;
        }
    }
#ifdef WriteUnsuccessfulPaths
    fclose(of);
#endif
    if(k==3)
        return 0.0;
    return alpha;
}



int FindApex(double* start,double* dir,double* result)
{
/* Given a point and a normal, finds the point along the normal which has the proper contour value.
    Do this by bisection.  A quasi-Newton method might be faster, but there are some continuity issues
    with the derivative which make that not so desirable.  I'm not sure what failed, but this is
    functional.*/
    double alpha=0.0;
    double alphamin=-0.015625,alphamax=0.015625;
    int cnvrg=0;
    double cv=0.0,cv2=0.0,cvt;
    int i;

    for(i=0;i<3;i++)
        result[i]=start[i]+alpha*dir[i];
    cv=InterpolateContourValue(result)-contval;
    while((cv*cv2>=0.0)&&(!cnvrg))
    {
        alpha=alphamin;
        for(i=0;i<3;i++)
            result[i]=start[i]+alpha*dir[i];
        cv2=InterpolateContourValue(result)-contval;
        if(cv*cv2>0.0)
        {
            alpha=alphamax;
            for(i=0;i<3;i++)
                result[i]=start[i]+alpha*dir[i];
            cv2=InterpolateContourValue(result)-contval;
            if(cv*cv2>0.0)
            {
                alphamin*=2.0;
                alphamax*=2.0;
            }
            else
            {
                alphamin=0.0;
            }
        }
        else
        {
            alphamax=0.0;
            cvt=cv2;
            cv2=cv;
            cv=cvt;
        }
//      printf("%g %g\n%g %g\n",alphamax,alphamin,cv,cv2);
        if(alphamax-alphamin>0.5)
        {
            cnvrg=2;
            if((alpha=WriteApexPath(start,dir,alphamin,alphamax))==0.0)
            {
                for(i=0;i<3;i++)
                    result[i]=start[i];
            }
            else
            {
                for(i=0;i<3;i++)
                    result[i]=start[i]+alpha*dir[i];
                if(alpha<0.0)
                {
                    alphamax=0.0;
                    alphamin=2.0*alpha;
                }
                else
                {
                    alphamax=2.0*alpha;
                    alphamin=0.0;
                }
                cnvrg=0;
            }
        }
    }
//  printf("%g %g\n",cv,cv2);
    while(!cnvrg)
    {
        alpha=0.5*(alphamin+alphamax);
        for(i=0;i<3;i++)
            result[i]=start[i]+alpha*dir[i];
        cvt=InterpolateContourValue(result)-contval;
        if(cvt*cv>0.0)
        {
            alphamin=alpha;
            cv=cvt;
        }
        else
        {
            alphamax=alpha;
            cv2=cvt;
        }
        cnvrg=(alphamax-alphamin>1.e-14);
    }
    if(fabs(alpha)>0.25)
    {
        printf("Large alpha= %g\n",alpha);
        gh=getchar();
    }
    return cnvrg-1;
//      printf("%g %g %g %g\n",result[0],result[1],result[2],cv);
}

void WriteTriangleVertex(double* vertex,double* normal)
{
/*  writes triangle vertex information to a text file.  For illustration and debugging.*/
    static int VertexWritten=0;
    FILE* of;

    if(!VertexWritten)
    {
        of=fopen("verfile.txt","w");
        VertexWritten=1;
    }
    else
    {
        of=fopen("verfile.txt","a");
    }
    if((vertex)&&(normal))
    {
        fprintf(of,"%g %g %g\n",vertex[0],vertex[1],vertex[2]);
        fprintf(of,"%g %g %g\n",normal[0],normal[1],normal[2]);
    }
    else
    {
        fprintf(of,"%g %g %g\n",-9.e99,-9.e99,-9.e99);
        fprintf(of,"%g %g %g\n",-9.e99,-9.e99,-9.e99);
    }
    fclose(of);
}

void WriteTriangleDetail(double** vertex,double** normals,double* midpt,double* mpnorm,double rad)
{
/*  Writes the details of a triangle to a text file.  Used for illustration and debugging.*/
    static int TriangleWritten=0;
    FILE* of;
    int i;
    double edge[3];

    if(!TriangleWritten)
    {
        of=fopen("trigfile.txt","w");
        TriangleWritten=1;
    }
    else
    {
        of=fopen("trigfile.txt","a");
    }
    for(i=0;i<3;i++)
    {
        fprintf(of,"%g %g %g\n",vertex[i][0],vertex[i][1],vertex[i][2]);
    }
    fprintf(of,"%g %g %g\n",midpt[0],midpt[1],midpt[2]);
    fprintf(of,"%g\n",rad);
    for(i=0;i<3;i++)
    {
        fprintf(of,"%g %g %g\n",normals[i][0],normals[i][1],normals[i][2]);
    }
    fprintf(of,"%g %g %g\n",mpnorm[0],mpnorm[1],mpnorm[2]);
    for(i=0;i<3;i++)
        edge[i]=vertex[1][i]-vertex[0][i];
    fprintf(of,"%g %g %g\n",edge[0],edge[1],edge[2]);
    for(i=0;i<3;i++)
        edge[i]=vertex[2][i]-vertex[1][i];
    fprintf(of,"%g %g %g\n",edge[0],edge[1],edge[2]);
    for(i=0;i<3;i++)
        edge[i]=vertex[0][i]-vertex[2][i];
    fprintf(of,"%g %g %g\n",edge[0],edge[1],edge[2]);
    fclose(of);
}

double SplitObtuse(double* v1,double* v2,double* v3)
{
/*Calculates the partitioning of the long edge of an obtuse triangle to make two right triangles.
    The long edge is between v1 and v2. */
    double M=0.0;
    double N=0.0;
    double L=0.0;
    double cph=0.0;
    double cga=0.0;
    double l1,l2;
    int i;

    for(i=0;i<3;i++)
    {
        M+=(v3[i]-v1[i])*(v3[i]-v1[i]);
        N+=(v3[i]-v2[i])*(v3[i]-v2[i]);
        L+=(v2[i]-v1[i])*(v2[i]-v1[i]);
        cph+=(v2[i]-v1[i])*(v3[i]-v1[i]);
        cga+=(v3[i]-v2[i])*(v1[i]-v2[i]);
    }
    M=sqrt(M);N=sqrt(N);L=sqrt(L);
    cph/=(M*L);
    cga/=(N*L);
    l2=N*cga;
    l1=M*cph;
//  printf("%g %g %g %g\n",l1,l2,l1+l2,L);
    return l1/(l1+l2);
}

double GetAverage(int Ref,double* apex,double** vertex)
{
    double omega[3]={0.125939180544827151595683945500,0.2250000000000000000000000000,
                        0.132394152788506180737649387833};
    int    nP[3]={3,1,3};
    double P[3][3][2]={{{A0,A1},{A1,A0},{A1,A1}},{{A2,A2},{0.0,0.0},{0.0,0.0}},{{A3,A3},{A3,A4},{A4,A3}}};
    double IntVal=0.0;
    double edge[3][3];
    double edgelength[3];
    double testpt[3];
    int i,j,k;

    if(Ref==CRUDE)
    {
        IntVal=InterpolateAveragedValue(apex);
    }
    else
    {
        for(i=0;i<3;i++)
        {
            edge[0][i]=vertex[1][i]-vertex[0][i];
            edge[1][i]=vertex[2][i]-vertex[1][i];
            edge[2][i]=vertex[0][i]-vertex[2][i];
            for(j=0;j<3;j++)
                edgelength[j]+=edge[j][i]*edge[j][i];
        }
        if(edgelength[0]>=edgelength[1])
        {
            if(edgelength[0]>=edgelength[2])
            {
                for(i=0;i<3;i++)
                {
                    edge[0][i]=vertex[1][i]-vertex[2][i];
                    edge[1][i]=vertex[0][i]-vertex[2][i];
                    edge[2][i]=vertex[2][i];
                }
            }
            else
            {
                for(i=0;i<3;i++)
                {
                    edge[0][i]=vertex[2][i]-vertex[1][i];
                    edge[1][i]=vertex[0][i]-vertex[1][i];
                    edge[2][i]=vertex[1][i];
                }
            }
        }
        else
        {
            if(edgelength[1]>=edgelength[2])
            {
                for(i=0;i<3;i++)
                {
                    edge[0][i]=vertex[1][i]-vertex[0][i];
                    edge[1][i]=vertex[2][i]-vertex[0][i];
                    edge[2][i]=vertex[0][i];
                }
            }
            else
            {
                for(i=0;i<3;i++)
                {
                    edge[0][i]=vertex[2][i]-vertex[1][i];
                    edge[1][i]=vertex[0][i]-vertex[1][i];
                    edge[2][i]=vertex[1][i];
                }
            }
        }
        IntVal=0.0;
        for(i=0;i<3;i++)
        {
            for(j=0;j<nP[i];j++)
            {
                for(k=0;k<3;k++)
                {
                    testpt[k]=edge[2][k]+P[i][j][0]*edge[0][k]+P[i][j][1]*edge[1][k];
                }
                IntVal+=omega[i]*InterpolateAveragedValue(testpt);
            }
        }
    }
    return IntVal;
}

void TriangleStuff(double **vertex,double **normals,char* retstr,int Avg)
{
    double midpt[3];
    double edge1[3],edge2[3],edge3[3];
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
    double delq,delr;
    double dA;
    double minr,maxr;
    double shv[3][3],shn[3][3];
    double **nn;
    double **nv;
    double nmidpt[3];
    double apex[3];
    double qval,qv;
    double pa,sfa;
    double dda[3];
    double nt[3];
    double dot12,dot23,dot31;
#ifdef spheretest
    float* transarray;
#endif
    double obfrac,bofrac;

    for(i=0;i<3;i++)
    {
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
        exptype=0;
        if(FindApex(vertex[i],normals[i],apex))
            WriteTriangleDetail(vertex,normals,apex,normals[0],-1.0);
        for(j=0;j<3;j++)
            vertex[i][j]=apex[j];
    }
    for(i=0;i<3;i++)
    {
        edge1[i]=vertex[1][i]-vertex[0][i];
        edge2[i]=vertex[2][i]-vertex[1][i];
        edge3[i]=vertex[0][i]-vertex[2][i];
    }
#ifdef WriteTriangleVertices
    if(!nestlevel)
    {
        for(i=0;i<3;i++)
            WriteTriangleVertex(vertex[i],normals[i]);
    }
#endif
    ++ntriangles;
    ++nestlevel;
    rmnorm[0]=ww[0]=edge1[1]*edge2[2]-edge1[2]*edge2[1];  //find the normal to the triangle
    rmnorm[1]=ww[1]=edge1[2]*edge2[0]-edge1[0]*edge2[2];
    rmnorm[2]=ww[2]=edge1[0]*edge2[1]-edge1[1]*edge2[0];
    parea=0.5*sqrt(ww[0]*ww[0]+ww[1]*ww[1]+ww[2]*ww[2]);
    if(nestlevel==1)
        P_AREA_CUTOFF=parea*P_CUT_SCALE;
    for(i=0;i<3;i++)
        rmnorm[i]*=0.5/parea;  //normalize the normal
/*  printf("%g %g %g\n%g %g %g\n",rmnorm[0],rmnorm[1],rmnorm[2],ww[0],ww[1],ww[2]);
    printf("%g\n",parea);
    gh=getch();*/
#ifdef spheretest
    transarray=gom_GetTranslateArray();
    printf("%g %g %g\n",transarray[0],transarray[1],transarray[2]);
    for(i=0;i<3;i++)
    {
        r[i]=0.0;
        dot12=dot23=0.0;
        for(j=0;j<3;j++)
        {
            r[i]+=vertex[i][j]*normals[i][j];
            dot12+=vertex[i][j]*vertex[i][j];
            dot23+=normals[i][j]*normals[i][j];
        }
        r[i]/=sqrt(dot12*dot23);
        printf("%g %g %g\n",r[i],dot12,dot23);
    }
    dot12=dot23=0.0;
    for(i=0;i<3;i++)
    {
        dot12+=rmnorm[i]*midpt[i];
        dot23+=midpt[i]*midpt[i];
    }
    dot12/=sqrt(dot23);
    printf("%g\n",dot12);
    if(fabs(dot12)<0.75)
    {
        dot12=0.0;
        for(i=0;i<3;i++)
        {
            rmnorm[i]=0.0;
            for(j=0;j<3;j++)
            {
                rmnorm[i]+=normals[j][i];
            }
            dot12+=rmnorm[i]*rmnorm[i];
        }
        dot12=1.0/sqrt(dot12);
        for(i=0;i<3;i++)
            rmnorm[i]*=dot12;
    }
    if((fabs(r[0])<0.95)||(fabs(r[1])<0.95)||(fabs(r[2])<0.95))
    {
        printf("%i\n",nestlevel);
        printf("%g %g %g\n",r[0],r[1],r[2]);
        printf("%g %g %g\n",vertex[0][0],vertex[0][1],vertex[0][2]);
        printf("%g %g %g\n",normals[0][0],normals[0][1],normals[0][2]);
        printf("%g %g %g\n",vertex[1][0],vertex[1][1],vertex[1][2]);
        printf("%g %g %g\n",normals[1][0],normals[1][1],normals[1][2]);
        printf("%g %g %g\n",vertex[2][0],vertex[2][1],vertex[2][2]);
        printf("%g %g %g\n",normals[2][0],normals[2][1],normals[2][2]);
//      gh=getch();
    }
#endif
    k=0;
    for(i=0;i<3;i++)
    {
        dotp=0.0;
        for(j=0;j<3;j++)
            dotp+=rmnorm[j]*normals[i][j];
        r[i]=dotp;
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
    dot12=dot23=dot31=0.0;
    for(i=0;i<3;i++)
    {
        dot12+=edge1[i]*edge2[i];
        dot23+=edge2[i]*edge3[i];
        dot31+=edge3[i]*edge1[i];
    }
    if(dot12>0.0)
    {
        if(parea<P_AREA_CUTOFF)
        {
            sarea=parea;
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,midpt,rmnorm,-1.0);
#endif
            if(Avg)
            {
                qval=sarea*GetAverage(Avg,midpt,vertex);
            }
            dA=0.0;
        }
        else if(nestlevel>=MAX_NEST)
        {
            sarea=parea;
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,midpt,rmnorm,-2.0);
#endif
            if(Avg)
            {
                qval=sarea*GetAverage(Avg,midpt,vertex);
            }
            dA=0.0;
        }
        else
        {
            obfrac=SplitObtuse(vertex[2],vertex[0],vertex[1]);
            bofrac=1.0-obfrac;
            qval=0.0;
            for(i=0;i<3;i++)
            {
                midpt[i]=obfrac*vertex[0][i]+bofrac*vertex[2][i];
                midptnorm[i]=obfrac*normals[0][i]+bofrac*normals[2][i];
            }
            exptype=11;
            FindApex(midpt,midptnorm,apex);
            nv=(double**)malloc(3*sizeof(double*));
            nn=(double**)malloc(3*sizeof(double*));
            for(i=0;i<3;i++)
            {
                nv[i]=(double*)malloc(3*sizeof(double));
                nn[i]=(double*)malloc(3*sizeof(double));
            }
            for(i=0;i<3;i++)
            {
                nv[0][i]=vertex[0][i];
                nn[0][i]=normals[0][i];
                nv[1][i]=vertex[1][i];
                nn[1][i]=normals[1][i];
                nv[2][i]=apex[i];
                nn[2][i]=midptnorm[i];
            }
#ifdef WriteTriangleVertices
            WriteTriangleVertex(nv[2],nn[2]);
#endif
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf %lf",&parea,&sarea,dda,&qv);
            qval=qv*sarea;
            for(i=0;i<3;i++)
            {
                nv[0][i]=vertex[1][i];
                nn[0][i]=normals[1][i];
                nv[1][i]=vertex[2][i];
                nn[1][i]=normals[2][i];
                nv[2][i]=apex[i];
                nn[2][i]=midptnorm[i];
            }
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf %lf",&pa,&sfa,dda+1,&qv);
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
    else if(dot23>0.0)
    {
        if(parea<P_AREA_CUTOFF)
        {
            sarea=parea;
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,midpt,rmnorm,-3.0);
#endif
            if(Avg)
            {
                qval=sarea*GetAverage(Avg,midpt,vertex);
            }
            dA=0.0;
        }
        else if(nestlevel>=MAX_NEST)
        {
            sarea=parea;
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,midpt,rmnorm,-4.0);
#endif
            if(Avg)
            {
                qval=sarea*GetAverage(Avg,midpt,vertex);
            }
            dA=0.0;
        }
        else
        {
            obfrac=SplitObtuse(vertex[0],vertex[1],vertex[2]);
            bofrac=1.0-obfrac;
            qval=0.0;
            for(i=0;i<3;i++)
            {
                midpt[i]=obfrac*vertex[1][i]+bofrac*vertex[0][i];
                midptnorm[i]=obfrac*normals[1][i]+bofrac*normals[0][i];
            }
            exptype=12;
            FindApex(midpt,midptnorm,apex);
            nv=(double**)malloc(3*sizeof(double*));
            nn=(double**)malloc(3*sizeof(double*));
            for(i=0;i<3;i++)
            {
                nv[i]=(double*)malloc(3*sizeof(double));
                nn[i]=(double*)malloc(3*sizeof(double));
            }
            for(i=0;i<3;i++)
            {
                nv[0][i]=vertex[1][i];
                nn[0][i]=normals[1][i];
                nv[1][i]=vertex[2][i];
                nn[1][i]=normals[2][i];
                nv[2][i]=apex[i];
                nn[2][i]=midptnorm[i];
            }
#ifdef WriteTriangleVertices
            WriteTriangleVertex(nv[2],nn[2]);
#endif
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf %lf",&parea,&sarea,dda,&qv);
            qval=qv*sarea;
            for(i=0;i<3;i++)
            {
                nv[0][i]=vertex[2][i];
                nn[0][i]=normals[2][i];
                nv[1][i]=vertex[0][i];
                nn[1][i]=normals[0][i];
                nv[2][i]=apex[i];
                nn[2][i]=midptnorm[i];
            }
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf %lf",&pa,&sfa,dda+1,&qv);
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
    else if(dot31>0.0)
    {
        if(parea<P_AREA_CUTOFF)
        {
            sarea=parea;
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,midpt,rmnorm,-5.0);
#endif
            if(Avg)
            {
                qval=sarea*GetAverage(Avg,midpt,vertex);
            }
            dA=0.0;
        }
        else if(nestlevel>=MAX_NEST)
        {
            sarea=parea;
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,midpt,rmnorm,-6.0);
#endif
            if(Avg)
            {
                qval=sarea*GetAverage(Avg,midpt,vertex);
            }
            dA=0.0;
        }
        else
        {
            obfrac=SplitObtuse(vertex[1],vertex[2],vertex[0]);
            bofrac=1.0-obfrac;
            qval=0.0;
            for(i=0;i<3;i++)
            {
                midpt[i]=obfrac*vertex[2][i]+bofrac*vertex[1][i];
                midptnorm[i]=obfrac*normals[2][i]+bofrac*normals[1][i];
            }
#ifdef WriteTriangleVertices
            WriteTriangleVertex(midpt,midptnorm);
#endif
            exptype=13;
            FindApex(midpt,midptnorm,apex);
            nv=(double**)malloc(3*sizeof(double*));
            nn=(double**)malloc(3*sizeof(double*));
            for(i=0;i<3;i++)
            {
                nv[i]=(double*)malloc(3*sizeof(double));
                nn[i]=(double*)malloc(3*sizeof(double));
            }
            for(i=0;i<3;i++)
            {
                nv[0][i]=vertex[0][i];
                nn[0][i]=normals[0][i];
                nv[1][i]=vertex[1][i];
                nn[1][i]=normals[1][i];
                nv[2][i]=apex[i];
                nn[2][i]=midptnorm[i];
            }
#ifdef WriteTriangleVertices
            WriteTriangleVertex(nv[2],nn[2]);
#endif
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf %lf",&parea,&sarea,dda,&qv);
            printf("%s\n",retstr);
            qval=qv*sarea;
            for(i=0;i<3;i++)
            {
                nv[0][i]=vertex[2][i];
                nn[0][i]=normals[2][i];
                nv[1][i]=vertex[0][i];
                nn[1][i]=normals[0][i];
                nv[2][i]=apex[i];
                nn[2][i]=midptnorm[i];
            }
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf %lf",&pa,&sfa,dda+1,&qv);
            printf("%s\n",retstr);
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
    else
    {
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
            nmidpt[i]=midpt[i]+(ravg-havg)*midptnorm[i];   //translate to put the center of the sphere at the origin
    }
    else
    {
        exptype=2;
        FindApex(midpt,midptnorm,h);
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
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,nmidpt,midptnorm,ravg);
#endif
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
        if(Avg)
        {
            qval=sarea*GetAverage(Avg,h,vertex);
        }
/*      printf("%g %g %g\n",vertex[0][0],vertex[0][1],vertex[0][2]);
        printf("%g %g %g\n",vertex[1][0],vertex[1][1],vertex[1][2]);
        printf("%g %g %g\n",vertex[2][0],vertex[2][1],vertex[2][2]);
        printf("%g %g %g\n",midpt[0],midpt[1],midpt[2]);
        printf("%g %g %g\n",normals[0][0],normals[0][1],normals[0][2]);
        printf("%g %g %g\n",normals[1][0],normals[1][1],normals[1][2]);
        printf("%g %g %g\n",normals[2][0],normals[2][1],normals[2][2]);
        printf("%g %g %g\n",midptnorm[0],midptnorm[1],midptnorm[2]);
        printf("%g %g %g %g\n",r[0],r[1],r[2],ravg);
        printf("%g %g %g %g\n",A,B,C,(A+B+C-pi));
        printf("%i Areas: %g %g %g %g\n",nestlevel,parea,sarea,dA,(sarea-parea)/sarea*100);*/
//      gh=getch();
    }
    else
    {
//if one of the angles is 0, then the "sphere" is really a plane.
        if(parea<P_AREA_CUTOFF)
        {
            sarea=parea;
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,midpt,rmnorm,-7.0);
#endif
            if(Avg)
            {
                qval=sarea*GetAverage(Avg,midpt,vertex);
            }
            dA=0.0;
//          printf("Small triangle threshold reached.  Area=%g\n",parea);
        }
        else if(nestlevel>=MAX_NEST)
        {
            sarea=parea;
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,midpt,rmnorm,-8.0);
#endif
            if(Avg)
            {
                qval=sarea*GetAverage(Avg,midpt,vertex);
            }
            dA=0.0;
//          printf("Maximum depth of mesh refinement reached.  Area=%g\n",parea);
        }
        else if((maxr-minr)/minr<=R_CUT_RATIO)
        {
/*          printf("Degenerate spherical triangle! ");
            if((1.0-ca*ca)>0.0)
                printf("Angle a=0\n");
            else if ((1.0-ca*ca)>0.0)
                printf("Angle b=0\n");
            else
                printf("Angle c=0\n");*/
            sarea=parea;
#ifdef WriteTriangleDetails
            WriteTriangleDetail(vertex,normals,midpt,rmnorm,-9.0);
#endif
            if(Avg)
            {
                qval=sarea*GetAverage(Avg,midpt,vertex);
            }
            dA=0.0;
        }
        else
        {
/*          printf("Increasing mesh Area=%g\n",parea);
            printf("Radii inconsistent by %g (%g %g)\n",(maxr-minr)/minr,maxr,minr);*/
            nv=(double**)malloc(3*sizeof(double*));
            nn=(double**)malloc(3*sizeof(double*));
            for(i=0;i<3;i++)
            {
                nv[i]=(double*)malloc(3*sizeof(double));
                nn[i]=(double*)malloc(3*sizeof(double));
            }
            exptype=3;
            if(!Avg)
                FindApex(midpt,midptnorm,apex);
            else
                FindApex(h,midptnorm,apex);
#ifdef WriteTriangleVertices
            WriteTriangleVertex(apex,midptnorm);
#endif
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
            qval=qv*sarea;
            for(i=0;i<3;i++)
            {
                nv[1][i]=vertex[0][i];
                nv[2][i]=vertex[2][i];
                nn[1][i]=normals[0][i];
                nn[2][i]=normals[2][i];
            }
            TriangleStuff(nv,nn,retstr,Avg);
            sscanf(retstr,"%lf %lf %*f %lf %lf",&pa,&sfa,dda+1,&qv);
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
            sscanf(retstr,"%lf %lf %*f %lf %lf",&pa,&sfa,dda+1,&qv);
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
    }
    --nestlevel;
//  printf("%i %i\n",nestlevel,ntriangles);
#ifdef WriteTriangleVertices
    if(!nestlevel)
    {
        WriteTriangleVertex(NULL,NULL);
    }
#endif
    sprintf(retstr,"%g %g %g %g %g",parea,sarea,(sarea-parea)/sarea*100,dA,qval);//return our result
#ifdef ReportIntermediates
    printf("%i %s\n",nestlevel,retstr);
#endif
}

int CalcTriangleAreas(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    char tstr[80];
    double **vertices;
    double **normals;
    int i;
    int Avg;
    static int ntlast=0;


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
    printf("%i %i \n",ntriangles,ntriangles-ntlast);
    ntlast=ntriangles;
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

int MakeEllipsoidalMesh(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    int nx,ny,nz;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    double xi,eta,zeta;
    double lambda;
    double omega;
    double xr,yr,zr;
    double dx,dy,dz;
    int i,j,k;
    float *data;

    nx=atoi(argv[1]);
    ny=atoi(argv[2]);
    nz=atoi(argv[3]);
    xmin=atof(argv[4]);
    xmax=atof(argv[5]);
    ymin=atof(argv[6]);
    ymax=atof(argv[7]);
    zmin=atof(argv[8]);
    zmax=atof(argv[9]);
    xi=atof(argv[10]);
    eta=atof(argv[11]);
    zeta=atof(argv[12]);
    lambda=atof(argv[13]);
    data=(float*)malloc(nx*ny*nz*sizeof(float));
    dx=(xmax-xmin)/(double)(nx-1);
    dy=(ymax-ymin)/(double)(ny-1);
    dz=(zmax-zmin)/(double)(nz-1);
    for(k=0;k<nz;k++)
    {
        zr=dz*(double)k+zmin;
        for(j=0;j<ny;j++)
        {
            yr=dy*(double)j+ymin;
            for(i=0;i<nx;i++)
            {
                xr=dx*(double)i+xmin;
                omega=sqrt((xr/xi)*(xr/xi)+(yr/eta)*(yr/eta)+(zr/zeta)*(zr/zeta));
                data[k*ny*nx+j*nx+i]=(float)exp(-lambda*omega);
            }
        }
    }
    gom_FillContourStructure("testa","",nx,ny,nz,(float)xmin,(float)xmax,(float)ymin,(float)ymax,(float)zmin,(float)zmax,data);
    return TCL_OK;
}

int SetTriangleParms(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    P_CUT_SCALE=atof(argv[1]);
    R_CUT_RATIO=atof(argv[2]);
    MAX_NEST=atoi(argv[3]);
    return TCL_OK;
}

int ResetTriangleCount(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    ntriangles=0;
    return TCL_OK;
}

int InterpolateMesh(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    return TCL_OK;
}

int SumMesh(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    return TCL_OK;
}

int DiffMesh(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    return TCL_OK;
}

int ProdMesh(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    return TCL_OK;
}

int QuotMesh(ClientData cd,Tcl_Interp *ti,int nv,const char **argv)
{
    return TCL_OK;
}

} // namespace Surfarea
} // namespace Plugin
