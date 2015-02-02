/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
/*Calculates wavefunction information under certain conditions */
#include "tclacc.h"

namespace Plugin {
namespace Fchk {
    
#define pi 3.1415926535897932384626433832795
#define sq2 1.4142135623730950488016887242097

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
extern int NumAlphaElectrons;
extern int NumBetaElectrons;

double EvalPsi(int Orbital,int beta,double *R)
{
    int i;
    double PsiRet=0.0;

    if(beta)
        Orbital+=NumIndepFun;
    for(i=0;i<NumGBasis;i++)
        PsiRet+=FullCoeffs[i][Orbital]*Basis[i].Evaluate(R);
    return PsiRet;
}

double EvalGrad(int Orbital,int beta,double *R)
{
    int i;
    double PsiRet=0.0;
    double Grad[3];

    if(beta)
        Orbital+=NumIndepFun;
    for(i=0;i<NumGBasis;i++)
        PsiRet+=FullCoeffs[i][Orbital]*Basis[i].EvaluateGradient(R,Grad);
    return PsiRet;
}


double GetFPlane(int atomi,int atomj,double frac,double* planeq,double* origin)
{
    int i;
    double dist=0.0;

    planeq[3]=0.0;
    for(i=0;i<3;i++)
    {
        planeq[i]=AtomList[atomj].Position[i]-AtomList[atomi].Position[i];
        origin[i]=frac*planeq[i]+AtomList[atomi].Position[i];
        planeq[3]+=planeq[i]*origin[i];
        dist+=(AtomList[atomi].Position[i]-origin[i])*(AtomList[atomi].Position[i]-origin[i]);
    }
    dist=sqrt(dist);
    printf("Plane distance: %g\n",dist);
    return dist;
}

double GetAPlane(int atomi,int atomj,double dist,double* planeq,double* origin)
{
    int i;
    double normfact=0.0;

    planeq[3]=0.0;
    for(i=0;i<3;i++)
    {
        planeq[i]=AtomList[atomi].Position[i]-AtomList[atomj].Position[i];
        normfact+=planeq[i]*planeq[i];
    }
    normfact=dist/sqrt(normfact);
    for(i=0;i<3;i++)
    {
        origin[i]=AtomList[atomi].Position[i]+normfact*planeq[i];
        planeq[3]+=planeq[i]*origin[i];
    }
    printf("Plane distance: %g\n",dist);
    return dist;
}

void GetBasisVecs(double* planeq,double *veca,double* vecb)
{
    double norm[3];
    double val;
    int i;

    val=0.0;
    for(i=0;i<3;i++)
        val+=planeq[i]*planeq[i];
    val=1.0/sqrt(val);
    for(i=0;i<3;i++)
        norm[i]=planeq[i]*val;
    if(fabs(norm[2])<1.e-10)
    {
        veca[0]=veca[1]=0.0;
        veca[2]=1.0;
    }
    else
    {
        veca[0]=norm[0];
        veca[1]=norm[1];
        veca[2]=-(norm[0]*norm[0]+norm[1]*norm[1])/norm[2];
        val=0.0;
        for(i=0;i<3;i++)
            val+=veca[i]*veca[i];
        val=1.0/sqrt(val);
        for(i=0;i<3;i++)
            veca[i]*=val;
    }
    vecb[0]=norm[1]*veca[2]-norm[2]*veca[1];
    vecb[1]=norm[2]*veca[0]-norm[0]*veca[2];
    vecb[2]=norm[0]*veca[1]-norm[1]*veca[0];
}

double CalcPlane(int argc,const char** argv)
{
    int atom1,atom2;
    int orbital,beta;
    int nxstep,nystep;
    double Rr[3];
    double frac;
    double xfrac,yfrac;
    double origin[3];
    double plane[4];
    double xaxis[3],yaxis[3];
    double xmin,xmax,ymin,ymax;
    double *PsiSec;
    int i,j,k,h;
    int error=-1;
    int ok=1;
    char filename[80];
    FILE *of;
    double g;
    int nslA,NumDegLevelsA;
    int nslB,NumDegLevelsB;
    int usefrac;
    double dist;


    if(argc<7)
        return error;
    atom1=atoi(argv[1])-1;
    atom2=atoi(argv[2])-1;
    frac=atof(argv[3]);
    if((argv[4][0]=='f')||(argv[4][0]=='F'))
        usefrac=1;
    else
        usefrac=0;
    sprintf(filename,"%s",argv[5]);
    orbital=atoi(argv[6])-1;
    beta=atoi(argv[7]);
    if(argc<10)
    {
        xmin=-10.0;
        xmax=10.0;
    }
    else
    {
        xmin=atof(argv[8]);
        xmax=atof(argv[9]);
    }
    if(argc<12)
    {
        ymin=-10.0;
        ymax=10.0;
    }
    else
    {
        ymin=atof(argv[10]);
        ymax=atof(argv[11]);
    }
    if(argc<13)
    {
        nxstep=31;
    }
    else
    {
        nxstep=atoi(argv[12]);
    }
    if(argc<14)
    {
        nystep=31;
    }
    else
    {
        nystep=atoi(argv[13]);
    }
    if(usefrac)
        dist=GetFPlane(atom1,atom2,frac,plane,origin);
    else
        dist=GetAPlane(atom1,atom2,frac,plane,origin);
    GetBasisVecs(plane,xaxis,yaxis);
    PsiSec=(double*)malloc(nxstep*nystep*sizeof(double));
    for(nslA=-3;nslA<=3;nslA++)
        printf("%i %g\n",NumAlphaElectrons+nslA,AlphaEnergies[NumAlphaElectrons+nslA]);
//  PauseForAction();
    NumDegLevelsA=1;
    while(fabs(AlphaEnergies[NumAlphaElectrons+NumDegLevelsA]-AlphaEnergies[NumAlphaElectrons])<1.e-4)
        NumDegLevelsA++;
    nslA=1;
    while(fabs(AlphaEnergies[NumAlphaElectrons-nslA]-AlphaEnergies[NumAlphaElectrons])<1.e-4)
        nslA++;
    NumDegLevelsA+=nslA-1;
    printf("%i %i\n",NumDegLevelsA,nslA);
    for(nslB=-3;nslB<=3;nslB++)
        printf("%i %g\n",NumBetaElectrons+nslB,BetaEnergies[NumBetaElectrons+nslB]);
//  PauseForAction();
    NumDegLevelsB=1;
    while(fabs(BetaEnergies[NumBetaElectrons+NumDegLevelsB]-BetaEnergies[NumBetaElectrons])<1.e-4)
        NumDegLevelsB++;
    nslB=1;
    while(fabs(BetaEnergies[NumBetaElectrons-nslB]-BetaEnergies[NumBetaElectrons])<1.e-4)
        nslB++;
    NumDegLevelsB+=nslB-1;
    printf("%i %i\n",NumDegLevelsB,nslB);
    for(i=0;i<nxstep;i++)
    {
        xfrac=xmin+(double)i*(xmax-xmin)/(double)(nxstep-1);
        for(j=0;j<nystep;j++)
        {
            yfrac=ymin+(double)j*(ymax-ymin)/(double)(nystep-1);
            for(k=0;k<3;k++)
                Rr[k]=xfrac*xaxis[k]+yfrac*yaxis[k]+origin[k];
            switch(beta)
            {
            case 0:
            case 1: g=EvalPsi(orbital,beta,Rr);
                    PsiSec[i*nystep+j]=g;
                    break;
            case 2: PsiSec[i*nystep+j]=1.0;
                    for(h=0;h<orbital;h++)
                    {
                        g=EvalPsi(h,0,Rr);
                        PsiSec[i*nystep+j]*=g;
                        g=EvalPsi(h,1,Rr);
                        PsiSec[i*nystep+j]*=g;
                    }
                    break;
            case 3: PsiSec[i*nystep+j]=0.0;
                    for(h=0;h<NumAlphaElectrons-nslA;h++)
                    {
                        g=EvalPsi(h,0,Rr);
                        PsiSec[i*nystep+j]+=g*g;
                    }
                    for(h=NumAlphaElectrons-nslA;h<NumAlphaElectrons-nslA+NumDegLevelsA;h++)
                    {
                        g=EvalPsi(h,0,Rr);
                        PsiSec[i*nystep+j]+=g*g*(double)nslA/(double)NumDegLevelsA;
                    }
                    for(h=0;h<NumBetaElectrons-nslB;h++)
                    {
                        g=EvalPsi(h,1,Rr);
                        PsiSec[i*nystep+j]+=g*g;
                    }
                    for(h=NumBetaElectrons-nslB;h<NumBetaElectrons-nslB+NumDegLevelsB;h++)
                    {
                        g=EvalPsi(h,1,Rr);
                        PsiSec[i*nystep+j]+=g*g*(double)nslB/(double)NumDegLevelsB;
                    }
                    break;
            default:g=EvalPsi(orbital,beta,Rr);
                    PsiSec[i*nystep+j]=g;
                    break;
            }
        }
    }
    of=fopen(filename,"w");
    for(i=0;i<nxstep;i++)
    {
        for(j=0;j<nystep;j++)
        {
            fprintf(of,"%g ",PsiSec[i*nystep+j]);
        }
        fprintf(of,"\n");
    }
    fclose(of);
    free((void*)PsiSec);
    return dist;
}

} // namespace Fchk
} // namespace Plugin
