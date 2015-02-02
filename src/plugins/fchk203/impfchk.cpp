#include "impfchk.h"

namespace Plugin {
namespace Fchk {

double** FullCoeffs;
GaussianOrbital* Basis;
GaussianOrbital* SBasisS;
GaussianOrbital* SBasisM;
GaussAtom* AtomList=NULL;
int NumGBasis;
int fun1,fun2;

int* ShellStart;

int* OrbitalToAtom;
int* ShellToAtom;

int NumContShells;
int NumUncontShells;
int NumIndepFun;
int NumAtoms;

double* AlphaEnergies;
double* BetaEnergies;
double* AlphaCoeffs;
double* BetaCoeffs;
ContractedShell* ContShells;
char** AlfOrbID;
char** BetaOrbID;
double DipoleMoment[3]={9.9e99,9.9e99,9.9e99};

int norb;

int NumAlphaElectrons;
int NumBetaElectrons;

double* Mulliken;

int ArraysClear=1;

int IntFlag(char* d)
{
    if((d[0]=='I')&&(d[1]=='\000'))
        return 1;
    else
        return 0;
}

int FloatFlag(char* d)
{
    if((d[0]=='R')&&(d[1]=='\000'))
        return 1;
    else
        return 0;
}

int CountFlag(char* d)
{
    if((d[0]=='N')&&(d[1]=='=')&&(d[2]=='\000'))
        return 1;
    else
        return 0;
}

int SaysNumber(char* d)
{
    return !cmp(d,"Number");
}

int SaysShells(char* d)
{
    return !cmp(d,"shells");
}

int SaysContracted(char* d)
{
    return !cmp(d,"contracted");
}

int SaysUncontracted(char* d)
{
    return !cmp(d,"primitive");
}

int SaysShell(char* d)
{
    return !cmp(d,"Shell");
}

int SaysTypes(char* d)
{
    return !cmp(d,"types");
}

int ReadFCHK(int argc,const char* argv[])
{
    char linein[132];
    char* temp[10];
    FILE* infile;
    int k=0;
    int s=0;
    int sold=0;
    int NumStr;
    int i;
    int j;
    int NumVal;
    int IntVal;
    int ng=0;
    double FVal;
    int betafound=0;
    int nshell;
    int n1=0;
    double tval[100];
    int no;

    if(!ArraysClear)
    {
        printf("Error in ImpFCHK.  Arrays not yet clear.");
    }
    for(i=0;i<10;i++)
        temp[i]=(char*)malloc(30*sizeof(char));
    infile=fopen(argv[1],"r");
    if(infile==NULL)
    {
        printf("Error in ImpFChk.  Unable to open file %s.\n",argv[1]);
        return -1;
    }
    while(k>=0)
    {
//      printf("%i\n",k=ReadToEOL(infile,linein));
//      printf("%s\n",linein);*/
        k=ReadToEOL(infile,linein);
        s=sold=0;
        NumStr=0;
        while((s=SplitString(linein,temp[NumStr],s))!=abs(k))
        {
            if(s>(sold+1))
                NumStr+=1;
            sold=s;
        }
/*      for(i=0;i<NumStr;i++)
        {
            if(IntFlag(temp[i]))
            {
                for(j=0;j<i;j++)
                {
                    printf("%s ",temp[j]);
                }
                if(CountFlag(temp[++i]))
                {
                    NumVal=atoi(temp[++i]);
                    printf("Has %i integer values\n",NumVal);
                }
                else
                {
                    NumVal=atoi(temp[i]);
                    printf("Has value %i\n",NumVal);
                }

            }
            if(FloatFlag(temp[i]))
            {
                for(j=0;j<i;j++)
                {
                    printf("%s ",temp[j]);
                }
                if(CountFlag(temp[++i]))
                {
                    NumVal=atoi(temp[++i]);
                    printf("Has %i real values",NumVal);
                }
                else
                {
                    FVal=atof(temp[i]);
                    printf("Has value %g\n",FVal);
                }
            }
        }*/
        if(SaysNumber(temp[0]))
        {
            if(SaysShells(temp[3]))
            {
                if(SaysContracted(temp[2]))
                {
                    NumContShells=atoi(temp[5]);
                    ContShells=(ContractedShell*)malloc(NumContShells*sizeof(ContractedShell));
                }
                if(SaysUncontracted(temp[2]))
                {
                    NumUncontShells=atoi(temp[5]);
                }
            }
            if(!cmp(temp[3],"electrons"))
            {
                if(!cmp(temp[2],"alpha"))
                    NumAlphaElectrons=atoi(temp[5]);
                else
                    NumBetaElectrons=atoi(temp[5]);
            }
            if(!cmp(temp[2],"independant"))
            {
                NumIndepFun=atoi(temp[5]);
                printf("%i independent functions\n",NumIndepFun);
                AlphaEnergies=(double*)malloc(NumIndepFun*sizeof(double));
                BetaEnergies=(double*)malloc(NumIndepFun*sizeof(double));
                AlphaCoeffs=(double*)malloc(NumIndepFun*NumIndepFun*sizeof(double));
                BetaCoeffs=(double*)malloc(NumIndepFun*NumIndepFun*sizeof(double));
                AlfOrbID=(char**)malloc(NumIndepFun*sizeof(char*));
                BetaOrbID=(char**)malloc(NumIndepFun*sizeof(char*));
                for(i=0;i<NumIndepFun;i++)
                {
                    AlfOrbID[i]=(char*)malloc(10*sizeof(char));
                    BetaOrbID[i]=(char*)malloc(10*sizeof(char));
                }
            }
        }
        if(SaysShell(temp[0]))
        {
            if(SaysTypes(temp[1]))
            {
                NumVal=atoi(temp[4]);
                for(i=0;i<NumVal;i++)
                {
                    fscanf(infile,"%i",&IntVal);
                    ContShells[i].type=IntVal;
//                  printf("Shell %i is of type %i\n",i,ContShells[i].type);
                }
            }
            if(!cmp(temp[2],"atom"))
            {
                nshell=atoi(temp[6]);
                ShellToAtom=(int*)malloc(nshell*sizeof(int));
                for(i=0;i<nshell;i++)
                {
                    fscanf(infile,"%i",&IntVal);
                    ShellToAtom[i]=IntVal;
                }
            }
        }
        if(SaysNumber(temp[0]))
        {
            if(!cmp(temp[2],"primitives"))
            {
                if(!cmp(temp[4],"shell"))
                {
                    NumVal=atoi(temp[7]);
                    for(i=0;i<NumVal;i++)
                    {
                        fscanf(infile,"%i",&IntVal);
                        ContShells[i].NumberOfPrimitives=IntVal;
                        printf("Shell %i has %i primitives.\n",i,ContShells[i].NumberOfPrimitives);
                        ContShells[i].PrimitiveExponents=(double*)malloc(ContShells[i].NumberOfPrimitives*sizeof(double));
                        if(ContShells[i].type!=-1)
                        {
                            ContShells[i].ContractionCoeffs=(double*)malloc(ContShells[i].NumberOfPrimitives*sizeof(double));
                        }
                        else
                        {
                            ContShells[i].ContractionCoeffs=(double*)malloc(2*ContShells[i].NumberOfPrimitives*sizeof(double));
                        }
                    }
                }
            }
            if(!cmp(temp[2],"atoms"))
            {
                NumVal=atoi(temp[4]);
                NumAtoms=NumVal;
                AtomList=(GaussAtom*)malloc(NumAtoms*sizeof(GaussAtom));
            }
        }
        if(!cmp(temp[0],"Atomic"))
        {
            if(!cmp(temp[1],"numbers"))
            {
                for(i=0;i<NumAtoms;i++)
                {
                    fscanf(infile,"%i",&IntVal);
                    AtomList[i].AtNum=IntVal;
                }
            }
        }
        if(!cmp(temp[0],"Current"))
        {
            if(!cmp(temp[1],"cartesian"))
            {
                if(!cmp(temp[2],"coordinates"))
                {
                    for(i=0;i<NumAtoms;i++)
                    {
                        for(j=0;j<3;j++)
                        {
                            fscanf(infile,"%lf",&FVal);
                            AtomList[i].Position[j]=FVal;
                        }
                    }
                }
            }
        }
        if(!cmp(temp[0],"Primitive"))
        {
            if(!cmp(temp[1],"exponents"))
            {
                for(i=0;i<NumContShells;i++)
                {
                    for(j=0;j<ContShells[i].NumberOfPrimitives;j++)
                    {
                        fscanf(infile,"%lf",&FVal);
                        ContShells[i].PrimitiveExponents[j]=FVal;
//                      printf("Shell %i primitive %i exponent is %g\n",i,j,ContShells[i].PrimitiveExponents[j]);
                    }
                }
            }
        }
        if(!cmp(temp[0],"Contraction"))
        {
            if(!cmp(temp[1],"coefficients"))
            {
                for(i=0;i<NumContShells;i++)
                {
                    for(j=0;j<ContShells[i].NumberOfPrimitives;j++)
                    {
                        fscanf(infile,"%lf",&FVal);
                        ContShells[i].ContractionCoeffs[j]=FVal;
//                      printf("Shell %i primitive %i coefficient is %g\n",i,j,ContShells[i].ContractionCoeffs[j]);
                    }
                }
            }
        }
        if(!cmp(temp[0],"P"))
        {
            if(!cmp(temp[2],"Contraction"))
            {
                if(!cmp(temp[3],"coefficients"))
                {
                    for(i=0;i<NumContShells;i++)
                    {
                        for(j=0;j<ContShells[i].NumberOfPrimitives;j++)
                        {
                            fscanf(infile,"%lf",&FVal);
                            if(ContShells[i].type==-1)
                            {
                                ContShells[i].ContractionCoeffs[j+ContShells[i].NumberOfPrimitives]=FVal;
//                              printf("Shell %i p-primitive %i coefficient is %g\n",i,j,ContShells[i].ContractionCoeffs[j+ContShells[i].NumberOfPrimitives]);
                            }
                        }
                    }
                }
            }
        }
        if(!cmp(temp[0],"Coordinates"))
        {
            if(!cmp(temp[2],"each"))
            {
                if(!cmp(temp[3],"shell"))
                {
                    for(i=0;i<NumContShells;i++)
                    {
                        for(j=0;j<3;j++)
                        {
                            fscanf(infile,"%lf",&FVal);
                            ContShells[i].Center[j]=FVal;
                        }   
                    }
                }
            }
        }
        if(!cmp(temp[0],"Alpha"))
        {
            if(!cmp(temp[2],"Energies"))
            {
//              energyfile=fopen("MOEnergies.txt","w");
                NumVal=atoi(temp[5]);
                for(i=0;i<NumVal;i++)
                {
                    fscanf(infile,"%lf",&FVal);
                    AlphaEnergies[i]=FVal;
                    printf("Alpha energy %i=%g\n",i,AlphaEnergies[i]);
//                  fprintf(energyfile,"%g\n",AlphaEnergies[i]);
                }
//              fprintf(energyfile,"\n");
//              fclose(energyfile);
//              PauseForAction();
            }
            if(!cmp(temp[2],"coefficients"))
            {
                NumVal=atoi(temp[5]);
                for(i=0;i<NumVal;i++)
                {
                    fscanf(infile,"%lf",&FVal);
                    AlphaCoeffs[i]=FVal;
                    printf("Alpha coefficient %i=%g\n",i,AlphaCoeffs[i]);
                }
//              PauseForAction();
            }
        }
        if(!cmp(temp[0],"Beta"))
        {
            if(!cmp(temp[2],"Energies"))
            {
                betafound=1;
//              energyfile=fopen("MOEnergies.txt","a");
                NumVal=atoi(temp[5]);
                for(i=0;i<NumVal;i++)
                {
                    fscanf(infile,"%lf",&FVal);
                    BetaEnergies[i]=FVal;
                    printf("Beta energy %i=%g\n",i,BetaEnergies[i]);
//                  fprintf(energyfile,"%g\n",BetaEnergies[i]);
                }
//              fclose(energyfile);
//              PauseForAction();
            }
            if(!cmp(temp[2],"coefficients"))
            {
                NumVal=atoi(temp[5]);
                for(i=0;i<NumVal;i++)
                {
                    fscanf(infile,"%lf",&FVal);
                    BetaCoeffs[i]=FVal;
                    printf("Beta coefficient %i=%g\n",i,BetaCoeffs[i]);
                }
//              PauseForAction();
            }
        }
        if(!cmp(temp[0],"Mulliken"))
        {
        }
        if(!cmp(temp[0],"Dipole"))
        {
            for(i=0;i<3;i++)
                fscanf(infile,"%lf",DipoleMoment+i);
        }
//      PauseForAction();
    }
    if(!betafound)
    {
        for(i=0;i<NumIndepFun;i++)
        {
            BetaEnergies[i]=AlphaEnergies[i];
            for(j=0;j<NumIndepFun;j++)
                BetaCoeffs[i*NumIndepFun+j]=AlphaCoeffs[i*NumIndepFun+j];
        }
    }
//All the relevant information is now in memory.  Report atom list
    for(i=0;i<NumAtoms;i++)
    {
        AtomList[i].AssignType();
        printf("Atom %i is of type %s\n Location ",i,AtomList[i].Element);
        for(j=0;j<3;j++)
            printf("%g ",AtomList[i].Position[j]);
        printf("\n");
    }
//  PauseForAction();
//Build Gaussian orbitals.
    NumUncontShells=0;
    for(i=0;i<NumContShells;i++)
    {
        NumUncontShells+=ContShells[i].NumberOfPrimitives;
        switch(ContShells[i].type)
        {
        case -3:    ng+=7*ContShells[i].NumberOfPrimitives;
                    break;
        case -2:    ng+=5*ContShells[i].NumberOfPrimitives;
                    break;
        case -1:    ng+=4*ContShells[i].NumberOfPrimitives;
                    NumUncontShells+=ContShells[i].NumberOfPrimitives;
                    break;
        case 0:     ng+=ContShells[i].NumberOfPrimitives;
                    break;
        case 1:     ng+=3*ContShells[i].NumberOfPrimitives;
                    break;
        case 2:     ng+=6*ContShells[i].NumberOfPrimitives;
                    break;
        case 3:     ng+=10*ContShells[i].NumberOfPrimitives;
                    break;
        default: ng+=0;
        }
    }
    NumGBasis=ng;
    printf("%i Gaussian orbitals needed for calculation.\n",ng);
//  PauseForAction();
    Basis=(GaussianOrbital*)malloc(ng*sizeof(GaussianOrbital));
    FullCoeffs=(double**)malloc(ng*sizeof(double*));
    for(i=0;i<ng;i++)
    {
        FullCoeffs[i]=(double*)malloc(2*NumIndepFun*sizeof(double));
    }
    ShellStart=(int*)malloc((NumUncontShells+1)*sizeof(int));
    ShellStart[NumUncontShells]=ng;
    ng=0;
    norb=0;
    nshell=0;
    for(i=0;i<NumContShells;i++)
    {
        printf("%i %i\n",i,norb);
//      PauseForAction();
        switch(ContShells[i].type)
        {
        case -3:    for(k=0;k<ContShells[i].NumberOfPrimitives;k++)
                    {
                        ShellStart[nshell++]=ng;
                        Basis[ng].AngMom=-3;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-3;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-3;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=-1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-3;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=2;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-3;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=-2;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+4]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+4]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-3;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=3;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+5]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+5]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-3;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=-3;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+6]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+6]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                    }
                    norb+=7;
                    break;
        case -2:    for(k=0;k<ContShells[i].NumberOfPrimitives;k++)
                    {
                        ShellStart[nshell++]=ng;
                        Basis[ng].AngMom=-2;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=1;
                        Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+4]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+4]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-2;
                        Basis[ng].AngMomV[1]=Basis[ng].AngMomV[2]=1;
                        Basis[ng].AngMomV[0]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-2;
                        Basis[ng].AngMomV[2]=Basis[ng].AngMomV[0]=1;
                        Basis[ng].AngMomV[1]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-2;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=-2;
                        Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=-2;
                        Basis[ng].AngMomV[0]=Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=-2;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                    }
                    norb+=5;
                    break;
        case -1:    for(k=0;k<ContShells[i].NumberOfPrimitives;k++)
                    {
                        ShellStart[nshell++]=ng;
                        Basis[ng].AngMom=0;
                        for(j=0;j<3;j++)
                            Basis[ng].AngMomV[j]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        ShellStart[nshell++]=ng;
                        Basis[ng].AngMom=1;
                        for(j=0;j<3;j++)
                            Basis[ng].AngMomV[j]=0;
                        Basis[ng].AngMomV[0]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[ContShells[i].NumberOfPrimitives+k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[ContShells[i].NumberOfPrimitives+k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=1;
                        for(j=0;j<3;j++)
                            Basis[ng].AngMomV[j]=0;
                        Basis[ng].AngMomV[1]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[ContShells[i].NumberOfPrimitives+k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[ContShells[i].NumberOfPrimitives+k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=1;
                        for(j=0;j<3;j++)
                            Basis[ng].AngMomV[j]=0;
                        Basis[ng].AngMomV[2]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[ContShells[i].NumberOfPrimitives+k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[ContShells[i].NumberOfPrimitives+k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                    }                   
                    norb+=4;
                    break;
        case 0:     for(k=0;k<ContShells[i].NumberOfPrimitives;k++)
                    {
                        ShellStart[nshell++]=ng;
                        Basis[ng].AngMom=0;
                        for(j=0;j<3;j++)
                            Basis[ng].AngMomV[j]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                    }
                    norb+=1;
                    break;
        case 1:     for(k=0;k<ContShells[i].NumberOfPrimitives;k++)
                    {
                        ShellStart[nshell++]=ng;
                        Basis[ng].AngMom=1;
                        for(j=0;j<3;j++)
                            Basis[ng].AngMomV[j]=0;
                        Basis[ng].AngMomV[0]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=1;
                        for(j=0;j<3;j++)
                            Basis[ng].AngMomV[j]=0;
                        Basis[ng].AngMomV[1]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=1;
                        for(j=0;j<3;j++)
                            Basis[ng].AngMomV[j]=0;
                        Basis[ng].AngMomV[2]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                    }                   
                    norb+=3;
                    break;
        case 2:     for(k=0;k<ContShells[i].NumberOfPrimitives;k++)
                    {
                        ShellStart[nshell++]=ng;
                        Basis[ng].AngMom=2;
                        Basis[ng].AngMomV[0]=2;Basis[ng].AngMomV[1]=0;Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=2;
                        Basis[ng].AngMomV[0]=1;Basis[ng].AngMomV[1]=1;Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=2;
                        Basis[ng].AngMomV[0]=1;Basis[ng].AngMomV[1]=0;Basis[ng].AngMomV[2]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+4]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+4]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=2;
                        Basis[ng].AngMomV[0]=0;Basis[ng].AngMomV[1]=2;Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=2;
                        Basis[ng].AngMomV[0]=0;Basis[ng].AngMomV[1]=1;Basis[ng].AngMomV[2]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+5]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+5]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=2;
                        Basis[ng].AngMomV[0]=0;Basis[ng].AngMomV[1]=0;Basis[ng].AngMomV[2]=2;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                    }                   
                    norb+=6;
                    break;
        case 3:     for(k=0;k<ContShells[i].NumberOfPrimitives;k++)
                    {
                        ShellStart[nshell++]=ng;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=3;
                        Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=0;
                        Basis[ng].AngMomV[1]=3;
                        Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+1]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=0;
                        Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=3;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+2]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=2;
                        Basis[ng].AngMomV[1]=1;
                        Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+3]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=2;
                        Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+4]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+4]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=1;
                        Basis[ng].AngMomV[1]=2;
                        Basis[ng].AngMomV[2]=0;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+5]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+5]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=0;
                        Basis[ng].AngMomV[1]=2;
                        Basis[ng].AngMomV[2]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+6]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+6]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=1;
                        Basis[ng].AngMomV[1]=0;
                        Basis[ng].AngMomV[2]=2;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+7]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+7]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=0;
                        Basis[ng].AngMomV[1]=1;
                        Basis[ng].AngMomV[2]=2;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+8]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+8]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                        Basis[ng].AngMom=3;
                        Basis[ng].AngMomV[0]=1;
                        Basis[ng].AngMomV[1]=1;
                        Basis[ng].AngMomV[2]=1;
                        Basis[ng].OrbExponent=ContShells[i].PrimitiveExponents[k];
                        for(j=0;j<NumIndepFun;j++)
                        {
                            FullCoeffs[ng][j]=AlphaCoeffs[j*NumIndepFun+norb+9]*ContShells[i].ContractionCoeffs[k];
                            FullCoeffs[ng][NumIndepFun+j]=BetaCoeffs[j*NumIndepFun+norb+9]*ContShells[i].ContractionCoeffs[k];
                        }
                        for(j=0;j<3;j++)
                            Basis[ng].Center[j]=ContShells[i].Center[j];
                        ng+=1;
                    }
                    norb+=10;
                    break;
        default:    break;
        }
    }
    OrbitalToAtom=(int*)malloc(ng*sizeof(int));
    for(i=0;i<NumContShells;i++)
    {
        int nn=0;
        switch(ContShells[i].type)
        {
        case -2: nn=5*ContShells[i].NumberOfPrimitives;
                 for(k=0;k<nn;k++)
                     OrbitalToAtom[n1++]=ShellToAtom[i]-1;
                 break;
        case -1: nn=4*ContShells[i].NumberOfPrimitives;
                 for(k=0;k<nn;k++)
                     OrbitalToAtom[n1++]=ShellToAtom[i]-1;
                 break;
        case 0:  nn=ContShells[i].NumberOfPrimitives;
                 for(k=0;k<nn;k++)
                     OrbitalToAtom[n1++]=ShellToAtom[i]-1;
                 break;
        case 1:  nn=3*ContShells[i].NumberOfPrimitives;
                 for(k=0;k<nn;k++)
                     OrbitalToAtom[n1++]=ShellToAtom[i]-1;
                 break;
        case 2:  nn=6*ContShells[i].NumberOfPrimitives;
                 for(k=0;k<nn;k++)
                     OrbitalToAtom[n1++]=ShellToAtom[i]-1;
                 break;
        }
    }
    newrecurse::InitMatrices();
    for(i=0;i<NumUncontShells;i++)
    {
        newrecurse::OverlapIntegral(Basis+ShellStart[i],Basis+ShellStart[i],tval);
        no=ShellStart[i+1]-ShellStart[i];
        for(k=0;k<no;k++)
        {
            FVal=tval[k*no+k];
            Basis[ShellStart[i]+k].NormCoeff=sqrt(1.0/FVal);
        }
    }
    for(i=0;i<10;i++)
        free((void*)temp[i]);
    fclose(infile);
    ArraysClear=0;
    return 0;
}

int FreeFCHKArrays()
{
    int i;

    if(ArraysClear)
        return 0;
    for(i=0;i<NumContShells;i++)
    {
        free((void*)ContShells[i].ContractionCoeffs);
        free((void*)ContShells[i].PrimitiveExponents);
    }
    free((void*)ShellStart);
    free((void*)ShellToAtom);
    free((void*)OrbitalToAtom);
    free((void*)ContShells);
    free((void*)AtomList);
    free((void*)AlphaEnergies);
    free((void*)BetaEnergies);
    free((void*)AlphaCoeffs);
    free((void*)BetaCoeffs);
    free((void*)Basis);
    free((void*)SBasisS);
    free((void*)SBasisM);
    for(i=0;i<NumGBasis;i++)
        free((void*)FullCoeffs[i]);
    free((void*)FullCoeffs);
    for(i=0;i<NumIndepFun;i++)
    {
        free((void*)AlfOrbID[i]);
        free((void*)BetaOrbID[i]);
    }
    free((void*)AlfOrbID);
    free((void*)BetaOrbID);
    newrecurse::ClearMatrices();
    ArraysClear=1;
    return 0;
}

} // namespace Fchk
} // namespace Plugin
