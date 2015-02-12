#include <gopenmolext.h>
#include <tcl.h>
#include "gaussians.h"
#include "impfchk.h"
#include "matrices.h"
#include "parser.h"
#include "tclacc.h"

namespace Plugin {
namespace Fchk {

#if TCL_RELEASE_LEVEL==TCL_ALPHA_RELEASE
#define TCL_CONST 
#else
#define TCL_CONST const
#endif

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

extern int* ShellStart;

extern int* OrbitalToAtom;

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
extern double DipoleMoment[3];

extern int NumAlphaElectrons;
extern int NumBetaElectrons;

extern int nstep;

double Box[4][3];
int nx=31,ny=31,nz=31;

double* Pmn;
double* Pmnd;
double* Smn;
extern double* Mulliken;
double* Lowdin;

static int InitBox(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int SetBox(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int ScaleBox(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int SetGrain(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int GetParms(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int GetMolecule(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int FreeArrays(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int ExtractAtomNumbers(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int NewFillGrid(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int SaveGrid(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int ListEnergies(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int GraphPlane(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int ResetBox(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int ReportBox(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int CompFileNames(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int AdjustFileName(ClientData ,Tcl_Interp *,int ,TCL_CONST char **);
static int CalcCharges(ClientData, Tcl_Interp *,int ,TCL_CONST char **);
static int GetMullikenCharge(ClientData, Tcl_Interp *,int ,TCL_CONST char **);
static int GetLowdinCharge(ClientData, Tcl_Interp *,int ,TCL_CONST char **);
static int ReportDipoleMoment(ClientData, Tcl_Interp *,int ,TCL_CONST char **);
static void SetInitialBox();


DYNEXPORT_C int Fchk_Init(Tcl_Interp *Interp)
//declare exported command
{
    printf("Creating Tcl Extensions.\n");
Tcl_CreateCommand(Interp,"initbox",InitBox,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"setbox",SetBox,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"scalebox",ScaleBox,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"repbox",ReportBox,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"setgrain",SetGrain,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"getparms",GetParms,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"readfchk",GetMolecule,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"clearmem",FreeArrays,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"getatomid",ExtractAtomNumbers,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"newfillgrid",NewFillGrid,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"listfchkenergies",ListEnergies,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"cc2dgrid",GraphPlane,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"resetbox",ResetBox,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"fncomp",CompFileNames,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"adjfn",AdjustFileName,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"calccharges",CalcCharges,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"mulliken",GetMullikenCharge,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"lowdin",GetLowdinCharge,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"savegrid",SaveGrid,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"getdipole",ReportDipoleMoment,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
SetInitialBox();
Tcl_PkgProvide(Interp,"Fchk","1.0");
    printf("Finished.\n");
return(TCL_OK);
}

void SetInitialBox()
{
    int ns,na;
    float xm=1.e6,xx=-1.e6,ym=1.e6,yx=-1.e6,zm=1.e6,zx=-1.e6;
    int i,j;
    const float *xp,*yp,*zp;


    ns=gom_GetNumMolecStructs();
    for(i=0;i<ns;i++)
    {
        na=gom_GetNumAtomsInMolecStruct(i);
        xp=gom_GetAtomXCoordPointer(i);
        yp=gom_GetAtomYCoordPointer(i);
        zp=gom_GetAtomZCoordPointer(i);
        for(j=0;j<na;j++)
        {
            if(xp[j]>xx)
                xx=xp[j];
            if(xp[j]<xm)
                xm=xp[j];
            if(yp[j]>yx)
                yx=yp[j];
            if(yp[j]<ym)
                ym=yp[j];
            if(zp[j]>zx)
                zx=zp[j];
            if(zp[j]<zm)
                zm=zp[j];
        }
    }
    Box[0][0]=2.0*xm-xx;
    Box[0][1]=2.0*ym-yx;
    Box[0][2]=2.0*zm-zx;
    Box[1][0]=2.0*xx-xm;
    Box[1][1]=2.0*ym-yx;
    Box[1][2]=2.0*zm-zx;
    Box[2][0]=2.0*xm-xx;
    Box[2][1]=2.0*yx-ym;
    Box[2][2]=2.0*zm-zx;
    Box[3][0]=2.0*xm-xx;
    Box[3][1]=2.0*ym-yx;
    Box[3][2]=2.0*zx-zm;
    for(i=0;i<3;i++)
    {
        if(fabs(Box[i+1][i]-Box[0][i])<2.0)
        {
            Box[i+1][i]+=1.0;
            Box[0][i]-=1.0;
            for(j=1;j<4;j++)
            {
                if((i+1)!=j)
                    Box[j][i]-=1.0;
            }
        }
    }
    for(i=0;i<4;i++)
        for(j=0;j<3;j++)
            Box[i][j]/=0.529;
}

int InitBox(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    SetInitialBox();
    return TCL_OK;
}

int SetBox(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    int i,j;

    if(argc<13)
        return TCL_ERROR;
    for(i=0;i<4;i++)
    {
        for(j=0;j<3;j++)
        {
            Box[i][j]=atof(argv[3*i+j+1]);
        }
    }
    return TCL_OK;
}

int ScaleBox(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    double ww[2];
    double scalf;

    if(argc<2)
        scalf=2.0;
    else
    {
        if(!cmp(argv[1],""))
            scalf=2.0;
        else
            scalf=atof(argv[1]);
    }
    ww[0]=0.5*(Box[0][0]+Box[1][0]);
    ww[1]=ww[0]-Box[0][0];
    Box[0][0]=Box[2][0]=Box[3][0]=ww[0]-scalf*ww[1];
    Box[1][0]=ww[0]+scalf*ww[1];
    ww[0]=0.5*(Box[0][1]+Box[2][1]);
    ww[1]=ww[0]-Box[0][1];
    Box[0][1]=Box[1][1]=Box[3][1]=ww[0]-scalf*ww[1];
    Box[2][1]=ww[0]+scalf*ww[1];
    ww[0]=0.5*(Box[0][2]+Box[3][2]);
    ww[1]=ww[0]-Box[0][2];
    Box[0][2]=Box[1][2]=Box[2][2]=ww[0]-scalf*ww[1];
    Box[3][2]=ww[0]+scalf*ww[1];
    return TCL_OK;
}

int SetGrain(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    if(argc<4)
        return TCL_ERROR;
    nx=atoi(argv[1]);
    ny=atoi(argv[2]);
    nz=atoi(argv[3]);
    return TCL_OK;
}

int GetParms(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    int id=atoi(argv[1]);
    char tcl_ret[512];

    switch(id)
    {
    case 0: sprintf(tcl_ret,"%i",nx);
            break;
    case 1: sprintf(tcl_ret,"%i",ny);
            break;
    case 2: sprintf(tcl_ret,"%i",nz);
            break;
    case 3: sprintf(tcl_ret,"%g %g %g",Box[0][0],Box[0][1],Box[0][2]);
            break;
    case 4: sprintf(tcl_ret,"%g %g %g",Box[1][0],Box[1][1],Box[1][2]);
            break;
    case 5: sprintf(tcl_ret,"%g %g %g",Box[2][0],Box[2][1],Box[2][2]);
            break;
    case 6: sprintf(tcl_ret,"%g %g %g",Box[3][0],Box[3][1],Box[3][2]);
            break;
    default: return TCL_ERROR;
    }
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int GetMolecule(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    int nslA,nslB;
    int NumDegLevelsA,NumDegLevelsB;
    int h,m,k;
    char tcl_ret[512];
    int nh,nm,ch,cm;
    int c1,c2;
    double tval[100];

    ReadFCHK(argc,argv);
    Pmn=(double*)malloc(NumGBasis*(NumGBasis)*sizeof(double));
    if(Pmn==NULL)
    {
      printf(Tcl_GetStringResult(ti),"Memory allocation error in density matrix.");
        return TCL_ERROR;
    }
    Smn=(double*)malloc(NumGBasis*(NumGBasis)*sizeof(double));
    if(Smn==NULL)
    {
      printf(Tcl_GetStringResult(ti),"Memory allocation error in overlap matrix.");
        free((void*)Pmn);
        return TCL_ERROR;
    }
    Pmnd=(double*)malloc(NumGBasis*(NumGBasis)*sizeof(double));
    if(Pmnd==NULL)
    {
      printf(Tcl_GetStringResult(ti),"Memory allocation error in density matrix.");
        free((void*)Pmnd);
        free((void*)Smn);
        return TCL_ERROR;
    }
    NumDegLevelsA=1;
    while(fabs(AlphaEnergies[NumAlphaElectrons+NumDegLevelsA]-AlphaEnergies[NumAlphaElectrons])<1.e-4)
        NumDegLevelsA++;
    nslA=1;
    while(fabs(AlphaEnergies[NumAlphaElectrons-nslA]-AlphaEnergies[NumAlphaElectrons])<1.e-4)
        nslA++;
    NumDegLevelsA+=nslA-1;
    printf("%i %i\n",NumDegLevelsA,nslA);
    NumDegLevelsB=1;
    while(fabs(BetaEnergies[NumBetaElectrons+NumDegLevelsB]-BetaEnergies[NumBetaElectrons])<1.e-4)
        NumDegLevelsB++;
    nslB=1;
    while(fabs(BetaEnergies[NumBetaElectrons-nslB]-BetaEnergies[NumBetaElectrons])<1.e-4)
        nslB++;
    NumDegLevelsB+=nslB-1;
    for(h=0;h<NumGBasis;h++)
    {
        for(m=0;m<NumGBasis;m++)
        {
            Pmn[h*NumGBasis+m]=0.0;
            Pmnd[h*NumGBasis+m]=0.0;
            for(k=0;k<NumAlphaElectrons-nslA;k++)
            {
                Pmn[h*NumGBasis+m]+=FullCoeffs[h][k]*FullCoeffs[m][k];
                Pmnd[h*NumGBasis+m]+=FullCoeffs[h][k]*FullCoeffs[m][k];
            }
            for(k=NumAlphaElectrons-nslA;k<NumAlphaElectrons-nslA+NumDegLevelsA;k++)
            {
                Pmn[h*NumGBasis+m]+=FullCoeffs[h][k]*FullCoeffs[m][k]*(double)nslA/(double)NumDegLevelsA;
                Pmnd[h*NumGBasis+m]+=FullCoeffs[h][k]*FullCoeffs[m][k]*(double)nslA/(double)NumDegLevelsA;
            }
            for(k=0;k<NumBetaElectrons-nslB;k++)
            {
                Pmn[h*NumGBasis+m]+=FullCoeffs[h][NumIndepFun+k]*FullCoeffs[m][NumIndepFun+k];
                Pmnd[h*NumGBasis+m]-=FullCoeffs[h][NumIndepFun+k]*FullCoeffs[m][NumIndepFun+k];
            }
            for(k=NumBetaElectrons-nslB;k<NumBetaElectrons-nslB+NumDegLevelsB;k++)
            {
                Pmn[h*NumGBasis+m]+=FullCoeffs[h][NumIndepFun+k]*FullCoeffs[m][NumIndepFun+k]*(double)nslB/(double)NumDegLevelsB;
                Pmnd[h*NumGBasis+m]-=FullCoeffs[h][NumIndepFun+k]*FullCoeffs[m][NumIndepFun+k]*(double)nslB/(double)NumDegLevelsB;
            }
//          Smn[h*NumGBasis+m]=Basis[h].NormCoeff*Basis[m].NormCoeff*OverlapIntegral(Basis[h],Basis[m]);
        }
    }
    for(h=0;h<NumUncontShells;h++)
    {
        nh=ShellStart[h+1]-ShellStart[h];
        for(m=0;m<NumUncontShells;m++)
        {
            nm=ShellStart[m+1]-ShellStart[m];
            newrecurse::OverlapIntegral(Basis+ShellStart[h],Basis+ShellStart[m],tval);
            for(c1=0,ch=ShellStart[h];c1<nh;c1++,ch++)
            {
                for(c2=0,cm=ShellStart[m];c2<nm;c2++,cm++)
                {
                    Smn[ch*NumGBasis+cm]=Basis[ch].NormCoeff*Basis[cm].NormCoeff*tval[c1*nm+cm];
                }
            }
        }
    }
    sprintf(tcl_ret,"%i",NumAtoms);
    SetInitialBox();
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int InList(char s,char* p)
{
    int i;

    i=0;
    while((p[i]!='\000')&&(s!=p[i]))
        i+=1;
    if(p[i]==s)
        return i+1;
    return 0;
}

void ParseAndDiscard(char* stbp,char* pl)
{
    int i,len,indx;

    len=0;
    while(stbp[len]!='\000')
        len+=1;
    indx=0;
    while(!InList(stbp[indx],pl))
    {
        indx+=1;
    }
    indx+=1;
    for(i=indx;i<=len;i++)
        stbp[i-indx]=stbp[i];
}

int ParseAndKeep(char* stbp,char* pl,char* ks)
{
    int i,len,indx,k;

    len=0;
    while(stbp[len]!='\000')
        len+=1;
    indx=0;
    k=InList(stbp[indx],pl);
    while(!k)
    {
        indx+=1;
        k=InList(stbp[indx],pl);
    }
    for(i=0;i<indx;i++)
        ks[i]=stbp[i];
    ks[indx]='\000';
    indx+=1;
    for(i=indx;i<=len;i++)
        stbp[i-indx]=stbp[i];
    return k-1;
}


int ExtractAtomNumbers(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** arggv)
{
    int indx,atom1,atom2;
    char ar1[10];
    char** argv;
    char tcl_ret[512];

    argv=(char**)malloc((argc-1)*sizeof(char*));
    for(indx=0;indx<argc-1;indx++)
    {
        argv[indx]=(char*)malloc(strlen(arggv[indx+1])+1);
        strcpy(argv[indx],arggv[indx+1]);
    }
    ParseAndDiscard(argv[1],"}");
    ParseAndDiscard(argv[1],"}");
    ParseAndDiscard(argv[1],"{");
    indx=ParseAndKeep(argv[1],"-,",ar1);
    atom1=atoi(ar1);
    if(indx==0)
    {
        atom2=atom1+1;
    }
    else
    {
        indx=ParseAndKeep(argv[1],",-}",ar1);
        atom2=atoi(ar1);
    }
    sprintf(tcl_ret,"%i %i",atom1,atom2);
    for(indx=0;indx<argc-1;indx++)
        free((void*)argv[indx]);
    free((void*)argv);
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int FreeArrays(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    FreeFCHKArrays();
    free((void*)Pmn);
    free((void*)Pmnd);
    free((void*)Smn);
    if(Mulliken!=NULL)
        free((void*)Mulliken);
    if(Lowdin!=NULL)
        free((void*)Lowdin);
    return TCL_OK;
}

int SaveGrid(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    float* data;
    float tval;
    int tival;
    FILE* of;

    data=(float*)atoi(argv[1]);
    of=fopen(argv[2],"wb");
    if(of==NULL)
    {
        printf("Error--Unable to open file %s\n",argv[2]);
        return TCL_ERROR;
    }
    tival=3;
    fwrite(&tival,1,sizeof(int),of);
    tival=25;
    fwrite(&tival,1,sizeof(int),of);
    tival=nz;
    fwrite(&tival,1,sizeof(int),of);
    tival=ny;
    fwrite(&tival,1,sizeof(int),of);
    tival=nz;
    fwrite(&tival,1,sizeof(int),of);
    tval=(float)(Box[0][2]*0.529);
    fwrite(&tval,1,sizeof(float),of);
    tval=(float)(Box[3][2]*0.529);
    fwrite(&tval,1,sizeof(float),of);
    tval=(float)(Box[0][1]*0.529);
    fwrite(&tval,1,sizeof(float),of);
    tval=(float)(Box[2][1]*0.529);
    fwrite(&tval,1,sizeof(float),of);
    tval=(float)(Box[0][0]*0.529);
    fwrite(&tval,1,sizeof(float),of);
    tval=(float)(Box[1][0]*0.529);
    fwrite(&tval,1,sizeof(float),of);
    fwrite(data,nx*ny*nz,sizeof(float),of);
    fclose(of);
    return TCL_OK;
}
    


int NewFillGrid(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    double CurrentPoint[3];
    int h,i,j,k,m,l;
    int orbid;
    double l1,l2,dp;
    double val,g1,g2;
    double gval[3],gv1[3],gv2[3];
    float* GridPtr;
    char ContID[512]="";
    int rv;
    int beta=0;
    int grad=0;
    int spin=0;
    int lplc=0;
    int dirid=0;
    int diff=0;
    int mesp=0;
    AtomicCore TestC;
    double rr;
    char tcl_ret[512];
    double tval[100];
    int nh,nm,c1,c2;
    int ch,cm;

    orbid=atoi(argv[1]);
    if(argc>2)
        sprintf(ContID,"%s",argv[2]);
    if(argc>3)
    {
        beta=((argv[3][0]=='b')||(argv[3][0]=='B'));
        diff=((argv[3][0]=='d')||(argv[3][0]=='D'));
    }
    if(argc>4)
    {
        grad=((argv[4][0]=='g')||(argv[4][0]=='G'));
        spin=3*((argv[4][0]=='s')||(argv[4][0]=='S'));
        lplc=((argv[4][0]=='l')||(argv[4][0]=='L'));
        mesp=((argv[4][0]=='m')||(argv[4][0]=='M'));
    }
    if((spin)&&(!diff))
    {
        spin=beta+1;
    }
    if(grad)
    {
        if(argc>5)
        {
            dirid=atoi(argv[5]);
        }
    }
    GridPtr=(float*)malloc(nx*ny*nz*sizeof(float));
    for(k=0;k<nz;k++)
    {
        printf("%i\n",k);
        for(j=0;j<ny;j++)
        {
            if(mesp)
                printf("%i\n",j);
            for(i=0;i<nx;i++)
            {
                CurrentPoint[0]=(double)i/(double)(nx-1)*(Box[1][0]-Box[0][0])+
                                    (double)j/(double)(ny-1)*(Box[2][0]-Box[0][0])+
                                    (double)k/(double)(nz-1)*(Box[3][0]-Box[0][0])
                                    +Box[0][0];
                CurrentPoint[1]=(double)i/(double)(nx-1)*(Box[1][1]-Box[0][1])+
                                    (double)j/(double)(ny-1)*(Box[2][1]-Box[0][1])+
                                    (double)k/(double)(nz-1)*(Box[3][1]-Box[0][1])
                                    +Box[0][1];
                CurrentPoint[2]=(double)i/(double)(nx-1)*(Box[1][2]-Box[0][2])+
                                    (double)j/(double)(ny-1)*(Box[2][2]-Box[0][2])+
                                    (double)k/(double)(nz-1)*(Box[3][2]-Box[0][2])
                                    +Box[0][2];
                if(orbid>=0)
                {
                    val=0.0;
                    for(h=0;h<NumGBasis;h++)
                    {
                        if(!beta)
                            val+=FullCoeffs[h][orbid]*Basis[h].Evaluate(CurrentPoint);
                        else
                            val+=FullCoeffs[h][orbid+NumIndepFun]*Basis[h].Evaluate(CurrentPoint);
                    }
                    *(GridPtr+(k*(nx*ny)+j*nx+i))=(float)val;
                }
                else
                {
                    if((!grad)&&(!spin)&&(!lplc)&&(!mesp))
                    {
                        val=0.0;
                        for(h=0;h<NumGBasis;h++)
                        {
                            g1=Basis[h].Evaluate(CurrentPoint);
                            for(m=0;m<NumGBasis;m++)
                            {
                                val+=Pmn[h*NumGBasis+m]*g1*Basis[m].Evaluate(CurrentPoint);
                            }
                        }
                        *(GridPtr+(k*(nx*ny)+j*nx+i))=(float)val;
                    }
                    else
                    {
                        if(spin)
                        {
                            if(spin==3)
                            {
                                val=0.0;
                                for(h=0;h<NumGBasis;h++)
                                {
                                    g1=Basis[h].Evaluate(CurrentPoint);
                                    for(m=0;m<NumGBasis;m++)
                                    {
                                        val+=Pmnd[h*NumGBasis+m]*g1*Basis[m].Evaluate(CurrentPoint);
                                    }
                                }   
                                *(GridPtr+(k*(nx*ny)+j*nx+i))=(float)val;
                            }
                            else if (spin==1)
                            {
                                val=0.0;
                                for(h=0;h<NumGBasis;h++)
                                {
                                    g1=Basis[h].Evaluate(CurrentPoint);
                                    for(m=0;m<NumGBasis;m++)
                                    {
                                        val+=(Pmn[h*NumGBasis+m]+Pmnd[h*NumGBasis+m])*g1*Basis[m].Evaluate(CurrentPoint);
                                    }
                                }   
                                *(GridPtr+(k*(nx*ny)+j*nx+i))=0.5f*(float)val;
                            }
                            else 
                            {
                                val=0.0;
                                for(h=0;h<NumGBasis;h++)
                                {
                                    g1=Basis[h].Evaluate(CurrentPoint);
                                    for(m=0;m<NumGBasis;m++)
                                    {
                                        val+=(Pmn[h*NumGBasis+m]-Pmnd[h*NumGBasis+m])*g1*Basis[m].Evaluate(CurrentPoint);
                                    }
                                }   
                                *(GridPtr+(k*(nx*ny)+j*nx+i))=0.5f*(float)val;
                            }
                        }
                        else
                        {
                            if(grad)
                            {
                                for(l=0;l<3;l++)
                                    gval[l]=0.0;
                                val=1.0;
                                for(h=0;h<NumGBasis;h++)
                                {
                                    g1=Basis[h].EvaluateGradient(CurrentPoint,gv1);
                                    for(m=0;m<NumGBasis;m++)
                                    {
                                        g2=Basis[m].EvaluateGradient(CurrentPoint,gv2);
                                        for(l=0;l<3;l++)
                                            gval[l]+=Pmn[h*NumGBasis+m]*(g1*gv2[l]+g2*gv1[l]);
                                    }
                                }
                                val=gval[dirid];
                                *(GridPtr+(k*(nx*ny)+j*nx+i))=(float)(val);
                            }
                            else if (lplc)
                            {
                                val=0.0;
                                for(h=0;h<NumGBasis;h++)
                                {
                                    g1=Basis[h].EvaluateGradient(CurrentPoint,gv1);
                                    l1=Basis[h].EvaluateLaplacian(CurrentPoint);
                                    for(m=0;m<NumGBasis;m++)
                                    {
                                        g2=Basis[m].EvaluateGradient(CurrentPoint,gv2);
                                        l2=Basis[m].EvaluateLaplacian(CurrentPoint);
                                        dp=0.0;
                                        for(l=0;l<3;l++)
                                            dp+=gv1[l]*gv2[l];
                                        val+=Pmn[h*NumGBasis+m]*(g2*l1+g1*l2+2.0*dp);
                                    }
                                }
                                *(GridPtr+(k*(nx*ny)+j*nx+i))=(float)(val);
                            }
                            else 
                            {
//                              printf("%i\n",i);
                                nstep=100;
                                TestC.NuclCharge=1.0;
                                val=0.0;
                                for(h=0;h<NumAtoms;h++)
                                {
                                    rr=0.0;
                                    for(m=0;m<3;m++)
                                    {
                                        rr+=(AtomList[h].Position[m]-CurrentPoint[m])*(AtomList[h].Position[m]-CurrentPoint[m]);
                                    }
                                    val+=(double)AtomList[h].AtNum/sqrt(rr);
                                }
                                for(h=0;h<3;h++)
                                    TestC.Position[h]=CurrentPoint[h];
/*                              for(h=0;h<NumGBasis;h++)
                                {
                                    for(m=0;m<NumGBasis;m++)
                                    {
                                        val-=Pmn[h*NumGBasis+m]*Basis[h].NormCoeff*Basis[m].NormCoeff*oldrecurse::NuclearIntegral(Basis[h],Basis[m],TestC);
                                    }
                                }*/
                                for(h=0;h<NumUncontShells;h++)
                                {
                                    nh=ShellStart[h+1]-ShellStart[h];
                                    for(m=0;m<NumUncontShells;m++)
                                    {
                                        nm=ShellStart[m+1]-ShellStart[m];
                                        newrecurse::NuclearIntegral(Basis+ShellStart[h],Basis+ShellStart[m],&TestC,tval);
                                        for(c1=0,ch=ShellStart[h];c1<nh;c1++,ch++)
                                        {
                                            for(c2=0,cm=ShellStart[m];c2<nm;c2++,cm++)
                                            {
                                                val-=Pmn[ch*NumGBasis+cm]*Basis[ch].NormCoeff*Basis[cm].NormCoeff*tval[c1*nm+c2];
                                            }
                                        }
                                    }
                                }
                                *(GridPtr+(k*(nx*ny)+j*nx+i))=(float)(val);
                            }
                        }
                    }
                }
            }
        }
    }
    rv=1;
    rv=gom_FillContourStructure("test1",ContID,nx,ny,nz,(float)(Box[0][0]*0.529),(float)(Box[1][0]*0.529),(float)(Box[0][1]*0.529),(float)(Box[2][1]*0.529),(float)(Box[0][2]*0.529),(float)(Box[3][2]*0.529),GridPtr);
    sprintf(tcl_ret,"%i",rv);
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int CalcCharges(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    double Qelec=0.0;
    int h,m;
    double *Smn12,*Smnp;
    double *Pmm;
    double TotElCountL;
    double TotElCountM;
    int CalcLowdin=0;
    int CalcMulliken=0;
    char tcl_ret[512];

    if(argc<2)
    {
        CalcMulliken=1;
    }
    else
    {
        if((argv[1][0]=='l')||(argv[1][0]=='L')||(argv[1][0]=='b')||(argv[1][0]=='B'))
        {
            CalcLowdin=1;
        }
        if((argv[1][0]=='m')||(argv[1][0]=='M')||(argv[1][0]=='b')||(argv[1][0]=='B'))
        {
            CalcMulliken=1;
        }
    }
    Smn12=(double*)malloc(NumGBasis*NumGBasis*sizeof(double));
    if(Smn12==NULL)
    {
      printf(Tcl_GetStringResult(ti),"Memory allocation error in transform matrix");
        return TCL_ERROR;
    }
    Smnp=(double*)malloc(NumGBasis*NumGBasis*sizeof(double));
    if(Smnp==NULL)
    {
      printf(Tcl_GetStringResult(ti),"Memory allocation error in transformed matrix");
        free((void*)Smn12);
        return TCL_ERROR;
    }
    Pmm=(double*)malloc(NumGBasis*sizeof(double));
    if(Pmm==NULL)
    {
      printf(Tcl_GetStringResult(ti),"Memory allocation error in Mulliken vector.");
        free((void*)Smn12);
        free((void*)Smnp);
        return TCL_ERROR;
    }
    Mulliken=(double*)malloc(NumAtoms*sizeof(double));
    Lowdin=(double*)malloc(NumAtoms*sizeof(double));
    if(CalcMulliken)
    {
        for(h=0;h<NumAtoms;h++)
            Mulliken[h]=(double)AtomList[h].AtNum;
        TotElCountM=0.0;
        for(h=0;h<NumGBasis;h++)
        {
            Pmm[h]=0.0;
            for(m=0;m<NumGBasis;m++)
            {
                Pmm[h]+=Pmn[h*NumGBasis+m]*Smn[h*NumGBasis+m];
            }
            Mulliken[OrbitalToAtom[h]]-=Pmm[h];
            TotElCountM+=Pmm[h];
        }
    }
    if(CalcLowdin)
    {
        MatSqrtS(Smn,Smn12,NumGBasis);
        MatMul(Smn12,Pmn,NumGBasis,NumGBasis,NumGBasis,NumGBasis,Smnp);
        for(h=0;h<NumAtoms;h++)
            Lowdin[h]=(double)AtomList[h].AtNum;
        TotElCountL=0.0;
        for(h=0;h<NumGBasis;h++)
        {
            Pmm[h]=0.0;
            for(m=0;m<NumGBasis;m++)
            {
                Pmm[h]+=Smnp[h*NumGBasis+m]*Smn12[h*NumGBasis+m];
            }
            Lowdin[OrbitalToAtom[h]]-=Pmm[h];
            TotElCountL+=Pmm[h];
        }
    }
    free((void*)Pmm);
    free((void*)Smn12);
    free((void*)Smnp);
    if(CalcLowdin)
    {
        if(CalcMulliken)
            sprintf(tcl_ret,"%lf %lf",TotElCountM,TotElCountL);
        else
            sprintf(tcl_ret,"%lf",TotElCountL);
    }
    else
    {
        if(CalcMulliken)
            sprintf(tcl_ret,"%lf",TotElCountM);
        else
            sprintf(tcl_ret,"No charge method specified");
    }
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int GetMullikenCharge(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    int atomid;
    char tcl_ret[512];

    atomid=atoi(argv[1]);
    sprintf(tcl_ret,"%lf",Mulliken[atomid]);
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int GetLowdinCharge(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    int atomid;
    char tcl_ret[512];

    atomid=atoi(argv[1]);
    sprintf(tcl_ret,"%lf",Lowdin[atomid]);
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int ListEnergies(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    char *retstr;
    int i;
    int ul;
    int ll=0;

    ul=NumIndepFun;
    retstr=(char*)malloc(32*ul);
    if(argc<2)
        return TCL_ERROR;
    int spin=atoi(argv[1]);
    if(argc>=3)
        ul=atoi(argv[2]);
    if(argc>=4)
        ll=atoi(argv[3]);
    if(ul<0)
    {
        if(spin)
        {
            ul=NumBetaElectrons;
        }
        else
        {
            ul=NumAlphaElectrons;
        }
    }
    else
    {
        ul=NumIndepFun;
    }
    retstr[0]='\000';
    if(spin)
    {
        for(i=ll;i<ul;i++)
        {
            sprintf(retstr,"%s %g",retstr,BetaEnergies[i]);
            if(i<NumBetaElectrons)
                sprintf(retstr,"%s(*)",retstr);
        }
    }
    else
    {
        for(i=ll;i<ul;i++)
        {
            sprintf(retstr,"%s %g",retstr,AlphaEnergies[i]);
            if(i<NumAlphaElectrons)
                sprintf(retstr,"%s(*)",retstr);
        }
    }
    Tcl_SetResult(ti,retstr,TCL_VOLATILE);
//  sprintf(ti->result,"%s",retstr);
    free((void*)retstr);
    return TCL_OK;
}

int GraphPlane(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    double result;
    char tcl_ret[512];

    result=CalcPlane(argc,argv);
    if(result<0.0)
        return TCL_ERROR;
    sprintf(tcl_ret,"%g",result);
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int ResetBox(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    SetInitialBox();
    return TCL_OK;
}

int ReportBox(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<3;j++)
            printf("%g ",Box[i][j]);
        printf("\n");
    }
    return TCL_OK;
}

int CompFileNames(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    char fn1[80],fn2[80];
    int i,result=0;
    char tcl_ret[512];

    sprintf(fn1,"%s",argv[1]);
    sprintf(fn2,"%s",argv[2]);
    i=0;
    while(fn1[i]!='\000')
    {
        if(fn1[i]!=fn2[i])
            result=1;
        i+=1;
    }
    if(!result)
    {
        if(fn2[i]!='\000')
            result=1;
    }
    sprintf(tcl_ret,"%i",result);
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int AdjustFileName(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    char fn1[80],ext[10];
    int i;
    char tcl_ret[512];

    sprintf(fn1,"%s",argv[1]);
    sprintf(ext,"%s",argv[2]);
    i=0;
    while((fn1[i]!='.')&&(fn1[i]!='\000'))
        i+=1;
    fn1[i]='\000';
    sprintf(tcl_ret,"%s.%s",fn1,ext);
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

int ReportDipoleMoment(ClientData cd,Tcl_Interp* ti,int argc,TCL_CONST char** argv)
{
    char tcl_ret[512];

    sprintf(tcl_ret,"%g %g %g",DipoleMoment[0],DipoleMoment[1],DipoleMoment[2]);
    Tcl_SetResult(ti,tcl_ret,TCL_VOLATILE);
    return TCL_OK;
}

} // namespace Fchk
} // namespace Plugin
