// symmetry.cpp : Defines the entry point for the DLL application.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <tcl.h>
#include <gopenmolext.h>

namespace Plugin {
namespace Filter {

static int ReadChem3DCoords(ClientData ,Tcl_Interp *,int ,const char **);
static int GetChem3DName(ClientData ,Tcl_Interp *,int ,const char **);
static int GetChem3DVec(ClientData ,Tcl_Interp *,int ,const char **);
static int FreeC3DPointers(ClientData ,Tcl_Interp *,int ,const char **);
static int ReadSpartanCoords(ClientData ,Tcl_Interp *,int ,const char **);


DYNEXPORT_C int Filters_Init(Tcl_Interp *Interp)
//declare exported command
{
    printf("Creating Tcl Extensions.\n");
Tcl_CreateCommand(Interp,"Chem3DImport",ReadChem3DCoords,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"releasec3d",FreeC3DPointers,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"getChem3DName",GetChem3DName,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"getChem3DVec",GetChem3DVec,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"SpartanImport",ReadSpartanCoords,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_PkgProvide(Interp,"symmetry","1.0");
    printf("Finished.\n");
return(TCL_OK);
}

#define DELTA 1.e-2

int lencmp(const char* c1,const char* c2,int minlen)
{
    int i;

    i=0;
    while(c1[i]!='\000')
    {
        if(c1[i]!=c2[i])
            return i+1;
        i++;
    }
    if(i<minlen)
        return minlen;
    return 0;
}

int cmp(const char* c1,const char* c2)
{
    int i;

    i=0;
    while(c1[i]!='\000')
    {
        if(c1[i]!=c2[i])
            return i+1;
        i++;
    }
    return 0;
}

class SpartanAtom
{
public:
    int AtomNumber;
    double pos[3];
    int IsValidAtom(char* instr);
};

char AtomLabel[103][3]={"H","He",
                        "Li","Be","B","C","N","O","F","Ne",
                        "Na","Mg","Al","Si","P","S","Cl","Ar",
                        "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
                        "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe"};


int SpartanAtom::IsValidAtom(char* instr)
{
    char ParseStr[47];
    int i;
    int nr;

    for(i=0;i<46;i++)
    {
        ParseStr[i]=instr[i];
    }
    ParseStr[46]='\000';
    nr=sscanf(ParseStr,"%i %lf %lf %lf",&AtomNumber,pos,pos+1,pos+2);
    return (nr==4);
}

int ExtractSpartanAtoms(char*** El,double*** rvec,FILE* infile)
{
    char buffer[15000];
    char* bufend;
    int i;
    int p;
    int NumAtom;
    SpartanAtom Atom[300];

    for(i=0;i<15000;i++)
        buffer[i]='\000';
    fread(buffer,7,1,infile);
    i=7;
    bufend=buffer;
    while(lencmp(bufend,"ENDCART",7)&&(i<15000))
    {
        fread(buffer+i,1,1,infile);
        bufend+=1;
        i+=1;
    }
    if(cmp(bufend,"ENDCART"))
    {
        while(!cmp(bufend,"ENDCART"))
        {
            for(i=0;i<14999;i++)
                buffer[i]=buffer[i-1];
            fread(buffer+15000,1,1,infile);
        }
    }
    bufend-=46;
    NumAtom=0;
    while(Atom[NumAtom].IsValidAtom(bufend))
    {
        NumAtom+=1;
        bufend-=46;
        printf("Atom type %i\n",Atom[NumAtom].AtomNumber);
    }
    El[0]=(char**)malloc(NumAtom*sizeof(char*));
    rvec[0]=(double**)malloc(NumAtom*sizeof(double*));
    for(i=0;i<NumAtom;i++)
    {
        El[0][i]=(char*)malloc(3*sizeof(char));
        rvec[0][i]=(double*)malloc(3*sizeof(double));
    }
    for(i=0;i<NumAtom;i++)
    {
        sprintf(El[0][i],"%s",AtomLabel[Atom[NumAtom-1-i].AtomNumber-1]);
        for(p=0;p<3;p++)
            rvec[0][i][p]=Atom[NumAtom-1-i].pos[p];
    }
    return NumAtom;
}

int ExtractAtoms(char*** El,double*** rvec,FILE* infile)
{
    char* c1;
    char c2;
    char nkrt[5];
    int NumAtoms;
    int i,j;
    int tvi;

    c1=&c2;
    for(i=0;i<=6;i++)
        fread(c1,1,1,infile);
    NumAtoms=(int)c2;
    printf("%i\n",NumAtoms);
    El[0]=(char**)malloc(NumAtoms*sizeof(char*));
    rvec[0]=(double**)malloc(NumAtoms*sizeof(double*));
    for(i=0;i<NumAtoms;i++)
    {
        El[0][i]=(char*)malloc(3*sizeof(char));
        rvec[0][i]=(double*)malloc(3*sizeof(double));
    }
    nkrt[4]='\000';
    nkrt[0]=nkrt[1]=nkrt[2]='P';
    c1=nkrt+3;
    fread(c1,1,1,infile);
    while((cmp(nkrt,"trac"))||(!cmp(nkrt,"")))
    {
        for(i=0;i<3;i++)
            nkrt[i]=nkrt[i+1];
        fread(c1,1,1,infile);
        printf("%s\n",nkrt);
    };
    for(i=0;i<4;i++)
        fread(c1,1,1,infile);
    for(j=0;j<3*NumAtoms;j++)
    {
        fread(&tvi,4,1,infile);
        rvec[0][j/3][j%3]=(double)(tvi)*1.526e-5;  //related to C-C bond in cyclohexane!!!
        printf("%i %i %i %g\n",j/3,j%3,tvi,rvec[0][j/3][j%3]);
    }
    while((cmp(nkrt,"bmys"))||(!cmp(nkrt,"")))
    {
        for(i=0;i<3;i++)
            nkrt[i]=nkrt[i+1];
        fread(c1,1,1,infile);
    };
    printf("%s\n",nkrt);
    fread(nkrt,3,1,infile);
    for(i=0;i<NumAtoms;i++)
    {
        fread(nkrt,4,1,infile);
        printf("%i %s\n",i,nkrt+2);
        El[0][i][0]=nkrt[2];
        El[0][i][1]=nkrt[3];
        El[0][i][2]='\000';
    }
    return NumAtoms;
}

double **C3DRvec;
char **C3DEl;


int ReadChem3DCoords(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    FILE* infile;
    int NAtom;
    char *result;
    
    infile=fopen(argv[1],"rb");
    printf("%s\n",argv[1]);
    NAtom=ExtractAtoms(&C3DEl,&C3DRvec,infile);
    fclose(infile);
    result = Tcl_Alloc(16);
    sprintf(result, "%i",NAtom);
    Tcl_SetResult(ti, result, TCL_DYNAMIC);
    return TCL_OK;
}

int ReadSpartanCoords(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    FILE* infile;
    int NAtom;
    char *result;

    infile=fopen(argv[1],"rb");
    printf("%s\n",argv[1]);
    NAtom=ExtractSpartanAtoms(&C3DEl,&C3DRvec,infile);
    fclose(infile);
    result = Tcl_Alloc(16);
    sprintf(result, "%i",NAtom);
    Tcl_SetResult(ti, result, TCL_DYNAMIC);
    return TCL_OK;
}

int GetChem3DName(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    int s,i;
    char *result;

    s=atoi(argv[1]);
    i=atoi(argv[2]);
    gom_PutAtomAtype(s-1,C3DEl[i-1],i-1);
    result = Tcl_Alloc(64);  // What is the expected max length of C3DEl strings?
    sprintf(result,"%s",C3DEl[i-1]);
    Tcl_SetResult(ti, result, TCL_DYNAMIC);
    printf("%s\n",C3DEl[i-1]);
    return TCL_OK;
}

int GetChem3DVec(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    int i;
    int coord;
    char *result;

    if(argc<3)
        return TCL_ERROR;
    i=atoi(argv[1]);
    coord=atoi(argv[2]);
    printf("%i %i\n",i,coord);
    result=Tcl_Alloc(16);
    sprintf(result,"%f7.5",C3DRvec[i-1][coord]);
    Tcl_SetResult(ti, result, TCL_DYNAMIC);
    return TCL_OK;
}

int FreeC3DPointers(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    int i;
    int NAtom;

    if(argc<2)
        return TCL_ERROR;
    NAtom=atoi(argv[1]);
    for(i=0;i<NAtom;i++)
    {
        free((void*)C3DEl[i]);
        free((void*)C3DRvec[i]);
    }
    free((void*)C3DEl);
    free((void*)C3DRvec);
    return TCL_OK;
}

} // namespace Filter
} // namespace Plugin
