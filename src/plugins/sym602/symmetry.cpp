// symmetry.cpp : Defines the entry point for the DLL application.
//

#include "symmetry.h"

namespace Plugin {
namespace Symmetry {

static int TellAtoms(ClientData ,Tcl_Interp *,int ,const char **);
static int ReportBackBone(ClientData ,Tcl_Interp *,int ,const char **);
static int VRMLBackBone(ClientData ,Tcl_Interp *,int ,const char **);
static int AlignStructures(ClientData ,Tcl_Interp *,int ,const char **);
static int ReadChem3DCoords(ClientData ,Tcl_Interp *,int ,const char **);
static int GetChem3DName(ClientData ,Tcl_Interp *,int ,const char **);
static int GetChem3DVec(ClientData ,Tcl_Interp *,int ,const char **);
static int FreeC3DPointers(ClientData ,Tcl_Interp *,int ,const char **);
static int ReportVectorRelationships(ClientData ,Tcl_Interp *,int ,const char **);

//declare exported command
DYNEXPORT_C int Symmetry_Init(Tcl_Interp *Interp)
{
    printf("Creating Tcl Extensions.\n");
Tcl_CreateCommand(Interp,"CheckAddresses",TellAtoms,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"getbackbone",ReportBackBone,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"vrmlbb",VRMLBackBone,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"alignmols",AlignStructures,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"Chem3DImport",ReadChem3DCoords,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"releasec3d",FreeC3DPointers,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"getChem3DName",GetChem3DName,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"getChem3DVec",GetChem3DVec,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_CreateCommand(Interp,"repvecs",ReportVectorRelationships,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
Tcl_PkgProvide(Interp,"symmetry","1.0");
    printf("Finished.\n");
return(TCL_OK);
}

#define DELTA 1.e-2

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
    int i;
    int NAtom;

    infile=fopen(argv[1],"rb");
    printf("%s\n",argv[1]);
    NAtom=ExtractAtoms(&C3DEl,&C3DRvec,infile);
    fclose(infile);
    sprintf(ti->result,"%i",NAtom);
    return TCL_OK;
}

int GetChem3DName(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    int s,i;

    s=atoi(argv[1]);
    i=atoi(argv[2]);
    gom_PutAtomAtype(s-1,C3DEl[i-1],i-1);
/*  switch(C3DEl[i-1][0])
    {
    case 'C':   gom_PutAtomType(s-1,1,i-1);
                gom_PutAtomBndRad(s-1,1.6,i-1);
                gom_PutAtomVdwRad(s-1,1.7,i-1);
                gom_PutAtomRmin(s-1,1.0,i-1);
                break;
    case 'N':   gom_PutAtomType(s-1,1,i-1);
                gom_PutAtomBndRad(s-1,1.6,i-1);
                gom_PutAtomVdwRad(s-1,1.7,i-1);
                gom_PutAtomRmin(s-1,1.0,i-1);
                break;
    case 'S':   gom_PutAtomType(s-1,1,i-1);
                gom_PutAtomBndRad(s-1,1.6,i-1);
                gom_PutAtomVdwRad(s-1,1.7,i-1);
                gom_PutAtomRmin(s-1,1.0,i-1);
                break;
    case 'H':   gom_PutAtomType(s-1,1,i-1);
                gom_PutAtomBndRad(s-1,1.3,i-1);
                gom_PutAtomVdwRad(s-1,1.5,i-1);
                gom_PutAtomRmin(s-1,0.7,i-1);
                break;
    case 'F':   gom_PutAtomType(s-1,1,i-1);
                gom_PutAtomBndRad(s-1,1.6,i-1);
                gom_PutAtomVdwRad(s-1,1.7,i-1);
                gom_PutAtomRmin(s-1,1.0,i-1);
                break;
    case 'O':   gom_PutAtomType(s-1,1,i-1);
                gom_PutAtomBndRad(s-1,1.6,i-1);
                gom_PutAtomVdwRad(s-1,1.7,i-1);
                gom_PutAtomRmin(s-1,1.0,i-1);
                break;
    default:    printf("Unrecognized element %s of atom %i\n",C3DEl[i-1],i);
                gom_PutAtomType(s-1,1,i-1);
                gom_PutAtomBndRad(s-1,1.6,i-1);
                gom_PutAtomVdwRad(s-1,1.7,i-1);
                gom_PutAtomRmin(s-1,1.0,i-1);
                break;
    }*/
    sprintf(ti->result,"%s",C3DEl[i-1]);
    printf("%s\n",C3DEl[i-1]);
    return TCL_OK;
}

int GetChem3DVec(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    int i;
    int coord;

    if(argc<3)
        return TCL_ERROR;
    i=atoi(argv[1]);
    coord=atoi(argv[2]);
    printf("%i %i\n",i,coord);
    sprintf(ti->result,"%f7.5",C3DRvec[i-1][coord]);
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


int NBBC;
float **BackBone=NULL;

int ReportBackBone(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    int i;

    if((argc<2)||(BackBone==NULL))
        return TCL_ERROR;
    i=atoi(argv[1]);
    printf("Atom %i of %i\n",i,NBBC);
    if(i>=NBBC)
        return TCL_ERROR;
    if(i<0)
    {
        for(i=0;i<NBBC;i++)
            free((void*)BackBone[i]);
        free((void*)BackBone);
    }
    else
    {
        sprintf(ti->result,"%g %g %g",BackBone[i][0],BackBone[i][1],BackBone[i][2]);
    }
    return TCL_OK;
}

int VRMLBackBone(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    int i;
    float rad;
    FILE* VRMLOut;
    double theta;
    float XSect[18][2];

    if((argc<3)||(BackBone==NULL))
        return TCL_ERROR;
    rad=(float)atof(argv[1]);
    for(i=0;i<18;i++)
    {
        theta=3.14159265359/9.0*i;
        XSect[i][0]=rad*(float)cos(theta);
        XSect[i][1]=rad*(float)sin(theta);
    }
    VRMLOut=fopen(argv[2],"w");
    fprintf(VRMLOut,"#VRML V2.0 utf8\n");
    fprintf(VRMLOut,"Transform {\n\tchildren [\n\t\tShape {\n\t\t\tappearance Appearance {\n");
    fprintf(VRMLOut,"\t\t\t\tmaterial Material {diffuseColor 0.1 0.1 1.0}\n\t\t\t}\n\t\t\t");
    fprintf(VRMLOut,"geometry Extrusion {\n\t\t\t\tcreaseAngle 2.8\n\t\t\t\tcrossSection [\n");
    for(i=0;i<18;i++)
    {
        fprintf(VRMLOut,"\t\t\t\t\t%g %g,\n",XSect[i][0],XSect[i][1]);
    }
    fprintf(VRMLOut,"\t\t\t\t]\n\t\t\t\tspine [\n");
    for(i=0;i<NBBC;i++)
    {
        fprintf(VRMLOut,"\t\t\t\t\t%g %g %g,\n",BackBone[i][0],BackBone[i][1],BackBone[i][2]);
    }
    fprintf(VRMLOut,"\t\t\t\t]\n\t\t\t}\n\t\t}\n\t]\n}\n");
    fclose(VRMLOut);
    return TCL_OK;
}

int AlignStructures(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    atomrec ar;
    char ts1[80];
    char ts2[80];
    char* tok;
    int i;

    if(argc<6)
        return TCL_ERROR;
    ar.NAtom=atoi(argv[1]);
    printf("%i\n%i\n",argc,ar.NAtom);
    ar.Str1=atoi(argv[2]);
    printf("%i\n",ar.Str1);
    ar.Str2=atoi(argv[3]);
    printf("%i\n",ar.Str2);
    ar.AtomList1=(int*)malloc(ar.NAtom*sizeof(int));
    ar.AtomList2=(int*)malloc(ar.NAtom*sizeof(int));
    sprintf(ts1,"%s",argv[4]);
    printf("%s\n",ts1);
    sprintf(ts2,"%s",argv[5]);
    printf("%s\n",ts2);
    i=0;
    tok=strtok(ts1," ,");
    if(tok!=NULL)
    {
        printf("%s\n",tok);
        ar.AtomList1[i++] = atoi(tok);
        printf("%i %i\n",i,ar.AtomList1[i-1]);
    }
    while ((i<ar.NAtom)&&(tok!=NULL))
    {
        tok=strtok(NULL," ,");
        if(tok!=NULL)
        {
            printf("%s\n",tok);
            ar.AtomList1[i++] = atoi(tok);
            printf("%i %i\n",i,ar.AtomList1[i-1]);
        }
    }
    i=0;
    tok=strtok(ts2," ,");
    if(tok!=NULL)
    {
        printf("%s\n",tok);
        ar.AtomList2[i++] = atoi(tok);
        printf("%i %i\n",i,ar.AtomList2[i-1]);
    }
    while ((i<ar.NAtom)&&(tok!=NULL))
    {
        tok=strtok(NULL," ,");
        if(tok!=NULL)
        {
            printf("%s\n",tok);
            ar.AtomList2[i++] = atoi(tok);
            printf("%i %i\n",i,ar.AtomList2[i-1]);
        }
    }
    if(argc<7)
    {
        ar.FixedAtom1=-1;
    }
    else
    {
        ar.FixedAtom1=atoi(argv[6]);
    }
    if(argc<8)
    {
        ar.FixedAtom2=-1;
    }
    else
    {
        ar.FixedAtom2=atoi(argv[7]);
    }
    printf("%i\t%i\tStructures\n",ar.Str1,ar.Str2);
    for(i=0;i<ar.NAtom;i++)
    {
        printf("%i\t%i\n",ar.AtomList1[i],ar.AtomList2[i]);
    }
    if((ar.FixedAtom1>=0)&&(ar.FixedAtom2>=0))
        printf("%i\t%i\n",ar.FixedAtom1,ar.FixedAtom2);
    AlignStructs(&ar);
    free((void*)ar.AtomList1);
    free((void*)ar.AtomList2);
    return TCL_OK;
}

int TellAtoms(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    int i,j,k,l;
    int NStruct;
    int NAt;
    float *xc,*yc,*zc,*rc;
    int *Z;
    int *ccl;
    int ncl=0;
    int xcl;
    int first;
    const int *Res;
    int MaxRes=0;
    char *ResName;
    char *AtomName;
    int ResNum;
    int lfl;
    float w,mw;
    float wx,wy,wz;
    float xcm,ycm,zcm;
    const float *CPKRadScale;

    NStruct=gom_GetNumMolecStructs();
    printf("Found %i Structures\n",NStruct);
    for(i=0;i<NStruct;i++)
    {
        CPKRadScale=gom_GetAtomCPKScalePointer(i);
        NAt=gom_GetNumAtomsInMolecStruct(i);
        Res=gom_GetAtomResNum1Pointer(i);
        xc=(float*)malloc(NAt*sizeof(float));
        yc=(float*)malloc(NAt*sizeof(float));
        zc=(float*)malloc(NAt*sizeof(float));
        rc=(float*)malloc(NAt*sizeof(float));
        Z=(int*)malloc(NAt*sizeof(int));
        ccl=(int*)malloc(NAt*sizeof(int));
        printf("Structure %i has %i atoms\n",NStruct,NAt);
        mw=0.0;
        wx=wy=wz=0.0;
        for(j=0;j<NAt;j++)
        {
            xc[j]=gom_GetAtomXCoord(i,j);
            yc[j]=gom_GetAtomYCoord(i,j);
            zc[j]=gom_GetAtomZCoord(i,j);
            Z[j]=(int)gom_GetAtomNucCharge(i,j);
            ccl[j]=-1;
            w=gom_GetAtomMass(i,j);
            mw+=w;
            wx+=xc[j]*w;
            wy+=yc[j]*w;
            wz+=zc[j]*w;
            if(Res[j]>MaxRes)
                MaxRes=Res[j];
/*          printf("Atom %i at %g,%g,%g has atomic number %i\n",j,xc[j],yc[j],zc[j],Z[j]);
            printf("Atom %i belongs to Residue %i \n",j,Res[j]);*/
        }
        printf("Structure %i has %i residues\n",i,MaxRes);
//      PauseForAction();
/*      for(j=0;j<MaxRes/18;j++)
        {
            lfl=0;
            for(k=0;k<18;k++)
            {
                ResNum=(j*18)+k;
                l=0;
                while((l<NAt)&&(Res[l]!=ResNum))
                    l++;
                if(Res[l]==ResNum)
                {
                    ResName=gom_GetAtomResName(i,l);
                    printf("%s ",ResName);
                    lfl=1;
                }
            }
            if(lfl)
                printf("\n");
        }
        PauseForAction();
        NBBC=0;
        for(j=0;j<NAt;j++)
        {
            AtomName=gom_GetAtomAtmName(i,j);
            if(!cmp(AtomName,"CA"))
            {
                printf("Backbone carbon %i\n",j);
                NBBC+=1;
            }
        }
        printf("%i BackBone carbon atoms.\n",NBBC);
        BackBone=(float**)malloc(NBBC*sizeof(float*));
        for(j=0;j<NBBC;j++)
            BackBone[j]=(float*)malloc(3*sizeof(float));
        lfl=0;
        for(j=0;j<NAt;j++)
        {
            AtomName=gom_GetAtomAtmName(i,j);
            if(!cmp(AtomName,"CA"))
            {
                BackBone[lfl][0]=gom_GetAtomXCoord(i,j);
                BackBone[lfl][1]=gom_GetAtomYCoord(i,j);
                BackBone[lfl++][2]=gom_GetAtomZCoord(i,j);
                printf("BackBone Carbon %i at %g %g %g\n",lfl,BackBone[lfl-1][0],BackBone[lfl-1][1],BackBone[lfl-1][2]);
            }
        }*/
/*      xcm=wx/mw;
        ycm=wy/mw;
        zcm=wz/mw;
        printf("Center of mass is at %g,%g,%g\n",xcm,ycm,zcm);
        for(j=0;j<NAt;j++)
        {
            xc[j]-=xcm;
            yc[j]-=ycm;
            zc[j]-=zcm;
            rc[j]=(float)sqrt(xc[j]*xc[j]+yc[j]*yc[j]+zc[j]*zc[j]);
        }
        for(j=0;j<NAt-1;j++)
        {
            if(ccl[j]==-1)
            {
                ncl++;
                ccl[j]=j;
                for(k=j+1;k<NAt;k++)
                {
                    if((fabs(rc[j]-rc[k])<DELTA)&&(Z[j]==Z[k]))
                    {
                        printf("Possible pair %i %i\n",j,k);
                        ccl[k]=j;
                    }
                }
            }
        }
        xcl=1;
        for(j=0;j<NAt;j++)
        {
            first=0;
            for(k=0;k<NAt;k++)
            {
                if(ccl[k]==j)
                {
                    if(first)
                    {
                        printf("%i, ",k);
                        CPKRadScale[k]=0.1f*xcl;
                    }
                    else
                    {
                        printf("\nCluster %i contains atoms %i, ",j,k);
                        CPKRadScale[k]=0.1f*++xcl;
                        first=1;
                    }
                }
            }
            printf("\n");
        }
        printf("%i Clusters found\n",ncl);*/
        free((void*)xc);
        free((void*)yc);
        free((void*)zc);
        free((void*)rc);
        free((void*)Z);
    }
    sprintf(ti->result,"%i",NBBC);
    return TCL_OK;
}

int ReportVectorRelationships(ClientData cd,Tcl_Interp *ti,int argc,const char** argv)
{
    double r1n[3],r2n[3],r1c[3],r2c[3],r1d[3],r2d[3];
    double q11[3],q12[3],q21[3],q22[3];
    double l11,l12,l21,l22;
    double h;
    double thhinge,thnorm;
    double degrad=57.2957795130823208767981548141052;
    double dotp;
    int at2[7];
    int mol;
    int i;

    if(argc<9)
        return TCL_ERROR;
    mol=atoi(argv[1])-1;
    for(i=0;i<7;i++)
        at2[i]=atoi(argv[i+2])-1;
    r1c[0]=(gom_GetAtomXCoord(mol,at2[0])+gom_GetAtomXCoord(mol,at2[1])+gom_GetAtomXCoord(mol,at2[2]))/3.0;
    r2c[0]=(gom_GetAtomXCoord(mol,at2[4])+gom_GetAtomXCoord(mol,at2[5])+gom_GetAtomXCoord(mol,at2[6]))/3.0;
    q11[0]=gom_GetAtomXCoord(mol,at2[1])-gom_GetAtomXCoord(mol,at2[0]);
    q12[0]=gom_GetAtomXCoord(mol,at2[2])-gom_GetAtomXCoord(mol,at2[0]);
    q21[0]=gom_GetAtomXCoord(mol,at2[5])-gom_GetAtomXCoord(mol,at2[4]);
    q22[0]=gom_GetAtomXCoord(mol,at2[6])-gom_GetAtomXCoord(mol,at2[4]);
    r1d[0]=gom_GetAtomXCoord(mol,at2[3])-r1c[0];
    r2d[0]=r2c[0]-gom_GetAtomXCoord(mol,at2[3]);
    r1c[1]=(gom_GetAtomYCoord(mol,at2[0])+gom_GetAtomYCoord(mol,at2[1])+gom_GetAtomYCoord(mol,at2[2]))/3.0;
    r2c[1]=(gom_GetAtomYCoord(mol,at2[4])+gom_GetAtomYCoord(mol,at2[5])+gom_GetAtomYCoord(mol,at2[6]))/3.0;
    q11[1]=gom_GetAtomYCoord(mol,at2[1])-gom_GetAtomYCoord(mol,at2[0]);
    q12[1]=gom_GetAtomYCoord(mol,at2[2])-gom_GetAtomYCoord(mol,at2[0]);
    q21[1]=gom_GetAtomYCoord(mol,at2[5])-gom_GetAtomYCoord(mol,at2[4]);
    q22[1]=gom_GetAtomYCoord(mol,at2[6])-gom_GetAtomYCoord(mol,at2[4]);
    r1d[1]=gom_GetAtomYCoord(mol,at2[3])-r1c[1];
    r2d[1]=r2c[1]-gom_GetAtomYCoord(mol,at2[3]);
    r1c[2]=(gom_GetAtomZCoord(mol,at2[0])+gom_GetAtomZCoord(mol,at2[1])+gom_GetAtomZCoord(mol,at2[2]))/3.0;
    r2c[2]=(gom_GetAtomZCoord(mol,at2[4])+gom_GetAtomZCoord(mol,at2[5])+gom_GetAtomZCoord(mol,at2[6]))/3.0;
    q11[2]=gom_GetAtomZCoord(mol,at2[1])-gom_GetAtomZCoord(mol,at2[0]);
    q12[2]=gom_GetAtomZCoord(mol,at2[2])-gom_GetAtomZCoord(mol,at2[0]);
    q21[2]=gom_GetAtomZCoord(mol,at2[5])-gom_GetAtomZCoord(mol,at2[4]);
    q22[2]=gom_GetAtomZCoord(mol,at2[6])-gom_GetAtomZCoord(mol,at2[4]);
    r1d[2]=gom_GetAtomZCoord(mol,at2[3])-r1c[2];
    r2d[2]=r2c[2]-gom_GetAtomZCoord(mol,at2[3]);
    l11=l12=l21=l22=0.0;
    for(i=0;i<3;i++)
    {
        l11+=q11[i]*q11[i];
        l21+=q21[i]*q21[i];
        l12+=q12[i]*q12[i];
        l22+=q22[i]*q22[i];
    }
    l11=1.0/sqrt(l11);
    l21=1.0/sqrt(l21);
    l12=1.0/sqrt(l12);
    l22=1.0/sqrt(l22);
    for(i=0;i<3;i++)
    {
        q11[i]*=l11;
        q21[i]*=l21;
        q12[i]*=l12;
        q22[i]*=l22;
    }
    r1n[0]=q11[1]*q12[2]-q12[1]*q11[2];
    r1n[1]=q11[2]*q12[0]-q12[2]*q11[0];
    r1n[2]=q11[0]*q12[1]-q12[0]*q11[1];
    dotp=0.0;
    for(i=0;i<3;i++)
        dotp+=r1n[i]*r1d[i];
    if(dotp<0.0)
    {
        for(i=0;i<3;i++)
            r1n[i]*=-1.0;
    }
    r2n[0]=q21[1]*q22[2]-q22[1]*q21[2];
    r2n[1]=q21[2]*q22[0]-q22[2]*q21[0];
    r2n[2]=q21[0]*q22[1]-q22[0]*q21[1];
    dotp=0.0;
    for(i=0;i<3;i++)
        dotp+=r2n[i]*r2d[i];
    if(dotp<0.0)
    {
        for(i=0;i<3;i++)
            r2n[i]*=-1.0;
    }
    h=0.0;
    for(i=0;i<3;i++)
        h+=(r2c[i]-r1c[i])*r1n[i];
    dotp=0.0;
    for(i=0;i<3;i++)
        dotp+=r1n[i]*r2n[i];
    thnorm=atan2(sqrt(1.0-dotp*dotp),dotp);
    l11=0.0;
    l22=0.0;
    dotp=0.0;
    for(i=0;i<3;i++)
    {
        l11+=r1d[i]*r1d[i];
        l22+=r2d[i]*r2d[i];
        dotp+=r2d[i]*r1d[i];
    }
    dotp/=sqrt(l11*l22);
    thhinge=atan2(sqrt(1.0-dotp*dotp),dotp);
    printf("Height %g\nNormal tilt %g\nHinge %g\n",h,thnorm*degrad,thhinge*degrad);
    return TCL_OK;
}

} // namespace Symmetry
} // namespace Plugin
