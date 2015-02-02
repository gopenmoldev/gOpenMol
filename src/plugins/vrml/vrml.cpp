// vrml.cpp : Defines the entry point for the DLL application.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <tcl.h>
#include <gopenmolext.h>

namespace Plugin {
namespace Vrml {

#define pi 3.14159265359f

static int SetDisplay(ClientData ,Tcl_Interp *,int ,const char **);
static int CreateMovie(ClientData,Tcl_Interp *,int, const char **);
static int AddFrame(ClientData,Tcl_Interp *,int, const char **);
static int CloseMovie(ClientData,Tcl_Interp *,int, const char **);

DYNEXPORT_C int Vrml_Init(Tcl_Interp *Interp)
//declare exported command
{
    printf("Creating Tcl Extensions.\n");
    Tcl_CreateCommand(Interp,"VRMLStill",SetDisplay,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(Interp,"VRMLmovie",CreateMovie,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(Interp,"VRMLframe",AddFrame,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(Interp,"VRMLwrite",CloseMovie,(ClientData)NULL,(Tcl_CmdDeleteProc *)NULL);
    Tcl_PkgProvide(Interp,"VRMLMake","1.0");
    printf("Finished.\n");
return(TCL_OK);
}

void GetCoord(int WStr,int Wat,float* r)
{
    r[0]=gom_GetAtomXCoord(WStr,Wat);
    r[1]=gom_GetAtomYCoord(WStr,Wat);
    r[2]=gom_GetAtomZCoord(WStr,Wat);
}

void WriteVRMLHeader(FILE* of)
{
    fprintf(of,"#VRML V2.0 utf8\n");
    fprintf(of,"#created using VRMLtool gOpenMol Extension\n");
    fprintf(of,"#written by Kevin J. Boyd, University of New Orleans\n");
}

void AddSphere(FILE* of,float* q,float r,float* color)
{
    fprintf(of,"Transform {\n\ttranslation %5.3f %5.3f %5.3f\n\t",q[0],q[1],q[2]);
    fprintf(of,"children [\n\t\t");
    fprintf(of,"Shape {\n\t\t\tgeometry Sphere {radius %5.3f}\n",r);
    fprintf(of,"\t\t\t\tappearance Appearance {\n\t\t\t\t");
    fprintf(of,"material Material {diffuseColor %4.2f %4.2f %4.2f}\n\t\t\t}\n\t\t}\n\t]\n}\n",color[0],color[1],color[2]);
}

void AddCylinder(FILE* of,float* q1,float* q2,float r,float* c1,float* c2)
{
    float xy[21][2];
    int i;
    float q3[3];
    float theta;

    for(i=0;i<3;i++)
    {
        q3[i]=(q1[i]+q2[i])*0.5f;
    }
    for(i=0;i<21;i++)
    {
        theta=pi*(float)i*0.1f;
        xy[i][0]=r*(float)sin((double)theta);
        xy[i][1]=r*(float)cos((double)theta);
//      printf("%i %g %g %g\n",i,theta,xy[i][0],xy[i][1]);
    }
//  PauseForAction();
    fprintf(of,"Shape {\n",r);
    fprintf(of,"\tappearance Appearance {\n\t\t");
    fprintf(of,"material Material {diffuseColor %4.2f %4.2f %4.2f}\n\t}\n",c1[0],c1[1],c1[2]);
    fprintf(of,"\tgeometry Extrusion{\n\t\t");
    fprintf(of,"crossSection [\n");
    for(i=0;i<21;i++)
    {
        fprintf(of,"\t\t\t%6.4f %6.4f,\n",xy[i][0],xy[i][1]);
    }
    fprintf(of,"\t\t]\n\t\tspine [\n\t\t\t%6.4f %6.4f %6.4f,\n\t\t\t%6.4f %6.4f %6.4f",q1[0],q1[1],q1[2],q3[0],q3[1],q3[2]);
    fprintf(of,"\n\t\t]\n\t}\n}\n");
    fprintf(of,"Shape {\n",r);
    fprintf(of,"\tappearance Appearance {\n\t\t");
    fprintf(of,"material Material {diffuseColor %4.2f %4.2f %4.2f}\n\t}\n",c2[0],c2[1],c2[2]);
    fprintf(of,"\tgeometry Extrusion{\n\t\t");
    fprintf(of,"crossSection [\n");
    for(i=0;i<21;i++)
    {
        fprintf(of,"\t\t\t%6.4f %6.4f,\n",xy[i][0],xy[i][1]);
    }
    fprintf(of,"\t\t]\n\t\tspine [\n\t\t\t%6.4f %6.4f %6.4f,\n\t\t\t%6.4f %6.4f %6.4f",q3[0],q3[1],q3[2],q2[0],q2[1],q2[2]);
    fprintf(of,"\n\t\t]\n\t}\n}\n");
}

int PickFNums(int NNum,char* rstr,float* vals)
{
    int i=0;
    int j=0;
    int nv=0;
    char ts[32];
    double* dvals;

    dvals=(double*)malloc(NNum*sizeof(double));
    while((rstr[i]!='\000')&&(nv<NNum))
    {
        ts[j++]=rstr[i++];
        if((ts[j-1]=='\040')||(ts[j-1]=='\000')||(ts[j-1]=='\n')||(ts[j-1]=='\r'))
        {
            ts[j]='\000';
            dvals[nv++]=atof(ts);
            j=0;
        }
    }
    if(nv<NNum)
    {
        ts[j]=rstr[i];
        dvals[nv++]=atof(ts);
    }
    for(i=0;i<nv;i++)
        vals[i]=(float)dvals[i];
    free((void*)dvals);
    return nv;
}


int SetDisplay(ClientData cd,Tcl_Interp *ti,int argc,const char **argv)
{
    int numstruct;
    int numatom;
    int i,j,k;
    float r0[3],r1[3];
    float sr,cr;
    float cpkr;
    char q1[2];
    int nn;
    const int *nei;
    char filename[80];
    int bands;
    FILE* VRMLfile;
    float col[3],col2[3];
    char comstr[80];

    if(argc<3)
        return TCL_ERROR;
    bands=atoi(argv[2]);
    sprintf(filename,"%s",argv[1]);
    VRMLfile=fopen(filename,"w");
    WriteVRMLHeader(VRMLfile);
    q1[1]='\000';
    numstruct=gom_GetNumMolecStructs();
    for(i=0;i<numstruct;i++)
    {
        numatom=gom_GetNumAtomsInMolecStruct(i);
        cr=gom_GetAtomLicoRadC(i);
        sr=gom_GetAtomLicoRadS(i);
        if(bands)
        {
            cr=0.05f;
            printf("Using sticks.\n");
            if(bands<2)
            {
                sr=0.05f;
                printf("Not using balls.\n");
            }
        }
        for(j=0;j<numatom;j++)
        {
            sprintf(comstr,"show atom color %i %i",j+1,i+1);
            Tcl_Eval(ti,comstr);
            PickFNums(3,ti->result,col);
            nei=gom_GetAtomConnection(i,j);
            nn=nei[0];
//          printf("\n");
            if(gom_GetAtomDisplayState(i,j))
            {
                GetCoord(i,j,r0);
                for(k=1;k<=nn;k++)
                {
                    if(gom_GetAtomDisplayState(i,nei[k]))
                    {
                        GetCoord(i,nei[k],r1);
                        if(j<nei[k])
                        {
                            sprintf(comstr,"show atom color %i %i",nei[k]+1,i+1);
                            printf("%s %i %i\n",ti->result,nei[k]+1,i+1);
                            Tcl_Eval(ti,comstr);
                            PickFNums(3,ti->result,col2);
                            AddCylinder(VRMLfile,r0,r1,cr,col,col2);
                        }
                    }
                }
                if(gom_GetAtomCPKDisplayState(i,j))
                {
                    cpkr=gom_GetAtomVdwRad(i,j)*gom_GetAtomCPKScale(i,j);
                    AddSphere(VRMLfile,r0,cpkr,col);
                }
                else
                {
                        AddSphere(VRMLfile,r0,sr,col);
                }
            }
            printf("\n");
//          if((j%20==0)&&(j!=0))
//              PauseForAction();
        }
    }
    fclose(VRMLfile);
    return TCL_OK;
}

class MovieFrame 
{
public:
    int NumAtoms;
    int NumBonds;
    int **Bonds;
    float **Rat;
};

MovieFrame *Frames;
int NumFrames;
int BandS;
float maxX,minX,maxY,minY,maxZ,minZ;

int CreateMovie(ClientData cd,Tcl_Interp *ti,int argc,const char **argv)
{
    int i;

    if(argc<3)
        return TCL_ERROR;
    NumFrames=atoi(argv[1]);
    BandS=atoi(argv[2]);
    printf("%i\n",BandS);
    Frames=(MovieFrame*)malloc(NumFrames*sizeof(MovieFrame));
    for(i=0;i<NumFrames;i++)
    {
        Frames[i].Bonds=NULL;
        Frames[i].Rat=NULL;
    }
    maxX=maxY=maxZ=-1.e6;
    minX=minY=minZ=1.e6;
    return TCL_OK;
}

int AddFrame(ClientData cd,Tcl_Interp *ti,int argc,const char **argv)
{
    int i,j,k;
    const int *cnct;
    int FrameID;

    if(argc<2)
        return TCL_ERROR;
    FrameID=atoi(argv[1])-1;
    Frames[FrameID].NumAtoms=gom_GetNumAtomsInMolecStruct(0);
    printf("Adding frame %i, %i atoms\n",FrameID,Frames[FrameID].NumAtoms);
//  PauseForAction();
    Frames[FrameID].NumBonds=0;
    Frames[FrameID].Rat=(float**)malloc(Frames[FrameID].NumAtoms*sizeof(float*));
    for(i=0;i<Frames[FrameID].NumAtoms;i++)
    {
        Frames[FrameID].Rat[i]=(float*)malloc(3*sizeof(float));
        Frames[FrameID].Rat[i][0]=gom_GetAtomXCoord(0,i);
        Frames[FrameID].Rat[i][1]=gom_GetAtomYCoord(0,i);
        Frames[FrameID].Rat[i][2]=gom_GetAtomZCoord(0,i);
        if(Frames[FrameID].Rat[i][0]>maxX)
            maxX=Frames[FrameID].Rat[i][0];
        if(Frames[FrameID].Rat[i][0]<minX)
            minX=Frames[FrameID].Rat[i][0];
        if(Frames[FrameID].Rat[i][1]>maxY)
            maxY=Frames[FrameID].Rat[i][1];
        if(Frames[FrameID].Rat[i][1]<minY)
            minY=Frames[FrameID].Rat[i][1];
        if(Frames[FrameID].Rat[i][2]>maxZ)
            maxZ=Frames[FrameID].Rat[i][2];
        if(Frames[FrameID].Rat[i][2]<minZ)
            minZ=Frames[FrameID].Rat[i][2];
        cnct=gom_GetAtomConnection(0,i);
        Frames[FrameID].NumBonds+=cnct[0];
    }
    Frames[FrameID].NumBonds/=2;
    Frames[FrameID].Bonds=(int**)malloc(Frames[FrameID].NumBonds*sizeof(int*));
    for(i=0;i<Frames[FrameID].NumBonds;i++)
        Frames[FrameID].Bonds[i]=(int*)malloc(2*sizeof(int));
    k=0;
    for(i=0;i<Frames[FrameID].NumAtoms;i++)
    {
        cnct=gom_GetAtomConnection(0,i);
        for(j=1;j<=cnct[0];j++)
        {
            if(cnct[j]>i)
            {
                Frames[FrameID].Bonds[k][0]=i;
                Frames[FrameID].Bonds[k][1]=cnct[j];
                k++;
            }
        }
    }
    return TCL_OK;
}

void WriteVRMLTimerObject(FILE* of,float Delay,float RMin)
{
    float bcx,bcy,bcz;

    if(0.1f*(maxX-minX)>RMin)
        bcx=maxX+0.1f*(maxX-minX);
    else
        bcx=maxX+RMin;
    if(0.1f*(maxY-minY)>RMin)
        bcy=maxY+0.1f*(maxY-minY);
    else
        bcy=maxY+RMin;
    if(0.1f*(maxZ-minZ)>RMin)
        bcz=maxZ+0.1f*(maxZ-minZ);
    else
        bcz=maxZ+RMin;
    float br=0.02f*(float)sqrt((maxX-minX)*(maxX-minX)+(maxY-minY)*(maxY-minY)+(maxZ-minZ)*(maxZ-minZ));
    Delay/=1000.0;
    fprintf(of,"DEF Timer TimeSensor{\n\tloop TRUE\n\tstartTime 0.0\n\tstopTime 15.0\n\t");
    fprintf(of,"cycleInterval  %5.4f\n}\n",Delay);
    fprintf(of,"Transform {\n\ttranslation %5.3f %5.3f %5.3f\n\tchildren [\n",bcx,bcy,bcz);
    fprintf(of,"\tDEF RedBall TouchSensor{}\n");
    fprintf(of,"\tShape {\n");
    fprintf(of,"\t\tappearance Appearance {\n\t\t\tmaterial Material {diffuseColor 1.0 0.0 0.0}\n");
    fprintf(of,"\t\t}\n\t\tgeometry Sphere {\n\t\t\tradius %5.3f\n\t\t}\n\t}\n\t]\n}\n\n",br);
    if(0.1f*(maxX-minX)>RMin)
        bcx=minX-0.1f*(maxX-minX);
    else
        bcx=minX-RMin;
    if(0.1f*(maxY-minY)>RMin)
        bcy=minY-0.1f*(maxY-minY);
    else
        bcy=minY-RMin;
    if(0.1f*(maxZ-minZ)>RMin)
        bcz=minZ-0.1f*(maxZ-minZ);
    else
        bcz=minZ-RMin;
    fprintf(of,"Transform {\n\ttranslation %5.3f %5.3f %5.3f\n\tchildren [\n",bcx,bcy,bcz);
    fprintf(of,"\tDEF RedBall2 TouchSensor{}\n");
    fprintf(of,"\tShape {\n");
    fprintf(of,"\t\tappearance Appearance {\n\t\t\tmaterial Material {diffuseColor 1.0 0.0 0.0}\n");
    fprintf(of,"\t\t}\n\t\tgeometry Sphere {\n\t\t\tradius %5.3f\n\t\t}\n\t}\n\t]\n}\n\n",br);
    fprintf(of,"DEF ToggleScript Script{\n");
    fprintf(of,"\teventIn SFTime TTime\n");
    fprintf(of,"\tfield SFBool Running FALSE\n");
    fprintf(of,"\teventOut SFTime GoTime\n");
    fprintf(of,"\teventOut SFTime StopTime\n");
    fprintf(of,"\turl \"javascript:\n");
    fprintf(of,"\t\tfunction initialize() {\n");
    fprintf(of,"\t\t}\n");
    fprintf(of,"\t\tfunction TTime (event,value) {\n");
    fprintf(of,"\t\t\tif(Running) {\n");
    fprintf(of,"\t\t\t\tRunning=FALSE;\n\t\t\t\tStopTime=value;\n\t\t\t}\n");
    fprintf(of,"\t\t\telse {\n");
    fprintf(of,"\t\t\t\tRunning=TRUE;\n\t\t\t\tGoTime=value\n\t\t\t}\n");
    fprintf(of,"\t\t}\n\t\"\n}\n\n");
}

void WriteVRMLScriptObject(FILE* of)
{
    int i,j;
    char tstr[32],tstr2[32],tstr3[32];
    float x,y,z;

    fprintf(of,"DEF AnimScript Script{\n\t");
    fprintf(of,"eventIn SFTime TTime\n");   
    fprintf(of,"\tfield SFInt32 NumFrames %i\n",NumFrames-1);
    fprintf(of,"\tfield SFInt32 FrameID 0\n");
    printf("%i frames; %i atoms\n",NumFrames,Frames[0].NumAtoms);
    switch(BandS)
    {
    case 0: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPos",i);
                fprintf(of,"\tfield MFVec3f %s [",tstr);
                for(j=0;j<NumFrames-1;j++)
                {
                    fprintf(of,"%6.4f %6.4f %6.4f,\n\t\t\t",Frames[j].Rat[i][0],Frames[j].Rat[i][1],Frames[j].Rat[i][2]);
                }
                fprintf(of,"%6.4f %6.4f %6.4f\n\t\t\t",Frames[j].Rat[i][0],Frames[j].Rat[i][1],Frames[j].Rat[i][2]);
                fprintf(of,"\t\t]\n");
            }
            for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                fprintf(of,"\teventOut SFVec3f %s\n",tstr);
            }
            break;
    case 1: 
    case 2:
    case 3: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPos",i);
                fprintf(of,"\tfield MFVec3f %s [",tstr);
                for(j=0;j<NumFrames-1;j++)
                {
                    fprintf(of,"%6.4f %6.4f %6.4f,\n\t\t\t",Frames[j].Rat[i][0],Frames[j].Rat[i][1],Frames[j].Rat[i][2]);
                }
                fprintf(of,"%6.4f %6.4f %6.4f\n\t\t\t",Frames[j].Rat[i][0],Frames[j].Rat[i][1],Frames[j].Rat[i][2]);
                fprintf(of,"\t\t]\n");
            }
            for(i=0;i<Frames[0].NumBonds;i++)
            {
                sprintf(tstr,"Bond%0iPos0",i);
                fprintf(of,"\tfield MFVec3f %s [\n",tstr);
                for(j=0;j<NumFrames-1;j++)
                {
                    x=Frames[j].Rat[Frames[j].Bonds[i][0]][0];
                    y=Frames[j].Rat[Frames[j].Bonds[i][0]][1];
                    z=Frames[j].Rat[Frames[j].Bonds[i][0]][2];
                    fprintf(of,"%6.4f %6.4f %6.4f,\n\t\t\t",x,y,z);
                }
                x=Frames[j].Rat[Frames[j].Bonds[i][0]][0];
                y=Frames[j].Rat[Frames[j].Bonds[i][0]][1];
                z=Frames[j].Rat[Frames[j].Bonds[i][0]][2];
                fprintf(of,"%6.4f %6.4f %6.4f\n\t\t\t",x,y,z);
                fprintf(of,"\t\t]\n");
                sprintf(tstr,"Bond%0iPos1",i);
                fprintf(of,"\tfield MFVec3f %s [\n",tstr);
                for(j=0;j<NumFrames-1;j++)
                {
                    x=Frames[j].Rat[Frames[j].Bonds[i][1]][0];
                    y=Frames[j].Rat[Frames[j].Bonds[i][1]][1];
                    z=Frames[j].Rat[Frames[j].Bonds[i][1]][2];
                    fprintf(of,"%6.4f %6.4f %6.4f,\n\t\t\t",x,y,z);
                }
                x=Frames[j].Rat[Frames[j].Bonds[i][1]][0];
                y=Frames[j].Rat[Frames[j].Bonds[i][1]][1];
                z=Frames[j].Rat[Frames[j].Bonds[i][1]][2];
                fprintf(of,"%6.4f %6.4f %6.4f\n\t\t\t",x,y,z);
                fprintf(of,"\t\t]\n");
            }
            for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                fprintf(of,"\teventOut SFVec3f %s\n",tstr);
            }
            for(i=0;i<Frames[0].NumBonds;i++)
            {
                sprintf(tstr,"Bond%0iSpine1",i);
                fprintf(of,"\teventOut MFVec3f %s\n",tstr);
                sprintf(tstr,"Bond%0iSpine2",i);
                fprintf(of,"\teventOut MFVec3f %s\n",tstr);
            }
            break;
    default: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPos",i);
                fprintf(of,"\tfield MFVec3f %s (0.0,0.0,0.0)\n",tstr);
            }
            for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                fprintf(of,"\teventOut SFVec3f %s\n",tstr);
            }
            break;
    }
    fprintf(of,"\turl \"javascript:\n");
    fprintf(of,"\t\tfunction initialize (){\n");
    switch(BandS)
    {
    case 0: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                fprintf(of,"\t\t\t%s = new SFVec3f (%6.4f,%6.4f,%6.4f);\n",tstr,Frames[0].Rat[i][0],Frames[0].Rat[i][1],Frames[0].Rat[i][2]);
            }
            break;
    case 1:
    case 2:
    case 3: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                fprintf(of,"\t\t\t%s = new SFVec3f (%6.4f,%6.4f,%6.4f);\n",tstr,Frames[0].Rat[i][0],Frames[0].Rat[i][1],Frames[0].Rat[i][2]);
            }
            for(i=0;i<Frames[0].NumBonds;i++)
            {
                sprintf(tstr,"Bond%0iSpine1",i);
                x=0.5f*(Frames[0].Rat[Frames[0].Bonds[i][0]][0]+Frames[0].Rat[Frames[0].Bonds[i][1]][0]);
                y=0.5f*(Frames[0].Rat[Frames[0].Bonds[i][0]][1]+Frames[0].Rat[Frames[0].Bonds[i][1]][1]);
                z=0.5f*(Frames[0].Rat[Frames[0].Bonds[i][0]][2]+Frames[0].Rat[Frames[0].Bonds[i][1]][2]);
                fprintf(of,"\t\t\t%s = new MFVec3f (new SFVec3f(%6.4f,%6.4f,%6.4f),new SFVec3f(%6.4f,%6.4f,%6.4f));\n",
                        tstr,Frames[0].Rat[Frames[0].Bonds[i][0]][0],Frames[0].Rat[Frames[0].Bonds[i][0]][1],Frames[0].Rat[Frames[0].Bonds[i][0]][2],x,y,z);
                sprintf(tstr,"Bond%0iSpine2",i);
                fprintf(of,"\t\t\t%s = new MFVec3f (new SFVec3f(%6.4f,%6.4f,%6.4f),new SFVec3f(%6.4f,%6.4f,%6.4f));\n",
                        tstr,x,y,z,Frames[0].Rat[Frames[0].Bonds[i][1]][0],Frames[0].Rat[Frames[0].Bonds[i][1]][1],Frames[0].Rat[Frames[0].Bonds[i][1]][2]);
            }
            break;
    default:for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                fprintf(of,"\t\t\t%s = new SFVec3f (%6.4f,%6.4f,%6.4f);\n",tstr,Frames[0].Rat[i][0],Frames[0].Rat[i][1],Frames[0].Rat[i][2]);
            }
    }
    fprintf(of,"\t\t}\n");
    fprintf(of,"\t\tfunction TTime (event,value) {\n");
    switch(BandS)
    {
    case 0: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                sprintf(tstr2,"Atom%0iPos",i);
                fprintf(of,"\t\t\t%s.x = %s[FrameID].x;\n",tstr,tstr2);
                fprintf(of,"\t\t\t%s.y = %s[FrameID].y;\n",tstr,tstr2);
                fprintf(of,"\t\t\t%s.z = %s[FrameID].z;\n",tstr,tstr2);
            }
            break;
    case 1:
    case 2:
    case 3: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                sprintf(tstr2,"Atom%0iPos",i);
                fprintf(of,"\t\t\t%s.x = %s[FrameID].x;\n",tstr,tstr2);
                fprintf(of,"\t\t\t%s.y = %s[FrameID].y;\n",tstr,tstr2);
                fprintf(of,"\t\t\t%s.z = %s[FrameID].z;\n",tstr,tstr2);
            }
            for(i=0;i<Frames[0].NumBonds;i++)
            {
                sprintf(tstr,"Bond%0iSpine1",i);
                sprintf(tstr2,"Bond%0iPos0",i);
                sprintf(tstr3,"Bond%0iPos1",i);
                fprintf(of,"\t\t\t%s[0].x = %s[FrameID].x;\n",tstr,tstr2);
                fprintf(of,"\t\t\t%s[1].x = 0.5*(%s[FrameID].x + %s[FrameID].x);\n",tstr,tstr2,tstr3);
                fprintf(of,"\t\t\t%s[0].y = %s[FrameID].y;\n",tstr,tstr2);
                fprintf(of,"\t\t\t%s[1].y = 0.5*(%s[FrameID].y + %s[FrameID].y);\n",tstr,tstr2,tstr3);
                fprintf(of,"\t\t\t%s[0].z = %s[FrameID].z;\n",tstr,tstr2);
                fprintf(of,"\t\t\t%s[1].z = 0.5*(%s[FrameID].z + %s[FrameID].z);\n",tstr,tstr2,tstr3);
                sprintf(tstr,"Bond%0iSpine2",i);
                fprintf(of,"\t\t\t%s[1].x = %s[FrameID].x;\n",tstr,tstr3);
                fprintf(of,"\t\t\t%s[0].x = 0.5*(%s[FrameID].x + %s[FrameID].x);\n",tstr,tstr2,tstr3);
                fprintf(of,"\t\t\t%s[1].y = %s[FrameID].y;\n",tstr,tstr3);
                fprintf(of,"\t\t\t%s[0].y = 0.5*(%s[FrameID].y + %s[FrameID].y);\n",tstr,tstr2,tstr3);
                fprintf(of,"\t\t\t%s[1].z = %s[FrameID].z;\n",tstr,tstr3);
                fprintf(of,"\t\t\t%s[0].z = 0.5*(%s[FrameID].z + %s[FrameID].z);\n",tstr,tstr2,tstr3);
            }
            break;
    default:for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                sprintf(tstr2,"Atom%0iPos",i);
                fprintf(of,"\t\t\t%s.x = %s[FrameID].x;\n",tstr,tstr2);
                fprintf(of,"\t\t\t%s.y = %s[FrameID].y;\n",tstr,tstr2);
                fprintf(of,"\t\t\t%s.z = %s[FrameID].z;\n",tstr,tstr2);
            }
            break;
    }
    fprintf(of,"\t\t\tFrameID+=1;\n\t\t\tif (FrameID>NumFrames) {\n");
    fprintf(of,"\t\t\t\tFrameID=0;\n\t\t\t}\n");
    fprintf(of,"\t\t}\n\t\"\n}\n\n");
}

int PickNums(int NNum,char* rstr,double* vals)
{
    int i=0;
    int j=0;
    int nv=0;
    char ts[32];

    while((rstr[i]!='\000')&&(nv<NNum))
    {
        ts[j++]=rstr[i++];
        if((ts[j-1]=='\040')||(ts[j-1]=='\000')||(ts[j-1]=='\n')||(ts[j-1]=='\r'))
        {
            ts[j]='\000';
            vals[nv++]=atof(ts);
            j=0;
        }
    }
    if(nv<NNum)
    {
        ts[j]=rstr[i];
        vals[nv++]=atof(ts);
    }
    return nv;
}

void WriteVRMLAtoms(FILE* of,Tcl_Interp *ti)
{
    char tstr[32];
    int i;
    float vdW,scf;
    float r,b,g;
    float rad;
    double colors[3];
    char comstr[80];
    
    for(i=0;i<Frames[0].NumAtoms;i++)
    {
        if(BandS==0)
        {
            vdW=gom_GetAtomVdwRad(0,i);
            scf=gom_GetAtomCPKScale(0,i);
            rad=vdW*scf;
        }
        else
        {
            rad=gom_GetAtomLicoRadS(0);
        }
        sprintf(comstr,"show atom color %0i 1",i+1);
        Tcl_Eval(ti,comstr); 
        PickNums(3,ti->result,colors);
        r=(float)colors[0];g=(float)colors[1];b=(float)colors[2];
        sprintf(tstr,"Atom%0iSph",i);
        fprintf(of,"DEF %s Transform {\n\t",tstr);
        fprintf(of,"children [\n\t\t");
        fprintf(of,"Shape {\n\t\t\tgeometry Sphere {radius %6.4f}\n",rad);
        fprintf(of,"\t\t\t\tappearance Appearance {\n\t\t\t\t");
        fprintf(of,"material Material {diffuseColor %4.3f %4.3f %4.3f}\n\t\t\t}\n\t\t}\n\t]\n}\n",r,g,b);
    }
}

void WriteVRMLSAtoms(FILE* of,Tcl_Interp* ti)
{
    char tstr[32];
    int i;
    float ar,dar;
    float r,b,g;
    char comstr[80];
    double colors[3];
    
    if(BandS<3)
    {
        ar=gom_GetAtomLicoRadC(0);
    }
    else
    {
        ar=0.05f;
    }
    for(i=0;i<Frames[0].NumAtoms;i++)
    {
        if(gom_GetAtomCPKDisplayState(0,i))
        {
            dar=gom_GetAtomVdwRad(0,i)*gom_GetAtomCPKScale(0,i);
        }
        else
        {
            dar=ar;
        }
        sprintf(comstr,"show atom color %0i 1",i+1);
        Tcl_Eval(ti,comstr); 
        PickNums(3,ti->result,colors);
        r=(float)colors[0];g=(float)colors[1];b=(float)colors[2];
        sprintf(tstr,"Atom%0iSph",i);
        fprintf(of,"DEF %s Transform {\n\t",tstr);
        fprintf(of,"children [\n\t\t");
        fprintf(of,"Shape {\n\t\t\tgeometry Sphere {radius %6.4f}\n",dar);
        fprintf(of,"\t\t\t\tappearance Appearance {\n\t\t\t\t");
        fprintf(of,"material Material {diffuseColor %4.3f %4.3f %4.3f}\n\t\t\t}\n\t\t}\n\t]\n}\n",r,g,b);
    }
}


void WriteVRMLBonds(FILE* of,Tcl_Interp* ti)
{
    char tstr[32];
    int i,j;
    float cr;
    float r,b,g;
    int at1,at2;
    double theta;
    float xx,yy;
    char comstr[80];
    double colors[3];
    
    for(i=0;i<Frames[0].NumBonds;i++)
    {
        at1=Frames[0].Bonds[i][0];
        at2=Frames[0].Bonds[i][1];
        if(BandS<3)
        {
            cr=gom_GetAtomLicoRadC(0);
        }
        else
        {
            cr=0.05f;
        }
        sprintf(comstr,"show atom color %0i 1",at1+1);
        Tcl_Eval(ti,comstr); 
        PickNums(3,ti->result,colors);
        r=(float)colors[0];g=(float)colors[1];b=(float)colors[2];
        sprintf(tstr,"Bond%0iCyl0",i);
        fprintf(of,"Shape {\n\tgeometry DEF %s Extrusion {\n",tstr);
        fprintf(of,"\t\tsolid TRUE\n");
        fprintf(of,"\t\tcrossSection [\n");
        for(j=0;j<21;j++)
        {
            theta=pi*(double)j/10.0;
            xx=cr*(float)cos(theta);
            yy=cr*(float)sin(theta);
            fprintf(of,"\t\t\t%6.4f %6.4f\n",xx,yy);
        }
        fprintf(of,"\t\t]\n\t}\n");
        fprintf(of,"\t\tappearance Appearance {\n\t\t");
        fprintf(of,"material Material {diffuseColor %4.3f %4.3f %4.3f}\n\t}\n}\n\n",r,g,b);
        sprintf(comstr,"show atom color %0i 1",at2+1);
        Tcl_Eval(ti,comstr); 
        PickNums(3,ti->result,colors);
        r=(float)colors[0];g=(float)colors[1];b=(float)colors[2];
        sprintf(tstr,"Bond%0iCyl1",i);
        fprintf(of,"Shape {\n\tgeometry DEF %s Extrusion {\n",tstr);
        fprintf(of,"\t\tcrossSection [\n");
        for(j=0;j<21;j++)
        {
            theta=pi*(double)j/10.0;
            xx=cr*(float)cos(theta);
            yy=cr*(float)sin(theta);
            fprintf(of,"\t\t\t%6.4f %6.4f\n",xx,yy);
        }
        fprintf(of,"\t\t]\n\t}\n");
        fprintf(of,"\t\tappearance Appearance {\n\t\t");
        fprintf(of,"material Material {diffuseColor %4.3f %4.3f %4.3f}\n\t}\n}\n\n",r,g,b);
    }
}


void WriteVRMLRoutes(FILE* of)
{
    int i;
    char tstr[32],tstr2[32];

    fprintf(of,"ROUTE RedBall.touchTime TO ToggleScript.TTime\n");
    fprintf(of,"ROUTE RedBall2.touchTime TO ToggleScript.TTime\n");
    fprintf(of,"ROUTE ToggleScript.GoTime TO Timer.startTime\n");
    fprintf(of,"ROUTE ToggleScript.StopTime TO Timer.stopTime\n");
    fprintf(of,"ROUTE Timer.cycleTime TO AnimScript.TTime\n");
    switch(BandS)
    {
    case 0: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                sprintf(tstr2,"Atom%0iSph",i);
                fprintf(of,"ROUTE AnimScript.%s TO %s.translation\n",tstr,tstr2);
            }
            break;
    case 1:
    case 2:
    case 3: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                sprintf(tstr2,"Atom%0iSph",i);
                fprintf(of,"ROUTE AnimScript.%s TO %s.translation\n",tstr,tstr2);
            }
            for(i=0;i<Frames[0].NumBonds;i++)
            {
                sprintf(tstr,"Bond%0iSpine1",i);
                sprintf(tstr2,"Bond%0iCyl0",i);
                fprintf(of,"ROUTE AnimScript.%s TO %s.set_spine\n",tstr,tstr2);
                sprintf(tstr,"Bond%0iSpine2",i);
                sprintf(tstr2,"Bond%0iCyl1",i);
                fprintf(of,"ROUTE AnimScript.%s TO %s.set_spine\n",tstr,tstr2);
            }
            break;
    default:for(i=0;i<Frames[0].NumAtoms;i++)
            {
                sprintf(tstr,"Atom%0iPI",i);
                sprintf(tstr2,"Atom%0iSph",i);
                fprintf(of,"ROUTE AnimScript.%s TO %s.translation\n",tstr,tstr2);
            }
            break;
    }
}

int CloseMovie(ClientData cd,Tcl_Interp *ti,int argc,const char **argv)
{
    FILE* of;
    char filename[80];
    float AnimDelay;
    float MaxRad=-0.1f;
    float dar;
    int i,j;

    if(argc<3)
        return TCL_ERROR;
    AnimDelay=(float)atof(argv[2]);
    sprintf(filename,"%s",argv[1]);
    of=fopen(filename,"w");
    printf("%s %f\n",filename,AnimDelay);
    WriteVRMLHeader(of);
    printf("Header written\n");
    switch(BandS)
    {
    case 0: for(i=0;i<Frames[0].NumAtoms;i++)
            {
                dar=gom_GetAtomVdwRad(0,i)*gom_GetAtomCPKScale(0,i);
                if(dar>MaxRad)
                    MaxRad=dar;
            }
            break;
    case 1: MaxRad=gom_GetAtomLicoRadS(0);
            break;
    case 2: MaxRad=gom_GetAtomLicoRadC(0);
            break;
    case 3: MaxRad=0.05f;
            break;
    default: MaxRad=gom_GetAtomLicoRadS(0);
    }
    WriteVRMLTimerObject(of,AnimDelay,MaxRad);
    printf("Timer written\n");
    WriteVRMLScriptObject(of);
    printf("Script written\n");
    if(BandS<2)
    {
        WriteVRMLAtoms(of,ti);
    }
    else
    {
        WriteVRMLSAtoms(of,ti);
    }
    if(BandS>0)
        WriteVRMLBonds(of,ti);
    WriteVRMLRoutes(of);
    fclose(of);
    for(i=0;i<NumFrames;i++)
    {
        for(j=0;j<Frames[i].NumBonds;j++)
        {
            free((void*)Frames[i].Bonds[j]);
        }
        free((void*)Frames[i].Bonds);
        for(j=0;j<Frames[i].NumAtoms;j++)
        {
            free((void*)Frames[i].Rat[j]);
        }
        free((void*)Frames[i].Rat);
    }
    free((void*)Frames);
    return TCL_OK;
}

} // namespace Vrml
} // namespace Plugin

