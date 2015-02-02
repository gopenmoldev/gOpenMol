/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
/* Provides somewhat improved implementation of Saito and Obara recusrion formulas
(J. Chem. Phys. 1986).*/
#include "gaussians.h"

namespace Plugin {
namespace Fchk {
    
#define pi 3.1415926535897932384626433832795

int nstep=10000;

namespace newrecurse {

int vecindx(int* k)
{
    switch(k[0]+k[1]+k[2])
    {
    case 0: return 0;
            break;
    case 1: if(k[0])
                return 1;
            else if(k[1])
                return 2;
            else
                return 3;
            break;
    case 2: if(k[0])
            {
                if(k[0]==2)
                    return 4;
                else if (k[1])
                    return 7;
                else
                    return 8;
            }
            else if (k[1])
            {
                if(k[1]==2)
                    return 5;
                else
                    return 9;
            }
            else
            {
                return 6;
            }
            break;
    case 3: switch(k[0])
            {
            case 0: switch(k[1])
                    {
                    case 0: return 12;
                            break;
                    case 1: return 18;
                            break;
                    case 2: return 16;
                            break;
                    case 3: return 11;
                            break;
                    }
                    break;
            case 1: switch(k[1])
                    {
                    case 0: return 17;
                            break;
                    case 1: return 19;
                            break;
                    case 2: return 15;
                            break;
                    }
            case 2: if(k[1])
                        return 13;
                    else
                        return 14;
                    break;
            case 3: return 10;
            }
    }
}

double Smat[400];
double Tmat[400];
double Kmat[400];
double *Nmat[6];

double OverlapRecursion(int* k,int* l,double* P,double* A,double* B,double zeta,double xi)
{
    int indx=20*vecindx(k)+vecindx(l);
    double rr;
    int kp[3],kpp[3];
    double rval;
    int i;

    if(Smat[indx]==-9.e99)
    {
        rr=0.0;
        for(i=0;i<3;i++)
        {
            rr+=(A[i]-B[i])*(A[i]-B[i]);
        }
        if(k[0]+k[1]+k[2]>=l[0]+l[1]+l[2])
        {
            if(k[2]>0)
            {
                kp[0]=k[0];kp[1]=k[1];kp[2]=k[2]-1;
                rval=(P[2]-A[2])*OverlapRecursion(kp,l,P,A,B,zeta,xi);
                if(kp[2])
                {
                    kpp[0]=k[0];kpp[1]=k[1];kpp[2]=k[2]-2;
                    rval+=(double)kp[2]*0.5/zeta*OverlapRecursion(kpp,l,P,A,B,zeta,xi);
                }
                if(l[2])
                {
                    kpp[0]=l[0];kpp[1]=l[1];kpp[2]=l[2]-1;
                    rval+=(double)l[2]*0.5/zeta*OverlapRecursion(kp,kpp,P,A,B,zeta,xi);
                }
            } 
            else if (k[1]>0)
            {
                kp[0]=k[0];kp[1]=k[1]-1;kp[2]=0;
                rval=(P[1]-A[1])*OverlapRecursion(kp,l,P,A,B,zeta,xi);
                if(kp[1])
                {
                    kpp[0]=k[0];kpp[1]=k[1]-2;kpp[2]=0;
                    rval+=(double)kp[1]*0.5/zeta*OverlapRecursion(kpp,l,P,A,B,zeta,xi);
                }
                if(l[1])
                {
                    kpp[0]=l[0];kpp[1]=l[1]-1;kpp[2]=l[2];
                    rval+=(double)l[1]*0.5/zeta*OverlapRecursion(kp,kpp,P,A,B,zeta,xi);
                }
            }
            else
            {
                if(k[0]==0)
                {
                    rval=exp(1.5*log(pi/zeta)-xi*rr);
                }
                else
                {
                    kp[0]=k[0]-1;kp[1]=0;kp[2]=0;
                    rval=(P[0]-A[0])*OverlapRecursion(kp,l,P,A,B,zeta,xi);
                    if(kp[0])
                    {
                        kpp[0]=k[0]-2;kpp[1]=0;kpp[2]=0;
                        rval+=(double)kp[0]*0.5/zeta*OverlapRecursion(kpp,l,P,A,B,zeta,xi);
                    }
                    if(l[0])
                    {
                        kpp[0]=l[0]-1;kpp[1]=l[1];kpp[2]=l[2];
                        rval+=(double)l[0]*0.5/zeta*OverlapRecursion(kp,kpp,P,A,B,zeta,xi);
                    }
                }
            }
        }
        else
        {
            if(l[2]>0)
            {
                kp[0]=l[0];kp[1]=l[1];kp[2]=l[2]-1;
                rval=(P[2]-B[2])*OverlapRecursion(k,kp,P,A,B,zeta,xi);
                if(kp[2])
                {
                    kpp[0]=l[0];kpp[1]=l[1];kpp[2]=l[2]-2;
                    rval+=(double)kp[2]*0.5/zeta*OverlapRecursion(k,kpp,P,A,B,zeta,xi);
                }
                if(k[2])
                {
                    kpp[0]=k[0];kpp[1]=k[1];kpp[2]=k[2]-1;
                    rval+=(double)k[2]*0.5/zeta*OverlapRecursion(kpp,kp,P,A,B,zeta,xi);
                }
            } 
            else if (l[1]>0)
            {
                kp[0]=l[0];kp[1]=l[1]-1;kp[2]=0;
                rval=(P[1]-B[1])*OverlapRecursion(k,kp,P,A,B,zeta,xi);
                if(kp[1])
                {
                    kpp[0]=l[0];kpp[1]=l[1]-2;kpp[2]=0;
                    rval+=(double)kp[1]*0.5/zeta*OverlapRecursion(k,kpp,P,A,B,zeta,xi);
                }
                if(k[1])
                {
                    kpp[0]=k[0];kpp[1]=k[1]-1;kpp[2]=k[2];
                    rval+=(double)k[1]*0.5/zeta*OverlapRecursion(kpp,kp,P,A,B,zeta,xi);
                }
            }
            else
            {
                kp[0]=l[0]-1;kp[1]=0;kp[2]=0;
                rval=(P[0]-B[0])*OverlapRecursion(k,kp,P,A,B,zeta,xi);
                if(kp[0])
                {
                    kpp[0]=l[0]-2;kpp[1]=0;kpp[2]=0;
                    rval+=(double)kp[0]*0.5/zeta*OverlapRecursion(k,kpp,P,A,B,zeta,xi);
                }
                if(k[0])
                {
                    kpp[0]=k[0]-1;kpp[1]=k[1];kpp[2]=k[2];
                    rval+=(double)k[0]*0.5/zeta*OverlapRecursion(kpp,kp,P,A,B,zeta,xi);
                }
            }
        }
        Smat[indx]=rval;
    }
    return Smat[indx];
}

double KineticRecursion(int* k,int* l,double* P,double* A,double* B,double zeta,double xi,double za,double zb)
{
    int indx=20*vecindx(k)+vecindx(l);
    double rr;
    int kp[3],kpp[3];
    double rval;
    int i;

    if(Tmat[indx]==-9.e99)
    {
        rr=0.0;
        for(i=0;i<3;i++)
        {
            rr+=(A[i]-B[i])*(A[i]-B[i]);
        }
        if(k[0]+k[1]+k[2]>=l[0]+l[1]+l[2])
        {
            if(k[2]>0)
            {
                kp[0]=k[0];kp[1]=k[1];kp[2]=k[2]-1;
                rval=(P[2]-A[2])*KineticRecursion(kp,l,P,A,B,zeta,xi,za,zb);
                rval+=2.0*xi*OverlapRecursion(k,l,P,A,B,zeta,xi);
                if(kp[2])
                {
                    kpp[0]=k[0];kpp[1]=k[1];kpp[2]=k[2]-2;
                    rval+=(double)kp[2]*0.5/zeta*KineticRecursion(kpp,l,P,A,B,zeta,xi,za,zb);
                    rval-=xi*(double)kp[2]/za*OverlapRecursion(kpp,l,P,A,B,zeta,xi);
                }
                if(l[2])
                {
                    kpp[0]=l[0];kpp[1]=l[1];kpp[2]=l[2]-1;
                    rval+=(double)l[2]*0.5/zeta*KineticRecursion(kp,kpp,P,A,B,zeta,xi,za,zb);
                }
            } 
            else if (k[1]>0)
            {
                kp[0]=k[0];kp[1]=k[1]-1;kp[2]=0;
                rval=(P[1]-A[1])*KineticRecursion(kp,l,P,A,B,zeta,xi,za,zb);
                rval+=2.0*xi*OverlapRecursion(k,l,P,A,B,zeta,xi);
                if(kp[1])
                {
                    kpp[0]=k[0];kpp[1]=k[1]-2;kpp[2]=0;
                    rval+=(double)kp[1]*0.5/zeta*KineticRecursion(kpp,l,P,A,B,zeta,xi,za,zb);
                    rval-=xi*(double)kp[1]/za*OverlapRecursion(kpp,l,P,A,B,zeta,xi);
                }
                if(l[1])
                {
                    kpp[0]=l[0];kpp[1]=l[1]-1;kpp[2]=l[2];
                    rval+=(double)l[1]*0.5/zeta*KineticRecursion(kp,kpp,P,A,B,zeta,xi,za,zb);
                }
            }
            else
            {
                if(k[0]==0)
                {
                    rval=xi*(3.0-2.0*xi*rr)*OverlapRecursion(k,l,P,A,B,zeta,xi);
                }
                else
                {
                    kp[0]=k[0]-1;kp[1]=0;kp[2]=0;
                    rval=(P[0]-A[0])*KineticRecursion(kp,l,P,A,B,zeta,xi,za,zb);
                    rval+=2.0*xi*OverlapRecursion(k,l,P,A,B,zeta,xi);
                    if(kp[0])
                    {
                        kpp[0]=k[0]-2;kpp[1]=0;kpp[2]=0;
                        rval+=(double)kp[0]*0.5/zeta*KineticRecursion(kpp,l,P,A,B,zeta,xi,za,zb);
                        rval-=xi*(double)kp[0]/za*OverlapRecursion(kpp,l,P,A,B,zeta,xi);
                    }
                    if(l[0])
                    {
                        kpp[0]=l[0]-1;kpp[1]=l[1];kpp[2]=l[2];
                        rval+=(double)l[0]*0.5/zeta*KineticRecursion(kp,kpp,P,A,B,zeta,xi,za,zb);
                    }
                }
            }
        }
        else
        {
            if(l[2]>0)
            {
                kp[0]=l[0];kp[1]=l[1];kp[2]=l[2]-1;
                rval=(P[2]-B[2])*KineticRecursion(k,kp,P,A,B,zeta,xi,za,zb);
                rval+=2.0*xi*OverlapRecursion(k,l,P,A,B,zeta,xi);
                if(kp[2])
                {
                    kpp[0]=l[0];kpp[1]=l[1];kpp[2]=l[2]-2;
                    rval+=(double)kp[2]*0.5/zeta*KineticRecursion(k,kpp,P,A,B,zeta,xi,za,zb);
                    rval-=xi*(double)kp[2]/zb*OverlapRecursion(k,kpp,P,A,B,zeta,xi);
                }
                if(k[2])
                {
                    kpp[0]=k[0];kpp[1]=k[1];kpp[2]=k[2]-1;
                    rval+=(double)k[2]*0.5/zeta*KineticRecursion(kpp,kp,P,A,B,zeta,xi,za,zb);
                }
            } 
            else if (l[1]>0)
            {
                kp[0]=l[0];kp[1]=l[1]-1;kp[2]=0;
                rval=(P[1]-B[1])*KineticRecursion(k,kp,P,A,B,zeta,xi,za,zb);
                rval+=2.0*xi*OverlapRecursion(k,l,P,A,B,zeta,xi);
                if(kp[1])
                {
                    kpp[0]=l[0];kpp[1]=l[1]-2;kpp[2]=0;
                    rval+=(double)kp[1]*0.5/zeta*KineticRecursion(k,kpp,P,A,B,zeta,xi,za,zb);
                    rval-=xi*(double)kp[1]/zb*OverlapRecursion(k,kpp,P,A,B,zeta,xi);
                }
                if(k[1])
                {
                    kpp[0]=k[0];kpp[1]=k[1]-1;kpp[2]=k[2];
                    rval+=(double)k[1]*0.5/zeta*KineticRecursion(kpp,kp,P,A,B,zeta,xi,za,zb);
                }
            }
            else
            {
                kp[0]=l[0]-1;kp[1]=0;kp[2]=0;
                rval=(P[0]-B[0])*KineticRecursion(k,kp,P,A,B,zeta,xi,za,zb);
                rval+=2.0*xi*OverlapRecursion(k,l,P,A,B,zeta,xi);
                if(kp[0])
                {
                    kpp[0]=l[0]-2;kpp[1]=0;kpp[2]=0;
                    rval+=(double)kp[0]*0.5/zeta*KineticRecursion(k,kpp,P,A,B,zeta,xi,za,zb);
                    rval-=xi*(double)kp[0]/zb*OverlapRecursion(k,kpp,P,A,B,zeta,xi);
                }
                if(k[0])
                {
                    kpp[0]=k[0]-1;kpp[1]=k[1];kpp[2]=k[2];
                    rval+=(double)k[0]*0.5/zeta*KineticRecursion(kpp,kp,P,A,B,zeta,xi,za,zb);
                }
            }
        }
        Tmat[indx]=rval;
    }
    return Tmat[indx];
}

double Flast[6];
double zlast[6]={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};

double F(int k,double z)
{
    int i;
    double dz=1.0/(double)(nstep);
    double zl;
    double zu;
    double vl,vu;
    double rval=0.0;

    if(z==zlast[k])
        return Flast[k];
    zlast[k]=z;
    if(z==0.0)
    {
        if(k>0)
            return 1.0/(2.0*(double)k+1.0);
        else
            return 1.0;
    }
    for(i=0;i<nstep;i++)
    {
        zl=(double)i*dz;
        zu=zl+dz;
        vl=pow(zl,2*k)*exp(-z*zl*zl);
        vu=pow(zu,2*k)*exp(-z*zu*zu);
        rval+=0.5*(vl+vu);
    }
    Flast[k]=rval*dz;
    return Flast[k];
}

int mmax;

double NuclearRecursion(int* k,int* l,double* P,double* A,double* B,double zeta,double xi,double za,double zb,double* C,int m)
{
    int indx=20*vecindx(k)+vecindx(l);
    double rr;
    double uu;
    int kp[3],kpp[3];
    double rval;
    int i;

    if (m>mmax)
        mmax=m;
    if(Nmat[m][indx]==-9.e99)
    {
        rr=0.0;
        uu=0.0;
        for(i=0;i<3;i++)
        {
            rr+=(A[i]-B[i])*(A[i]-B[i]);
            uu+=(P[i]-C[i])*(P[i]-C[i]);
        }
        uu*=zeta;
        if(k[0]+k[1]+k[2]>=l[0]+l[1]+l[2])
        {
            if(k[2]>0)
            {
                kp[0]=k[0];kp[1]=k[1];kp[2]=k[2]-1;
                rval=(P[2]-A[2])*NuclearRecursion(kp,l,P,A,B,zeta,xi,za,zb,C,m);
                rval-=(P[2]-C[2])*NuclearRecursion(kp,l,P,A,B,zeta,xi,za,zb,C,m+1);
                if(kp[2])
                {
                    kpp[0]=k[0];kpp[1]=k[1];kpp[2]=k[2]-2;
                    rval+=(double)kp[2]*0.5/zeta*(NuclearRecursion(kpp,l,P,A,B,zeta,xi,za,zb,C,m)
                                -NuclearRecursion(kpp,l,P,A,B,zeta,xi,za,zb,C,m+1));
                }
                if(l[2])
                {
                    kpp[0]=l[0];kpp[1]=l[1];kpp[2]=l[2]-1;
                    rval+=(double)l[2]*0.5/zeta*(NuclearRecursion(kp,kpp,P,A,B,zeta,xi,za,zb,C,m)
                                -NuclearRecursion(kp,kpp,P,A,B,zeta,xi,za,zb,C,m+1));
                }
            } 
            else if (k[1]>0)
            {
                kp[0]=k[0];kp[1]=k[1]-1;kp[2]=0;
                rval=(P[1]-A[1])*NuclearRecursion(kp,l,P,A,B,zeta,xi,za,zb,C,m);
                rval-=(P[1]-C[1])*NuclearRecursion(kp,l,P,A,B,zeta,xi,za,zb,C,m+1);
                if(kp[1])
                {
                    kpp[0]=k[0];kpp[1]=k[1]-2;kpp[2]=0;
                    rval+=(double)kp[1]*0.5/zeta*(NuclearRecursion(kpp,l,P,A,B,zeta,xi,za,zb,C,m)
                                -NuclearRecursion(kpp,l,P,A,B,zeta,xi,za,zb,C,m+1));
                }
                if(l[1])
                {
                    kpp[0]=l[0];kpp[1]=l[1]-1;kpp[2]=l[2];
                    rval+=(double)l[1]*0.5/zeta*(NuclearRecursion(kp,kpp,P,A,B,zeta,xi,za,zb,C,m)
                                    -NuclearRecursion(kp,kpp,P,A,B,zeta,xi,za,zb,C,m+1));
                }
            }
            else
            {
                if(k[0]==0)
                {
                    rval=2.0*sqrt(zeta/pi)*F(m,uu)*exp(1.5*log(pi/zeta)-xi*rr);
                }
                else
                {
                    kp[0]=k[0]-1;kp[1]=0;kp[2]=0;
                    rval=(P[0]-A[0])*NuclearRecursion(kp,l,P,A,B,zeta,xi,za,zb,C,m);
                    rval-=(P[0]-C[0])*NuclearRecursion(kp,l,P,A,B,zeta,xi,za,zb,C,m+1);
                    if(kp[0])
                    {
                        kpp[0]=k[0]-2;kpp[1]=0;kpp[2]=0;
                        rval+=(double)kp[0]*0.5/zeta*(NuclearRecursion(kpp,l,P,A,B,zeta,xi,za,zb,C,m)
                                    -NuclearRecursion(kpp,l,P,A,B,zeta,xi,za,zb,C,m+1));
                    }
                    if(l[0])
                    {
                        kpp[0]=l[0]-1;kpp[1]=l[1];kpp[2]=l[2];
                        rval+=(double)l[0]*0.5/zeta*(NuclearRecursion(kp,kpp,P,A,B,zeta,xi,za,zb,C,m)
                                    -NuclearRecursion(kp,kpp,P,A,B,zeta,xi,za,zb,C,m+1));
                    }
                }
            }
        }
        else
        {
            if(l[2]>0)
            {
                kp[0]=l[0];kp[1]=l[1];kp[2]=l[2]-1;
                rval=(P[2]-B[2])*NuclearRecursion(k,kp,P,A,B,zeta,xi,za,zb,C,m);
                rval-=(P[2]-C[2])*NuclearRecursion(k,kp,P,A,B,zeta,xi,za,zb,C,m+1);
                if(kp[2])
                {
                    kpp[0]=l[0];kpp[1]=l[1];kpp[2]=l[2]-2;
                    rval+=(double)kp[2]*0.5/zeta*(NuclearRecursion(k,kpp,P,A,B,zeta,xi,za,zb,C,m)
                                -NuclearRecursion(k,kpp,P,A,B,zeta,xi,za,zb,C,m+1));
                }
                if(k[2])
                {
                    kpp[0]=k[0];kpp[1]=k[1];kpp[2]=k[2]-1;
                    rval+=(double)k[2]*0.5/zeta*(NuclearRecursion(kpp,kp,P,A,B,zeta,xi,za,zb,C,m)
                                -NuclearRecursion(kpp,kp,P,A,B,zeta,xi,za,zb,C,m+1));
                }
            } 
            else if (l[1]>0)
            {
                kp[0]=l[0];kp[1]=l[1]-1;kp[2]=0;
                rval=(P[1]-B[1])*NuclearRecursion(k,kp,P,A,B,zeta,xi,za,zb,C,m);
                rval-=(P[1]-C[1])*NuclearRecursion(k,kp,P,A,B,zeta,xi,za,zb,C,m+1);
                if(kp[1])
                {
                    kpp[0]=l[0];kpp[1]=l[1]-2;kpp[2]=0;
                    rval+=(double)kp[1]*0.5/zeta*(NuclearRecursion(k,kpp,P,A,B,zeta,xi,za,zb,C,m)
                                -NuclearRecursion(k,kpp,P,A,B,zeta,xi,za,zb,C,m+1));
                }
                if(k[1])
                {
                    kpp[0]=k[0];kpp[1]=k[1]-1;kpp[2]=k[2];
                    rval+=(double)k[1]*0.5/zeta*(NuclearRecursion(kpp,kp,P,A,B,zeta,xi,za,zb,C,m)
                                -NuclearRecursion(kpp,kp,P,A,B,zeta,xi,za,zb,C,m+1));
                }
            }
            else
            {
                kp[0]=l[0]-1;kp[1]=0;kp[2]=0;
                rval=(P[0]-B[0])*NuclearRecursion(k,kp,P,A,B,zeta,xi,za,zb,C,m);
                rval-=(P[0]-C[0])*NuclearRecursion(k,kp,P,A,B,zeta,xi,za,zb,C,m+1);
                if(kp[0])
                {
                    kpp[0]=l[0]-2;kpp[1]=0;kpp[2]=0;
                    rval+=(double)kp[0]*0.5/zeta*(NuclearRecursion(k,kpp,P,A,B,zeta,xi,za,zb,C,m)
                                -NuclearRecursion(k,kpp,P,A,B,zeta,xi,za,zb,C,m+1));
                }
                if(k[0])
                {
                    kpp[0]=k[0]-1;kpp[1]=k[1];kpp[2]=k[2];
                    rval+=(double)k[0]*0.5/zeta*(NuclearRecursion(kpp,kp,P,A,B,zeta,xi,za,zb,C,m)
                                -NuclearRecursion(kpp,kp,P,A,B,zeta,xi,za,zb,C,m+1));
                }
            }
        }
        Nmat[m][indx]=rval;
    }
    return Nmat[m][indx];
}

int FillAngMom(int** k,int angmom)
{
    int nk,i;

    switch(angmom)
    {
    case 0: nk=1;
            k[0]=(int*)malloc(3*sizeof(int));
            k[0][0]=k[0][1]=k[0][2]=0;
            break;
    case 1: nk=3;
            for(i=0;i<3;i++)
            {
                k[i]=(int*)malloc(3*sizeof(int));
            }
            k[0][0]=k[1][1]=k[2][2]=1;
            k[0][1]=k[0][2]=k[1][0]=k[1][2]=k[2][0]=k[2][1]=0;
            break;
    case 2: nk=6;
            for(i=0;i<6;i++)
            {
                k[i]=(int*)malloc(3*sizeof(int));
            }
            k[0][0]=2;k[0][1]=0;k[0][2]=0;
            k[1][0]=0;k[1][1]=2;k[1][2]=0;
            k[2][0]=0;k[2][1]=0;k[2][2]=2;
            k[3][0]=1;k[3][1]=1;k[3][2]=0;
            k[4][0]=1;k[4][1]=0;k[4][2]=1;
            k[5][0]=0;k[5][1]=1;k[5][2]=1;
            break;
    case 3: nk=10;
            for(i=0;i<10;i++)
            {
                k[i]=(int*)malloc(3*sizeof(int));
            }
            k[0][0]=3;k[0][1]=0;k[0][2]=0;
            k[1][0]=0;k[1][1]=3;k[1][2]=0;
            k[2][0]=0;k[2][1]=0;k[2][2]=3;
            k[3][0]=2;k[3][1]=1;k[3][2]=0;
            k[4][0]=2;k[4][1]=0;k[4][2]=1;
            k[5][0]=1;k[5][1]=2;k[5][2]=0;
            k[6][0]=0;k[6][1]=2;k[6][2]=1;
            k[7][0]=1;k[7][1]=0;k[7][2]=2;
            k[8][0]=0;k[8][1]=1;k[8][2]=2;
            k[9][0]=1;k[9][1]=1;k[9][2]=1;
            break;
    }
    return nk;
}

void FillResult(double* rspace,double* workspace,int AM1,int AM2,int nk,int nl)
{
    int ik,il;

    if(AM1>=0)
    {
        for(ik=0;ik<nk;ik++)
        {
            if(AM2>=0)
            {
                for(il=0;il<nl;il++)
                {
                    rspace[ik*nl+il]=workspace[ik*nl+il];
                }
            }
            else
            {
                if(AM2=-2)
                {
                    rspace[ik*5]=2.0*workspace[ik*nl+2]-workspace[ik*nl]-workspace[ik*nl+1];
                    rspace[ik*5+1]=workspace[ik*nl+4];
                    rspace[ik*5+2]=workspace[ik*nl+5];
                    rspace[ik*5+3]=workspace[ik*nl+3];
                    rspace[ik*5+4]=workspace[ik*nl]-workspace[ik*nl+1];
                }
                else
                {
                    rspace[ik*7]=2.0*workspace[ik*nl+2]-3.0*workspace[ik*nl+4]-3.0*workspace[ik*nl+6];
                    rspace[ik*7+1]=4.0*workspace[ik*nl+7]-workspace[ik*nl+5]-workspace[ik*nl];
                    rspace[ik*7+2]=4.0*workspace[ik*nl+8]-workspace[ik*nl+3]-workspace[ik*nl+1];
                    rspace[ik*7+3]=workspace[ik*nl+4]-workspace[ik*nl+6];
                    rspace[ik*7+4]=workspace[ik*nl+9];
                    rspace[ik*7+5]=workspace[ik*nl]-3.0*workspace[ik*nl+5];
                    rspace[ik*7+6]=3.0*workspace[ik*nl+3]-workspace[ik*nl+1];
                }
            }
        }
    }
    else
    {
        if(AM1==-2)
        {
            if(AM2>=0)
            {
                for(il=0;il<nl;il++)
                {
                    rspace[il]=2.0*workspace[2*nl+il]-workspace[il]-workspace[nl+il];
                    rspace[nl+il]=workspace[4*nl+il];
                    rspace[2*nl+il]=workspace[5*nl+il];
                    rspace[3*nl+il]=workspace[3*nl+il];
                    rspace[4*nl+il]=workspace[il]-workspace[nl+il];
                }
            }
            else
            {
                if(AM2==-2)
                {
                    rspace[0]=4.0*workspace[14]-2.0*(workspace[12]+workspace[13]+workspace[2]+workspace[8])
                        +workspace[1]+workspace[0]+workspace[6]+workspace[7];
                    rspace[1]=2.0*workspace[16]-workspace[4]-workspace[10];
                    rspace[2]=2.0*workspace[17]-workspace[5]-workspace[11];
                    rspace[3]=2.0*workspace[15]-workspace[3]-workspace[9];
                    rspace[4]=2.0*(workspace[12]-workspace[13])-
                                workspace[0]+workspace[1]-workspace[6]+workspace[7];
                    rspace[5]=2.0*workspace[26]-workspace[24]-workspace[25];
                    rspace[6]=workspace[28];
                    rspace[7]=workspace[29];
                    rspace[8]=workspace[27];
                    rspace[9]=workspace[24]-workspace[25];
                    rspace[10]=2.0*workspace[32]-workspace[30]-workspace[31];
                    rspace[11]=workspace[34];
                    rspace[12]=workspace[35];
                    rspace[13]=workspace[33];
                    rspace[14]=workspace[30]-workspace[31];
                    rspace[15]=2.0*workspace[20]-workspace[18]-workspace[19];
                    rspace[16]=workspace[22];
                    rspace[17]=workspace[23];
                    rspace[18]=workspace[21];
                    rspace[19]=workspace[18]-workspace[17];
                    rspace[20]=2.0*(workspace[2]-workspace[8])
                                -workspace[0]+workspace[1]-workspace[6]+workspace[7];
                    rspace[21]=workspace[4]-workspace[10];
                    rspace[22]=workspace[5]-workspace[11];
                    rspace[23]=workspace[3]-workspace[9];
                    rspace[24]=workspace[0]-workspace[1]-workspace[6]+workspace[7];
                }
                else
                {
                    rspace[0]=4.0*workspace[22]-6.0*workspace[24]-6.0*workspace[26]
                        -2.0*workspace[2]+3.0*workspace[4]+3.0*workspace[6]
                        -2.0*workspace[12]+3.0*workspace[14]+3.0*workspace[16];
                    rspace[1]=2.0*workspace[42]-3.0*workspace[44]-3.0*workspace[46];
                    rspace[2]=2.0*workspace[52]-3.0*workspace[54]-3.0*workspace[56];
                    rspace[3]=2.0*workspace[32]-3.0*workspace[34]-3.0*workspace[36];
                    rspace[4]=2.0*(workspace[2]-workspace[12])
                        -3.0*(workspace[4]-workspace[14])
                        -3.0*(workspace[6]-workspace[16]);
                    rspace[5]=8.0*workspace[27]-2.0*workspace[20]-2.0*workspace[25]
                        -4.0*workspace[7]+workspace[0]+workspace[5]
                        -4.0*workspace[17]+workspace[10]+workspace[15];
                    rspace[6]=4.0*workspace[47]-workspace[40]-workspace[45];
                    rspace[7]=4.0*workspace[57]-workspace[50]-workspace[55];
                    rspace[8]=4.0*workspace[37]-workspace[30]-workspace[35];
                    rspace[9]=4.0*(workspace[7]-workspace[17])
                        -(workspace[0]-workspace[10])
                        -(workspace[5]-workspace[15]);
                    rspace[10]=8.0*workspace[28]-2.0*workspace[21]-2.0*workspace[23]
                        -4.0*workspace[8]+workspace[1]+workspace[3]
                        -4.0*workspace[18]+workspace[11]+workspace[13];
                    rspace[11]=4.0*workspace[48]-workspace[41]-workspace[43];
                    rspace[12]=4.0*workspace[58]-workspace[51]-workspace[53];
                    rspace[13]=4.0*workspace[38]-workspace[31]-workspace[33];
                    rspace[14]=4.0*(workspace[8]-workspace[18])
                        -(workspace[1]-workspace[11])
                        -(workspace[3]-workspace[13]);
                    rspace[15]=2.0*(workspace[24]-workspace[26])
                        -(workspace[4]-workspace[6])
                        -(workspace[14]-workspace[16]);
                    rspace[16]=workspace[44]-workspace[46];
                    rspace[17]=workspace[54]-workspace[56];
                    rspace[18]=workspace[34]-workspace[36];
                    rspace[19]=(workspace[4]-workspace[6])
                        -(workspace[14]-workspace[16]);
                    rspace[20]=2.0*workspace[29]-workspace[9]-workspace[19];
                    rspace[21]=workspace[49];
                    rspace[22]=workspace[59];
                    rspace[23]=workspace[39];
                    rspace[24]=workspace[9]-workspace[19];
                    rspace[25]=2.0*(workspace[20]-3.0*workspace[25])
                        -(workspace[0]-3.0*workspace[5])
                        -(workspace[10]-3.0*workspace[15]);
                    rspace[26]=workspace[40]-3.0*workspace[45];
                    rspace[27]=workspace[50]-3.0*workspace[55];
                    rspace[28]=workspace[30]-3.0*workspace[35];
                    rspace[29]=(workspace[0]-workspace[10])
                        -3.0*(workspace[5]-workspace[15]);
                    rspace[30]=3.0*(2.0*workspace[23]-workspace[3]-workspace[13])
                        -(2.0*workspace[21]-workspace[1]-workspace[11]);
                    rspace[31]=3.0*workspace[43]-workspace[41];
                    rspace[32]=3.0*workspace[53]-workspace[51];
                    rspace[33]=3.0*workspace[33]-workspace[31];
                    rspace[34]=3.0*(workspace[3]-workspace[1])
                        -(workspace[13]-workspace[11]);
                }
            }
        }
        else
        {
            if(AM2>=0)
            {
                for(il=0;il<nl;il++)
                {
                    rspace[nl*7]=2.0*workspace[2*nl+il]-3.0*workspace[4*nl+il]-3.0*workspace[6*nl+il];
                    rspace[nl*7+1]=4.0*workspace[7*nl+il]-workspace[5*nl+il]-workspace[il];
                    rspace[nl*7+2]=4.0*workspace[8*nl+il]-workspace[3*nl+il]-workspace[1*nl+il];
                    rspace[nl*7+3]=workspace[4*nl+il]-workspace[6*nl+il];
                    rspace[nl*7+4]=workspace[9*nl+il];
                    rspace[nl*7+5]=workspace[il]-3.0*workspace[5*nl+il];
                    rspace[nl*7+6]=3.0*workspace[3*nl+il]-workspace[1*nl+il];
                }
            }
            else
            {
                if(AM2==-2)
                {
                    rspace[0]=2.0*(2.0*workspace[14]-workspace[12]-workspace[13])
                        -3.0*(2.0*workspace[26]-workspace[24]-workspace[25])
                        -3.0*(2.0*workspace[38]-workspace[36]-workspace[37]);
                    rspace[1]=2.0*workspace[16]-3.0*workspace[28]-3.0*workspace[40];
                    rspace[2]=2.0*workspace[17]-3.0*workspace[29]-3.0*workspace[41];
                    rspace[3]=2.0*workspace[15]-3.0*workspace[27]-3.0*workspace[39];
                    rspace[4]=2.0*(workspace[12]-workspace[13])
                        -3.0*(workspace[24]-workspace[25])
                        -3.0*(workspace[36]-workspace[37]);
                    rspace[5]=4.0*(2.0*workspace[44]-workspace[42]-workspace[43])
                        -(2.0*workspace[2]-workspace[0]-workspace[1])
                        -(2.0*workspace[32]-workspace[30]-workspace[31]);
                    rspace[6]=4.0*workspace[46]-workspace[4]-workspace[34];
                    rspace[7]=4.0*workspace[47]-workspace[5]-workspace[35];
                    rspace[8]=4.0*workspace[45]-workspace[3]-workspace[33];
                    rspace[9]=4.0*(workspace[42]-workspace[43])
                        -(workspace[0]-workspace[1])
                        -(workspace[30]-workspace[31]);
                    rspace[10]=4.0*(2.0*workspace[50]-workspace[48]-workspace[49])
                        -(2.0*workspace[8]-workspace[6]-workspace[7])
                        -(2.0*workspace[20]-workspace[18]-workspace[19]);
                    rspace[11]=4.0*workspace[52]-workspace[10]-workspace[22];
                    rspace[12]=4.0*workspace[53]-workspace[11]-workspace[23];
                    rspace[13]=4.0*workspace[51]-workspace[9]-workspace[21];
                    rspace[14]=4.0*(workspace[48]-workspace[49])
                        -(workspace[6]-workspace[7])
                        -(workspace[18]-workspace[19]);
                    rspace[15]=(2.0*workspace[26]-workspace[24]-workspace[25])
                        -(2.0*workspace[38]-workspace[36]-workspace[37]);
                    rspace[16]=workspace[28]-workspace[40];
                    rspace[17]=workspace[29]-workspace[41];
                    rspace[18]=workspace[27]-workspace[38];
                    rspace[19]=(workspace[24]-workspace[25])
                        -(workspace[36]-workspace[37]);
                    rspace[20]=2.0*workspace[56]-workspace[54]-workspace[55];
                    rspace[21]=workspace[58];
                    rspace[22]=workspace[59];
                    rspace[23]=workspace[57];
                    rspace[24]=workspace[54]-workspace[55];
                    rspace[25]=(2.0*workspace[2]-workspace[0]-workspace[1])
                        -3.0*(2.0*workspace[32]-workspace[30]-workspace[31]);
                    rspace[26]=workspace[4]-3.0*workspace[34];
                    rspace[27]=workspace[5]-3.0*workspace[35];
                    rspace[28]=workspace[3]-3.0*workspace[33];
                    rspace[29]=(workspace[0]-workspace[1])
                        -3.0*(workspace[30]-workspace[31]);
                    rspace[30]=3.0*(2.0*workspace[20]-workspace[18]-workspace[19])
                        -(2.0*workspace[8]-workspace[6]-workspace[7]);
                    rspace[31]=3.0*workspace[22]-workspace[10];
                    rspace[32]=3.0*workspace[23]-workspace[11];
                    rspace[33]=3.0*workspace[21]-workspace[9];
                    rspace[34]=3.0*(workspace[18]-workspace[19])
                        -(workspace[6]-workspace[7]);
                }
                else
                {
                    rspace[0]=2.0*(2.0*workspace[22]-3.0*workspace[24]-3.0*workspace[26])
                        -3.0*(2.0*workspace[42]-3.0*workspace[44]-3.0*workspace[46])
                        -3.0*(2.0*workspace[62]-3.0*workspace[64]-3.0*workspace[66]);
                    rspace[1]=2.0*(4.0*workspace[27]-workspace[20]-workspace[25])
                        -3.0*(4.0*workspace[47]-workspace[40]-workspace[45])
                        -3.0*(4.0*workspace[67]-workspace[60]-workspace[65]);
                    rspace[2]=2.0*(4.0*workspace[28]-workspace[21]-workspace[23])
                        -3.0*(4.0*workspace[48]-workspace[41]-workspace[43])
                        -3.0*(4.0*workspace[68]-workspace[61]-workspace[63]);
                    rspace[3]=2.0*(workspace[24]-workspace[26])
                        -3.0*(workspace[44]-workspace[46])
                        -3.0*(workspace[64]-workspace[66]);
                    rspace[4]=2.0*workspace[29]-3.0*workspace[49]-3.0*workspace[69];
                    rspace[5]=2.0*(workspace[20]-3.0*workspace[25])
                        -3.0*(workspace[40]-3.0*workspace[45])
                        -3.0*(workspace[60]-3.0*workspace[65]);
                    rspace[6]=2.0*(3.0*workspace[23]-workspace[21])
                        -3.0*(3.0*workspace[43]-workspace[41])
                        -3.0*(3.0*workspace[63]-workspace[61]);
                    rspace[7]=4.0*(2.0*workspace[72]-3.0*workspace[74]-3.0*workspace[76])
                        -(2.0*workspace[2]-3.0*workspace[4]-3.0*workspace[6])
                        -(2.0*workspace[52]-3.0*workspace[54]-3.0*workspace[56]);
                    rspace[8]=4.0*(4.0*workspace[77]-workspace[70]-workspace[75])
                        -(4.0*workspace[7]-workspace[0]-workspace[5])
                        -(4.0*workspace[57]-workspace[50]-workspace[55]);
                    rspace[9]=4.0*(4.0*workspace[78]-workspace[71]-workspace[73])
                        -(4.0*workspace[8]-workspace[1]-workspace[3])
                        -(4.0*workspace[58]-workspace[51]-workspace[53]);
                    rspace[10]=4.0*(workspace[74]-workspace[76])
                        -(workspace[4]-workspace[6])
                        -(workspace[54]-workspace[56]);
                    rspace[11]=4.0*workspace[79]-workspace[9]-workspace[59];
                    rspace[12]=4.0*(workspace[70]-3.0*workspace[75])
                        -(workspace[0]-3.0*workspace[5])
                        -(workspace[50]-3.0*workspace[55]);
                    rspace[13]=4.0*(3.0*workspace[73]-workspace[71])
                        -(3.0*workspace[3]-workspace[1])
                        -(3.0*workspace[53]-workspace[51]);
                    rspace[14]=4.0*(2.0*workspace[82]-3.0*workspace[84]-3.0*workspace[86])
                        -(2.0*workspace[12]-3.0*workspace[14]-3.0*workspace[16])
                        -(2.0*workspace[32]-3.0*workspace[34]-3.0*workspace[36]);
                    rspace[15]=4.0*(4.0*workspace[87]-workspace[80]-workspace[85])
                        -(4.0*workspace[17]-workspace[10]-workspace[15])
                        -(4.0*workspace[37]-workspace[30]-workspace[35]);
                    rspace[16]=4.0*(4.0*workspace[88]-workspace[81]-workspace[83])
                        -(4.0*workspace[18]-workspace[11]-workspace[13])
                        -(4.0*workspace[38]-workspace[31]-workspace[33]);
                    rspace[17]=4.0*(workspace[84]-workspace[86])
                        -(workspace[14]-workspace[16])
                        -(workspace[34]-workspace[36]);
                    rspace[18]=4.0*workspace[89]-workspace[19]-workspace[39];
                    rspace[19]=4.0*(workspace[80]-3.0*workspace[85])
                        -(workspace[10]-3.0*workspace[15])
                        -(workspace[30]-3.0*workspace[35]);
                    rspace[20]=4.0*(3.0*workspace[83]-workspace[81])
                        -(3.0*workspace[13]-workspace[11])
                        -(3.0*workspace[33]-workspace[31]);
                    rspace[21]=(2.0*workspace[42]-3.0*workspace[44]-3.0*workspace[46])
                        -(2.0*workspace[62]-3.0*workspace[64]-3.0*workspace[66]);
                    rspace[22]=(4.0*workspace[47]-workspace[40]-workspace[45])
                        -(4.0*workspace[67]-workspace[60]-workspace[65]);
                    rspace[23]=(4.0*workspace[48]-workspace[41]-workspace[43])
                        -(4.0*workspace[68]-workspace[61]-workspace[63]);
                    rspace[24]=(workspace[44]-workspace[46])
                        -(workspace[64]-workspace[66]);
                    rspace[25]=workspace[49]-workspace[69];
                    rspace[26]=(workspace[40]-3.0*workspace[45])
                        -(workspace[60]-3.0*workspace[65]);
                    rspace[27]=(3.0*workspace[43]-workspace[41])
                        -(3.0*workspace[63]-workspace[61]);
                    rspace[28]=2.0*workspace[92]-3.0*workspace[94]-3.0*workspace[96];
                    rspace[29]=4.0*workspace[97]-workspace[90]-workspace[95];
                    rspace[30]=4.0*workspace[98]-workspace[91]-workspace[93];
                    rspace[31]=workspace[94]-workspace[96];
                    rspace[32]=workspace[99];
                    rspace[33]=workspace[90]-3.0*workspace[95];
                    rspace[34]=3.0*workspace[93]-workspace[91];
                    rspace[35]=(2.0*workspace[2]-3.0*workspace[4]-3.0*workspace[6])
                        -3.0*(2.0*workspace[52]-3.0*workspace[54]-3.0*workspace[56]);
                    rspace[36]=(4.0*workspace[7]-workspace[0]-workspace[5])
                        -3.0*(4.0*workspace[57]-workspace[50]-workspace[55]);
                    rspace[37]=(4.0*workspace[8]-workspace[1]-workspace[3])
                        -3.0*(4.0*workspace[58]-workspace[51]-workspace[53]);
                    rspace[38]=(workspace[4]-workspace[6])
                        -3.0*(workspace[54]-workspace[56]);
                    rspace[39]=workspace[9]-3.0*workspace[59];
                    rspace[40]=(workspace[0]-3.0*workspace[5])
                        -3.0*(workspace[50]-3.0*workspace[55]);
                    rspace[41]=(3.0*workspace[3]-workspace[1])
                        -3.0*(3.0*workspace[53]-workspace[51]);
                    rspace[42]=3.0*(2.0*workspace[32]-3.0*workspace[34]-3.0*workspace[36])
                        -(2.0*workspace[12]-3.0*workspace[14]-3.0*workspace[16]);
                    rspace[43]=3.0*(4.0*workspace[37]-workspace[30]-workspace[35])
                        -(4.0*workspace[17]-workspace[10]-workspace[15]);
                    rspace[44]=3.0*(4.0*workspace[38]-workspace[31]-workspace[33])
                        -(4.0*workspace[18]-workspace[11]-workspace[13]);
                    rspace[45]=3.0*(rspace[34]-rspace[36])
                        -(rspace[14]-rspace[16]);
                    rspace[46]=3.0*workspace[39]-workspace[19];
                    rspace[47]=3.0*(workspace[30]-3.0*workspace[35])
                        -(workspace[10]-3.0*workspace[15]);
                    rspace[48]=3.0*(3.0*workspace[33]-workspace[31])
                        -(3.0*workspace[13]-workspace[11]);
                }
            }
        }
    }
}

int OverlapIntegral(GaussianOrbital *g1,GaussianOrbital *g2,double *rspace)
{
    double zeta,xi,P[3];
    int AngMom[2];
    int *k[10];
    int *l[10];
    int i;
    int nk,nl;
    double *workspace;
    double* A;
    double* B;
    int ik,il;

    memcpy(Smat,Kmat,400*sizeof(double));
    zeta=g1->OrbExponent+g2->OrbExponent;
    xi=g1->OrbExponent*g2->OrbExponent/zeta;
    A=g1->Center;
    B=g2->Center;
    for(i=0;i<3;i++)
    {
        P[i]=(g1->OrbExponent*A[i]+g2->OrbExponent*B[i])/zeta;
    }
    AngMom[0]=abs(g1->AngMom);
    AngMom[1]=abs(g2->AngMom);
    nk=FillAngMom(k,AngMom[0]);
    nl=FillAngMom(l,AngMom[1]);
    workspace=(double*)malloc(nk*nl*sizeof(double));
    for(ik=0;ik<nk;ik++)
    {
        for(il=0;il<nl;il++)
        {
            workspace[ik*nl+il]=OverlapRecursion(k[ik],l[il],P,A,B,zeta,xi);
        }
    }
    FillResult(rspace,workspace,g1->AngMom,g2->AngMom,nk,nl);
    free((void*)workspace);
    for(i=0;i<nk;i++)
        free((void*)k[i]);
    for(i=0;i<nl;i++)
        free((void*)l[i]);
    return 0;
}

int KineticIntegral(GaussianOrbital *g1,GaussianOrbital *g2,double *rspace)
{
    double zeta,xi,P[3];
    int AngMom[2];
    int *k[10];
    int *l[10];
    int i;
    int nk,nl;
    double *workspace;
    double* A;
    double* B;
    int ik,il;
    double za,zb;

    memcpy(Smat,Kmat,400*sizeof(double));
    memcpy(Tmat,Kmat,400*sizeof(double));
    zeta=g1->OrbExponent+g2->OrbExponent;
    xi=g1->OrbExponent*g2->OrbExponent/zeta;
    A=g1->Center;
    B=g2->Center;
    za=g1->OrbExponent;
    zb=g2->OrbExponent;
    for(i=0;i<3;i++)
    {
        P[i]=(g1->OrbExponent*A[i]+g2->OrbExponent*B[i])/zeta;
    }
    AngMom[0]=abs(g1->AngMom);
    AngMom[1]=abs(g2->AngMom);
    nk=FillAngMom(k,AngMom[0]);
    nl=FillAngMom(l,AngMom[1]);
    workspace=(double*)malloc(nk*nl*sizeof(double));
    for(ik=0;ik<nk;ik++)
    {
        for(il=0;il<nl;il++)
        {
            workspace[ik*nl+il]=KineticRecursion(k[ik],l[il],P,A,B,zeta,xi,za,zb);
        }
    }
    FillResult(rspace,workspace,g1->AngMom,g2->AngMom,nk,nl);
    free((void*)workspace);
    for(i=0;i<nk;i++)
        free((void*)k[i]);
    for(i=0;i<nl;i++)
        free((void*)l[i]);
    return 0;
}

int NuclearIntegral(GaussianOrbital *g1,GaussianOrbital *g2,AtomicCore *AC,double *rspace)
{
    double zeta,xi,P[3];
    int AngMom[2];
    int *k[10];
    int *l[10];
    int i;
    int nk,nl;
    double *workspace;
    double* A;
    double* B;
    double* C;
    int ik,il;
    double za,zb;
    double Q;

    memcpy(Smat,Kmat,400*sizeof(double));
    for(i=0;i<6;i++)
    {
        memcpy(Nmat[i],Kmat,400*sizeof(double));
        zlast[i]=-1.0;
    }
    zeta=g1->OrbExponent+g2->OrbExponent;
    xi=g1->OrbExponent*g2->OrbExponent/zeta;
    A=g1->Center;
    B=g2->Center;
    C=AC->Position;
    Q=(double)AC->NuclCharge;
    za=g1->OrbExponent;
    zb=g2->OrbExponent;
    for(i=0;i<3;i++)
    {
        P[i]=(g1->OrbExponent*A[i]+g2->OrbExponent*B[i])/zeta;
    }
    AngMom[0]=abs(g1->AngMom);
    AngMom[1]=abs(g2->AngMom);
    nk=FillAngMom(k,AngMom[0]);
    nl=FillAngMom(l,AngMom[1]);
    workspace=(double*)malloc(nk*nl*sizeof(double));
    for(ik=0;ik<nk;ik++)
    {
        for(il=0;il<nl;il++)
        {
            workspace[ik*nl+il]=NuclearRecursion(k[ik],l[il],P,A,B,zeta,xi,za,zb,C,0);
        }
    }
    FillResult(rspace,workspace,g1->AngMom,g2->AngMom,nk,nl);
    if(g1->AngMom>=0)
    {
        for(ik=0;ik<nk;ik++)
        {
            if(g2->AngMom>=0)
            {
                for(il=0;il<nl;il++)
                {
                    rspace[ik*nl+il]*=Q;
                }
            }
            else
            {
                if(g2->AngMom=-2)
                {
                    for(il=0;il<5;il++)
                        rspace[ik*5+il]*=Q;
                }
                else
                {
                    for(il=0;il<7;il++)
                        rspace[ik*7+il]*=Q;
                }
            }
        }
    }
    else
    {
        if(g1->AngMom==-2)
        {
            if(g2->AngMom>=0)
            {
                for(il=0;il<nl;il++)
                {
                    for(ik=0;ik<5;ik++)
                        rspace[ik*nl+il]*=Q;
                }
            }
            else
            {
                if(g2->AngMom==-2)
                {
                    for(ik=0;ik<25;ik++)
                        rspace[ik]*=Q;
                }
                else
                {
                    for(ik=0;ik<35;ik++)
                        rspace[ik]*=Q;
                }
            }
        }
        else
        {
            if(g2->AngMom>=0)
            {
                for(il=0;il<nl;il++)
                {
                    for(ik=0;ik<7;ik++)
                        rspace[ik*7+il]*=Q;
                }
            }
            else
            {
                if(g2->AngMom==-2)
                {
                    for(il=0;il<35;il++)
                        rspace[il]*=Q;
                }
                else
                {
                    for(il=0;il<49;il++)
                        rspace[il]*=Q;
                }
            }
        }
    }
    free((void*)workspace);
    for(i=0;i<nk;i++)
        free((void*)k[i]);
    for(i=0;i<nl;i++)
        free((void*)l[i]);
    return 0;
}

void InitMatrices()
{
    int i;

    for(i=0;i<6;i++)
        Nmat[i]=(double*)malloc(400*sizeof(double));
    for(i=0;i<400;i++)
        Kmat[i]=-9.e99;
}

void ClearMatrices()
{
    for(int i=0;i<6;i++)
        free((void*)Nmat[i]);
}
}

 
} // namespace Fchk
} // namespace Plugin
