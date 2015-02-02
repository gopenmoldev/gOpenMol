/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
/*Implements GaussianOrbital and GaussAtom classes.  Also contains crude and ineffective
implementation of OS recursion formulas.  In particular, Nuclear Attraction integrals are
fatal for non-atomic calculations.*/
#include "gaussians.h"

namespace Plugin {
namespace Fchk {
    
#define pi 3.1415926535897932384626433832795

void GaussAtom::AssignType()
{
    char pt[56][3]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl",
                "Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As",
                "Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In",
                "Sn","Sb","Te","I","Xe","Cs","Ba"};
    
    sprintf(Element,"%s",pt[AtNum-1]);
}



GaussianOrbital GaussianOrbital::operator=(GaussianOrbital g)
{
    GaussianOrbital gr;
    int i;

    gr.OrbExponent=OrbExponent=g.OrbExponent;
    gr.AngMom=AngMom=g.AngMom;
    gr.NormCoeff=NormCoeff=g.NormCoeff;
    for(i=0;i<3;i++)
    {
        gr.AngMomV[i]=AngMomV[i]=g.AngMomV[i];
        gr.Center[i]=Center[i]=g.Center[i];
    }
    return gr;
}

double GaussianOrbital::Evaluate(double *r)
{
    int i,j;
    double r0[3];
    double rr;
    double geval;

    for(i=0;i<3;i++)
    {
        r0[i]=r[i]-Center[i];
    }
    rr=r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2];
    geval=NormCoeff*exp(-OrbExponent*rr);
    if(AngMom>=0)
    {
        for(i=0;i<3;i++)
        {
            for(j=0;j<AngMomV[i];j++)
            {
                geval*=r0[i];
            }
        }
    }
    else
    {
        if(AngMom==-2)
        {
/*Now we have the native d-functions to deal with.  We adopt the following convention.
For the xy,xz, and yz orbitals, leave things as they are.  For z2, angmomV[2]=-2
for x2-y2, angmomv[0]=angmomv[1]=-2*/
            if((AngMomV[0]>=0)&&(AngMomV[1]>=0)&&(AngMomV[2]>=0))
            {
                for(i=0;i<3;i++)
                {
                    for(j=0;j<AngMomV[i];j++)
                    {
                        geval*=r0[i];
                    }
                }
            }
            else
            {
                if(AngMomV[2]==-2)
                {
                    geval*=3.0*r0[2]*r0[2]-rr;
                }   
                else
                {
                    geval*=r0[0]*r0[0]-r0[1]*r0[1];
                }
            }
        }
        else
        {
            switch(AngMomV[2])
            {
            case -3:    geval*=(3.0*r0[0]*r0[0]-r0[1]*r0[1])*r0[1];
                        break;
            case -2:    geval*=r0[0]*r0[1]*r0[2];
                        break;
            case -1:    geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[1];
                        break;
            case 0:     geval*=(2.0*r0[2]*r0[2]-3.0*r0[1]*r0[1]-3.0*r0[0]*r0[0])*r0[2];
                        break;
            case 1:     geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[0];
                        break;
            case 2:     geval*=(r0[0]*r0[0]-r0[1]*r0[1])*r0[2];
                        break;
            case 3:     geval*=(r0[0]*r0[0]-3.0*r0[1]*r0[1])*r0[0];
                        break;
            default:    break;
            }
        }
    }
    return geval;
}

double GaussianOrbital::Evaluate(double x,double y,double z)
{
    int i,j;
    double r0[3];
    double rr;
    double geval;
    
    r0[0]=x-Center[0];
    r0[1]=y-Center[1];
    r0[2]=z-Center[2];
    rr=r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2];
    geval=NormCoeff*exp(-OrbExponent*rr);
    if(AngMom>=0)
    {
        for(i=0;i<3;i++)
        {
            for(j=0;j<AngMomV[i];j++)
            {
                geval*=r0[i];
            }
        }
    }
    else
    {
        if(AngMom==-2)
        {
/*Now we have the native d-functions to deal with.  We adopt the following convention.
For the xy,xz, and yz orbitals, leave things as they are.  For z2, angmomV[2]=-2
for x2-y2, angmomv[0]=angmomv[1]=-2*/
            if((AngMomV[0]>=0)&&(AngMomV[1]>=0)&&(AngMomV[2]>=0))
            {
                for(i=0;i<3;i++)
                {
                    for(j=0;j<AngMomV[i];j++)
                    {
                        geval*=r0[i];
                    }
                }
            }
            else
            {
                if(AngMomV[2]==-2)
                {
                    geval*=3.0*r0[2]*r0[2]-rr;
                }   
                else
                {
                    geval*=r0[0]*r0[0]-r0[1]*r0[1];
                }
            }
        }
        else
        {
            switch(AngMomV[2])
            {
            case -3:    geval*=(3.0*r0[0]*r0[0]-r0[1]*r0[1])*r0[1];
                        break;
            case -2:    geval*=r0[0]*r0[1]*r0[2];
                        break;
            case -1:    geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[1];
                        break;
            case 0:     geval*=(2.0*r0[2]*r0[2]-3.0*r0[1]*r0[1]-3.0*r0[0]*r0[0])*r0[2];
                        break;
            case 1:     geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[0];
                        break;
            case 2:     geval*=(r0[0]*r0[0]-r0[1]*r0[1])*r0[2];
                        break;
            case 3:     geval*=(r0[0]*r0[0]-3.0*r0[1]*r0[1])*r0[0];
                        break;
            default:    break;
            }
        }
    }
    return geval;
}

double GaussianOrbital::EvaluateGradient(double* r,double *Grad)
{
    double geval;
    double g2,gx;
    double r0[3];
    double rr;
    int i,j,k;  
    
    for(i=0;i<3;i++)
        r0[i]=r[i]-Center[i];
    rr=r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2];
    gx=geval=NormCoeff*exp(-OrbExponent*rr);
    if(AngMom>=0)
    {
        for(i=0;i<3;i++)
        {
            for(j=0;j<AngMomV[i];j++)
            {
                geval*=r0[i];
            }
        }
    }
    else
    {
        if(AngMom==-2)
        {
/*Now we have the native d-functions to deal with.  We adopt the following convention.
For the xy,xz, and yz orbitals, leave things as they are.  For z2, angmomV[2]=-2
for x2-y2, angmomv[0]=angmomv[1]=-2*/
            if((AngMomV[0]>=0)&&(AngMomV[1]>=0)&&(AngMomV[2]>=0))
            {
                for(i=0;i<3;i++)
                {
                    for(j=0;j<AngMomV[i];j++)
                    {
                        geval*=r0[i];
                    }
                }
            }
            else
            {
                if(AngMomV[2]==-2)
                {
                    geval*=3.0*r0[2]*r0[2]-rr;
                }   
                else
                {
                    geval*=r0[0]*r0[0]-r0[1]*r0[1];
                }
            }
        }
        else
        {
            switch(AngMomV[2])
            {
            case -3:    geval*=(3.0*r0[0]*r0[0]-r0[1]*r0[1])*r0[1];
                        break;
            case -2:    geval*=r0[0]*r0[1]*r0[2];
                        break;
            case -1:    geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[1];
                        break;
            case 0:     geval*=(2.0*r0[2]*r0[2]-3.0*r0[1]*r0[1]-3.0*r0[0]*r0[0])*r0[2];
                        break;
            case 1:     geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[0];
                        break;
            case 2:     geval*=(r0[0]*r0[0]-r0[1]*r0[1])*r0[2];
                        break;
            case 3:     geval*=(r0[0]*r0[0]-3.0*r0[1]*r0[1])*r0[0];
                        break;
            default:    break;
            }
        }
    }
    for(i=0;i<3;i++)
    {
        Grad[i]=-2.0*OrbExponent*r0[i]*geval;
    }
    switch(AngMom)
    {
    case -3:    switch(AngMomV[2])
                {
                case -3:    Grad[0]+=3.0*(r0[0]*r0[0]-r0[1]*r0[1])*gx;
                            Grad[1]-=6.0*r0[0]*r0[1]*gx;
                            break;
                case -2:    Grad[0]+=r0[1]*r0[2]*gx;
                            Grad[1]+=r0[0]*r0[2]*gx;
                            Grad[2]+=r0[0]*r0[1]*gx;
                            break;
                case -1:    Grad[0]-=2.0*r0[1]*r0[0]*gx;
                            Grad[1]+=(4.0*r0[2]*r0[2]-3.0*r0[1]*r0[1]-r0[0]*r0[0])*gx;
                            Grad[2]+=8.0*r0[1]*r0[2]*gx;
                            break;
                case 0:     Grad[0]-=6.0*r0[0]*r0[2]*gx;
                            Grad[1]-=6.0*r0[1]*r0[2]*gx;
                            Grad[2]+=(2.0*r0[2]*r0[2]-3.0*(r0[0]*r0[0]+r0[1]*r0[1]))*gx;
                            break;
                case 1:     Grad[0]+=(4.0*r0[2]*r0[2]-3.0*r0[0]*r0[0]-r0[1]*r0[1])*gx;
                            Grad[1]-=2.0*r0[0]*r0[1]*gx;
                            Grad[2]+=8.0*r0[0]*r0[2]*gx;
                            break;
                case 2:     Grad[0]+=2.0*r0[0]*r0[2]*gx;
                            Grad[1]-=2.0*r0[1]*r0[2]*gx;
                            Grad[2]+=(r0[0]*r0[0]-r0[1]*r0[1])*gx;
                            break;
                case 3:     Grad[0]=3.0*(r0[0]*r0[0]-r0[1]*r0[1])*gx;
                            Grad[1]-=6.0*r0[0]*r0[1];
                            break;
                default:    break;
                }
                break;
    case -2:    if(AngMomV[0]<0)
                {
                    Grad[0]+=2.0*r0[0]*gx;
                    Grad[1]-=2.0*r0[1]*gx;
                }
                else
                {
                    if(AngMomV[2]<0)
                    {
                        Grad[0]-=2.0*r0[0]*gx;
                        Grad[1]-=2.0*r0[1]*gx;
                        Grad[2]+=4.0*r0[2]*gx;
                    }
                    else
                    {
                        for(i=0;i<3;i++)
                        {
                            g2=AngMomV[i];
                            if(g2>0)
                            {
                                for(j=0;j<3;j++)
                                {
                                    if(i==j)
                                    {
                                        for(k=0;k<AngMomV[j]-1;k++)
                                            g2*=r0[j];
                                    }
                                    else
                                    {
                                        for(k=0;k<AngMomV[j];k++)
                                            g2*=r0[j];
                                    }
                                }
                            }
                            Grad[i]+=g2*gx;
                        }
                    }
                }
                break;
    case 0:     break;
    case 1:     for(i=0;i<3;i++)
                {
                    g2=AngMomV[i]*gx;
                    Grad[i]+=g2;
                }
                break;
    case 2:
    case 3:     for(i=0;i<3;i++)
                {
                    g2=AngMomV[i];
                    if(g2>0)
                    {
                        for(j=0;j<3;j++)
                        {
                            if(i!=j)
                            {
                                for(k=0;k<AngMomV[j]-1;k++)
                                    g2*=r0[j];
                            }
                            else
                            {
                                for(k=0;k<AngMomV[j];k++)
                                    g2*=r0[j];
                            }
                        }
                    }
                    Grad[i]+=g2*gx;
                }
                break;
    default:    break;
    }
    return geval;
}

double GaussianOrbital::EvaluateGradient(double x,double y,double z,double *Grad)
{
    double geval;
    double g2,gx;
    double r0[3];
    double rr;
    int i,j,k;  
    
    r0[0]=x-Center[0];
    r0[1]=y-Center[1];
    r0[2]=z-Center[2];
    rr=r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2];
    gx=geval=NormCoeff*exp(-OrbExponent*rr);
    if(AngMom>=0)
    {
        for(i=0;i<3;i++)
        {
            for(j=0;j<AngMomV[i];j++)
            {
                geval*=r0[i];
            }
        }
    }
    else
    {
        if(AngMom==-2)
        {
/*Now we have the native d-functions to deal with.  We adopt the following convention.
For the xy,xz, and yz orbitals, leave things as they are.  For z2, angmomV[2]=-2
for x2-y2, angmomv[0]=angmomv[1]=-2*/
            if((AngMomV[0]>=0)&&(AngMomV[1]>=0)&&(AngMomV[2]>=0))
            {
                for(i=0;i<3;i++)
                {
                    for(j=0;j<AngMomV[i];j++)
                    {
                        geval*=r0[i];
                    }
                }
            }
            else
            {
                if(AngMomV[2]==-2)
                {
                    geval*=3.0*r0[2]*r0[2]-rr;
                }   
                else
                {
                    geval*=r0[0]*r0[0]-r0[1]*r0[1];
                }
            }
        }
        else
        {
            switch(AngMomV[2])
            {
            case -3:    geval*=(3.0*r0[0]*r0[0]-r0[1]*r0[1])*r0[1];
                        break;
            case -2:    geval*=r0[0]*r0[1]*r0[2];
                        break;
            case -1:    geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[1];
                        break;
            case 0:     geval*=(2.0*r0[2]*r0[2]-3.0*r0[1]*r0[1]-3.0*r0[0]*r0[0])*r0[2];
                        break;
            case 1:     geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[0];
                        break;
            case 2:     geval*=(r0[0]*r0[0]-r0[1]*r0[1])*r0[2];
                        break;
            case 3:     geval*=(r0[0]*r0[0]-3.0*r0[1]*r0[1])*r0[0];
                        break;
            default:    break;
            }
        }
    }
    for(i=0;i<3;i++)
    {
        Grad[i]=-2.0*OrbExponent*r0[i]*geval;
    }
    switch(AngMom)
    {
    case -3:    switch(AngMomV[2])
                {
                case -3:    Grad[0]+=3.0*(r0[0]*r0[0]-r0[1]*r0[1])*gx;
                            Grad[1]-=6.0*r0[0]*r0[1]*gx;
                            break;
                case -2:    Grad[0]+=r0[1]*r0[2]*gx;
                            Grad[1]+=r0[0]*r0[2]*gx;
                            Grad[2]+=r0[0]*r0[1]*gx;
                            break;
                case -1:    Grad[0]-=2.0*r0[1]*r0[0]*gx;
                            Grad[1]+=(4.0*r0[2]*r0[2]-3.0*r0[1]*r0[1]-r0[0]*r0[0])*gx;
                            Grad[2]+=8.0*r0[1]*r0[2]*gx;
                            break;
                case 0:     Grad[0]-=6.0*r0[0]*r0[2]*gx;
                            Grad[1]-=6.0*r0[1]*r0[2]*gx;
                            Grad[2]+=(2.0*r0[2]*r0[2]-3.0*(r0[0]*r0[0]+r0[1]*r0[1]))*gx;
                            break;
                case 1:     Grad[0]+=(4.0*r0[2]*r0[2]-3.0*r0[0]*r0[0]-r0[1]*r0[1])*gx;
                            Grad[1]-=2.0*r0[0]*r0[1]*gx;
                            Grad[2]+=8.0*r0[0]*r0[2]*gx;
                            break;
                case 2:     Grad[0]+=2.0*r0[0]*r0[2]*gx;
                            Grad[1]-=2.0*r0[1]*r0[2]*gx;
                            Grad[2]+=(r0[0]*r0[0]-r0[1]*r0[1])*gx;
                            break;
                case 3:     Grad[0]=3.0*(r0[0]*r0[0]-r0[1]*r0[1])*gx;
                            Grad[1]-=6.0*r0[0]*r0[1];
                            break;
                default:    break;
                }
                break;
    case -2:    if(AngMomV[0]<0)
                {
                    Grad[0]+=2.0*r0[0]*gx;
                    Grad[1]-=2.0*r0[1]*gx;
                }
                else
                {
                    if(AngMomV[2]<0)
                    {
                        Grad[0]-=2.0*r0[0]*gx;
                        Grad[1]-=2.0*r0[1]*gx;
                        Grad[2]+=4.0*r0[2]*gx;
                    }
                    else
                    {
                        for(i=0;i<3;i++)
                        {
                            g2=AngMomV[i];
                            if(g2>0)
                            {
                                for(j=0;j<3;j++)
                                {
                                    if(i==j)
                                    {
                                        for(k=0;k<AngMomV[j]-1;k++)
                                            g2*=r0[j];
                                    }
                                    else
                                    {
                                        for(k=0;k<AngMomV[j];k++)
                                            g2*=r0[j];
                                    }
                                }
                            }
                            Grad[i]+=g2*gx;
                        }
                    }
                }
                break;
    case 0:     break;
    case 1:     for(i=0;i<3;i++)
                {
                    g2=AngMomV[i]*gx;
                    Grad[i]+=g2;
                }
                break;
    case 2:
    case 3:     for(i=0;i<3;i++)
                {
                    g2=AngMomV[i];
                    if(g2>0)
                    {
                        for(j=0;j<3;j++)
                        {
                            if(i!=j)
                            {
                                for(k=0;k<AngMomV[j]-1;k++)
                                    g2*=r0[j];
                            }
                            else
                            {
                                for(k=0;k<AngMomV[j];k++)
                                    g2*=r0[j];
                            }
                        }
                    }
                    Grad[i]+=g2*gx;
                }
                break;
    default:    break;
    }
    return geval;
}

double GaussianOrbital::EvaluateLaplacian(double *r)
{
    int i,j,k;
    double r0[3];
    double rr;
    double geval,gx;
    double Grad[3],Gr2[3];
    double lplc;
    double g2;
    
    for(i=0;i<3;i++)
        r0[i]=r[i]-Center[i];
    rr=r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2];
    gx=geval=NormCoeff*exp(-OrbExponent*rr);
    if(AngMom>=0)
    {
        for(i=0;i<3;i++)
        {
            for(j=0;j<AngMomV[i];j++)
            {
                geval*=r0[i];
            }
        }
    }
    else
    {
        if(AngMom==-2)
        {
/*Now we have the native d-functions to deal with.  We adopt the following convention.
For the xy,xz, and yz orbitals, leave things as they are.  For z2, angmomV[2]=-2
for x2-y2, angmomv[0]=angmomv[1]=-2*/
            if((AngMomV[0]>=0)&&(AngMomV[1]>=0)&&(AngMomV[2]>=0))
            {
                for(i=0;i<3;i++)
                {
                    for(j=0;j<AngMomV[i];j++)
                    {
                        geval*=r0[i];
                    }
                }
            }
            else
            {
                if(AngMomV[2]==-2)
                {
                    geval*=3.0*r0[2]*r0[2]-rr;
                }   
                else
                {
                    geval*=r0[0]*r0[0]-r0[1]*r0[1];
                }
            }
        }
        else
        {
            switch(AngMomV[2])
            {
            case -3:    geval*=(3.0*r0[0]*r0[0]-r0[1]*r0[1])*r0[1];
                        break;
            case -2:    geval*=r0[0]*r0[1]*r0[2];
                        break;
            case -1:    geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[1];
                        break;
            case 0:     geval*=(2.0*r0[2]*r0[2]-3.0*r0[1]*r0[1]-3.0*r0[0]*r0[0])*r0[2];
                        break;
            case 1:     geval*=(4.0*r0[2]*r0[2]-r0[1]*r0[1]-r0[0]*r0[0])*r0[0];
                        break;
            case 2:     geval*=(r0[0]*r0[0]-r0[1]*r0[1])*r0[2];
                        break;
            case 3:     geval*=(r0[0]*r0[0]-3.0*r0[1]*r0[1])*r0[0];
                        break;
            default:    break;
            }
        }
    }
    for(i=0;i<3;i++)
    {
        Grad[i]=0.0;
        Gr2[i]=0.0;
    }
    switch(AngMom)
    {
    case -3:    switch(AngMomV[2])
                {
                case -3:    Grad[0]+=6.0*r0[0]*r0[1];
                            Grad[1]+=3.0*(r0[0]*r0[0]-r0[1]*r0[1]);
                            Gr2[0]+=6.0*r0[1];
                            Gr2[1]-=6.0*r0[1];
                            break;
                case -2:    Grad[0]+=r0[1]*r0[2];
                            Grad[1]+=r0[0]*r0[2];
                            Grad[2]+=r0[0]*r0[1];
                            break;
                case -1:    Grad[0]-=2.0*r0[1]*r0[0];
                            Grad[1]+=4.0*r0[2]*r0[2]-r0[0]*r0[0]-3.0*r0[1]*r0[1];
                            Grad[2]+=8.0*r0[2]*r0[1];
                            Gr2[0]-=2.0*r0[1];
                            Gr2[1]-=6.0*r0[1];
                            Gr2[2]+=8.0*r0[1];
                            break;
                case 0:     Grad[0]-=6.0*r0[0]*r0[2];
                            Grad[1]-=6.0*r0[1]*r0[2];
                            Grad[2]+=6.0*r0[2]*r0[2]-3.0*(r0[1]*r0[1]+r0[0]*r0[0]);
                            Gr2[0]-=6.0*r0[2];
                            Gr2[1]-=6.0*r0[2];
                            Gr2[2]+=12.0*r0[2];
                            break;
                case 1:     Grad[0]+=4.0*r0[2]*r0[2]-3.0*r0[0]*r0[0]-r0[1]*r0[1];
                            Grad[1]-=2.0*r0[1]*r0[0];
                            Grad[2]+=8.0*r0[0]*r0[2];
                            Gr2[0]-=6.0*r0[0];
                            Gr2[1]-=2.0*r0[0];
                            Gr2[2]+=8.0*r0[0];
                            break;
                case 2:     Grad[0]+=2.0*r0[0]*r0[2];
                            Grad[1]-=2.0*r0[0]*r0[2];
                            Grad[2]+=(r0[0]*r0[0]-r0[1]*r0[1]);
                            Gr2[0]+=2.0*r0[2];
                            Gr2[1]-=2.0*r0[2];
                            break;
                case 3:     Grad[0]+=3.0*r0[0]*r0[0]-3.0*r0[1]*r0[1];
                            Grad[1]-=6.0*r0[0]*r0[1];
                            Gr2[0]+=6.0*r0[0];
                            Gr2[1]-=6.0*r0[0];
                            break;
                default:    break;
                }
                break;
    case -2:    if(AngMomV[0]<0)
                {
                    Grad[0]+=2.0*r0[0];
                    Grad[1]-=2.0*r0[1];
                    Gr2[0]+=2.0;
                    Gr2[1]-=2.0;
                }
                else
                {
                    if(AngMomV[2]<0)
                    {
                        Grad[0]-=2.0*r0[0];
                        Grad[1]-=2.0*r0[1];
                        Grad[2]+=4.0*r0[2];
                        Gr2[0]+=2.0;
                        Gr2[1]-=2.0;
                    }
                    else
                    {
                        for(i=0;i<3;i++)
                        {
                            g2=AngMomV[i];
                            if(g2==2)
                                Gr2[i]+=2.0;
                            if(g2>0)
                            {
                                for(j=0;j<3;j++)
                                {
                                    if(i==j)
                                    {
                                        for(k=0;k<AngMomV[j]-1;k++)
                                            g2*=r0[j];
                                    }
                                    else
                                    {
                                        for(k=0;k<AngMomV[j];k++)
                                            g2*=r0[j];
                                    }
                                }
                            }
                            Grad[i]+=g2;
                        }
                    }
                }
                break;
    case 0:     break;
    case 1:     for(i=0;i<3;i++)
                {
                    Grad[i]+=AngMomV[i];
                }
                break;
    case 2:
    case 3:     for(i=0;i<3;i++)
                {
                    g2=AngMomV[i];
                    if(g2==2)
                        Gr2[i]+=2.0;
                    if(g2==3)
                        Gr2[i]+=6.0*r0[i];
                    if(g2>0)
                    {
                        for(j=0;j<3;j++)
                        {
                            if(i==j)
                            {
                                for(k=0;k<AngMomV[j]-1;k++)
                                    g2*=r0[j];
                            }
                            else
                            {
                                for(k=0;k<AngMomV[j];k++)
                                    g2*=r0[j];
                            }
                        }
                    }
                    Grad[i]+=g2;
                }
                break;
    default:    break;
    }
    lplc=0.0;
    for(i=0;i<3;i++)
    {
        lplc+=gx*(Gr2[i]-4.0*OrbExponent*r0[i]*Grad[i])+geval*2.0*OrbExponent*(2.0*OrbExponent*r0[i]*r0[i]-1.0);
    }
    return lplc;
}

double GaussianOrbital::EvaluateLaplacian(double x,double y,double z)
{
    int i,j,k;
    double r0[3];
    double rr;
    double geval,gx;
    double Grad[3],Gr2[3];
    double lplc;
    double g2;
    
    r0[0]=x-Center[0];
    r0[1]=y-Center[1];
    r0[2]=z-Center[2];
    rr=r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2];
    gx=geval=NormCoeff*exp(-OrbExponent*rr);
    if((AngMomV[0]>=0)&&(AngMomV[1]>=0)&&(AngMomV[2]>=0))
    {
        for(i=0;i<3;i++)
        {
            for(j=0;j<AngMomV[i];j++)
            {
                geval*=r0[i];
            }
        }
    }
    else
    {
/*Now we have the native d-functions to deal with.  We adopt the following convention.
For the xy,xz, and yz orbitals, leave things as they are.  For z2, angmomV[2]=-2
for x2-y2, angmomv[0]=angmomv[1]=-2*/
        if(AngMomV[2]==-2)
        {
            geval*=3.0*r0[2]*r0[2]-rr;
        }
        else
        {
            geval*=r0[0]*r0[0]-r0[1]*r0[1];
        }
    }
    for(i=0;i<3;i++)
    {
        Grad[i]=0.0;
        Gr2[i]=0.0;
    }
    switch(AngMom)
    {
    case -2:    if(AngMomV[0]<0)
                {
                    Grad[0]+=2.0*r0[0];
                    Grad[1]-=2.0*r0[1];
                    Gr2[0]+=2.0;
                    Gr2[1]-=2.0;
                }
                else
                {
                    if(AngMomV[2]<0)
                    {
                        Grad[0]-=2.0*r0[0];
                        Grad[1]-=2.0*r0[1];
                        Grad[2]+=4.0*r0[2];
                        Gr2[0]+=2.0;
                        Gr2[1]-=2.0;
                    }
                    else
                    {
                        for(i=0;i<3;i++)
                        {
                            g2=AngMomV[i];
                            if(g2==2)
                                Gr2[i]+=2.0;
                            if(g2>0)
                            {
                                for(j=0;j<3;j++)
                                {
                                    if(i==j)
                                    {
                                        for(k=0;k<AngMomV[j]-1;k++)
                                            g2*=r0[j];
                                    }
                                    else
                                    {
                                        for(k=0;k<AngMomV[j];k++)
                                            g2*=r0[j];
                                    }
                                }
                            }
                            Grad[i]+=g2;
                        }
                    }
                }
                break;
    case 0:     break;
    case 1:     for(i=0;i<3;i++)
                {
                    Grad[i]+=AngMomV[i];
                }
                break;
    case 2:     for(i=0;i<3;i++)
                {
                    g2=AngMomV[i];
                    if(g2==2)
                        Gr2[i]+=2.0;
                    if(g2>0)
                    {
                        for(j=0;j<3;j++)
                        {
                            if(i==j)
                            {
                                for(k=0;k<AngMomV[j]-1;k++)
                                    g2*=r0[j];
                            }
                            else
                            {
                                for(k=0;k<AngMomV[j];k++)
                                    g2*=r0[j];
                            }
                        }
                    }
                    Grad[i]+=g2;
                }
                break;
    default:    break;
    }
    lplc=0.0;
    for(i=0;i<3;i++)
    {
        lplc+=gx*(Gr2[i]-4.0*OrbExponent*r0[i]*Grad[i])+geval*2.0*OrbExponent*(2.0*OrbExponent*r0[i]*r0[i]-1.0);
    }
    return lplc;
}

namespace oldrecurse {

double OIss(double zeta,double xi,double rr)
{
    return exp(1.5*log(pi/zeta)-xi*rr);
}

double OIps(double zeta,double xi,double rr,double dx)
{
    return dx*OIss(zeta,xi,rr);
}

double OIpp(double zeta,double xi,double rr,double dx,double dy,int i,int j)
{
    double rval;

    if(i==j)
        rval=OIss(zeta,xi,rr)*0.5/zeta;
    else
        rval=0.0;
    return dy*OIps(zeta,xi,rr,dx)+rval;
}

double OIds(double zeta,double xi,double rr,double dx,double dy,int i,int j)
{
    double rval;

    if(i==j)
        rval=OIss(zeta,xi,rr)*0.5/zeta;
    else
        rval=0.0;
    return dy*OIps(zeta,xi,rr,dx)+rval;
}

double OIdp(double zeta,double xi,double rr,double dx,double dy,double dz,int i,int j,int k)
{
    double rval;

    rval=dz*OIds(zeta,xi,rr,dx,dy,i,j);
    if(i==k)
        rval+=OIps(zeta,xi,rr,dy)*0.5/zeta;
    if(j==k)
        rval+=OIps(zeta,xi,rr,dx)*0.5/zeta;
    return rval;
}

double OIdd(double zeta,double xi,double rr,double dx,double dy,double dz,double dw,int i,int j,int k,int l)
{
    double rval;

    rval=dw*OIdp(zeta,xi,rr,dx,dy,dz,i,j,k);
    if(i==l)
        rval+=OIpp(zeta,xi,rr,dy,dz,j,k)*0.5/zeta;
    if(j==l)
        rval+=OIpp(zeta,xi,rr,dx,dz,i,k)*0.5/zeta;
    if(k==l)
        rval+=OIds(zeta,xi,rr,dx,dy,i,j)*0.5/zeta;
    return rval;
}

double OIfs(double zeta,double xi,double rr,double dx,double dy,double dz,int i,int j,int k)
{
    double rval;

    rval=dz*OIds(zeta,xi,rr,dx,dy,i,j);
    if(i==k)
        rval+=OIps(zeta,xi,rr,dy)*0.5/zeta;
    if(j==k)
        rval+=OIps(zeta,xi,rr,dx)*0.5/zeta;
    return rval;
}

double OIfp(double zeta,double xi,double rr,double dx,double dy,double dz,double dp,int i,int j,int k,int r)
{
    double rval;

    rval=dz*OIdp(zeta,xi,rr,dx,dy,dp,i,j,r);
    if(i==k)
        rval+=OIpp(zeta,xi,rr,dy,dp,j,r)*0.5/zeta;
    if(j==k)
        rval+=OIpp(zeta,xi,rr,dx,dp,i,r)*0.5/zeta;
    if(r==k)
        rval+=OIds(zeta,xi,rr,dx,dy,i,j)*0.5/zeta;
    return rval;
}

double OIfd(double zeta,double xi,double rr,double dx,double dy,double dz,double dp,double dq,int i,int j,int k,int r,int s)
{
    double rval;

    rval=dz*OIdd(zeta,xi,rr,dx,dy,dp,dq,i,j,r,s);
    if(i==k)
        rval+=OIdp(zeta,xi,rr,dp,dq,dy,r,s,j)*0.5/zeta;
    if(j==k)
        rval+=OIdp(zeta,xi,rr,dp,dq,dx,r,s,i)*0.5/zeta;
    if(r==k)
        rval+=OIdp(zeta,xi,rr,dx,dy,dq,i,j,s)*0.5/zeta;
    if(s==k)
        rval+=OIdp(zeta,xi,rr,dx,dy,dp,i,j,r)*0.5/zeta;
    return rval;
}

double OIff(double zeta,double xi,double rr,double dx,double dy,double dz,double dp,double dq,double dr,int i,int j,int k,int r,int s,int t)
{
    double rval;

    rval=dz*OIfd(zeta,xi,rr,dp,dq,dr,dx,dy,r,s,t,i,j);
    if(i==k)
        rval+=OIfp(zeta,xi,rr,dp,dq,dr,dy,r,s,t,j)*0.5/zeta;
    if(j==k)
        rval+=OIfp(zeta,xi,rr,dp,dq,dr,dx,r,s,t,i)*0.5/zeta;
    if(r==k)
        rval+=OIdd(zeta,xi,rr,dx,dy,dq,dr,i,j,s,t)*0.5/zeta;
    if(s==k)
        rval+=OIdd(zeta,xi,rr,dx,dy,dp,dr,i,j,r,t)*0.5/zeta;
    if(t==k)
        rval+=OIdd(zeta,xi,rr,dx,dy,dp,dq,i,j,r,s)*0.5/zeta;
    return rval;
}

double KIss(double zeta,double xi,double rr)
{
    return xi*(3.0-2.0*xi*rr)*OIss(zeta,xi,rr);
}

double KIps(double zeta,double xi,double rr,double dx)
{
    return dx*KIss(zeta,xi,rr)+2*xi*OIps(zeta,xi,rr,dx);
}

double KIpp(double zeta,double xi,double rr,double dx,double dy,int i,int j)
{
    double rval;

    rval=dy*KIps(zeta,xi,rr,dx)+2.0*xi*OIpp(zeta,xi,rr,dx,dy,i,j);
    if(i==j)
        rval+=KIss(zeta,xi,rr)*0.5/zeta;
    return rval;
}

double KIds(double zeta,double xi,double rr,double dx,double dy,int i,int j,double za)
{
    double rval;

    rval=dy*KIps(zeta,xi,rr,dx)+2.0*xi*OIds(zeta,xi,rr,dx,dy,i,j);
    if(i==j)
        rval+=KIss(zeta,xi,rr)*0.5/zeta-xi/za*OIss(zeta,xi,rr);
    return rval;
}

double KIdp(double zeta,double xi,double rr,double dx,double dy,double dz,int i,int j,int k,double za)
{
    double rval;

    rval=dz*KIds(zeta,xi,rr,dx,dy,i,j,za)+2.0*xi*OIdp(zeta,xi,rr,dx,dy,dz,i,j,k);
    if(i==k)
        rval+=KIps(zeta,xi,rr,dy)*0.5/zeta;
    if(j==k)
        rval+=KIps(zeta,xi,rr,dx)*0.5/zeta;
    return rval;
}

double KIdd(double zeta,double xi,double rr,double dx,double dy,double dz,double dw,int i,int j,int k,int l,double za,double zb)
{
    double rval;

    rval=dw*KIdp(zeta,xi,rr,dx,dy,dz,i,j,k,za)+2.0*xi*OIdd(zeta,xi,rr,dx,dy,dz,dw,i,j,k,l);
    if(i==l)
        rval+=KIpp(zeta,xi,rr,dy,dz,j,k)*0.5/zeta;
    if(j==l)
        rval+=KIpp(zeta,xi,rr,dx,dz,i,k)*0.5/zeta;
    if(k==l)
        rval+=KIds(zeta,xi,rr,dx,dy,i,j,za)*0.5/zeta-xi/zb*OIds(zeta,xi,rr,dx,dy,i,j);
    return rval;
}

double KIfs(double zeta,double xi,double rr,double dx,double dy,double dz,int i,int j,int k,double za)
{
    double rval,tval;

    rval=dz*KIds(zeta,xi,rr,dx,dy,i,j,za);
    if(i==k)
        rval+=KIps(zeta,xi,rr,dy)*0.5/zeta;
    if(j==k)
        rval+=KIps(zeta,xi,rr,dx)*0.5/zeta;
    tval=OIfs(zeta,xi,rr,dx,dy,dz,i,j,k);
    if(i==k)
        tval-=OIps(zeta,xi,rr,dy)*0.5/za;
    if(j==k)
        tval-=OIps(zeta,xi,rr,dx)*0.5/za;
    tval*=2.0*xi;
    return rval+tval;
}

double KIfp(double zeta,double xi,double rr,double dx,double dy,double dz,double dp,int i,int j,int k,int r,double za)
{
    double rval,tval;

    rval=dz*KIdp(zeta,xi,rr,dx,dy,dp,i,j,r,za);
    if(i==k)
        rval+=KIpp(zeta,xi,rr,dy,dp,j,r)*0.5/zeta;
    if(j==k)
        rval+=KIpp(zeta,xi,rr,dx,dp,i,r)*0.5/zeta;
    if(r==k)
        rval+=KIds(zeta,xi,rr,dx,dy,i,j,za);
    tval=OIfp(zeta,xi,rr,dx,dy,dz,dp,i,j,k,r);
    if(i==k)
        tval-=OIpp(zeta,xi,rr,dy,dp,j,r)*0.5/za;
    if(j==k)
        tval-=OIpp(zeta,xi,rr,dy,dp,i,r)*0.5/za;
    tval*=2.0*xi;
    return rval+tval;
}

double KIfd(double zeta,double xi,double rr,double dx,double dy,double dz,double dp,double dq,int i,int j,int k,int r,int s,double za,double zb)
{
    double rval,tval;

    rval=dz*KIdd(zeta,xi,rr,dx,dy,dp,dq,i,j,r,s,za,zb);
    if(i==k)
        rval+=KIdp(zeta,xi,rr,dp,dq,dy,r,s,j,zb)*0.5/zeta;
    if(j==k)
        rval+=KIdp(zeta,xi,rr,dp,dq,dx,r,s,i,zb)*0.5/zeta;
    if(r==k)
        rval+=KIdp(zeta,xi,rr,dx,dy,dq,i,j,s,za)*0.5/zeta;
    if(s==k)
        rval+=KIdp(zeta,xi,rr,dx,dy,dp,i,j,r,za)*0.5/zeta;
    tval=OIfd(zeta,xi,rr,dx,dy,dz,dp,dq,i,j,k,r,s);
    if(i==k)
        tval-=OIdp(zeta,xi,rr,dp,dq,dy,r,s,j)*0.5/za;
    if(j==k)
        tval-=OIdp(zeta,xi,rr,dp,dq,dx,r,s,i)*0.5/za;
    tval*=2.0*xi;
    return rval+tval;
}

double KIff(double zeta,double xi,double rr,double dx,double dy,double dz,double dp,double dq,double dr,int i,int j,int k,int r,int s,int t,double za,double zb)
{
    double rval,tval;

    rval=dz*KIfd(zeta,xi,rr,dp,dq,dr,dx,dy,r,s,t,i,j,zb,za);
    if(i==k)
        rval+=KIfp(zeta,xi,rr,dp,dq,dr,dy,r,s,t,j,zb)*0.5/zeta;
    if(j==k)
        rval+=KIfp(zeta,xi,rr,dp,dq,dr,dx,r,s,t,i,zb)*0.5/zeta;
    if(r==k)
        rval+=KIdd(zeta,xi,rr,dx,dy,dq,dr,i,j,s,t,za,zb)*0.5/zeta;
    if(s==k)
        rval+=KIdd(zeta,xi,rr,dx,dy,dp,dr,i,j,r,t,za,zb)*0.5/zeta;
    if(t==k)
        rval+=KIdd(zeta,xi,rr,dx,dy,dp,dq,i,j,r,s,za,zb)*0.5/zeta;
    tval=OIff(zeta,xi,rr,dx,dy,dz,dp,dq,dr,i,j,k,r,s,t);
    if(i==k)
        tval-=OIfp(zeta,xi,rr,dp,dq,dr,dy,r,s,t,j)*0.5/za;
    if(j==k)
        tval-=OIfp(zeta,xi,rr,dp,dq,dr,dx,r,s,t,i)*0.5/za;
    tval*=2.0*xi;
    return rval+tval;
}


double EvalOverlapIntegral(GaussianOrbital g1,GaussianOrbital g2)
{
    int i,j,k,l,c1,m,n;
    int val=10*g1.AngMom+g2.AngMom;
    double zeta,xi;
    double rr,P[3];
    double dx,dy,dz,dw,du,dt;
    
    zeta=g1.OrbExponent+g2.OrbExponent;
    xi=g1.OrbExponent*g2.OrbExponent/zeta;
    rr=0.0;
    for(c1=0;c1<3;c1++)
    {
        rr+=(g1.Center[c1]-g2.Center[c1])*(g1.Center[c1]-g2.Center[c1]);
        P[c1]=(g1.OrbExponent*g1.Center[c1]+g2.OrbExponent*g2.Center[c1])/zeta;
    }
    switch(val)
    {
    case 0: return OIss(zeta,xi,rr);
            break;
    case 1: for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    i=c1;
            }
            dx=P[i]-g2.Center[i];
            return OIps(zeta,xi,rr,dx);
            break;
    case 2: j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            return OIds(zeta,xi,rr,dx,dy,i,j);
            break;
    case 3: k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g2.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            dz=P[k]-g2.Center[k];
            return OIfs(zeta,xi,rr,dx,dy,dz,i,j,k);
            break;
    case 10:for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    i=c1;
            }
            dx=P[i]-g1.Center[i];
            return OIps(zeta,xi,rr,dx);
            break;
    case 11:for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    i=c1;
            }
            dx=P[i]-g1.Center[i];
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    j=c1;
            }
            dy=P[j]-g2.Center[j];
            return OIpp(zeta,xi,rr,dx,dy,i,j);
            break;
    case 12:j=5;
            for(c1=0;c1<3;c1++)
            {   
                if(g2.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    k=c1;
            }
            dz=P[k]-g1.Center[k];
            return OIdp(zeta,xi,rr,dx,dy,dz,i,j,k);
            break;
    case 13:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g2.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            dz=P[k]-g2.Center[k];
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    l=c1;
            }
            dw=P[l]-g1.Center[l];
            return OIfp(zeta,xi,rr,dx,dy,dz,dw,i,j,k,l);
            break;
    case 20:j=5;
            for(c1=0;c1<3;c1++)
            {   
                if(g1.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            return OIds(zeta,xi,rr,dx,dy,i,j);
            break;
    case 21:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    k=c1;
            }
            dz=P[k]-g2.Center[k];
            return OIdp(zeta,xi,rr,dx,dy,dz,i,j,k);
            break;
    case 22:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            l=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(l==5)
                        l=c1;
                    else
                        k=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    k=c1;
                    l=c1;
                }
            }
            dz=P[k]-g2.Center[k];
            dw=P[l]-g2.Center[l];
            return OIdd(zeta,xi,rr,dx,dy,dz,dw,i,j,k,l);
            break;
    case 23:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g2.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            dz=P[k]-g2.Center[k];
            l=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(l==5)
                        l=c1;
                    else
                        m=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    l=c1;
                    m=c1;
                }
            }
            du=P[m]-g1.Center[m];
            dw=P[l]-g1.Center[l];
            return OIfd(zeta,xi,rr,dx,dy,dz,dw,du,i,j,k,l,m);           
            break;
    case 30:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g1.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g1.Center[k];
            return OIfs(zeta,xi,rr,dx,dy,dz,i,j,k);
            break;
    case 31:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g1.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g1.Center[k];
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    l=c1;
            }
            dw=P[l]-g2.Center[l];
            return OIfp(zeta,xi,rr,dx,dy,dz,dw,i,j,k,l);
            break;
    case 32:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g1.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g1.Center[k];
            l=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(l==5)
                        l=c1;
                    else
                        m=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    l=c1;
                    m=c1;
                }
            }
            du=P[m]-g2.Center[m];
            dw=P[l]-g2.Center[l];
            return OIfd(zeta,xi,rr,dx,dy,dz,dw,du,i,j,k,l,m);           
            break;
    case 33:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g1.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g1.Center[k];
            m=5;
            n=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(m==5)
                        m=c1;
                    else if (n==5)
                        n=c1;
                    else
                        l=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    if(m==5)
                    {
                        m=c1;
                        n=c1;
                    }
                    else
                    {
                        n=c1;
                        l=c1;
                    }
                }
                if(g2.AngMomV[c1]==3)
                {
                    m=c1;
                    n=c1;
                    l=c1;
                }
            }
            dw=P[l]-g2.Center[l];
            du=P[n]-g2.Center[n];
            dt=P[m]-g2.Center[m];
            return OIff(zeta,xi,rr,dx,dy,dz,dw,du,dt,i,j,k,l,n,m);
            break;
    default:return 0.0;
    }
}   

double EvalKineticIntegral(GaussianOrbital g1,GaussianOrbital g2)
{
    int i,j,k,l,c1,m,n;
    int val=10*g1.AngMom+g2.AngMom;
    double zeta,xi;
    double rr,P[3];
    double dx,dy,dz,dw,du,dt;
    double za,zb;

    zeta=g1.OrbExponent+g2.OrbExponent;
    za=g1.OrbExponent;
    zb=g2.OrbExponent;
    xi=g1.OrbExponent*g2.OrbExponent/zeta;
    rr=0.0;
    for(c1=0;c1<3;c1++)
    {
        rr+=(g1.Center[c1]-g2.Center[c1])*(g1.Center[c1]-g2.Center[c1]);
        P[c1]=(g1.OrbExponent*g1.Center[c1]+g2.OrbExponent*g2.Center[c1])/zeta;
    }
    switch(val)
    {
    case 0: return KIss(zeta,xi,rr);
            break;
    case 1: for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    i=c1;
            }
            dx=P[i]-g2.Center[i];
            return KIps(zeta,xi,rr,dx);
            break;
    case 2: j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            return KIds(zeta,xi,rr,dx,dy,i,j,zb);
            break;
    case 3: k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g2.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            dz=P[k]-g2.Center[k];
            return KIfs(zeta,xi,rr,dx,dy,dz,i,j,k,zb);
            break;
    case 10:for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    i=c1;
            }
            dx=P[i]-g1.Center[i];
            return KIps(zeta,xi,rr,dx);
            break;
    case 11:for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    i=c1;
            }
            dx=P[i]-g1.Center[i];
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    j=c1;
            }
            dy=P[j]-g2.Center[j];
            return KIpp(zeta,xi,rr,dx,dy,i,j);
            break;
    case 12:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    k=c1;
            }
            dz=P[k]-g1.Center[k];
            return KIdp(zeta,xi,rr,dx,dy,dz,i,j,k,zb);
            break;
    case 13:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g2.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            dz=P[k]-g2.Center[k];
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    l=c1;
            }
            dw=P[l]-g1.Center[l];
            return KIfp(zeta,xi,rr,dx,dy,dz,dw,i,j,k,l,zb);
            break;
    case 20:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            return KIds(zeta,xi,rr,dx,dy,i,j,za);
            break;
    case 21:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    k=c1;
            }
            dz=P[k]-g2.Center[k];
            return KIdp(zeta,xi,rr,dx,dy,dz,i,j,k,za);
            break;
    case 22:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            l=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(l==5)
                        l=c1;
                    else
                        k=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    k=c1;
                    l=c1;
                }
            }
            dz=P[k]-g2.Center[k];
            dw=P[l]-g2.Center[l];
            return KIdd(zeta,xi,rr,dx,dy,dz,dw,i,j,k,l,za,zb);
            break;
    case 23:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g2.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            dz=P[k]-g2.Center[k];
            l=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(l==5)
                        l=c1;
                    else
                        m=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    l=c1;
                    m=c1;
                }
            }
            du=P[m]-g1.Center[m];
            dw=P[l]-g1.Center[l];
            return KIfd(zeta,xi,rr,dx,dy,dz,dw,du,i,j,k,l,m,zb,za);         
            break;
    case 30:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g1.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g1.Center[k];
            return KIfs(zeta,xi,rr,dx,dy,dz,i,j,k,za);
            break;
    case 31:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g1.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g1.Center[k];
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    l=c1;
            }
            dw=P[l]-g2.Center[l];
            return KIfp(zeta,xi,rr,dx,dy,dz,dw,i,j,k,l,za);
            break;
    case 32:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g1.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g1.Center[k];
            l=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(l==5)
                        l=c1;
                    else
                        m=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    l=c1;
                    m=c1;
                }
            }
            du=P[m]-g2.Center[m];
            dw=P[l]-g2.Center[l];
            return KIfd(zeta,xi,rr,dx,dy,dz,dw,du,i,j,k,l,m,za,zb);         
            break;
    case 33:k=5;
            j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(k==5)
                        k=c1;
                    else if (j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    if(k==5)
                    {
                        k=c1;
                        j=c1;
                    }
                    else
                    {
                        j=c1;
                        i=c1;
                    }
                }
                if(g1.AngMomV[c1]==3)
                {
                    k=c1;
                    j=c1;
                    i=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g1.Center[k];
            m=5;
            n=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(m==5)
                        m=c1;
                    else if (n==5)
                        n=c1;
                    else
                        l=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    if(m==5)
                    {
                        m=c1;
                        n=c1;
                    }
                    else
                    {
                        n=c1;
                        l=c1;
                    }
                }
                if(g2.AngMomV[c1]==3)
                {
                    m=c1;
                    n=c1;
                    l=c1;
                }
            }
            dw=P[l]-g2.Center[l];
            du=P[n]-g2.Center[n];
            dt=P[m]-g2.Center[m];
            return KIff(zeta,xi,rr,dx,dy,dz,dw,du,dt,i,j,k,l,n,m,za,zb);
            break;
    default:return 0.0;
    }
}

double dOverlap(GaussianOrbital g1,GaussianOrbital g2)
{
    double val;
    GaussianOrbital g3;

    if(g1.AngMomV[2]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[0]=g3.AngMomV[1]=0;
        g3.AngMomV[2]=2;
        val=2.0*EvalOverlapIntegral(g3,g2);
        g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val-=EvalOverlapIntegral(g3,g2);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=EvalOverlapIntegral(g3,g2);
        return val;
    }
    else if (g1.AngMomV[0]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[1]=g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val=EvalOverlapIntegral(g3,g2);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=EvalOverlapIntegral(g3,g2);
        return val;
    }
    else
    {
        g3=g1;
        g3.AngMom=2;
        return EvalOverlapIntegral(g3,g2);
    }
}

double fOverlap(GaussianOrbital g1,GaussianOrbital g2)
{
    GaussianOrbital g3;
    double val;

    g3=g1;
    g3.AngMom=3;
    switch(g1.AngMomV[2])
    {
    case -3:    g3.AngMomV[0]=2;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=0;
                val=3.0*EvalOverlapIntegral(g3,g2);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=EvalOverlapIntegral(g3,g2);
                break;
    case -2:    g3.AngMomV[0]=1;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=1;
                val=EvalOverlapIntegral(g3,g2);
                break;
    case -1:    g3.AngMomV[0]=0;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=2;
                val=4.0*EvalOverlapIntegral(g3,g2);
                g3.AngMomV[0]=2;
                g3.AngMomV[2]=0;
                val-=EvalOverlapIntegral(g3,g2);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=EvalOverlapIntegral(g3,g2);
                break;
    case 0:     g3.AngMomV[0]=0;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=3;
                val=2.0*EvalOverlapIntegral(g3,g2);
                g3.AngMomV[2]=1;
                g3.AngMomV[0]=2;
                val-=3.0*EvalOverlapIntegral(g3,g2);
                g3.AngMomV[1]=2;
                g3.AngMomV[0]=0;
                val-=3.0*EvalOverlapIntegral(g3,g2);
                break;
    case 1: g3.AngMomV[1]=0;
                g3.AngMomV[0]=1;
                g3.AngMomV[2]=2;
                val=4.0*EvalOverlapIntegral(g3,g2);
                g3.AngMomV[1]=2;
                g3.AngMomV[2]=0;
                val-=EvalOverlapIntegral(g3,g2);
                g3.AngMomV[1]=0;
                g3.AngMomV[0]=3;
                val-=EvalOverlapIntegral(g3,g2);
                break;
    case 2:     g3.AngMomV[0]=2;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=1;
                val=EvalOverlapIntegral(g3,g2);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=2;
                val-=EvalOverlapIntegral(g3,g2);
                break;
    case 3:     g3.AngMomV[0]=3;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=0;
                val=EvalOverlapIntegral(g3,g2);
                g3.AngMomV[0]=1;
                g3.AngMomV[1]=2;
                val-=3.0*EvalOverlapIntegral(g3,g2);
                break;
    default:    break;
    }
    return val;
}

double ddOverlap(GaussianOrbital g1,GaussianOrbital g2)
{
    GaussianOrbital g3;
    double val;

    if(g1.AngMomV[2]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[0]=g3.AngMomV[1]=0;
        g3.AngMomV[2]=2;
        val=2.0*dOverlap(g2,g3);
        g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val-=dOverlap(g2,g3);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=dOverlap(g2,g3);
        return val;
    }
    else if (g1.AngMomV[0]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[1]=g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val=dOverlap(g2,g3);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=dOverlap(g2,g3);
        return val;
    }
    else
    {
        g3=g1;
        g3.AngMom=2;
        return dOverlap(g2,g3);
    }
}

double fdOverlap(GaussianOrbital g1,GaussianOrbital g2)
{
    GaussianOrbital g3;
    double val;

    if(g2.AngMomV[2]==-2)
    {
        g3=g2;
        g3.AngMom=2;
        g3.AngMomV[0]=g3.AngMomV[1]=0;
        g3.AngMomV[2]=2;
        val=2.0*fOverlap(g1,g3);
        g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val-=fOverlap(g1,g3);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=fOverlap(g1,g3);
        return val;
    }
    else if (g2.AngMomV[0]==-2)
    {
        g3=g2;
        g3.AngMom=2;
        g3.AngMomV[1]=g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val=fOverlap(g1,g3);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=fOverlap(g1,g3);
        return val;
    }
    else
    {
        g3=g2;
        g3.AngMom=2;
        return fOverlap(g1,g3);
    }
}

double ffOverlap(GaussianOrbital g1,GaussianOrbital g2)
{
    GaussianOrbital g3;
    double val;

    g3=g2;
    g3.AngMom=3;
    switch(g2.AngMomV[2])
    {
    case -3:    g3.AngMomV[0]=2;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=0;
                val=3.0*fOverlap(g1,g3);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=fOverlap(g1,g3);
                break;
    case -2:    g3.AngMomV[0]=1;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=1;
                val=fOverlap(g1,g3);
                break;
    case -1:    g3.AngMomV[0]=0;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=2;
                val=4.0*fOverlap(g1,g3);
                g3.AngMomV[0]=2;
                g3.AngMomV[2]=0;
                val-=fOverlap(g1,g3);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=fOverlap(g1,g3);
                break;
    case 0:     g3.AngMomV[0]=0;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=3;
                val=2.0*fOverlap(g1,g3);
                g3.AngMomV[2]=1;
                g3.AngMomV[0]=2;
                val-=3.0*fOverlap(g1,g3);
                g3.AngMomV[1]=2;
                g3.AngMomV[0]=0;
                val-=3.0*fOverlap(g1,g3);
                break;
    case 1: g3.AngMomV[1]=0;
                g3.AngMomV[0]=1;
                g3.AngMomV[2]=2;
                val=4.0*fOverlap(g1,g3);
                g3.AngMomV[1]=2;
                g3.AngMomV[2]=0;
                val-=fOverlap(g1,g3);
                g3.AngMomV[1]=0;
                g3.AngMomV[0]=3;
                val-=fOverlap(g1,g3);
                break;
    case 2:     g3.AngMomV[0]=2;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=1;
                val=fOverlap(g1,g3);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=2;
                val-=fOverlap(g1,g3);
                break;
    case 3:     g3.AngMomV[0]=3;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=0;
                val=fOverlap(g1,g3);
                g3.AngMomV[0]=1;
                g3.AngMomV[1]=2;
                val-=3.0*fOverlap(g1,g3);
                break;
    default:    break;
    }
    return val;
}

double OverlapIntegral(GaussianOrbital g1,GaussianOrbital g2)
{
    if((g1.AngMom>=0)&&(g2.AngMom>=0))
    {
        return EvalOverlapIntegral(g1,g2);
    }
    else
    {
        if(g1.AngMom==-2)
        {   //g1 is a pure d function
            if(g2.AngMom==-2)
            {
                return ddOverlap(g1,g2);
            }
            else if(g2.AngMom==-3)
            {
                return fdOverlap(g2,g1);
            }
            else
            {
                return dOverlap(g1,g2);
            }
        }
        else //if (g1.AngMom==-3)
        {
            if(g2.AngMom==-2)
            {
                return fdOverlap(g1,g2);
            }
            else if (g2.AngMom==-3)
            {
                return ffOverlap(g1,g2);
            }
            else
            {
                return fOverlap(g1,g2);
            }
        }
    }
}

double dKinetic(GaussianOrbital g1,GaussianOrbital g2)
{
    double val;
    GaussianOrbital g3;

    if(g1.AngMomV[2]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[0]=g3.AngMomV[1]=0;
        g3.AngMomV[2]=2;
        val=2.0*EvalKineticIntegral(g3,g2);
        g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val-=EvalKineticIntegral(g3,g2);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=EvalKineticIntegral(g3,g2);
        return val;
    }
    else if (g1.AngMomV[0]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[1]=g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val=EvalKineticIntegral(g3,g2);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=EvalKineticIntegral(g3,g2);
        return val;
    }
    else
    {
        g3=g1;
        g3.AngMom=2;
        return EvalKineticIntegral(g3,g2);
    }
}

double fKinetic(GaussianOrbital g1,GaussianOrbital g2)
{
    GaussianOrbital g3;
    double val;

    g3=g1;
    g3.AngMom=3;
    switch(g1.AngMomV[2])
    {
    case -3:    g3.AngMomV[0]=2;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=0;
                val=3.0*EvalKineticIntegral(g3,g2);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=EvalKineticIntegral(g3,g2);
                break;
    case -2:    g3.AngMomV[0]=1;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=1;
                val=EvalKineticIntegral(g3,g2);
                break;
    case -1:    g3.AngMomV[0]=0;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=2;
                val=4.0*EvalKineticIntegral(g3,g2);
                g3.AngMomV[0]=2;
                g3.AngMomV[2]=0;
                val-=EvalKineticIntegral(g3,g2);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=EvalKineticIntegral(g3,g2);
                break;
    case 0:     g3.AngMomV[0]=0;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=3;
                val=2.0*EvalKineticIntegral(g3,g2);
                g3.AngMomV[2]=1;
                g3.AngMomV[0]=2;
                val-=3.0*EvalKineticIntegral(g3,g2);
                g3.AngMomV[1]=2;
                g3.AngMomV[0]=0;
                val-=3.0*EvalKineticIntegral(g3,g2);
                break;
    case 1: g3.AngMomV[1]=0;
                g3.AngMomV[0]=1;
                g3.AngMomV[2]=2;
                val=4.0*EvalKineticIntegral(g3,g2);
                g3.AngMomV[1]=2;
                g3.AngMomV[2]=0;
                val-=EvalKineticIntegral(g3,g2);
                g3.AngMomV[1]=0;
                g3.AngMomV[0]=3;
                val-=EvalKineticIntegral(g3,g2);
                break;
    case 2:     g3.AngMomV[0]=2;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=1;
                val=EvalKineticIntegral(g3,g2);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=2;
                val-=EvalKineticIntegral(g3,g2);
                break;
    case 3:     g3.AngMomV[0]=3;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=0;
                val=EvalKineticIntegral(g3,g2);
                g3.AngMomV[0]=1;
                g3.AngMomV[1]=2;
                val-=3.0*EvalKineticIntegral(g3,g2);
                break;
    default:    break;
    }
    return val;
}

double ddKinetic(GaussianOrbital g1,GaussianOrbital g2)
{
    GaussianOrbital g3;
    double val;

    if(g1.AngMomV[2]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[0]=g3.AngMomV[1]=0;
        g3.AngMomV[2]=2;
        val=2.0*dKinetic(g2,g3);
        g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val-=dKinetic(g2,g3);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=dKinetic(g2,g3);
        return val;
    }
    else if (g1.AngMomV[0]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[1]=g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val=dKinetic(g2,g3);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=dKinetic(g2,g3);
        return val;
    }
    else
    {
        g3=g1;
        g3.AngMom=2;
        return dKinetic(g2,g3);
    }
}

double fdKinetic(GaussianOrbital g1,GaussianOrbital g2)
{
    GaussianOrbital g3;
    double val;

    if(g2.AngMomV[2]==-2)
    {
        g3=g2;
        g3.AngMom=2;
        g3.AngMomV[0]=g3.AngMomV[1]=0;
        g3.AngMomV[2]=2;
        val=2.0*fKinetic(g1,g3);
        g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val-=fKinetic(g1,g3);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=fKinetic(g1,g3);
        return val;
    }
    else if (g2.AngMomV[0]==-2)
    {
        g3=g2;
        g3.AngMom=2;
        g3.AngMomV[1]=g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val=fKinetic(g1,g3);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=fKinetic(g1,g3);
        return val;
    }
    else
    {
        g3=g2;
        g3.AngMom=2;
        return fKinetic(g1,g3);
    }
}

double ffKinetic(GaussianOrbital g1,GaussianOrbital g2)
{
    GaussianOrbital g3;
    double val;

    g3=g2;
    g3.AngMom=3;
    switch(g2.AngMomV[2])
    {
    case -3:    g3.AngMomV[0]=2;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=0;
                val=3.0*fKinetic(g1,g3);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=fKinetic(g1,g3);
                break;
    case -2:    g3.AngMomV[0]=1;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=1;
                val=fKinetic(g1,g3);
                break;
    case -1:    g3.AngMomV[0]=0;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=2;
                val=4.0*fKinetic(g1,g3);
                g3.AngMomV[0]=2;
                g3.AngMomV[2]=0;
                val-=fKinetic(g1,g3);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=fKinetic(g1,g3);
                break;
    case 0:     g3.AngMomV[0]=0;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=3;
                val=2.0*fKinetic(g1,g3);
                g3.AngMomV[2]=1;
                g3.AngMomV[0]=2;
                val-=3.0*fKinetic(g1,g3);
                g3.AngMomV[1]=2;
                g3.AngMomV[0]=0;
                val-=3.0*fKinetic(g1,g3);
                break;
    case 1: g3.AngMomV[1]=0;
                g3.AngMomV[0]=1;
                g3.AngMomV[2]=2;
                val=4.0*fKinetic(g1,g3);
                g3.AngMomV[1]=2;
                g3.AngMomV[2]=0;
                val-=fKinetic(g1,g3);
                g3.AngMomV[1]=0;
                g3.AngMomV[0]=3;
                val-=fKinetic(g1,g3);
                break;
    case 2:     g3.AngMomV[0]=2;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=1;
                val=fKinetic(g1,g3);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=2;
                val-=fKinetic(g1,g3);
                break;
    case 3:     g3.AngMomV[0]=3;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=0;
                val=fKinetic(g1,g3);
                g3.AngMomV[0]=1;
                g3.AngMomV[1]=2;
                val-=3.0*fKinetic(g1,g3);
                break;
    default:    break;
    }
    return val;
}

double KineticIntegral(GaussianOrbital g1,GaussianOrbital g2)
{
    if((g1.AngMom>=0)&&(g2.AngMom>=0))
    {
        return EvalKineticIntegral(g1,g2);
    }
    else
    {
        if(g1.AngMom==-2)
        {   //g1 is a pure d function
            if(g2.AngMom==-2)
            {
                return ddKinetic(g1,g2);
            }
            else if(g2.AngMom==-3)
            {
                return fdKinetic(g2,g1);
            }
            else
            {
                return dKinetic(g1,g2);
            }
        }
        else //if (g1.AngMom==-3)
        {
            if(g2.AngMom==-2)
            {
                return fdKinetic(g1,g2);
            }
            else if (g2.AngMom==-3)
            {
                return ffKinetic(g1,g2);
            }
            else
            {
                return fKinetic(g1,g2);
            }
        }
    }
}

double FF(int k,double z)
{
    int i;
    double dz=0.0001;
    double zl;
    double zu;
    double vl,vu;
    double rval=0.0;
    static double zmax=0.0;
    static double zmin=1.0e10;
    static int kmax=0;

/*/ printf("FF called parameter %i argument %g\n",k,z);
    if(z<zmin)
        zmin=z;
    if(z>zmax)
        zmax=z;
    if(k>kmax)
        kmax=k;
    printf("Max:%g Min:%g MaxK:%i\n",zmin,zmax,k);*/
    if(z==0.0)
    {
        if(k>0)
            return 1.0/(2.0*(double)k+1.0);
        else
            return 1.0;
    }
    for(i=0;i<10000;i++)
    {
        zl=(double)i*dz;
        zu=zl+dz;
        vl=pow(zl,2*k)*exp(-z*zl*zl);
        vu=pow(zu,2*k)*exp(-z*zu*zu);
        rval+=0.5*(vl+vu);
    }
    return rval*dz;
}

double NAIss(double zeta,double xi,double rr,double U,int m)
{
    return 2.0*sqrt(zeta/pi)*OIss(zeta,xi,rr)*FF(m,U);
}

double NAIps(double zeta,double xi,double rr,double U,double dx,double dUx,int m)
{
    double rval;

    rval=dx*NAIss(zeta,xi,rr,U,m)-dUx*NAIss(zeta,xi,rr,U,m+1);
    return rval;
}

double NAIpp(double zeta,double xi,double rr,double U,double dx,double dy,double dUx,double dUy,int i,int j,int m)
{
    double rval;

    rval=dy*NAIps(zeta,xi,rr,U,dx,dUx,m)-dUy*NAIps(zeta,xi,rr,U,dx,dUx,m+1);
    if(i==j)
        rval+=(NAIss(zeta,xi,rr,U,m)-NAIss(zeta,xi,rr,U,m+1))*0.5/zeta;
    return rval;
}

double NAIds(double zeta,double xi,double rr,double U,double dx,double dy,double dUx,double dUy,int i,int j,int m)
{
    double rval;

    rval=dy*NAIps(zeta,xi,rr,U,dx,dUx,m)-dUy*NAIps(zeta,xi,rr,U,dx,dUx,m+1);
    if(i==j)
        rval+=(NAIss(zeta,xi,rr,U,m)-NAIss(zeta,xi,rr,U,m+1))*0.5/zeta;
    return rval;
}

double NAIdp(double zeta,double xi,double rr,double U,double dx,double dy,double dz,double dUx,double dUy,double dUz,int i,int j,int k,int m)
{
    double rval;

    rval=dz*NAIds(zeta,xi,rr,U,dx,dy,dUx,dUy,i,j,m)-dUz*NAIds(zeta,xi,rr,U,dx,dy,dUx,dUy,i,j,m+1);
    if(i==k)
        rval+=(NAIps(zeta,xi,rr,U,dy,dUy,m)-NAIps(zeta,xi,rr,U,dy,dUy,m))*0.5/zeta;
    if(j==k)
        rval+=(NAIps(zeta,xi,rr,U,dx,dUx,m)-NAIps(zeta,xi,rr,U,dx,dUx,m))*0.5/zeta;
    return rval;
}

double NAIdd(double zeta,double xi,double rr,double U,double dx,double dy,double dz,double dw,double dUx,double dUy,double dUz,double dUw,int i,int j,int k,int l,int m)
{
    double rval;

    rval=dw*NAIdp(zeta,xi,rr,U,dx,dy,dz,dUx,dUy,dUz,i,j,k,m)-dUw*NAIdp(zeta,xi,rr,U,dx,dy,dz,dUx,dUy,dUz,i,j,k,m+1);
    if(i==l)
        rval+=(NAIpp(zeta,xi,rr,U,dy,dz,dUy,dUz,j,k,m)-NAIpp(zeta,xi,rr,U,dy,dz,dUy,dUz,j,k,m+1))*0.5/zeta;
    if(j==l)
        rval+=(NAIpp(zeta,xi,rr,U,dx,dz,dUx,dUz,i,k,m)-NAIpp(zeta,xi,rr,U,dx,dz,dUx,dUz,i,k,m+1))*0.5/zeta;
    if(k==l)
        rval+=(NAIds(zeta,xi,rr,U,dx,dy,dUx,dUy,i,j,m)-NAIds(zeta,xi,rr,U,dx,dy,dUx,dUy,i,j,m+1))*0.5/zeta;
    return rval;
}

double EvalNuclearIntegral(GaussianOrbital g1,GaussianOrbital g2,AtomicCore A)
{
    int i,j,k,l,c1;
    int val=10*g1.AngMom+g2.AngMom;
    double zeta,xi;
    double rr,P[3],U;
    double dx,dy,dz,dw;
    double dUx,dUy,dUz,dUw;
    double za,zb;

    zeta=g1.OrbExponent+g2.OrbExponent;
    za=g1.OrbExponent;
    zb=g2.OrbExponent;
    xi=g1.OrbExponent*g2.OrbExponent/zeta;
    rr=0.0;
    U=0.0;
    for(c1=0;c1<3;c1++)
    {
        rr+=(g1.Center[c1]-g2.Center[c1])*(g1.Center[c1]-g2.Center[c1]);
        P[c1]=(g1.OrbExponent*g1.Center[c1]+g2.OrbExponent*g2.Center[c1])/zeta;
        U+=(P[c1]-A.Position[c1])*(P[c1]-A.Position[c1]);
    }
    U*=zeta;
    switch(val)
    {
    case 0: return NAIss(zeta,xi,rr,U,0);
            break;
    case 1: for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    i=c1;
            }
            dx=P[i]-g2.Center[i];
            dUx=P[i]-A.Position[i];
            return NAIps(zeta,xi,rr,U,dx,dUx,0);
            break;
    case 2: j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            dUx=P[i]-A.Position[i];
            dUy=P[j]-A.Position[j];
            return NAIds(zeta,xi,rr,U,dx,dy,dUx,dUy,i,j,0);
            break;
    case 10:for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    i=c1;
            }
            dx=P[i]-g1.Center[i];
            dUx=P[i]-A.Position[i];
            return NAIps(zeta,xi,rr,U,dx,dUx,0);
            break;
    case 11:for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    i=c1;
            }
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    j=c1;
            }
            dx=P[i]-g1.Center[i];
            dUx=P[i]-A.Position[i];
            dy=P[j]-g2.Center[j];
            dUy=P[j]-A.Position[j];
            return NAIpp(zeta,xi,rr,U,dx,dy,dUx,dUy,i,j,0);
            break;
    case 12:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                    k=c1;
            }
            dx=P[i]-g2.Center[i];
            dy=P[j]-g2.Center[j];
            dz=P[k]-g1.Center[k];
            dUx=P[i]-A.Position[i];
            dUy=P[j]-A.Position[j];
            dUz=P[k]-A.Position[k];
            return NAIdp(zeta,xi,rr,U,dx,dy,dz,dUx,dUy,dUz,i,j,k,0);
            break;
    case 20:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dUx=P[i]-A.Position[i];
            dUy=P[j]-A.Position[j];
            return NAIds(zeta,xi,rr,U,dx,dy,dUx,dUy,i,j,0);
            break;
    case 21:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                    k=c1;
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g2.Center[k];
            dUx=P[i]-A.Position[i];
            dUy=P[j]-A.Position[j];
            dUz=P[k]-A.Position[k];
            return NAIdp(zeta,xi,rr,U,dx,dy,dz,dUx,dUy,dUz,i,j,k,0);
            break;
    case 22:j=5;
            for(c1=0;c1<3;c1++)
            {
                if(g1.AngMomV[c1]==1)
                {
                    if(j==5)
                        j=c1;
                    else
                        i=c1;
                }
                if(g1.AngMomV[c1]==2)
                {
                    i=c1;
                    j=c1;
                }
            }
            l=5;
            for(c1=0;c1<3;c1++)
            {
                if(g2.AngMomV[c1]==1)
                {
                    if(l==5)
                        l=c1;
                    else
                        k=c1;
                }
                if(g2.AngMomV[c1]==2)
                {
                    k=c1;
                    l=c1;
                }
            }
            dx=P[i]-g1.Center[i];
            dy=P[j]-g1.Center[j];
            dz=P[k]-g2.Center[k];
            dw=P[l]-g2.Center[l];
            dUx=P[i]-A.Position[i];
            dUy=P[j]-A.Position[j];
            dUz=P[k]-A.Position[k];
            dUw=P[l]-A.Position[l];
            return NAIdd(zeta,xi,rr,U,dx,dy,dz,dw,dUx,dUy,dUz,dUw,i,j,k,l,0);
            break;
    default:return 0.0;
    }
}

double dNuclear(GaussianOrbital g1,GaussianOrbital g2,AtomicCore A)
{
    double val;
    GaussianOrbital g3;

    if(g1.AngMomV[2]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[0]=g3.AngMomV[1]=0;
        g3.AngMomV[2]=2;
        val=2.0*EvalNuclearIntegral(g3,g2,A);
        g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val-=EvalNuclearIntegral(g3,g2,A);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=EvalNuclearIntegral(g3,g2,A);
        return val;
    }
    else if (g1.AngMomV[0]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[1]=g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val=EvalNuclearIntegral(g3,g2,A);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=EvalNuclearIntegral(g3,g2,A);
        return val;
    }
    else
    {
        g3=g1;
        g3.AngMom=2;
        return EvalNuclearIntegral(g3,g2,A);
    }
}

double fNuclear(GaussianOrbital g1,GaussianOrbital g2,AtomicCore A)
{
    GaussianOrbital g3;
    double val;

    g3=g1;
    g3.AngMom=3;
    switch(g1.AngMomV[2])
    {
    case -3:    g3.AngMomV[0]=2;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=0;
                val=3.0*EvalNuclearIntegral(g3,g2,A);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=EvalNuclearIntegral(g3,g2,A);
                break;
    case -2:    g3.AngMomV[0]=1;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=1;
                val=EvalNuclearIntegral(g3,g2,A);
                break;
    case -1:    g3.AngMomV[0]=0;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=2;
                val=4.0*EvalNuclearIntegral(g3,g2,A);
                g3.AngMomV[0]=2;
                g3.AngMomV[2]=0;
                val-=EvalNuclearIntegral(g3,g2,A);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=EvalNuclearIntegral(g3,g2,A);
                break;
    case 0:     g3.AngMomV[0]=0;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=3;
                val=2.0*EvalNuclearIntegral(g3,g2,A);
                g3.AngMomV[2]=1;
                g3.AngMomV[0]=2;
                val-=3.0*EvalNuclearIntegral(g3,g2,A);
                g3.AngMomV[1]=2;
                g3.AngMomV[0]=0;
                val-=3.0*EvalNuclearIntegral(g3,g2,A);
                break;
    case 1: g3.AngMomV[1]=0;
                g3.AngMomV[0]=1;
                g3.AngMomV[2]=2;
                val=4.0*EvalNuclearIntegral(g3,g2,A);
                g3.AngMomV[1]=2;
                g3.AngMomV[2]=0;
                val-=EvalNuclearIntegral(g3,g2,A);
                g3.AngMomV[1]=0;
                g3.AngMomV[0]=3;
                val-=EvalNuclearIntegral(g3,g2,A);
                break;
    case 2:     g3.AngMomV[0]=2;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=1;
                val=EvalNuclearIntegral(g3,g2,A);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=2;
                val-=EvalNuclearIntegral(g3,g2,A);
                break;
    case 3:     g3.AngMomV[0]=3;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=0;
                val=EvalNuclearIntegral(g3,g2,A);
                g3.AngMomV[0]=1;
                g3.AngMomV[1]=2;
                val-=3.0*EvalNuclearIntegral(g3,g2,A);
                break;
    default:    break;
    }
    return val;
}

double ddNuclear(GaussianOrbital g1,GaussianOrbital g2,AtomicCore A)
{
    GaussianOrbital g3;
    double val;

    if(g1.AngMomV[2]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[0]=g3.AngMomV[1]=0;
        g3.AngMomV[2]=2;
        val=2.0*dNuclear(g2,g3,A);
        g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val-=dNuclear(g2,g3,A);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=dNuclear(g2,g3,A);
        return val;
    }
    else if (g1.AngMomV[0]==-2)
    {
        g3=g1;
        g3.AngMom=2;
        g3.AngMomV[1]=g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val=dNuclear(g2,g3,A);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=dNuclear(g2,g3,A);
        return val;
    }
    else
    {
        g3=g1;
        g3.AngMom=2;
        return dNuclear(g2,g3,A);
    }
}

double fdNuclear(GaussianOrbital g1,GaussianOrbital g2,AtomicCore A)
{
    GaussianOrbital g3;
    double val;

    if(g2.AngMomV[2]==-2)
    {
        g3=g2;
        g3.AngMom=2;
        g3.AngMomV[0]=g3.AngMomV[1]=0;
        g3.AngMomV[2]=2;
        val=2.0*fNuclear(g1,g3,A);
        g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val-=fNuclear(g1,g3,A);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=fNuclear(g1,g3,A);
        return val;
    }
    else if (g2.AngMomV[0]==-2)
    {
        g3=g2;
        g3.AngMom=2;
        g3.AngMomV[1]=g3.AngMomV[2]=0;
        g3.AngMomV[0]=2;
        val=fNuclear(g1,g3,A);
        g3.AngMomV[0]=0;
        g3.AngMomV[1]=2;
        val-=fNuclear(g1,g3,A);
        return val;
    }
    else
    {
        g3=g2;
        g3.AngMom=2;
        return fNuclear(g1,g3,A);
    }
}

double ffNuclear(GaussianOrbital g1,GaussianOrbital g2,AtomicCore A)
{
    GaussianOrbital g3;
    double val;

    g3=g2;
    g3.AngMom=3;
    switch(g2.AngMomV[2])
    {
    case -3:    g3.AngMomV[0]=2;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=0;
                val=3.0*fNuclear(g1,g3,A);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=fNuclear(g1,g3,A);
                break;
    case -2:    g3.AngMomV[0]=1;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=1;
                val=fNuclear(g1,g3,A);
                break;
    case -1:    g3.AngMomV[0]=0;
                g3.AngMomV[1]=1;
                g3.AngMomV[2]=2;
                val=4.0*fNuclear(g1,g3,A);
                g3.AngMomV[0]=2;
                g3.AngMomV[2]=0;
                val-=fNuclear(g1,g3,A);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=3;
                val-=fNuclear(g1,g3,A);
                break;
    case 0:     g3.AngMomV[0]=0;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=3;
                val=2.0*fNuclear(g1,g3,A);
                g3.AngMomV[2]=1;
                g3.AngMomV[0]=2;
                val-=3.0*fNuclear(g1,g3,A);
                g3.AngMomV[1]=2;
                g3.AngMomV[0]=0;
                val-=3.0*fNuclear(g1,g3,A);
                break;
    case 1: g3.AngMomV[1]=0;
                g3.AngMomV[0]=1;
                g3.AngMomV[2]=2;
                val=4.0*fNuclear(g1,g3,A);
                g3.AngMomV[1]=2;
                g3.AngMomV[2]=0;
                val-=fNuclear(g1,g3,A);
                g3.AngMomV[1]=0;
                g3.AngMomV[0]=3;
                val-=fNuclear(g1,g3,A);
                break;
    case 2:     g3.AngMomV[0]=2;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=1;
                val=fNuclear(g1,g3,A);
                g3.AngMomV[0]=0;
                g3.AngMomV[1]=2;
                val-=fNuclear(g1,g3,A);
                break;
    case 3:     g3.AngMomV[0]=3;
                g3.AngMomV[1]=0;
                g3.AngMomV[2]=0;
                val=fNuclear(g1,g3,A);
                g3.AngMomV[0]=1;
                g3.AngMomV[1]=2;
                val-=3.0*fNuclear(g1,g3,A);
                break;
    default:    break;
    }
    return val;
}

double NuclearIntegral(GaussianOrbital g1,GaussianOrbital g2,AtomicCore A)
{
    if((g1.AngMom>=0)&&(g2.AngMom>=0))
    {
        return EvalNuclearIntegral(g1,g2,A);
    }
    else
    {
        if(g1.AngMom==-2)
        {   //g1 is a pure d function
            if(g2.AngMom==-2)
            {
                return ddNuclear(g1,g2,A);
            }
            else if(g2.AngMom==-3)
            {
                return fdNuclear(g2,g1,A);
            }
            else
            {
                return dNuclear(g1,g2,A);
            }
        }
        else //if (g1.AngMom==-3)
        {
            if(g2.AngMom==-2)
            {
                return fdNuclear(g1,g2,A);
            }
            else if (g2.AngMom==-3)
            {
                return ffNuclear(g1,g2,A);
            }
            else
            {
                return fNuclear(g1,g2,A);
            }
        }
    }
}

}

} // namespace Fchk
} // namespace Plugin
