/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

namespace Plugin {
namespace Fchk {
    
#ifndef gaussians_included

class GaussianOrbital
{
public:
    int AngMom;
    int AngMomV[3];
    double NormCoeff;
    double OrbExponent;
    double Center[3];
    GaussianOrbital operator=(GaussianOrbital);
    double Evaluate(double,double,double);
    double Evaluate(double*);
    double EvaluateGradient(double,double,double,double*);
    double EvaluateGradient(double*,double*);
    double EvaluateLaplacian(double,double,double);
    double EvaluateLaplacian(double*);
};

class ContractedShell
{
public:
    int NumberOfPrimitives;
    int type;
    double* PrimitiveExponents;
    double* ContractionCoeffs;  
    double Center[3];
};

class GaussAtom
{
public:
    int AtNum;
    char Element[3];
    double Position[3];
    void AssignType();
};


class AtomicCore
{
public:
    double NuclCharge;
    double Position[3];
};

namespace newrecurse {
int OverlapIntegral(GaussianOrbital*,GaussianOrbital*,double*);
int KineticIntegral(GaussianOrbital*,GaussianOrbital*,double*);
int NuclearIntegral(GaussianOrbital*,GaussianOrbital*,AtomicCore*,double*);
void InitMatrices();
void ClearMatrices();
}
namespace oldrecurse {
double OverlapIntegral(GaussianOrbital,GaussianOrbital);
double KineticIntegral(GaussianOrbital,GaussianOrbital);
double NuclearIntegral(GaussianOrbital,GaussianOrbital,AtomicCore);
}
#define gaussians_included 1
#endif

 
} // namespace Fchk
} // namespace Plugin
