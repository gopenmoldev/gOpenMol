#include <stdio.h>
#include <stdlib.h>
#include "gaussians.h"
#include "keys.h"
#include "rannos.hpp"

#ifndef Fchkrdr_Included

namespace Plugin {
namespace Fchk {
    
class FitLaguerre
{
public:
    int NPoly;
    double coeffs[25];
    double rscale;
    double maxerror;
    double averror;
    double Eval(double x);
};


double SHConstrain(double *Parms);
void dSHConstrain(double *x,double *Parms,double *dFdP);
double SHFit(double *x,double *Parms);
void dSHFit(double *x,double *Parms,double *dFdP);
void GaussBackSubst(double **M,double *x,double *b,int *p,int N);
void MSolve(double **M,double *x,double *b,int Order);
double levmar(double **X,double *y,double *weight,int NPts,
            double *Parms,int *vParms,int NParms,
            double (*Func)(double *x,double *Parms),
            void (*dFunc)(double *x,double *Parms,double *dFdP),
            double (*Constraints)(double* Parms)=NULL,
            void (*dConstraints)(double*,double*,double*)=NULL);
double gammaln(double xx);
double aslg(int l,int m,double x);
double fact(int N);
double SphHarm(int l,int m,double cth,double phi);
double EvalLaguerre(int order,double alpha,double xval,double* dlag=NULL);
void SetGauLagPar(double *a,double *w,int order,double alpha);
double func(double cl);
double LagFit(int NMin,int NMax,double aver,double maxer,double* xx,double* yy,int numpts,double xscal,FitLaguerre* FL=NULL);
FitLaguerre DetermineFit(int NP,double* r,double* iv,double tol1,double tol2);
char ExtractLetter(char* orb);
void Rotate(int axis,double* vec);
void Reflect(int ax1,int ax2,double* vec);
void Invert(double* vec);
void Reflect45(int ax1,int ax2,double* vec);
int GetSymmetryDesignation(int orbid);
void InsertID(char* orb,char* msg);
void GetProperpDesignation(char* orb,int orbid);
void GetProperdDesignation(char* orb,int orbid);
void WritePlot(char* fn,int orbid);

} // namespace Fchk
} // namespace Plugin

#define Fchkrdr_Included 1
#endif
