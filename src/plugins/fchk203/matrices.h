/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "keys.h"

namespace Plugin {
namespace Fchk {
    
void MatMul(double*,double *,int,int,int,int,double*);
void MatPr(double*,int,int);
void MatTr(double*,double*,int);
int BiDiag(double*,int,int,double *U=NULL,double *V=NULL);
void Givens(double,double,double[]);
void GolubKahan(double *,int,int,double *U=NULL,double *V=NULL);
int Cholesky(double*,int,double *G=NULL);
void MatSqrtSVD(double*,double*,int);
void MatSqrtS(double*,double*,int);

} // namespace Fchk
} // namespace Plugin
