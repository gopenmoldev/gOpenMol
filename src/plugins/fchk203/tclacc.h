/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaussians.h"
#include "impfchk.h"

#ifndef tclacc_included
namespace Plugin {
namespace Fchk {
    

double EvalPsi(int,int,double*);
double EvalGrad(int,int,double*);
double CalcPlane(int,const char**);

} // namespace Fchk
} // namespace Plugin
#define tclacc_included 1
#endif
