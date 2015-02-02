#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef BunchKaufmanIncluded

namespace Plugin {
namespace Surfarea {
    
void SISolver(int*,double*,double*,double*,double*,double*,int);
void BunchKaufman(double **Mat,int Order,int Inertia[],double *d=NULL,double *e=NULL,int *P=NULL,double *L=NULL,int Mag=9);
void BunchKaufman(double *Mat,int Order,int Inertia[],double *d=NULL,double *e=NULL,int *P=NULL,double *L=NULL,int Mag=9);

} // namespace Surfarea
} // namespace Plugin

#define BunchKaufmanIncluded 1
#endif
