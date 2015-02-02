#ifndef Matrix_Included
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "complex.h"

namespace Plugin {
namespace Symmetry {

double dotp(double*,double*,int);
complex dotp(complex*,complex*,int);
void tred2(double **a,int n,double *d,double *e);
void tred2ne(double **a,int n,double *d,double *e);
int tqli(double *d,double *e,int n,double **z);
int tqline(double *d,double *e,int n,double **z);
int choldc(double **a,int n,double p[]);
int ccholdc(complex **a,int n,double p[]);
int ccholsl(complex **a,int n,double p[],complex b[],complex x[]);
void geneigset(double **a,double **b,double **aeff,double **lt,int n);
void reseteigs(double **evc,double **lt,int n);
void eigsrt(double d[],double **v,int n);
void aeigsrt(double d[],double **v,int n);
void ceigsrt(double d[],complex **v,int n);
void veigsrt(double d[],double **v,int n,complex **w,complex *r,int o);
void SortCRots(double**,double*,int);

} // namespace Symmetry
} // namespace Plugin
#define Matrix_Included 1
#endif
