#ifndef RannosLoaded
#include <stdio.h>

namespace Plugin {
namespace Fchk {
    
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define MPI 3.14159265359
#define SPI 6.0*MPI
#define SP2 1.253314137316
#define ONET 0.083333333333333
#define MAXBIT 30
#define MAXDIM 6

#define IMIN(a,b) (a<b ? a:b)


static long rangen;  //seed value for the ranx function  Define outside program to give program scope
static long nrans;


void raninit(long seed);
long getranseed();
long numrans();
void setranseed(long nseed);
double ranx();          //ran2 function from numerical recipes
void rmarin(int IJ,int KL);
double ranmar(void);
void outrmarset(FILE *of);
void getrmarset(FILE* inf,int nt);
double granx(double sigma);
void funranx(double(*Fn)(int,double*),double Fnmax,int ndim,double* xmax,double* x);
void sobseq(int *n,double x[]);

} // namespace Fchk
} // namespace Plugin

#define RannosLoaded 1
#endif
