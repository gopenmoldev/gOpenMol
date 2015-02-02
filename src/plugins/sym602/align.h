#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tba/complex.h"
#include "tba/matrix.h"

namespace Plugin {
namespace Symmetry {

void PrintMatrix(complex**,int,char*);
void MatMul(complex**,complex**,int,complex**);
void MatHMul(complex**,complex**,int,complex**);
void MatMulH(complex**,complex**,int,complex**);
int BiDiag(double**,int,int,double **U=NULL,double **V=NULL);
void Givens(complex,complex,complex[]);
void GolubKahan(double **,int,int,double **U=NULL,double **V=NULL);
int SVD(double**,int,int,double **U=NULL,double **V=NULL);
void Align(double**,double**,int,double*,double*,double[]);
void UpperTriang(complex**,int,complex **Q=NULL);
void TriDiag(complex**,int,complex **Q=NULL);
void TQRImp(complex**,int,complex **Q=NULL);

} // namespace Symmetry
} // namespace Plugin
