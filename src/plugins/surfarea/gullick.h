#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* Copyright 2002-2003 Kevin J. Boyd and the University of New Orleans.  
Permission granted to distribute and modify this code for personal, educational, and
research use. */


void MatPr(double *Ma,int ma,int na);
void MatTr(double *Ma,double *Mat,int ma);
void MatVec(double* Ma,double *V,int N,int M,double *RVec);
int RPR(double *R,double *mu,double *qvec,int* pcol,int m,int n,double *Q=NULL);
int RPRSolve(double *A,double *x,double *b,double *mu,int* pcol,double *qvec,double *lambda,int m,int n,int p);
void LUDecomp(double*,int,int*);
void LUSolve(double*,double*,double*,int,int*);
void Regularize(double*,double*,double*,int,int,double h=0.1);
