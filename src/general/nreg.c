/*
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     PROGRAM FOR POLYNOMIAL REGRESSION
C
      DIMENSION X(250),Y(250),A(0:10,0:11)
C
C     READ INPUT DATA
C
      WRITE (6,*) 'GIVE THE END OF DATA INDICATORS XEND AND YEND'
      READ (5,*) XEND,YEND
      N = 1
C
      WRITE (6,*) 'GIVE THE DATA POINTS X(N) AND Y(N)'
    1 WRITE (6,100) N
  100 FORMAT (5X,'N=',I3,' X(N),Y(N)=  '$)
      READ (5,*) X(N),Y(N)
      IF (X(N).EQ.XEND .AND. Y(N).EQ.YEND)   GO TO 2
      N = N + 1
      GO TO 1
C
    2 N = N - 1
C
C     ASK DEGREE OF THE POLYNOMIAL
C
      WRITE (6,101)
  101 FORMAT (5X,'GIVE DEGREE OF THE POLYNOMIAL  '$)
      READ (5,*) M
C
C     FORM COEFFICIENT MATRIX
C
      DO 10 J=0,M
      DO 10 K=0,M+1
      A(J,K) = 0.0
   10 CONTINUE
C
      DO 40 I=1,N
C
      A(0,M+1) = A(0,M+1) + Y(I)
      A(0,0) = A(0,0) + 1
      DO 15 K=1,M
      A(0,K) = A(0,K) + X(I)**K
   15 CONTINUE
C
      DO 30 J=1,M
      A(J,M+1) = A(J,M+1) + Y(I)*X(I)**J
      DO 20 K=0,M
      A(J,K) = A(J,K) + X(I)**(J+K)
   20 CONTINUE
   30 CONTINUE
   40 CONTINUE
C
C     DIAGONALIZE
C
      CALL GAUSS (A,M)
C
C     PRINT RESULT
C
      WRITE (6,*) 'THE COEFFICIENTS ARE:'
      WRITE (6,*) (A(I,M+1),I=0,M)
C      R = -A(1,3)/A(2,3)/2.0
C      WRITE (6,*) 'R = ',R
C      E = A(0,3)+A(1,3)*R+A(2,3)*R*R
C      WRITE (6,*) 'E = ',E
C
      WRITE (6,*) 'THE OBSERVED CALCULATED AND DIFFERENCE'
      DO 60 I=1,N
      CALC = A(0,M+1)
      DO 50 J=1,M
      CALC = CALC + A(J,M+1)*X(I)**J
   50 CONTINUE
      WRITE (6,102) X(I),Y(I),CALC,(Y(I)-CALC)
  102 FORMAT (5X,4F12.5)
   60 CONTINUE
C
      STOP
      END
      SUBROUTINE GAUSS (A,M)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     GAUSS ELIMINATION PROCEDURE
C
      DIMENSION A(0:10,0:11)
C
      DO 10 I=0,M
      PIVOT = A(I,I)
      IF (ABS(PIVOT) .LE. 1.0E-10)   GO TO 9
C
      DO 20 K=0,M+1
      A(I,K) = A(I,K)/PIVOT
   20 CONTINUE
C
      DO 30 J=0,M
      IF (J .EQ. I)   GO TO 30
      AJI = A(J,I)
      DO 40 K=0,M+1
      A(J,K) = A(J,K) - AJI*A(I,K)
   40 CONTINUE
   30 CONTINUE
C
   10 CONTINUE
C
      RETURN
    9 WRITE (6,*) 'SINGULARITY IN GAUSS'
      STOP
      END
*/

/*
C
C     PROGRAM FOR POLYNOMIAL REGRESSION
C
      Original version Matti hotokka
      modified to c by
      Leif Laaksonen   1989
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "math_oper.h"
/* #include "printmsg.h" */

#include "stdafx.h"

#define MAXdeg   10        /* maximum degree of polynom      */

#define Rabs(x)   ( (x) > 0.0 ? (x) : -(x) )

static void GaussElim(float a[][MAXdeg+1] , int m);

int gomp_nreg(const float *xobs, const float *yobs, int n, int m, float *RegCoeff)
/*
  const float *xobs;
  const float *yobs;
  int m;                 polynomial degree                  
  int n;                 number of observations             
  const float *RegCoeff;       regression coefficients           
*/
{
    float a[MAXdeg][MAXdeg+1];
/*  float ssqr;*/

/*  float calc=0.0;*/
    int i,j,k;

/*  char OutText[BUFF_LEN];*/
/*
  C
  C     FORM COEFFICIENT MATRIX
  C
*/

    for(j = 0 ; j <= m ; j++) {
        for(k = 0 ; k <= m+1 ; k++ ) {
            a[j][k] = 0.0;
        }
    }

    for( i = 1 ; i <= n ; i++) {
        a[0][m+1] = a[0][m+1] + yobs[i];
        a[0][0] = a[0][0] + 1.;

        for(k = 1 ; k <= m ; k++) {
            a[0][k] = a[0][k] + pow(xobs[i],(float) k);
        }

        for(j = 1 ; j <= m ; j++) {
            a[j][m+1] = a[j][m+1] + yobs[i]*pow(xobs[i],(float) j);
            for(k = 0 ; k <= m ; k++) {
                a[j][k] = a[j][k] + pow(xobs[i],(float) (j+k));
            }
        }
    }

/*
  C
  C     DIAGONALIZE
  C
*/
    GaussElim(a , m);
/*
  C
  C     PRINT RESULT
  C
*/
    for(i = 0 ; i <= m ; i++) {
        RegCoeff[i] = a[i][m+1];
/*          sprintf(OutText,"Coefficient a%d : %f\n",i,a[i][m+1]);
            gomp_PrintMessage(OutText); */}
/*
  PrintMessage("\n x(i)      y(i)      calc    (y(i)-calc) ");

  ssqr=0.0;
  for( i = 1 ; i <= n ; i++) {
  calc=a[0][m+1];
  for( j = 1 ; j <= m ; j++) {
  calc = calc + a[j][m+1]*pow(xobs[i],(float) j); }
  ssqr+=pow((yobs[i]-calc),2);
  sprintf(OutText," %f %f %f %f \n",xobs[i],yobs[i],calc,(yobs[i]-calc));
  PrintMessage(OutText);
  }
  sprintf(OutText,"\nSum of squares : %f\n",ssqr);
  PrintMessage(OutText);
*/
    return(0);
}
void GaussElim (float a[][MAXdeg+1] , int m)
/*
  float a[][MAXdeg+1];
  int m;
*/
{
    float aji,pivot;
    int i,j,k;
/*
  C
  C     GAUSS ELIMINATION PROCEDURE
  C
*/

    for( i = 0 ; i <= m ; i++ ) {
        pivot = a[i][i];

        if ( Rabs(pivot) < 1.0e-10)  {
            printf(" Singularity in Gauss \n");
            exit(1); }

        for( k = 0 ; k <= m+1 ; k++ ) {
            a[i][k] = a[i][k]/pivot; }

        for( j = 0 ; j <= m ; j++ ) {
            if( j == i ) continue;
            aji = a[j][i];
            for( k = 0 ; k <= m+1 ; k++ ) {
                a[j][k] = a[j][k] - aji*a[i][k]; }
        }
    }
}

