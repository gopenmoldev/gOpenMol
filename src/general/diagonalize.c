#include "maindefs.h"
 
#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>

#include "math_oper.h"

#include "stdafx.h"

/*
  main()
  {

  double a[4] = {1.0 , 0.0 , 0.0 , 1.0};
  double b[4] = {1.0 , 0.0 , 0.0 , 1.0};
  double x[4];

  double eigv[2];
  double d[2];
  int n = 2;
  double rtol = 10.e-12;
  int i;

  gomp_jacobi(a,b,x,eigv,d,n,rtol);

  for(i = 0 ; i < 2 ; i++) {
  printf("%f (%f %f %f %f)\n",eigv[i],x[0],x[1],x[2],x[3]);
  printf("%f (%f %f %f %f)\n",eigv[i],a[0],a[1],a[2],a[3]);
  printf("%f (%f %f %f %f)\n",eigv[i],b[0],b[1],b[2],b[3]);
  }
  }

*/

/************************************************************************
  P R O G R A M 
  To solve the generalized eigenproblem using the
  generalized Jacobi iteration

  INPUT

  a(n,n)    = Stiffness matrix
  b(n,n)    = Mass matrix
  x(n,n)    = Matrix storing eigenvectors on solution exit
  eigv(n)   = Vector storing eigenvalues on solution exit
  d(n)      = Working vector
  n         = Order of the matrices a and b
  rtol      = Convergence tolerance (usually set to 10.**-12)
  nsmax     = Maximum number of sweeps allowed (usually set to 15)

  OUTPUT

  a(n,n)    = Diagonalized stiffness matrix
  b(n,n)    = Diagonalized mass matrix
  x(n,n)    = Eigenvectors stored columnwise
  eigv(n)   = Eigenvalues

*************************************************************************/

int gomp_jprJacobi(double  *a,double  *b,double  *x,double  *eigv,double  *d,int n,double rtol)
{
    register int i,j,k,ii,jj;
    int    nsmax=50,        /* Max number of sweeps */
        nsweep,          /* Current sweep number */
        nr,           
        jp1,jm1,kp1,km1,
        convergence;

    double eps,
        eptola,eptolb,
        akk,ajj,ab,
        check,sqch,
        d1,d2,den,
        ca,cg,
        aj,bj,ak,bk,
        xj,xk,
        tol,dif,
        epsa,epsb,
        bb;              /* Scale mass matrix */

/************************************************************************
  Initialize eigenvalue and eigenvector matrices 
*************************************************************************/
    for( i=0 ; i<n ; i++)
    {
        ii = i*n+i;            /* address of diagonal element */
/*    if (a[ii]<=0.0 || b[ii]<=0.0) return 0;*/
        eigv[i] = d[i] = a[ii]/b[ii];
        x[ii] = 1.0;           /* make an unit matrix */
    }
    if(n==1) return 1;      /* Return if single degree of freedom system */

/************************************************************************
  Initialize sweep counter and begin iteration
*************************************************************************/
    nr=n - 1;
    for ( nsweep = 0; nsweep < nsmax; nsweep++)
    {
/************************************************************************
    Check if present off-diagonal element is large enough to require zeroing
*************************************************************************/

        eps = pow(0.01, 2.0 * (nsweep + 1));
        for(j=0 ; j<nr ; j++)
        {
            jj=j+1;
            for(k=jj ; k<n ; k++)
            {
                eptola=( a[j*n+k]*a[j*n+k] / ( a[j*n+j]*a[k*n+k] ) );
                eptolb=( b[j*n+k]*b[j*n+k] / ( b[j*n+j]*b[k*n+k] ) );
                if ( eptola >= eps || eptolb >=eps )
                {

/************************************************************************
          if zeroing is required, calculate the rotation matrix elements
*************************************************************************/

                    akk=a[k*n+k]*b[j*n+k] - b[k*n+k]*a[j*n+k];
                    ajj=a[j*n+j]*b[j*n+k] - b[j*n+j]*a[j*n+k];
                    ab =a[j*n+j]*b[k*n+k] - a[k*n+k]*b[j*n+j];
                    check=(ab*ab+4.0*akk*ajj)/4.0;

                    if (check<=0.0)
                    {
                        printf("***Error   solution stop in *gomp_jacobi*\n");
                        printf("        check = %20.14e\n", check);
                        return 1;
                    }
                    sqch= sqrt(check);
                    d1  = ab/2.0+sqch;
                    d2  = ab/2.0-sqch;
                    den = d1;
                    if ( fabs(d2) > fabs(d1) ) den=d2;

                    if (den >= 0.0 && den <= 0.0)
                    {
                        ca=0.0;
                        cg= -a[j*n+k] / a[k*n+k];
                    }
                    else
                    {
                        ca= akk/den;
                        cg= -ajj/den;
                    }
/************************************************************************
  Perform the generalized rotation to zero the present off-diagonal element
*************************************************************************/
                    if (n != 2)
                    {
                        jp1=j+1;
                        jm1=j-1;
                        kp1=k+1;
                        km1=k-1;
/**************************************/
                        if ( jm1 >= 0)
                            for (i=0 ; i<=jm1 ; i++)
                            {
                                aj = a[i*n+j];
                                bj = b[i*n+j];
                                ak = a[i*n+k];
                                bk = b[i*n+k];
                                a[i*n+j] = aj+cg*ak;
                                b[i*n+j] = bj+cg*bk;
                                a[i*n+k] = ak+ca*aj;
                                b[i*n+k] = bk+ca*bj;
                            };
/**************************************/
                        if ((kp1-n+1) <= 0)
                            for (i=kp1 ; i<n ; i++)
                            {
                                aj = a[j*n+i];
                                bj = b[j*n+i];
                                ak = a[k*n+i];
                                bk = b[k*n+i];
                                a[j*n+i] = aj+cg*ak;
                                b[j*n+i] = bj+cg*bk;
                                a[k*n+i] = ak+ca*aj;
                                b[k*n+i] = bk+ca*bj;
                            };
/**************************************/
                        if ((jp1-km1) <= 0.0)
                            for (i=jp1 ; i<=km1 ; i++)
                            {
                                aj = a[j*n+i];
                                bj = b[j*n+i];
                                ak = a[i*n+k];
                                bk = b[i*n+k];
                                a[j*n+i] = aj+cg*ak;
                                b[j*n+i] = bj+cg*bk;
                                a[i*n+k] = ak+ca*aj;
                                b[i*n+k] = bk+ca*bj;
                            };
                    };
   
                    ak = a[k*n+k];
                    bk = b[k*n+k];
                    a[k*n+k] = ak+2.0*ca*a[j*n+k]+ca*ca*a[j*n+j];
                    b[k*n+k] = bk+2.0*ca*b[j*n+k]+ca*ca*b[j*n+j];
                    a[j*n+j] = a[j*n+j]+2.0*cg*a[j*n+k]+cg*cg*ak;
                    b[j*n+j] = b[j*n+j]+2.0*cg*b[j*n+k]+cg*cg*bk;
                    a[j*n+k] = 0.0;
                    b[j*n+k] = 0.0;

/************************************************************************
          Update the eigenvector matrix after each rotation
*************************************************************************/
                    for (i=0 ; i<n ; i++)
                    {
                        xj = x[i*n+j];
                        xk = x[i*n+k];
                        x[i*n+j] = xj+cg*xk;
                        x[i*n+k] = xk+ca*xj;
                    };
                };
            };
        };
/************************************************************************
    Update the eigenvalues after each sweep
*************************************************************************/
        for (i=0 ; i<n ; i++)
        {
            ii=i*n+i;
            if (a[ii]<=0.0 || b[ii]<=0.0)
            {
                printf( "*** Error  solution stop in *gomp_jacobi*\n Matrix not positive definite.");
                return(1);
            };
            eigv[i] = a[ii] / b[ii];
        };

/************************************************************************
    check for convergence
*************************************************************************/
        convergence = 1;
        for(i=0 ; i<n ; i++)
        {
            tol = rtol*d[i];
            dif = fabs(eigv[i]-d[i]);
            if (dif > tol) convergence = 0;
            if (! convergence) break;
        };
/************************************************************************
    check if all off-diag elements are satisfactorily small
*************************************************************************/
        if (convergence)
        {
            eps=rtol*rtol;
            for (j=0 ; j<nr ; j++)
            {
                jj=j+1;
                for (k=jj ; k<n ; k++)
                {
                    epsa = ( a[j*n+k]*a[j*n+k] / ( a[j*n+j]*a[k*n+k] ) );
                    epsb = ( b[j*n+k]*b[j*n+k] / ( b[j*n+j]*b[k*n+k] ) );
                    if ( epsa >= eps || epsb >=eps ) convergence = 0;
                    if (! convergence) break;
                };
                if (! convergence) break;
            };
        };
        if (! convergence)
        {
            for (i=0 ; i<n ; i++)
                d[i] = eigv[i];
        };
        if (convergence) break;
    };
/************************************************************************
  Fill out bottom triangle of resultant matrices and scale eigenvectors
*************************************************************************/
    for (i=0 ; i<n ; i++)
        for (j=i ; j<n ; j++)
        {
            b[j*n+i] = b[i*n+j];
            a[j*n+i] = a[i*n+j];
        };

    for (j=0 ; j<n ; j++)
    {
        bb = sqrt( b[j*n+j] );
        for (k=0 ; k<n ; k++)
            x[k*n+j] /= bb;
    };
    printf( "gomp_jacobi: nsweeps %d\n", nsweep );
    return 1 ;
}

