/*
     interface to a c routine to perform a dynamical
     memory allocation

     Leif Laaksonen (1992)

     Centre for Scientific Computing, Espoo, FINLAND
*/

#include <stdio.h>
#include <sys/types.h>
#include <malloc.h>

/*
      bv:     is an integer to contain the address to the reserved
              memory

      reqlen: contains the number of bytes to be reserved
*/
#if defined(WIN32)
void  FMALLOC(int * , int *);
#else
void  fmalloc_(int *bv , int *reqlen);
#endif

/************************************************************************/
#if defined(WIN32)
void  FMALLOC(int *bv , int *reqlen)
#else
void  fmalloc_(int *bv , int *reqlen)
#endif
/************************************************************************/
{
      float *BV;

      BV = (float *)malloc(*reqlen);
      if(BV == NULL) {
       perror("Can't allocate memory");
       exit(1);}

      *bv = (int)BV;
}
