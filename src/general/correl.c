
/*     These subroutines are based on ideas from the "Numerical Recipes in C"
       subroutine library by
       W.H. Press, B.P. Flannery, S.A. Teukolsky and W.T. Vetterling

       Leif Laaksonen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <stdlib.h>

#include "correl.h"

#include "stdafx.h"

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

static void four1(float data[], int nn, int isign);
static void twofft(const float data1[],const float data2[],float fft1[],float fft2[],int n);
static void realft(float data[], int n, int isign);

/**************************************************************************/
void gomp_Correl(const float *data1, const float *data2, float *fft, int n, float *ans)
/*
  const float *data1,*data2,*fft,*ans;
  int n;
*/
/**************************************************************************/
{

/*  fft has to be defined 2*length of data1 and data2 ( = 2*n) */
/*  and n has to be a power of 2 ( => n = 2**x)                */

    int no2,i;
    float dum;

    twofft(data1,data2,fft,ans,n);
    no2=n/2;
    for (i=2;i<=n+2;i+=2) {
        ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
        ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
    }
    ans[2]=ans[n+1];
    realft(ans,no2,-1);

}
void twofft(const float data1[],const float data2[],float fft1[],float fft2[],int n)
/*
  float data1[],data2[],fft1[],fft2[];
  int n;
*/
{
    int nn3,nn2,jj,j;
    float rep,rem,aip,aim;

    nn3=1+(nn2=2+n+n);
    for (j=1,jj=2;j<=n;j++,jj+=2) {
        fft1[jj-1]=data1[j];
        fft1[jj]=data2[j];
    }
    four1(fft1,n,1);
    fft2[1]=fft1[2];
    fft1[2]=fft2[2]=0.0;
    for (j=3;j<=n+1;j+=2) {
        rep=0.5*(fft1[j]+fft1[nn2-j]);
        rem=0.5*(fft1[j]-fft1[nn2-j]);
        aip=0.5*(fft1[j+1]+fft1[nn3-j]);
        aim=0.5*(fft1[j+1]-fft1[nn3-j]);
        fft1[j]=rep;
        fft1[j+1]=aim;
        fft1[nn2-j]=rep;
        fft1[nn3-j] = -aim;
        fft2[j]=aip;
        fft2[j+1] = -rem;
        fft2[nn2-j]=aip;
        fft2[nn3-j]=rem;
    }
}

void realft(float data[], int n, int isign)
/*
  float data[];
  int n,isign;
*/
{
    int i,i1,i2,i3,i4,n2p3;
    float c1=0.5,c2,h1r,h1i,h2r,h2i;
    double wr,wi,wpr,wpi,wtemp,theta;

    theta=3.141592653589793/(double) n;
    if (isign == 1) {
        c2 = -0.5;
        four1(data,n,1);
    } else {
        c2=0.5;
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    n2p3=2*n+3;
    for (i=2;i<=n/2;i++) {
        i4=1+(i3=n2p3-(i2=1+(i1=i+i-1)));
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[1] = (h1r=data[1])+data[2];
        data[2] = h1r-data[2];
    } else {
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        four1(data,n,-1);
    }
}


void four1(float data[], int nn, int isign)
/*
  float data[];
  int nn,isign;
*/
{
    int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    float tempr,tempi;

    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=2*mmax;
        theta=6.28318530717959/(isign*mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

int gomp_check_pow_2(int n)
{

    int i;

    i = 1;
    for(;;) {
        if(i >= n)
            break;
        i *= 2;
    }

    return(i);
}

#undef SWAP


