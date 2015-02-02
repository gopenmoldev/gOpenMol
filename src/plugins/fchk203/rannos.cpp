#include "rannos.hpp"
#include <stdio.h>
/*I doubt I can assert copyright on this stuff- it's all pretty simple implementations of
other people's code.*/
//static long rangen;

namespace Plugin {
namespace Fchk {
    
void raninit(long seed){

if(rangen>0)
    rangen=-seed;
else
    rangen=seed;
ranx();
}

long getranseed()
{
    return rangen;
}

void setranseed(long nseed)
{
    rangen=nseed;
}

long numrans()
{
    return nrans;
}

double ranx(){
//returns uniform random deviate in (0,1).  From Numerical Recipes
int j;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
double temp;

if(rangen<=0){
 nrans=0;
 if (-(rangen)<1)rangen=1;
 else rangen=-(rangen);
 idum2=(rangen);
 for(j=NTAB+7;j>=0;j--){
   k=(rangen)/IQ1;
   rangen=IA1*(rangen-k*IQ1)-k*IR1;
   if (rangen<0) rangen+=IM1;
   if (j<NTAB) iv[j]=(rangen);}
 iy=iv[0];}
nrans++;
k=(rangen/IQ1);
rangen=IA1*(rangen-k*IQ1)-k*IR1;
if(rangen<0)rangen+=IM1;
k=idum2/IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2;
if(idum2<0) idum2+=IM2;
j=iy/NDIV;
iy=iv[j]-idum2;
iv[j]= (rangen);
if(iy<1)iy+=IMM1;
if((temp=AM*iy)>RNMX)return RNMX;
else return temp;}

/*This is a pretty good generator based on Fibonacci sequences.  My wife gave me
a fortran file which I converted. I think its the same as this:
This Random Number Generator is based on the algorithm in a FORTRAN version published 
by George Marsaglia and Arif Zaman, Florida State University. 

THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE. 
It passes ALL of the tests for random number generators and has a period of 2^144, is 
completely portable (gives bit identical results on all machines with at least 24-bit 
mantissas in the floating point representation). 

The algorithm is a combination of a Fibonacci sequence (with lags of 97 and 33, and 
operation "subtraction plus one, modulo one") and an "arithmetic sequence" (using 
subtraction). */


static double U[97];
static double C,CD,CM;
static int I97,J97;
static int RanMarInit=0;

void rmarin(int IJ,int KL)
{
    double s,t;
    int i,j,k,l,ii,jj,m;

    if((IJ<0)||(IJ>31328)||(KL<0)||(KL>30081))
    {
        printf("First seed must be in [0,31328]\n");
        printf("Second seed must be in [0,30081]\n");
    }
    i=((IJ/177)%177)+2;
    j=IJ%177 +2;
    k=((KL/169)%178)+1;
    l=KL%169;
    for(ii=0;ii<97;ii++)
    {
        s=0.0;
        t=0.5;
        for(jj=1;jj<=24;jj++)
        {
            m=(((i*j)%179)*k)%179;
            i=j;
            j=k;
            k=m;
            l=(53*l+1)%169;
            if((l*m)%64>=32)
                s+=t;
            t*=0.5;
        }
        U[ii]=s;
    }
    C=362436.0/16777216.0;
    CD=7654321.0/16777216.0;
    CM=16777213.0/16777216.0;
    I97=96;
    J97=32;
    RanMarInit=1;
}

double ranmar(void)
{
    double uni;

    if(RanMarInit==0)
    {
        printf("Random number generator not initialized\n");
    }
    uni=U[I97]-U[J97--];
    if(uni<0.0)
        uni+=1.0;
    U[I97--]=uni;
    if(I97<0)I97=96;
    if(J97<0)J97=96;
    C-=CD;
    if(C<0.0)
        C+=CM;
    uni-=C;
    if(uni<0.0)
        uni+=1.0;
    return uni;
}

void outrmarset(FILE* of)
{
    fwrite(&C,sizeof(double),1,of);
    fwrite(&CD,sizeof(double),1,of);
    fwrite(&CM,sizeof(double),1,of);
    fwrite(&I97,sizeof(int),1,of);
    fwrite(&J97,sizeof(int),1,of);
    fwrite(U,sizeof(double),97,of);
}

void getrmarset(FILE* inf,int nt)
{
    int i;

    for(i=0;i<=nt;i++)
    {
        fread(&C,sizeof(double),1,inf);
        fread(&CD,sizeof(double),1,inf);
        fread(&CM,sizeof(double),1,inf);
        fread(&I97,sizeof(int),1,inf);
        fread(&J97,sizeof(int),1,inf);
        fread(U,sizeof(double),97,inf);
    }
}

/*End of FSU generator-related code*/

double granx(double sigma){
//returns (approximate) Gaussian deviate with sigma input
double q,z,f,rr;

sigma/=1.48;
q=ranmar();
if(q==0.5){
 return 0.0;}
if(q>0.5){
 f=1.0;
 z=2*(q-0.5);}
else{
 f=-1.0;
 z=2*q;}
rr=SP2*z+ONET*(1.0/(1-0.95*z)-1);
return rr*f*sigma;}

void funranx(double(*Fn)(int,double*),double Fnmax,int ndim,double* xmax,double* x){
//returns a random vector x with weight function Fn.  Takes maximal values and number of dimensions as additional inputs

double q,qq;
int flag=0,i;

while(flag==0){
 q=Fnmax*ranx();
 for(i=0;i<ndim;i++){
     x[i]=xmax[i]*ranx();}
 qq=Fn(ndim,x);
 if(q<qq){
    flag=1;}}}

void sobseq(int *n,double x[]){
//implements a low-dimensional Sobol' sequence.  From Numerical Recipes
int j,k,l;
unsigned long i,im,ipp;
static double fac;
static unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4};
static unsigned long iv[MAXDIM*MAXBIT+1]={
    0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

if(*n<0){
    for(k=1;k<=MAXDIM;k++)ix[k]=0;
    in=0;
    if(iv[1]!=1) return;
    fac=1.0/(1L<<MAXBIT);
    for(j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM)iu[j]=&iv[k];
    for(k=1;k<=MAXDIM;k++){
        for(j=1;j<=mdeg[k];j++)iu[j][k]<<=(MAXBIT-j);
        for(j=mdeg[k]+1;j<=MAXBIT;j++){
            ipp=ip[k];
            i=iu[j-mdeg[k]][k];
            i^=(i>>mdeg[k]);
            for(l=mdeg[k]-1;l>=1;l--){
                if(ipp&1)i^=iu[j-l][k];
                ipp>>=1;
            }
            iu[j][k]=i;
        }
    }
}
else {
    im=in++;
    for(j=1;j<=MAXBIT;j++){
        if(!(im&1))break;
        im>>=1;
    }
    if(j>MAXBIT)
        printf("Problem here! %i\n",j);
    im=(j-1)*MAXDIM;
    for(k=1;k<=IMIN(*n,MAXDIM);k++){
        ix[k]^=iv[im+k];
        x[k-1]=ix[k]*fac;
    }
}
}                                  

} // namespace Fchk
} // namespace Plugin
