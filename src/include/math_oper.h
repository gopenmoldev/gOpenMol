extern int gomp_CalcPlaneFrom3Points(
    const float * , const float * , const float * , 
    float * , float * , float * , float *);

extern float  gomp_VecSum(const float * , int);
extern float  gomp_FMaxi(const float * , int);
extern float  gomp_FMini(const float * , int);

extern float gomp_Simpson38(const float * , int , float);
extern float gomp_Simpson(const float * , int , float);
extern float gomp_Trapez(const float * , int , float);
extern int gomp_nreg(const float *, const float *, int, int, float *);
extern int gomp_jprJacobi(
    double  *,double  *,double  *,double  *,double  *,int,double);
