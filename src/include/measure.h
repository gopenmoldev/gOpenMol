extern float gomp_CalculateDistance(const char *, const char *, const char *,
                                    const char *, const char *, const char *);
extern float gomp_CalculateAngle(const char *, const char *, const char *,
                                 const char *, const char *, const char *,
                                 const char *, const char *, const char *);
extern float gomp_CalculateTorsionAngle(
    const char *, const char *, const char *, 
    const char *, const char *, const char *,
    const char *, const char *, const char *, 
    const char *, const char *, const char *);
extern int gomp_CalculateCenterofMass(
    const char *, const char *, const char *, float *,float *,float *);
extern int gomp_CalcCoordinateCenter(
    const char *, const char *, const char *, float *,float *,float *);

extern void gomp_floDihedAngle(float, float, float ,
                               float, float, float ,
                               float, float, float ,
                               float, float, float ,
                               float *);
extern void gomp_BondAngle(float , float , float ,
                           float , float , float ,
                           float , float , float ,
                           float *);

extern float gomp_QuatFit(const float *, const float *, const float *, int,
                          float *, float *, float *, int,
                          const int *, int,
                          const int *, int,
                          int);
extern float gomp_CalculateQuatfit(
    const char *, const char *, const char *, const char *,
    const char *, const char *, const char *, const char *,
    int);
