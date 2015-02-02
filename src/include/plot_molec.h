extern int   gomp_PlotMoleculeStick(int);
extern int   gomp_PlotMoleculeCPK(int);
extern int   gomp_PlotMoleculeLicorice(int);

extern int   gomp_SetMoleculeLineWidth(int);
extern int   gomp_GetMoleculeLineWidth(void);

extern int   gomp_GetSphereQuality(void);
extern int   gomp_SetSphereQuality(int);

extern int   gomp_GetCylinderQuality(void);
extern int   gomp_SetCylinderQuality(int);

extern int   gomp_SetCrossLen(float);
extern float gomp_GetCrossLen(void);

extern int   gomp_GetBondDisplayStyle(void);
extern int   gomp_SetBondDisplayStyle(int);

extern void  gomp_PlotCylinder(float , float , float ,
                             float , float , float ,
                             float , float , float ,
                             float , float , float ,
                             float , float , float);

extern int   gomp_SetArrowAD2CD(float);
extern int   gomp_SetArrowL2H(float);

extern void gomp_Arrow( float , float , float , float , float , float , float);
extern void gomp_Vector2Angles(float,float,float,float *,float *,float *);

extern int gomp_SetSphereQuad(void);
