#define COLOR_MAPPING_TYPE_TEXTURE 0
#define COLOR_MAPPING_TYPE_RAINBOW 1

extern int gomp_SetColorMappingType(int);
extern int gomp_GetColorMappingType(void);

extern int gomp_ColourName2RGB(const char *, float * , float * , float *);
extern int gomp_ColourName2RGBSilent(const char *, float *, float *, float *);

extern int gomp_ColorByCharge(int , int , const int *, double *, double *);
extern int gomp_ColorByFourth(int , int , const int *);
extern int gomp_ColorByResidueNumber(int , int , const int *);
extern int gomp_ColorByVector(int, int , const int *, double *, double *);

extern int gomp_ReadColourTable(const char *);
extern int gomp_SetDisplayColourType(int);
extern int gomp_GetDisplayColourType(void);

extern int gomp_Prepare1DTexture(void);
extern int gomp_ResetPrepare1DTexture(void); 

extern int gomp_RGB2Grayscale(float *, float *, float *);

extern int gomp_SetBGColor(float , float , float);
extern int gomp_GetBGColor(float *, float *, float *);

extern void gomp_PreRainbow(double ,float *,float *,float *);

extern int gomp_GetColourTableLength(void);
