#define STEREO_DISPLAY_OFF 0
#define STEREO_DISPLAY_ON  1

extern int   gomp_SetStereoPlotState(int);
extern int   gomp_GetStereoPlotState(void);
extern int   gomp_SetStereoPlotTranslate(float);
extern float gomp_GetStereoPlotTranslate(void);
extern int   gomp_SetStereoPlotAngle(float);
extern float gomp_GetStereoPlotAngle(void);

extern void  gomp_QuadStereoOff(void);
extern void  gomp_QuadStereoOn(void);
extern int   gomp_QuadStereoIsOn(void);

extern void  gomp_SetQuadStereoHalfAngle(float a);
extern float gomp_GetQuadStereoHalfAngle(void);

extern int   gomp_CheckHardwareStereo(void);

extern int   gomp_SetStereoDisplayState(int);
extern int   gomp_GetStereoDisplayState(void);
