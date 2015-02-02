extern int  gomp_PushText2AnnotateStack(int, const char *, const char *,
                                        float , float  , float ,
                                        float , float  , float );
extern int  gomp_DeleteTextStack(void);
extern int  gomp_ReadText2ModelFile(FILE *);
extern int  gomp_WriteText2ModelFile(FILE *);
extern void gomp_PrintString(const char *,const char *);
extern int  gomp_GetEntriesInTextStack(void);
extern int  gomp_PlotTextStack(void);
extern int  gomp_ResetNameStack(void);

