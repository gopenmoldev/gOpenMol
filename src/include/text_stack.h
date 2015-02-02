extern int    gomp_PushText2Stack(const char *,const char *);
extern int    gomp_PushText2StackTop(const char *);
extern int    gomp_PrepareTextStack(void);
extern int    gomp_IsTextStackSet(void);
extern int    gomp_EntriesInTextStack(void);
extern int    gomp_TextStackLineLength(void);
extern const char *const*gomp_ReturnTextInStack(void);
extern int    gomp_SetOutput2Widget(void);


extern const char *gomp_GetGlobalTextString(void);
extern int         gomp_SetGlobalTextString(const char *);

extern int    gomp_PushNode2Stack(const char *NodeName);
