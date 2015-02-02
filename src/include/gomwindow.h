#define SINGLE_WINDOWING      0
#define MULTI_WINDOWING       1

/* type of windows            */
#define STRUCTURE_WINDOW      0
#define LDP_WINDOW            1

extern int gomp_SetWindowingStyle(int);
extern int gomp_GetWindowingStyle(void);

extern int gomp_GetWindowID(void);
extern int gomp_SetWindow(int);
extern int gomp_IconifyWindow(void);
extern int gomp_DeIconifyWindow(void);
extern int gomp_MoveWindow(int, int);
extern int gomp_ResizeWindow(int, int);
extern int gomp_FullScreenWindow(void);
extern int gomp_GetWindowIDFromStack(int);

extern int gomp_PopWindow(void);
extern int gomp_GetNumDefinedWindows(void);
extern int gomp_GetWindowParameters(int *, int *, int *, int *);
extern int gomp_GetWindowInfo(int *, int *, int *, int *);

extern int gomp_GetWindowTypeFromStack(int);
extern int gomp_PushToWindowStack(int, int, const char *);
