/*
                           Copyright (c) 2002 - 2004 by:
        Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
                      Confidential unpublished property of 
                              Leif Laaksonen  
                            All rights reserved

*/

extern int   gomp_DeletePlotAxisList(void);

extern int   gomp_PlotCoordAxis(void *, int, int);

extern int   gomp_ParsePlotAxisList(const char *, const char *,
                                    const char *, const char *,
                                    const char *, const char *);

extern int   gomp_GetPlotAxisListLength(void);
extern int   gomp_GetPlotAxisCoordPoints(void);
extern int   gomp_GetPlotAxisLength(float * , float * , float *);
extern const float *gomp_GetPlotAxisXCoord(void);
extern const float *gomp_GetPlotAxisYCoord(void);
extern const float *gomp_GetPlotAxisZCoord(void);

extern const int *gomp_GetPlotAxisLengthP(void);
extern const int *gomp_GetPlotAxisStructureP(void);
extern const int *gomp_GetPlotAxisListP(void);
