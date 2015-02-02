
/*
                           Copyright (c) 2002 - 2004 by:
        Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
                      Confidential unpublished property of 
                              Leif Laaksonen  
                            All rights reserved

*/

#include "gomstdio.h"

extern int   gomp_FileNameIsURL(char *);
extern char *gomp_Fgets(char *, int, FILE *);
extern int   gomp_Check_if_file_exists(const char *);
extern int   gomp_GetFileSize(const char *);

extern int   gomp_OpenPipe2Program(const char *, const char *);
extern int   gomp_ClosePipe2Program(void);
extern int   gomp_SendPipe2Program(const char *);
extern int   gomp_CheckTextFileLineEnding(const char *);

#ifdef WIN32
#    define DIR_SEP "\\"
#else
#    define DIR_SEP "/"
#endif
