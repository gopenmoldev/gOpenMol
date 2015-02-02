/*  

Copyright (c) 2004 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Coded by:
Eero Häkkinen
*/

#include "maindefs.h"

#ifdef WIN32
#    define WIN32_LEAN_AND_MEAN
#    include <windows.h>
#else
#    include <sys/types.h>
#    include <sys/stat.h>
#    include <fcntl.h>
#    include <unistd.h>
#endif

#include <stdlib.h>
#include <string.h>

#define KEEP_STDIO_MACROS

#include "gomstdio.h"
#include "gomenv.h"
#include "gomfile.h"

#include "stdafx.h"

/**
 * MSVRCT creates temporary files to the `C:\' directory.
 * That is insane (and may be denied by directory permissions).
 * Here is a working solution.
 */
FILE *gomp_tmpfile()
{
    const char *tmpdir;
    char  filename[BUFF_LEN];
    FILE *file;
    int   i,j;

    tmpdir = gomp_ShowTempDir();
    if ( strlen(tmpdir) < BUFF_LEN - 16 ) {
        for ( i = rand() , j = 0 ; j < 0x100 ; j++ ) {
#ifdef WIN32
            HANDLE hdl;
#else
            int    fd;
#endif
            sprintf(
                filename, "%s" DIR_SEP "gom%x.tmp",
                tmpdir, ( i + j ) & 0xffff);
#ifdef WIN32
            /* Create a new file. */
            hdl = CreateFile(
                filename,
                GENERIC_READ | GENERIC_WRITE,
                0,
                NULL,
                CREATE_NEW,
                FILE_ATTRIBUTE_NORMAL | FILE_ATTRIBUTE_TEMPORARY,
                NULL);
            if ( hdl == INVALID_HANDLE_VALUE )
                continue;
            CloseHandle(hdl);
            /* Open a newly created file and delete it on close. */
            file = fopen(filename, "w+bD");
#else
            /* Create a new file. */
            fd = open(filename, O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
            if ( fd == -1 )
                continue;
            /* Detach it from a file system. */
            remove(filename);
            /* Get a FILE pointer to a file descriptor. */
            file = fdopen(fd, "w+b");
            if ( ! file )
                close(fd);
#endif
            if ( file )
                return file;
        }
    }
    return tmpfile();
}
