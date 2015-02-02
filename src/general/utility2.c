/*  

Copyright (c) 1990 - 2004 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved

  
Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>

#if defined(WIN32)
#include <windows.h>
#include <time.h>
#else
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <unistd.h>
#endif

#include <tcl.h>

#include "gomfile.h"
#include "gomproc.h"
#include "printmsg.h"

#include "stdafx.h"

#if 0
/**************************************************************************/
void gomp_Rest_sec(float seconds) 
    /* rest in this routine for the 'second' number of secs */
/**************************************************************************/
{

#if !defined(WIN32)
    long   ret_in;
    long   clk_tck;
    long   ret_out;
    float  retf_in;
    float  retf_out;

    struct tms buffer;

    clk_tck = sysconf(_SC_CLK_TCK);

    ret_in = times(&buffer);  /*start time */
    retf_in = (float)ret_in / (float)(clk_tck);
    while(1) {
        ret_out = times(&buffer);
        retf_out = (float)ret_out / (float)(clk_tck);
        if( (retf_out - retf_in) > seconds )
            return;
    }
#else
    gomp_PrintMessage("Dummy in rest seconds");
#endif
}
#endif
/**************************************************************************/
void gomp_Get_cpu_secs(float *Msecs , float *Csecs)  /* get used cpu seconds */
/**************************************************************************/
{
/* return used cpu seconds. This routine can to be called  
   with a NULL value to set the cpu time = 0 */

    static float start       = 0.0;
    static float start_child = 0.0;

#if defined(WIN32)
    static float temp;
    clock_t Itemp;

    if(!strcmp(gomp_GetOS(), WINDOWS_NT)) {
        FILETIME Ctime;
        FILETIME Etime;
        FILETIME Ktime;
        FILETIME Utime;

        GetProcessTimes(
            GetCurrentProcess(), /* specifies the process of interest        */
            &Ctime,              /* when the process was created             */
            &Etime,              /* when the process exited                  */
            &Ktime,              /* time the process has spent in kernel mode*/
            &Utime               /* time the process has spent in user mode  */
            );

        /* Ctime.dwLowDateTime,Etime.dwLowDateTime,Ktime.dwLowDateTime,Utime.dwLowDateTime
         * Ctime.dwHighDateTime,Etime.dwHighDateTime,Ktime.dwHighDateTime,Utime.dwHighDateTime
         */
   
        if( ! Msecs ) {
            start       = (float)(Utime.dwHighDateTime * pow(2.0 , 32.0)) + /* high 32 bits */
                (float)(Utime.dwLowDateTime);                     /* low  32 bits */
            start      *= 100.0e-09f;
            start_child = 0.0;
        }
        else {
            temp        = (float)(Utime.dwHighDateTime * pow(2.0 , 32.0)) + 
                (float)(Utime.dwLowDateTime);
            *Msecs      = 100.0e-09f * temp - start;
            *Csecs      = 0.0;
        }
    }
    else { /* non Windows NT  */

        Itemp = clock();
        if(Itemp < 0) {
            gomp_PrintERROR("can't get used cpu time");
            start       = 0.0;
            start_child = 0.0;
        }

        if( ! Msecs ) {
            start       = (float)Itemp / (float)CLOCKS_PER_SEC;
            start_child = (float)0.0;
        }
        else {
            temp        = (float)Itemp / (float)CLOCKS_PER_SEC;
            *Msecs      = temp - start;
            *Csecs      = 0.0;
        }
    }

#else /* non Windows */

    long ret_in;
    long clk_tck;

    struct tms buffer;

    clk_tck = sysconf(_SC_CLK_TCK);

    ret_in = times(&buffer);  /*start time */

    if( ! Msecs ) {
        start       = (float)(buffer.tms_utime  + buffer.tms_stime)  / 
            (float)(clk_tck);
        start_child = (float)(buffer.tms_cutime + buffer.tms_cstime) / 
            (float)(clk_tck);
    }
    else {
        *Msecs      = (float)(buffer.tms_utime  + buffer.tms_stime)  / 
            (float)(clk_tck) - start;
        *Csecs      = (float)(buffer.tms_cutime + buffer.tms_cstime) / 
            (float)(clk_tck) - start_child;
    }

#endif

    return;
}
/***********************************************************************/
int gomp_Check_if_file_exists(const char *file_name)  
    /* on return:   0 ok, file is there */
    /*           != 0 can't access file */
/***********************************************************************/
{
    static int retv;
    char OutText[BUFF_LEN];

#if defined(WIN32)
    FILE *File_p;
    File_p = fopen(file_name , "r");
    if(File_p != NULL) {
        fclose(File_p);
        return(0);
    }
    else {
        sprintf(OutText,"File '%s' does not exist ",file_name);
        gomp_PrintERROR(OutText);
        return(1);
    }
#else   
/* first check existence of file */
    retv = access(file_name,F_OK);

    if(retv == -1) {
        sprintf(OutText,"?ERROR - File '%s' does not exist ",file_name);
        gomp_PrintMessage(OutText);
        return(1);
    }    /* problems with the file */
    else { /* check read permission */
        retv = access(file_name,R_OK);
        if(retv == -1) {
            sprintf(OutText,"?ERROR - no read permission to file '%s' ",file_name);
            gomp_PrintMessage(OutText);
            return(1);
        }
    }
#endif

#if defined(JUNK)
    if(DebugS.DebugL > 0)
        (void)file_info(file_name);
#endif

    return(0);
}
#if 0
/***********************************************************************/
int CheckIfFileExists(const char *file_name)  
    /* on return:   0 ok, file is there */
    /*           != 0 can't access file */
/***********************************************************************/
{
    static int retv;
    char OutText[BUFF_LEN];

#if defined(WIN32)
    FILE *File_p;
    File_p = fopen(file_name , "r");
    if(File_p != NULL) {
        fclose(File_p);
        return(0);
    }
    else {
        sprintf(OutText,"File '%s' does not exist ",file_name);
        gomp_PrintERROR(OutText);
        return(1);
    }

#else
/* first check existence of file */
    retv = access(file_name,F_OK);

    if(retv == -1) {
        sprintf(OutText,"?ERROR - File '%s' does not exist ",file_name);
        gomp_PrintMessage(OutText);
        return(1);
    }    /* problems with the file */
    else { /* check read permission */
        retv = access(file_name,R_OK);
        if(retv == -1) {
            sprintf(OutText,"?ERROR - no read permission to file '%s' ",file_name);
            gomp_PrintMessage(OutText);
            return(1);
        }
    }

#endif

#if defined(JUNK)
    if(DebugS.DebugL > 0)
        (void)file_info(file_name);
#endif

    return(0);
}
#endif
#if 0
/**************************************************************************/
float gomp_GetCPUtime()  /* get used cpu seconds */
/**************************************************************************/
{
/* return used cpu seconds. This routine can to be called  
   with a negative value to set the cpu time = 0 */

    static float temp;

#if defined(WIN32)

    if(!strcmp(gomp_GetOS(), WINDOWS_NT)) {
        FILETIME Ctime;
        FILETIME Etime;
        FILETIME Ktime;
        FILETIME Utime;

        GetProcessTimes(
            GetCurrentProcess(),/* specifies the process of interest         */
            &Ctime,             /* when the process was created              */
            &Etime,             /* when the process exited                   */
            &Ktime,             /* time the process has spent in kernel mode */
            &Utime              /* time the process has spent in user mode   */
            );

        /* Ctime.dwLowDateTime,Etime.dwLowDateTime,Ktime.dwLowDateTime,Utime.dwLowDateTime
         * Ctime.dwHighDateTime,Etime.dwHighDateTime,Ktime.dwHighDateTime,Utime.dwHighDateTime
         */
        temp        = 1.e-07f * ((float)(Utime.dwHighDateTime * pow(2.0 , 33.0)) + 
                                 (float)(Utime.dwLowDateTime));
    }
    else {
        clock_t Itemp;

        Itemp = clock();
        if(Itemp < 0) {
            gomp_PrintERROR("can't get used cpu time");
            Itemp = 0;
        }
        temp  = (float)Itemp / (float)CLOCKS_PER_SEC;
    }

#else /* non windows */

    long ret_in;
    long clk_tck;

    struct tms buffer;

    clk_tck = sysconf(_SC_CLK_TCK);

    ret_in = times(&buffer);  /*start time */

    temp        = (float)(buffer.tms_utime  + buffer.tms_stime)  / 
        (float)(clk_tck);

#endif

    return(temp);
}
#endif
