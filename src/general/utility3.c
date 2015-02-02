/*

System statistics for the running process

Leif Laaksonen 1992 - 2004

Centre for Scientific Computing, Espoo, Finland

*/

#include "maindefs.h"

#include "gomstdio.h"

#ifndef WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/signal.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>
#endif

#if defined(IRIX)
#include <sys/fault.h>
#include <sys/syscall.h>
#include <sys/procfs.h>
#endif

/*#include "gomfile.h"*/
#include "gomproc.h"
#include "printmsg.h"

#include "stdafx.h"

#define String_Len  20

#if HAVE_GETRUSAGE
struct rusage gomp_RUsage;
#elif defined(IRIX)
prusage_t gomp_ProcessInfo;
#endif

/***************************************************************************/
int gomp_RunStatistics()
/***************************************************************************/
{
#if HAVE_GETRUSAGE

    if(getrusage(RUSAGE_SELF , &gomp_RUsage)) {
        gomp_PrintMessage("?ERROR - can't get process info");
        return(-1);
    }

#elif defined(IRIX)

    pid_t CurrentPID;
    int Fdes;
    static char PidString[String_Len];
    CurrentPID = getpid();              /* get the pid */
    sprintf(PidString,"/debug/%5.5d",CurrentPID);  

    if((Fdes = open(PidString , O_RDONLY)) == -1) {
        gomp_PrintMessage("?ERROR - can't open process file");
        return(-1);
    }

    if(ioctl(Fdes , PIOCUSAGE , &gomp_ProcessInfo) == -1) {
        gomp_PrintMessage("?ERROR - can't get process info");
        close(Fdes);
        return(-1);
    }

    close(Fdes);

#else

    gomp_PrintMessage("Dummy in run statistics");

#endif

    return(0);
}
