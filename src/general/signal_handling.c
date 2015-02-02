/*

Copyright (c) 1991 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <setjmp.h>
#include <errno.h>
#include <tcl.h>

#include "gommain.h"
#include "printmsg.h"

#include "stdafx.h"

static void  Sig_handler(int);
static void  SigSIGUSR2(int);

/**************************************************************************/
/* define the signal handler */
void Sig_handler(int i)
/**************************************************************************/
{
    int  err;

    err = errno;
    fprintf(stderr,"**** Signal was caught ****\n");
    errno = err;
    perror("=> ");
    fprintf(stderr,"Signal code is: %d, errno : %d.\n",i,err);
    fprintf(stderr,"Type of error: %s\n",strerror(err));
    exit(i);
}
#if 0
/**************************************************************************/
/* define the signal handler */
void gomp_Sig_die(int i)
/**************************************************************************/
{
    exit(i);
}
#endif
/**************************************************************************/
/* define the signal handler */
void SigSIGUSR2(int i)
/**************************************************************************/
{
    gomp_PrintMessage("Received SIGUSR2 signal, will stop the redisplay");

    (void)gomp_SetDisplayInterrupt(ON);
    return;
}

/**************************************************************************/
int gomp_DefineSignals()
/**************************************************************************/
{
/*  handle signals   */

#if defined(WIN32)
/* 
   No signal handling for Windows because I had severe problems
   on Windows 95/98 but strangely not on NT and 2000.

   I decided to disable the Windows signals until I understand
   better the procedure.
*/
/*   signal(SIGTERM,Sig_handler);  *//*  illegal instruction         */

#else
    signal(SIGBUS,Sig_handler);  /*  bus error                   */
    signal(SIGUSR2,SigSIGUSR2); 
/*  signal(SIGINT,SIG_IGN);            *  control C handle            */
    signal(SIGSEGV,Sig_handler); /*  segmentation violation      */
    signal(SIGILL,Sig_handler);  /*  illegal instruction         */
    signal(SIGFPE,SIG_IGN);           /*  don't listen to floating    */
                                      /*  point exceptions            */
#endif

/*********************/

    return(0);
}
