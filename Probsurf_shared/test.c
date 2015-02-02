/*
This code is based on the tcl-command language
Leif Laaksonen Center for Scientific Computing 1995, 1996
*/
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <tcl.h>


/* commands */
//int lulTest(ClientData , Tcl_Interp *, int , char **);
int Test_Init(void);

/* external */

//extern Tcl_Interp *lulGetTclInterp(void); /* get the tcl interp pointer */

/* ........ */


/*********************************************************************/
int Test_Init()
/*********************************************************************/
{
int code;
//Tcl_Interp *Interp;

//Interp = lulGetTclInterp();

/* test */
//Tcl_CreateCommand(Interp,"test",lulTest,(ClientData)NULL,
//(Tcl_CmdDeleteProc *)NULL);
printf("huuhaa \n");
return(TCL_OK);
}

/*********************************************************************/
int lulTest(ClientData clientdata , Tcl_Interp *interp,
int argc , char *argv[])
/*********************************************************************/
{
printf("Hello World\n");
return(0);
}

