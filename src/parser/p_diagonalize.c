/*

Copyright (c) 1996 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing , ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <sys/types.h>

#include "gomstring.h"
#include "gomtcl.h"
#include "math_oper.h"
#include "memalloc.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_DiagonalizeCommand(ClientData clientdata, Tcl_Interp *interp,
                          int argc, const char **argv)
/*********************************************************************/
{
    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    static char  Text2[BUFF_LEN];
    static const char *Value;
    static int   Code;
    static int   Dimension;

/* #1   copy ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "diag$onalize")) {

        gomp_CopyString(Text ,gomp_GetNextFromParserStack(argc,(const char **)NULL),BUFF_LEN);

        if(Text[0] == (char)NULL) {
            gomp_PrintERROR("matrix name is missing");
            return(TCL_ERROR);
        }

/* check if there is a matrix(0,0) element available */
        sprintf(Text1,"set gomTemp [info exists %s(0,0)]",Text);

        Code = Tcl_GlobalEval(gomp_GetTclInterp(), Text1);
        if(Code != TCL_OK) {
            gomp_PrintERROR("can't eval command in 'diagonalize (1)'");
            return(TCL_ERROR);
        }

        Value = Tcl_GetVar(gomp_GetTclInterp(),"gomTemp",0);

        if(!Value) {
            gomp_PrintERROR("matrix element (0,0) is missing");
            return(TCL_ERROR);
        }


/* get the dimension for the matrix (array) */
        sprintf(Text1,"set gomTemp [array size %s]",Text);

        Code = Tcl_GlobalEval(gomp_GetTclInterp(), Text1);
        if(Code != TCL_OK) {
            gomp_PrintERROR("can't eval command in 'diagonalize (3)'");
            return(TCL_ERROR);
        }

        Value = Tcl_GetVar(gomp_GetTclInterp(),"gomTemp",0);

        if(Value) {
            Dimension = atoi(Value);
        } else {
            gomp_PrintERROR("can't get matrix dimension");
            return(TCL_ERROR);
        }

/* check some things before we do the heavy job ...*/
        if(Dimension%2) {
            gomp_PrintERROR("the defined matrix has to be a n x n matrix");
            return(TCL_ERROR);
        }

        Dimension /= 2;
        {
            int     iLoop;
            int     jLoop;
            int     ij;

            double  *a;
            double  *b;
            double  *x;

            double  *eigv;
            double  *d;
            double  rtol = 10.e-12;

            a    = gomp_AllocateDoubleVector(Dimension * Dimension);
            b    = gomp_AllocateDoubleVector(Dimension * Dimension);
            x    = gomp_AllocateDoubleVector(Dimension * Dimension);

            eigv = gomp_AllocateDoubleVector(Dimension);
            d    = gomp_AllocateDoubleVector(Dimension);

            for (jLoop = 0 ; jLoop < Dimension  ; jLoop++) {

                eigv[jLoop] = 0.0;
                d[jLoop]    = 0.0;

                ij = Dimension * jLoop;

                for (iLoop = 0 ; iLoop < Dimension ; iLoop++) {

                    a[iLoop + ij] = 0.0;
                    b[iLoop + ij] = 0.0;
                    x[iLoop + ij] = 0.0;

                    if(iLoop == jLoop)
                        b[iLoop + ij] = 1.0;

/* get the values from the Tcl side ... */
                    sprintf(Text1,"%s(%d,%d)",Text,jLoop,iLoop);
                    Value = Tcl_GetVar(gomp_GetTclInterp(),Text1,0);

                    if(Value) {
                        a[iLoop+ij] = atof(Value);
                    } else {
                        free(a);free(b);free(x);free(eigv);free(d);
                        gomp_PrintERROR("can't get matrix element");
                        return(TCL_ERROR);
                    }
                }
            }
/*
  for(jLoop = 0 ; jLoop < Dimension ; jLoop++) {
  for(iLoop = 0 ; iLoop < Dimension ; iLoop++) {
  printf("%d %d %f %f\n",iLoop,jLoop,a[iLoop + Dimension*jLoop],
  b[iLoop + Dimension*jLoop]);
  }
  }
*/
            iLoop = gomp_jprJacobi(a,b,x,eigv,d,Dimension,rtol);
                
            printf("\n");
            for(jLoop = 0; jLoop < Dimension ; jLoop++) {
                printf("%f : ",eigv[jLoop]);
                for(iLoop = 0; iLoop < Dimension ; iLoop++) {
                    printf("%f ",x[iLoop + Dimension * jLoop]);
                }
                printf("\n");
            }

            for(jLoop = 0 ; jLoop < Dimension ; jLoop++) {
                sprintf(Text1,"eigval(%d)",jLoop);
                sprintf(Text2,"%f",eigv[jLoop]);
                Value = Tcl_SetVar(gomp_GetTclInterp(),Text1,Text2,0);

                if(!Value) {
                    free(a);free(b);free(x);free(eigv);free(d);
                    gomp_PrintERROR("can't save eigenvalue");
                    return(TCL_ERROR);
                }

                for(iLoop = 0 ; iLoop < Dimension ; iLoop++) {
                    sprintf(Text1,"eigvec(%d,%d)",jLoop,iLoop);
                    sprintf(Text2,"%f",x[iLoop + Dimension * jLoop]);
                    Value = Tcl_SetVar(gomp_GetTclInterp(),Text1,Text2,0);

                    if(!Value) {
                        free(a);free(b);free(x);free(eigv);free(d);
                        gomp_PrintERROR("can't save eigenvalue");
                        return(TCL_ERROR);
                    }
                }
            }

            free(a);free(b);free(x);free(eigv);free(d);
        }

        return(TCL_OK);
    }

/*  E R R O R command not recognized         */
    gomp_PrintERROR("'diagonalize' command not recognized");
    return(TCL_ERROR);
}

