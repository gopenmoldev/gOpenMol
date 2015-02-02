/*

Copyright (c) 1993 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved


Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <tcl.h>

#include "memalloc.h"
#include "printmsg.h"
#include "projview.h"
#include "text_stack.h"
#include "tclutils.h"

#include "stdafx.h"

/****************************************************************************/
int  gomp_SetOutput2Widget()
/****************************************************************************/
{
    int    i;
    int    Total;
    int    TextLen;
    const char *const*Text;
    char *WText;
    char   TempText[BUFF_LEN];
    int    Code;
     
    if(gomp_IsTextStackSet()) {

        Total = gomp_EntriesInTextStack() * 
#if defined(WIN32)
            (gomp_TextStackLineLength() + 2)  +
#else
            (gomp_TextStackLineLength() + 1)  +
#endif
            strlen("lulDisplayTextText {}") + 1;

        WText = gomp_AllocateCharVector(Total);
        strcpy(WText,"lulDisplayTextText {");

        Text = gomp_ReturnTextInStack();

        Total = strlen("lulDisplayTextText {");

        for(i  = gomp_EntriesInTextStack() - 1 ; 
            i >= 0 ; 
            i--) {
            TextLen = strlen(Text[i]);
            strncpy(TempText , Text[i] , BUFF_LEN-1);
            TempText[TextLen] = '\n';
            TextLen++;
            strncpy(&WText[Total],TempText,TextLen);
            Total += TextLen;        
        }

        strcpy(&WText[Total],"}\0");

/* check if graphics is available */
        if(gomp_GetTermType() == GRAPHICS_AVAILABLE) {
            Code = Tcl_GlobalEval(gomp_GetTclInterp(), WText);
            if(Code != TCL_OK) {
                gomp_PrintERROR("can't put text to text output widget");
                return(1);
            }
        }

        free(WText);
        return(0);
    }
    else
        return(1);
}
