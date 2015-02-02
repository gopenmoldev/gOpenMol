/*  

Copyright (c) 1994 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Enhancements 2003 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include <assert.h>
#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>

#if defined(WIN32)
#include <time.h>
#else
#include <sys/types.h>
#include <sys/time.h>
#endif

#include <stdlib.h>

#if defined(IRIX)
#include <sys/procfs.h>
#endif

#include <tcl.h>

/*#include "gomenv.h"*/
#include "gomfile.h"
#include "gomstring.h"
#include "model_file.h"
#include "printmsg.h"
#include "text_stack.h"
#include "tclutils.h"

#include "stdafx.h"

#define RABS(a)    ( ( a ) > 0.0 ? (a) : -(a))

/* Text stack to contain information */
#define  TXT_LINE   256
#define  TXT_DEEP   256

static struct {
    char **Text;               /* contains the stack      */
    int    StackDeep;          /* total number of entries */
    int    StackNow;           /* entries now             */
    int    StackSet;           /* stack control info      */
} OutputStack = { NULL , 0 , 0 , 0};

#if 0
static int   GetInputStackSize(void);
static const char *GetLineFromInputStack(int);
static int   SendInputStack2Tcl(void);
#endif
static int   PushIntoInputStack(const char *);
static int   DeleteInputStack(void);

static struct {
    char **Text;               /* contains the stack      */
    int    StackDeep;          /* total number of entries */
} InputStack = { NULL , 0 };

/*   Utility function to convert character string to lower case  */
/***********************************************************************/
int gomp_String2Lower(char *string )
/***********************************************************************/
{
    int i,j;

    i=0;

    while(string[i] != '\0')
        ++i;
    for( j = 0 ; j < i ; j++)
        string[j]=tolower(string[j]);

    return(0);
}

/***********************************************************************/
int  gomp_Indexo(const char *text ,const char *look)  
    /* imitates the FORTRAN index function */
/***********************************************************************/
{
    int i , i1 , i2;

    i1=strlen(text);
    i2=strlen(look);

    if( i2 > i1) {
        return (0); }

    for(i = 0 ; i < (i1-i2+1)  ; i++) {
        if(strncmp(text+(i),look,i2) == 0) return (i+1);
    }
    return (0);
}

/***********************************************************************/
int gomp_IsStringAFloat(const char *InString)
/* test if input is a float value
   On return:
   0 it is not a float
   >0 it is a const float */
/***********************************************************************/
{

    register int  i,j;
    static   int  Numbers;
    static   int  Characters;
    static   char SaveChar;
    
    j = strlen(InString);

    if(j < 1) return(0);

    Characters = 0;
    Numbers    = 0;

/* -1- */
    if(InString[0] == '+') return(1); /* yes */
    if(InString[0] == '-') return(1); /* yes */
    if(InString[0] == '.') return(1); /* yes */
    

/* -2- */
    for(i = 0 ; i < j ; i++) {
        if(InString[i] == '.') return(1); /* yes */
        if(isalpha(InString[i])) {
            Characters++;
            SaveChar = InString[i];
        }
        if(isdigit(InString[i])) Numbers++;
/* some special cases: '*' and '?' */
        if(InString[i] == '*')  return(0);
        if(InString[i] == '?') return(0);
    }

    if(Numbers && (Characters == 1) ) {
        if(SaveChar == 'e' || SaveChar == 'E') return(1);
    }
/* it's integer */
    if(Numbers && !Characters) return(0);
   
    return(0);
}
/***********************************************************************/
char* gomp_CopyString(char *to, const char *from, size_t size)
/***********************************************************************/
{
    char *const start = to;

    assert(size > 0);

    /* Copy until we reach end of the destination or source buffer. */
    while ( size > 0 && (*to++ = *from++) )
        size--;

    /* There is two possible cases:
     * - We have reached the end of the source buffer.
     *   We have assigned nul byte to *(to - 1).
     *   It makes no harm to do it again.
     * - We have reached the end of the destination buffer.
     *   Assigning something to *to would be a buffer overflow.
     *   We have to assign nul byte to *(to - 1).
     *   That will truncate the text but there is nothing else to do.
     */
    *(to - 1) = '\0';

    return start;
}
/***********************************************************************/
int gomp_PushText2Stack(const char *Prefix,const char *Text)
/***********************************************************************/
{
    int   i;
    int   count;
    char  Temp[TXT_LINE];

    if ( ! gomp_IsTextStackSet() )
        gomp_PrepareTextStack();

    if(Text[0] == '\0') {
        if(gomp_PushText2StackTop(Text))
            return(1);
        return(0);
    }

    strcpy(Temp,Prefix);
    count = strlen(Temp);

    for ( i = 0 ; i <= (int)strlen(Text) ; i++ ) {
        Temp[count++] = Text[i];

        if ( ( Text[i] == '\r' ) || ( Text[i] == '\n' ) ) {
            /* End of line. */
            Temp[count-1] = '\0';
            if ( gomp_PushText2StackTop(Temp) )
                return(1);
            count = 0;
        }
        else if ( Text[i] == '\0' ) {
            /* End of text. */
            Temp[count] = '\0';
            if ( gomp_PushText2StackTop(Temp) )
                return(1);
            return(0);
        }
        else if ( count == TXT_LINE - 1 ) {
            /* Too long line. Wrap. */
            Temp[count] = '\0';
            if ( gomp_PushText2StackTop(Temp) )
                return(1);
            count = 0;
        }
    }

    return(0);
}
/***********************************************************************/
int gomp_PushText2StackTop(const char *Text)
/***********************************************************************/
{
    int i;

    for ( i = OutputStack.StackDeep - 2 ; i >= 0 ; i-- ) 
        strcpy(OutputStack.Text[i+1],OutputStack.Text[i]);

    strcpy(OutputStack.Text[0],Text);

    OutputStack.StackNow++;

    if ( OutputStack.StackNow > OutputStack.StackDeep )
        OutputStack.StackNow = OutputStack.StackDeep;

    return(0);
}
/***********************************************************************/
int gomp_PrepareTextStack()
/***********************************************************************/
{
    int i;

    OutputStack.StackDeep = TXT_DEEP;

    OutputStack.Text = malloc(OutputStack.StackDeep * sizeof(const char *));

    if(OutputStack.Text == NULL) return(1);

    for(i = 0 ; i < OutputStack.StackDeep ; i++) {
        OutputStack.Text[i]    = malloc(TXT_LINE);
        if ( OutputStack.Text[i] == NULL )
            return(1);
        OutputStack.Text[i][0] = '\0';
    }

    OutputStack.StackSet = 1;
    OutputStack.StackNow = 0;
    return(0);
}
/***********************************************************************/
int gomp_IsTextStackSet()
/***********************************************************************/
{

    return(OutputStack.StackSet);

}
/***********************************************************************/
int gomp_EntriesInTextStack()
/***********************************************************************/
{
    return(OutputStack.StackNow);
}

/***********************************************************************/
const char *const*gomp_ReturnTextInStack()
/***********************************************************************/
{
    char **text = OutputStack.Text;
    return((const char*const*)text);
}
/***********************************************************************/
int    gomp_TextStackLineLength()
/***********************************************************************/
{
    return((int)TXT_LINE);
}

/*********************************************************************/
int gomp_StringMatch(const char *Text1 , const char *Text2)
/*********************************************************************/
{
    int i,j;
    int Long1;
    int Long2;
    int Equal;

    Long1 = strlen(Text1);
    Long2 = strlen(Text2);

    Equal = 0;
    for(i = 0 ; i < Long2 ; i++) {
        if(Text2[i] == '$') {
            Equal = i;
            break;
        }
        else
            Equal = i+1;
    }

    j = 0 ;
    for(i = 0 ; i < Long1 ; i++) {

        if(Text2[j] == '\0') {  /* end of string */
            return(0);
        }

        if(Text2[j] != '$') {
            if(Text1[i] != Text2[j])   return(0);
        }
        else {
            j++;
            if(Text2[j] == '\0')       return(0);
            if(Text1[i] != Text2[j])   return(0);
        }
        j++;
    }

    if(Long1 < Equal) 
        return(0);
    else
        return(1);
}

/***********************************************************************/
int PushIntoInputStack(const char *Text)
/***********************************************************************/
{
    if(!InputStack.StackDeep) {
        InputStack.Text    = malloc(sizeof(const char *));
        if(InputStack.Text == NULL) return(1);
        InputStack.Text[0] = malloc(strlen(Text) + 1);
        if(InputStack.Text[InputStack.StackDeep] == NULL) return(1);
        strncpy(InputStack.Text[0] , Text , strlen(Text));
        InputStack.StackDeep++;
    }
    else {
        InputStack.Text = realloc(InputStack.Text , 
                                           (InputStack.StackDeep + 1) * sizeof(const char *));
        if(InputStack.Text == NULL) return(1);
        InputStack.Text[InputStack.StackDeep] = 
            malloc(strlen(Text) + 1);
        if(InputStack.Text[InputStack.StackDeep] == NULL) return(1);
        strncpy(InputStack.Text[InputStack.StackDeep] , Text , strlen(Text));
        InputStack.StackDeep++;
    }

    return(0);
}
/***********************************************************************/
int DeleteInputStack()
/***********************************************************************/
{
    int i;

    if(InputStack.StackDeep) {

        for(i = 0 ; i < InputStack.StackDeep ; i++) {
            free(InputStack.Text[i]);
        }

        free(InputStack.Text);

        InputStack.StackDeep = 0;
    }

    return(0);
}

#if 0
/***********************************************************************/
int GetInputStackSize()
/***********************************************************************/
{
    return(InputStack.StackDeep);
}
/***********************************************************************/
const char *GetLineFromInputStack(int Number)
/***********************************************************************/
{

    if(!InputStack.StackDeep) {
        gomp_PrintERROR("no input lines in input stack");
        return((const char *)NULL);
    }
    if(Number < 1 || Number > InputStack.StackDeep) {
        gomp_PrintERROR("line index out of allowed range");
        return((const char *)NULL);
    }

    return(InputStack.Text[Number - 1]);
}

/***********************************************************************/
int   SendInputStack2Tcl()
/***********************************************************************/
{
    int i;

    for(i = 0 ; i < InputStack.StackDeep ; i++) {
        if(gomp_SendCommand2Parser(InputStack.Text[i])) return(1);
    }

    (void)DeleteInputStack();

    return(0);
}
#endif
/***********************************************************************/
int    gomp_ReadTclInfo2FromFile(FILE *Model_f)
/***********************************************************************/
{
    int   lb;
    char  InputText[BUFF_LEN];


/* TCL - tag */
    lb = 0;
    while(gomp_Fgets(InputText,BUFF_LEN,Model_f) != NULL) {
        if(!strncmp(InputText,"[tcl ** tag ** end]",
                    strlen("[tcl ** tag ** end]"))) {
            sprintf(InputText,"Number of tcl command lines: %d",lb);
            gomp_PrintMessage(InputText);
            return(0);
        }

        if(PushIntoInputStack(InputText)) {
            gomp_PrintERROR("can't push tcl commands into stack");
            (void)DeleteInputStack();
            return(1);
        }
        lb++;
    }

    return(1);
}

/***********************************************************************/
int  gomp_StringTrim(char *Text)
/***********************************************************************/
{
    char *temp;

    if( ! *Text )
        return(0);

    temp = Text + strlen(Text) - 1;
    for( ; temp >= Text ; temp-- ) {
        if( *temp == ' ' ) {
            *temp = '\0';
        } else {
            break;
        }
    }

    return(0);
}
