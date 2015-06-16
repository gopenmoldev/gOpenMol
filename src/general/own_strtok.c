/*

Copyright (c) 1990 - 2005 by:
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

#include <sys/types.h>
#include <stdlib.h>

#include "gomstring.h"

#include "stdafx.h"

#define LBR  '['
#define RBR  ']'
#define LPA  '('
#define RPA  ')'
#define LCB  '{'
#define RCB  '}'
#define QUO  '\"'

static int TokenHook(char , const char *);
static int TokenPause(char);

/*
  This is my own implementation of the "strtok" function.
  Look at the man pages to see how it works.

  Leif Laaksonen
*/

/****************************************************************************/
const char *gomp_StrTok(char *String1,const char *String2)
/*      const char *String1;              string to separate into tokens */
/*      const char *String2;              separtator characters          */
/****************************************************************************/
{
    static int i;               /* loop counter               */
    static int FromLoop = 0;    /* index into the string      */
    static int SaveLen  = 0;    /* length of string           */
    static int StartLoop= 0;
    static int Hit      = 0;             
    register char TempChar;
    static   char *OutLoop;

    Hit = 0;

    if(String1 != NULL) {
        if(*String1 == '\0') return(NULL);
        SaveLen  = strlen(String1);
        if(SaveLen >= BUFF_LEN) 
            return(NULL);

        OutLoop = String1;
        FromLoop = 0;
        StartLoop= 0;
    }
    else {
        if(FromLoop >= SaveLen) 
            return(NULL);
        StartLoop = FromLoop;
    }

    for(i = FromLoop ; i < SaveLen ; i++) 
    {

        if((TempChar = OutLoop[i]) == '\0') 
            return(OutLoop+StartLoop);

        if(TokenPause(TempChar)) 
        {
            Hit = 1;
            continue;
        }

        if(TokenHook(TempChar,String2)) 
        {
            if(Hit) 
            {
                OutLoop[i] = '\0';
                FromLoop   = i + 1;
                return(OutLoop+StartLoop);
            }
            else
                StartLoop++;
        }
        else 
        {
            Hit = 1;
        }
    }
    if(OutLoop[StartLoop] == '\0') 
        return(NULL);
    FromLoop = i + 1;
    return(OutLoop+StartLoop);
}
/****************************************************************************/
int TokenHook(char CheckChar, const char *StringHook)
/****************************************************************************/
{
    register int i;
    register int j;

    j = strlen(StringHook);

    for(i = 0 ; i < j ; i++) 
    {
        if(CheckChar == StringHook[i]) return(1);
    }
    return(0);
}
/****************************************************************************/
int TokenPause(char CheckChar)
/****************************************************************************/
{
    static int InBr  = 0;
    static int InPa  = 0;
    static int InCBr = 0;
    static int InQu  = 0;

    if(CheckChar == LBR) InBr++;    
    if(CheckChar == RBR) InBr--;

    if(CheckChar == LPA) InPa++;    
    if(CheckChar == RPA) InPa--;

    if(CheckChar == LCB) InCBr++;    
    if(CheckChar == RCB) InCBr--;

    if(CheckChar == QUO) InQu = (InQu == 0) ? 1 : 0;    

    return(InBr | InPa | InCBr | InQu);
}
