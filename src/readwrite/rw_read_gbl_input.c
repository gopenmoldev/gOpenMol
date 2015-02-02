/*

Copyright (c) 1991 - 2004 by:
Leif Laaksonen , Centre for Scientific Computing , ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved

Enhancements 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <stdlib.h>

#if defined(IRIX)
#include <bstring.h>
#endif

#include "gaussian.h"
#include "gomfile.h"
#include "gomstring.h"
#include "printmsg.h"

#include "stdafx.h"

#define GAUSSIAN_BASIS ": gaussian_basis"
#define TAG            ":: tag"

static struct {
    int    Bases;
    int *StartBase;
    char **LineText;
    char *TagText;
} GBasisText;

static int GlobalGBasisCounter;

/*********************************************************************/
int gomp_ReadGaussianBasisLineInput(const char *FileName)
/*********************************************************************/
{

    char InputText[BUFF_LEN];
    char OutText[BUFF_LEN];

    FILE  *Gaussian_f;

    Gaussian_f = fopen(FileName,"r");
    if(Gaussian_f == NULL) {
        sprintf(OutText,"Can't open input file : %s",FileName);
        gomp_PrintMessage(OutText);
        return(1);
    }

/* delete/free all reserved space if this has already been run once */

    (void)gomp_DeleteGaussianBasisLineInput();

/* ................................................................. */

    GlobalGBasisCounter = 0;
    GBasisText.Bases    = 0;

    while(gomp_Fgets(InputText,BUFF_LEN,Gaussian_f) != NULL) {
        gomp_PushBasis2Stack(InputText);
        GlobalGBasisCounter++;

        if(!strncmp(InputText,GAUSSIAN_BASIS,strlen(GAUSSIAN_BASIS))) {
            (void)gomp_PushGBasisStartingPoint(GlobalGBasisCounter - 1);
            GBasisText.Bases++;
            gomp_ReadBasisSetRecords(Gaussian_f);
        }
    }

    fclose(Gaussian_f);

/* push one extra for the end of the labels */
    (void)gomp_PushGBasisStartingPoint(GlobalGBasisCounter);

    return(0);
}

/*********************************************************************/
int gomp_ReadBasisSetRecords(FILE *File_p)
/*********************************************************************/
{
    char InputText[BUFF_LEN];

    while(gomp_Fgets(InputText,BUFF_LEN,File_p) != NULL) {
        gomp_PushBasis2Stack(InputText);
        GlobalGBasisCounter++;

        if(!strncmp(InputText,GAUSSIAN_BASIS,strlen(GAUSSIAN_BASIS))) {
            (void)gomp_PushGBasisStartingPoint(GlobalGBasisCounter - 1);
            GBasisText.Bases++;
            gomp_ReadBasisSetRecords(File_p);
        }
    }

    return(0);
}

/*********************************************************************/
int gomp_PushBasis2Stack(const char * Text)
/*********************************************************************/
{
    if(!GlobalGBasisCounter) {
        GBasisText.LineText    = malloc(sizeof(const char *));
        GBasisText.LineText[0] = malloc(BUFF_LEN);
    }
    else {
        GBasisText.LineText    = realloc(GBasisText.LineText ,
                                                  sizeof(const char *) * 
                                                  (GlobalGBasisCounter + 1));
        GBasisText.LineText[GlobalGBasisCounter] =
            malloc(BUFF_LEN);
    }

    strncpy(GBasisText.LineText[GlobalGBasisCounter] , Text , BUFF_LEN-1);

    return(0);
}

/*********************************************************************/
int gomp_PushGBasisStartingPoint(int Position)
/*********************************************************************/
{
    if(!GBasisText.Bases) {
        GBasisText.StartBase  = malloc(sizeof(int));
    }
    else {
        GBasisText.StartBase  = realloc(GBasisText.StartBase ,
                                                sizeof(int) * 
                                                (GBasisText.Bases + 1));
    }

    GBasisText.StartBase[GBasisText.Bases] = Position;

    return(0);

}
/*********************************************************************/
int gomp_DeleteGaussianBasisLineInput()
/*********************************************************************/
{
    int i;

    if(!GBasisText.Bases && !GlobalGBasisCounter) return(0); 
    /* nothing to free/delete */

    for(i = 0 ; i < GlobalGBasisCounter ; i++)
        free(GBasisText.LineText[i]);
    
    free(GBasisText.LineText);

    GlobalGBasisCounter = 0;

    free(GBasisText.StartBase);

    GBasisText.Bases = 0;

    return(0);
}
/*********************************************************************/
int   gomp_GetNumberOfGaussianBasisSets()
/*********************************************************************/
{
    return(GBasisText.Bases);
}
/*********************************************************************/
const char *gomp_GetGaussianBasisTag(int Position)
/*********************************************************************/
{
    static int  i;
    static char TempText[BUFF_LEN];

    for(i = GBasisText.StartBase[Position] ; 
        i < GBasisText.StartBase[Position + 1] ;
        i++) {
        if(!strncmp(GBasisText.LineText[i],TAG,strlen(TAG))) {
            gomp_CopyString(TempText,GBasisText.LineText[i + 1],BUFF_LEN); 
            return(TempText);
        }
    }

    return(NULL);
}
/*********************************************************************/
const char *gomp_GetGaussianBasisSetEntry(int Position)
/*********************************************************************/
{
    static int   i;
    static int   StatusValue = 0;
    static int   iTemp;
    static char *TempText;

    if(StatusValue) {
        free(TempText);
        StatusValue   = 0;
        TempText      = (char *)NULL;
    }

    for(i = GBasisText.StartBase[Position] ; 
        i < GBasisText.StartBase[Position + 1] ;
        i++) {
        if(TempText == NULL) {
            TempText = malloc(strlen(GBasisText.LineText[i]) + 2);
            sprintf(TempText,"%s\n",GBasisText.LineText[i]);
        }
        else {
            iTemp = strlen(TempText);
            TempText = realloc (TempText , 
                                         (iTemp + strlen(GBasisText.LineText[i]) + 2));
            sprintf(&TempText[iTemp],"%s\n",GBasisText.LineText[i]);
        }
        StatusValue++;
    }

    return(TempText);
}
