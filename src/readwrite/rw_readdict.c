/*  

Copyright (c) 1990 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2005 by:
Eero Häkkinen

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "coord_file.h"
#include "gomenv.h"
#include "gomstring.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

#define MAX(a,b)  ( ( a ) > (b) ? (a) : (b))
#define MIN(a,b)  ( ( a ) < (b) ? (a) : (b))
#define Rabs(a)    ( ( a ) > 0 ? (a) : -(a))


static int DictionaryAtoms = 0;     
struct dtype {
    char resnam[MAX_RES_NAME_LEN];
    char atnam[MAX_ATM_NAME_LEN];
    int  type;
    float charge;
} ;

static struct dtype *DictionaryEntries;

static int SetAtomTypeFromDictionary(int);
static int ReadDictionaryEntries(const char *);

/***********************************************************************/
int ReadDictionaryEntries(const char *type_file)
    /* name of the type file to be read */
/***********************************************************************/
{

    int    i;
    char   input[BUFF_LEN];
    char   chelp[BUFF_LEN];
    char   OutText[BUFF_LEN];
    char   TempRes[BUFF_LEN];
    char   TempAtm[BUFF_LEN];
    int    TempType;
    float  TempCharge;
    FILE  *par_p;

/* look only into the data directory */
    if((par_p = fopen(type_file , "r")) == NULL) {

#if defined(WIN32)
        sprintf(chelp,"%s\\%s",gomp_ShowDataDir(),type_file);
#else
        sprintf(chelp,"%s/%s",gomp_ShowDataDir(),type_file);
#endif

        par_p = fopen(chelp,"r");

        if(par_p == NULL) {
            sprintf(OutText,"?Can't open input file : '%s'",chelp);
            (void)gomp_PrintERROR(OutText);
            return(1); }
    }

/*   We are ready now to start reading.
     First line HAS to start with a "*". After that it is free to have or
     not to have a coment. Last line has to be just a "*"  */

    sprintf(OutText,"Reading dictionary file: %s ...",type_file);
    gomp_PrintMessage(OutText);

/* 1 line */
    fgets(input,BUFF_LEN,par_p);   /*    read first line and check for "*" */
    if(input[0] != '*') {
        sprintf(OutText,"\n>>>>>> First line in a dictionary file has to start with a '*'-sign <<<<<<");
        gomp_PrintMessage(OutText);
        fclose(par_p);
        return (1);
    }

    if(DictionaryAtoms) {
        free(DictionaryEntries);
        DictionaryAtoms = 0;
    }

/*  read rest of the file (main loop)    */
    i=0;
    while(fgets(input,BUFF_LEN,par_p) != NULL) {

        if(input[strlen(input) - 1] == '\n')
            input[strlen(input) - 1] = (char)NULL;

        if(!sscanf(input,"%s %s %d %f",TempRes,TempAtm,&TempType,&TempCharge))
            continue;

/* allow for a comment starting with '#' at first or second position */
        if(input[0] == '#' || input[1] == '#')
            continue;

/* memory handling */
        if(!i) {
            DictionaryEntries = malloc(sizeof(*DictionaryEntries));
        }
        else {
            DictionaryEntries = realloc(DictionaryEntries,
                (i + 1) * sizeof(*DictionaryEntries));
        }

        if(strlen(TempRes) > MAX_RES_NAME_LEN) 
            gomp_PrintWARNING("Length of residue name in dictionary file is > MAX_RES_NAME_LEN");
        strncpy(DictionaryEntries[i].resnam,TempRes,MAX_RES_NAME_LEN);

        if(strlen(TempAtm) > MAX_ATM_NAME_LEN) 
            gomp_PrintWARNING("Length of atom name in dictionary file is > MAX_ATM_NAME_LEN");
        strncpy(DictionaryEntries[i].atnam,TempAtm,MAX_ATM_NAME_LEN);

        DictionaryEntries[i].type   = TempType;
        DictionaryEntries[i].charge = TempCharge;

        ++i; 

    }

    DictionaryAtoms = i;
    fclose(par_p);

    return (0);
}

/****************************************************************************/
int SetAtomTypeFromDictionary(int Wstr)  
    /* set atom type from a dictionary file */
/****************************************************************************/
{
    float tot_charge;
    float min_charge;
    int   min_charge_id;
    float max_charge;
    int   max_charge_id;
    float tempcharge;

    register int i,j;
    char OutText[BUFF_LEN];

    static char qresnam[MAX_RES_NAME_LEN + 1];
    static char wresnam[MAX_RES_NAME_LEN + 1];
    static char qatnam[MAX_ATM_NAME_LEN  + 1];
    static char watnam[MAX_ATM_NAME_LEN  + 1];

    for( i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
        strncpy(qatnam,gomp_GetAtomAtmName(Wstr , i) ,MAX_ATM_NAME_LEN);
        gomp_Sign_char(qatnam , MAX_ATM_NAME_LEN);
        strncpy(qresnam,gomp_GetAtomResName(Wstr , i),MAX_RES_NAME_LEN);
        gomp_Sign_char(qresnam , MAX_RES_NAME_LEN);

        for( j = 0 ; j < DictionaryAtoms ; j++) {

            strncpy(wresnam,DictionaryEntries[j].resnam,MAX_ATM_NAME_LEN);
            gomp_Sign_char(wresnam , MAX_ATM_NAME_LEN);
            strncpy(watnam,DictionaryEntries[j].atnam,MAX_RES_NAME_LEN);
            gomp_Sign_char(watnam  , MAX_RES_NAME_LEN);

            if( gomp_CompareStrings(qatnam,watnam  ,MAX_ATM_NAME_LEN) && 
                gomp_CompareStrings(qresnam,wresnam,MAX_RES_NAME_LEN)) {
                gomp_PutAtomType(Wstr , DictionaryEntries[j].type   , i);
                gomp_PutAtomCharge(  Wstr , DictionaryEntries[j].charge , i);
                (void)gomp_AssignAtomInfoFromDictionary( Wstr , i );
                break;
            }
        }
    }

    tot_charge    = 0.0;
    min_charge    = 1.e+20f;
    max_charge    = -1.e+20f;
    min_charge_id = -1;
    max_charge_id = -1;

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {

        tempcharge = gomp_GetAtomCharge(Wstr , i);

        tot_charge += tempcharge;

        if(tempcharge < min_charge) {
            min_charge = tempcharge;
            min_charge_id = i;
        }

        if(tempcharge > max_charge) {
            max_charge = tempcharge;
            max_charge_id = i;
        }
    }

    sprintf(OutText,"***************> Charge analysis (%d)<***************",(Wstr+1));
    gomp_PrintMessage(OutText);
    sprintf(OutText,"Total charge of system : %f\n",tot_charge);
    gomp_PrintMessage(OutText);
    i = max_charge_id;
/*
  sprintf(OutText,"Max charge: >%.4s<>%.4s<>%.4s< %f\n",segment+4*i,
  resnam+4*i,
  atnam+4*i,max_charge);
  PrintMessage(OutText);
  i = min_charge_id;
  sprintf(OutText,"Min charge: >%.4s<>%.4s<>%.4s< %f\n",segment+4*i,
  resnam+4*i,
  atnam+4*i,min_charge);
  PrintMessage(OutText);
*/

    return(0);
}
/****************************************************************************/
int gomp_ImportDictionaryAndApply(const char *FileName)
/****************************************************************************/
{
    int Ret;
    int Wstr;

    if(!gomp_GetNumMolecStructs())
        return(0);

    Wstr = 0;

    if((Ret = ReadDictionaryEntries(FileName)))
        return(Ret);

    if((Ret = SetAtomTypeFromDictionary(Wstr)))
        return(Ret);

    return(0);
}
