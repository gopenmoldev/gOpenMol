/*

Copyright (c) 1993 - 2004 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved



Utility tools

Leif Laaksonen 1993, 1996

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include  <ctype.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <sys/types.h>
#include  <time.h>

#if defined(WIN32)
#include  <winsock.h>
#else
#include  <unistd.h>
#endif

#include <tcl.h>

#include "gomenv.h"
#include "gommain.h"
#include "gomstring.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "tclutils.h"
#include "text_stack.h"

#include "stdafx.h"

#if 0
static int GetTerminalType();
/*int     PushIntoAtomTable(int , const char *, const char * , const char *,
                          const char *, int , const char *, 
                          float , float , float, float , int);
int     RawPushIntoAtomTable(const char * , const char * , const char *,
                             const char * , const char * , 
                             float , float , float, float);
*/
static int DelMolStruct();
#endif

#if 0
static int   FillPerspectiveStructure(void);
static float GetPerspectiveFovy();
static float GetPerspectiveAspect();
static float GetPerspectiveNear();
static float GetPerspectiveFar();
static float GetPerspectiveTranslateX();
static float GetPerspectiveTranslateY();
static float GetPerspectiveTranslateZ();
static int   CheckAtomName(const char *);
#endif

#if 0
static struct {
    float  Fovy;   /* Fovy angle          */
    float  Aspect; /* Aspect ratio        */
    float  Near;   /* Near clipping plane */
    float  Far;    /* Far clipping plane  */
    float  Lookx;
    float  Looky;
    float  Lookz;
} Perspective;
#endif

static struct STARTUP_info {
    char HomeDir[BUFF_LEN];          /* path to home dir */
    char StartUpFile[BUFF_LEN];      /* Start up file    */
    char HostName[BUFF_LEN];         /* Host name        */
} StartUp = {"\0" , STARTUP , "\0" };

static struct NODE_TABLE {
    int Entries;
    char NodeName[BUFF_LEN];
} *NodeTable = { NULL };

/*
  int OpenOpenMolLogFile();
  int CloseOpenMolLogFile();
*/

#if 0
/*****************************************************************************/
int GetTerminalType()
/*****************************************************************************/
{
    char TermString[BUFF_LEN];

    strncpy(TermString,getenv("TERM"),BUFF_LEN-1);

#ifdef DEBUG
    printf("TermString: '%s'\n",TermString);
#endif

    if(TermString[0] == (char)NULL) {
        gomp_PrintMessage("?ERROR - can't get type of terminal (will die now)");
        exit(1);
    }

    return(1);
}

/***************************************************************************/
int  PushIntoAtomTable(int Position , const char *Symbol , const char *Label ,
                       const char *Segment , const char *Residue, int ResidueNumber ,
                       const char *BasisSetTag,
                       float Xc , float Yc , float Zc, float Charge ,
                       int AtomIndex)
/***************************************************************************/
{
    int i;

    if(gomp_GetNumAtomsInMolecStruct(0) == 0) { /* no atoms so far */
        gomp_GetSpaceFor1Atom( 0 );
        i = 0;
    }

    else if(Position == gomp_GetNumAtomsInMolecStruct(0)) {

        gomp_GetSpaceFor1Atom( 0 );

        i = gomp_GetNumAtomsInMolecStruct(0) - 1;
    }
    else 
        i = Position;

/*     (void)gomp_PutAtomAtmName( 0 , Symbol , i); */
    (void)gomp_PutAtomAtmName( 0 , Label , i);
    (void)gomp_PutAtomBasisSetTag( 0 , BasisSetTag , i);

    (void)gomp_PutAtomXCoord( 0 , Xc , i);
    (void)gomp_PutAtomYCoord( 0 , Yc , i);
    (void)gomp_PutAtomZCoord( 0 , Zc , i);

    (void)gomp_PutAtomResNum1(0 , ResidueNumber , i);
    (void)gomp_PutAtomResNum2(0 , ResidueNumber , i);

    (void)gomp_PutAtomBValue(0 , 0.0 , i);

    (void)gomp_PutAtomNucCharge( 0 , Charge , i);

    if(Residue[0] == '\0')
        (void)gomp_PutAtomResName( 0 ,DEFAULT_RESIDUE_NAME , i);
    else
        (void)gomp_PutAtomResName( 0 ,Residue , i);

    if(Segment[0] == '\0')
        (void)gomp_PutAtomSegName( 0 , DEFAULT_SEGMENT_NAME, i);
    else 
        (void)gomp_PutAtomSegName( 0 , Segment, i);


    if(AtomIndex > 0) {
        (void)SaveAtomLabelIndex(0 , i , (AtomIndex - 1));
        return(0);
    }
    else {
        (void)SaveAtomLabelIndex(0 , i ,  -1);
        return(1);
    }

    return(gomp_GetNumAtomsInMolecStruct(0));

}

/***************************************************************************/
int  RawPushIntoAtomTable(const char *Symbol , const char *Label ,
                          const char *Segment , const char *Residue,
                          const char *BasisSetTag,
                          float Xc , float Yc , float Zc, float Charge)
/***************************************************************************/
{
    int i;


    gomp_GetSpaceFor1Atom( 0 );

    i = gomp_GetNumAtomsInMolecStruct(0) - 1;

    (void)gomp_PutAtomAtmName( 0 , Symbol , i);
    (void)gomp_PutAtomBasisSetTag( 0 , BasisSetTag , i);

    (void)gomp_PutAtomXCoord( 0 , Xc , i);
    (void)gomp_PutAtomYCoord( 0 , Yc , i);
    (void)gomp_PutAtomZCoord( 0 , Zc , i);
    (void)gomp_PutAtomNucCharge( 0 , Charge , i);

    if(Residue[0] == '\0')
        (void)gomp_PutAtomResName( 0 ,DEFAULT_RESIDUE_NAME , i);
    else
        (void)gomp_PutAtomResName( 0 ,Residue , i);

    if(Segment[0] == '\0')
        (void)gomp_PutAtomSegName( 0 , DEFAULT_SEGMENT_NAME, i);
    else 
        (void)gomp_PutAtomSegName( 0 , Segment, i);

    return(gomp_GetNumAtomsInMolecStruct(0));

}
#endif
/************************************************************************/
int gomp_GetEnv(int argc , const char *argv[])   
    /* go and hunt for environment parameters */
/************************************************************************/
{
    (void)gomp_StartUpStuff(argc, argv);

    strncpy(StartUp.HomeDir,gomp_ShowHomeDir(),BUFF_LEN-1);

    strncpy(StartUp.HostName,gomp_LongHostName(),BUFF_LEN-1);

    return(0);
}

/*
  There are two "startup" files used by gOpenMol:

  (1) User startup file, located in the login directory
  (2) System startup file, located in the gopenmol data directory

  The name of the files is the same.

*/
/************************************************************************/
int gomp_OpenStartFile()
/************************************************************************/
{
    char  OutText[BUFF_LEN];
    FILE  *StartFP;
    int   iSystem;
    int   iUser;

    if(StartUp.StartUpFile[0] == '\0') {
        gomp_PrintMessage("?WARNING - Start up file not defined");
        return(1);
    }

/* read first the system startup file ... */
    iSystem = 1;

#if defined(WIN32)
    sprintf(OutText,"%s\\%s",gomp_ShowDataDir(),StartUp.StartUpFile);
#else
    sprintf(OutText,"%s/%s",gomp_ShowDataDir(),StartUp.StartUpFile);
#endif

    StartFP = fopen(OutText,"r");

    if( StartFP == NULL) {
        sprintf(OutText,"?WARNING - can't open system start up file: '%s'",StartUp.StartUpFile);
        gomp_PrintMessage(OutText);
        iSystem = 1;
    }
    else {
        fclose(StartFP);
/* send file to Tcl parser ... */
        if(gomp_SendFile2TclParser(OutText) != TCL_OK) {
            gomp_PrintMessage(Tcl_GetStringResult(gomp_GetTclInterp()));
            gomp_PrintMessage("can't parse the input tcl file(s)");
        }
        iSystem = 0;
    }

/* read then user startup file ... */
    iUser = 1;

/* the windows (win32) user don't really have a home directory ... */

#if !defined(WIN32)
    sprintf(OutText,"%s/%s",gomp_ShowHomeDir(),StartUp.StartUpFile);

    StartFP = fopen(OutText,"r");

    if( StartFP == NULL) {
        iUser = 1;
    }
    else {
        fclose(StartFP);
/* send file to Tcl parser ... */
        (void)gomp_SendFile2TclParser(OutText);
        iUser = 0;
    }
#endif

    return(0);
}


/***************************************************************************/
int gomp_SplitString(char *InString , const char *Separate , const char *OutVector[BUFF_LEN])
/***************************************************************************/
{
    int Loop;
    int i;

/* the '#' marks the beginning of a comment */ 
    for(i = 0 ; i < (int)strlen(InString) ; i++) {
        if(InString[i] == '#') {
            InString[i] = '\0';
            break;
        }
    }

    Loop = 0;

    OutVector[Loop] = STRTOK(InString , Separate);

    if(OutVector[0] == NULL) return(0);

    Loop++;

    while ((OutVector[Loop] = STRTOK(NULL , Separate)) != '\0') Loop++;

    return(Loop);
}
/*   Utility function to convert character string to lower case  */
#if 0
#ifndef SCARECROW
/***********************************************************************/
void gomp_toller( char string[] )
/*
  char string[];
*/
/***********************************************************************/
{
    int i,j;

    i=0;

    while(string[i] != '\0')
        ++i;
    for( j = 0 ; j < i ; j++)
        string[j]=tolower(string[j]);
}
#endif
#endif
/***************************************************************************/
int gomp_PushNode2Stack(const char *NodeName)
/***************************************************************************/
{
    static int Loop = 0;
    char OutText[BUFF_LEN];

    if(!Loop) {
        NodeTable = (struct NODE_TABLE *) 
            malloc(sizeof(*NodeTable));
        gomp_CopyString(NodeTable[0].NodeName,NodeName,BUFF_LEN);
    }
    else {
        NodeTable = (struct NODE_TABLE *) 
            realloc(NodeTable, 
                    (Loop + 1) * sizeof(*NodeTable));
        gomp_CopyString(NodeTable[Loop].NodeName,NodeName,BUFF_LEN);
    }

    if(NodeTable == NULL) {
        gomp_PrintMessage("?ERROR - can't allocate space in gomp_PushNode2Stack");
        return(1);
    }

    sprintf(OutText,"Node '%s' defined",NodeTable[Loop].NodeName);
    gomp_PrintMessage(OutText);

    Loop++;

    NodeTable[0].Entries = Loop;

    return(0);
}
#if 0
/***************************************************************************/
static int PushNodes2Structure(const char *NodeName)
/***************************************************************************/
{
    int   Loop = 0;
    char  OutText[BUFF_LEN];
    const char *OutVec[BUFF_LEN];
    int   GotIt;

    if(NodeTable != NULL) free(NodeTable);

    GotIt = gomp_SplitString(OutText," \n,",OutVec);

    if(!GotIt) {
        NodeTable = NULL;
        return(0);
    }

    for(Loop = 0 ; Loop < GotIt ; Loop++) {

        if(!Loop) {
            NodeTable = (struct NODE_TABLE *) 
                malloc(sizeof(*NodeTable));
            gomp_CopyString(NodeTable[0].NodeName,NodeName,BUFF_LEN);
        }
        else {
            NodeTable = (struct NODE_TABLE *) 
                realloc(NodeTable, 
                        (Loop + 1) * sizeof(*NodeTable));
            gomp_CopyString(NodeTable[Loop].NodeName,NodeName,BUFF_LEN);
        }

        if(NodeTable == NULL) {
            gomp_PrintMessage("?ERROR - can't allocate space in gomp_PushNode2Stack");
            return(1);
        }

        sprintf(OutText,"Node '%s' defined",NodeTable[Loop].NodeName);
        gomp_PrintMessage(OutText);

    }

    NodeTable[0].Entries = GotIt;

    return(0);
}
/***************************************************************************/
static const char *GetNodesFromStructure()
/***************************************************************************/
{
    int   Loop = 0;
    static char *OutVec;

    if(NodeTable == NULL) return(NULL);

    for(Loop = 0 ; Loop < NodeTable[0].Entries ; Loop++) {

        if(!Loop) {
            OutVec = malloc(sizeof(NodeTable[0].NodeName) + 1);
            strcpy(OutVec,NodeTable[0].NodeName);
            strcat(OutVec,"\n");
        }
        else {
            OutVec = realloc(OutVec, 
                (sizeof(OutVec) + sizeof(NodeTable[Loop].NodeName) + 1));
            strcat(OutVec,NodeTable[Loop].NodeName);
            strcat(OutVec,"\n");
        }


        if(OutVec == NULL) {
            gomp_PrintMessage("?ERROR - can't allocate space in GetNodesFromStructure");
            return(NULL);
        }

    }

    return(OutVec);
}
/*
  Check the atom name.
  On return:    0 OK
  != 0 Not an allowed name
*/
/***************************************************************************/
int CheckAtomName(const char *InName)
/***************************************************************************/
{
    int i;
    char OutText[BUFF_LEN];

    if((i = strlen(InName)) > 2) {
        sprintf(OutText,"?Atom name '%s' is longer than 2 characters",InName);
        gomp_PrintERROR(OutText);
        return(-1);
    }

    if(i < 1) {
        gomp_PrintMessage("?ERROR - NULL string supplied to atomic symbol test");
        return(-1);
    }

    i = gomp_MatchAtom(InName);

    return(i);
}
/***************************************************************************/
int DelMolStruct()
/***************************************************************************/
{

    (void)gomp_DeleteMolecStructures();

    return(0);
}
/***************************************************************************/
int   FillPerspectiveStructure()
/***************************************************************************/
{
    return(0);
}
/***************************************************************************/
float GetPerspectiveFovy()
/***************************************************************************/
{
    return(0.0);
}
/***************************************************************************/
float GetPerspectiveAspect()
/***************************************************************************/
{
    return(0.0);
}
/***************************************************************************/
float GetPerspectiveNear()
/***************************************************************************/
{
    return(0.0);
}
/***************************************************************************/
float GetPerspectiveFar()
/***************************************************************************/
{
    return(0.0);
}
/***************************************************************************/
float GetPerspectiveTranslateX()
/***************************************************************************/
{
    return(0.0);
}
/***************************************************************************/
float GetPerspectiveTranslateY()
/***************************************************************************/
{
    return(0.0);
}
/***************************************************************************/
float GetPerspectiveTranslateZ()
/***************************************************************************/
{
    return(0.0);
}
#endif
