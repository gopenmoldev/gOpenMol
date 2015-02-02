/*  

Copyright (c) 1994 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  
Enhancements 2003 - 2005 by:
Eero Häkkinen

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <math.h>
#include <string.h>
#include <ctype.h>

#if !defined(WIN32)
#include <sys/types.h>
#include <sys/time.h>
/*#include <sys/resource.h>*/
#endif

#if defined(IRIX)
#include <sys/procfs.h>
#endif

#include <stdlib.h>

#include "atom_param.h"
#include "gomenv.h"
#include "gomstring.h"
#include "molecule.h"
#include "memalloc.h"
#include "printmsg.h"

#include "stdafx.h"

AtomTypes_t gomp_AtomTypes;

/*
  Format of the atom parameter file. Most values are the same as in
  the parameter file in QUANTA (Polygen) or the equivalent parameter
  file in CHARMM.
 
  Leif Laaksonen  1989 - 2000

  I will also include some own parameters when I have the time.

  All this data is taken from APPENDIX B of the QUANTA Reference Manual 1990.
  ===========================================================================

  The titles of the columns in the file are:

  (1)    (2)     (3)     (4)     (5)     (6)    (7)   (8)   (9)    (10)  (11)  (12)
  Number bndrad  vdwrad  plurad  global  emin   rmin  patom hbond  atype mass   ncharge
  1     0.4000  0.9500  0.100     F    -0.0498  0.8  0.044  D     H     1.008  1
  .
  .
  .

  (1):    Atom type number (type).
  (2):    Bonding radius (bndrad).
  (3):    van der Waals spheres radius (vdwrad).
  (4):    Sphere radius (plurad) in plots. (Not used in gOpenMol).
  (5):    Value (global) is used in global search for bonds.
  (6,7):  Value (emin and rmin) used in the calculation of van der
  Waals and electrostatic energy. 
  (8):    Atom polarizabilities (patom).
  (9):    The atom is either a hydrogen bond acceptor (A), donor (D or E), or
  not hydrogen bonded (N)
  (10):   The CHARMm atom type (atype).
  (11):   Atom mass (mass).

  The parameters in QUANTA 3.0 are similar to those used in Vesion 21 of CHARMm.
  The emin rmin values differ for HA and all metals; bndrad differs for 
  hydrogens. Atoms with Mxx are metals; atoms with Xxx are halogens.

  If you wish to extend the parameter file to include your own atom types, it is
  recommended to use atom numbers greater than 300.

*/

/* F U N C T I O N S  T O  H A N D L E  T H E  A T O M S */

/***********************************************************************/
int gomp_ReadAtomParams(const char *File)
/*    const char *File;                      name of the par file to be read */
/***********************************************************************/
{


    int  i;
    char input[BUFF_LEN];
    char temp1[BUFF_LEN];
    FILE *File_p;

    char chelp[BUFF_LEN];
    char OutText[BUFF_LEN];

    if((File_p = fopen(File,"r")) == NULL) {

#if defined(WIN32)
        sprintf(chelp,"%s\\%s",gomp_ShowDataDir(),File);
#else
        sprintf(chelp,"%s/%s",gomp_ShowDataDir(),File);
#endif

        File_p = fopen(chelp,"r");
        if(File_p == NULL) {
            sprintf(OutText,"**** ERROR. Can't open input file : '%s'",chelp);
            (void)gomp_PrintMessage(OutText);
            return(1); }
    }

/*   Check if there is already a structure */
    if(gomp_AtomTypes.params) {
        free(gomp_AtomTypes.AtomParams);
    }

/*   We are ready now to start reading.
     First line HAS to start with a "*". After that it is free to have or
     have not a coment. Last line has to be just a "*"  */

#ifdef DEBUG
    sprintf(OutText,"Reading atom parameter file: %s ...",chelp);
    gomp_PrintMessage(OutText);
#endif

/* 1 line */
    fgets(input,BUFF_LEN,File_p);   /*    read first line and check for "*" */
    if(input[0] != '*') {
        gomp_PrintERROR("?Atom parameter file - ERROR");
        gomp_PrintMessage("$ First line in a par-file has to start with a '*'-sign");
        return (1); }

/*  read rest of the file (main loop)    */
    i=0;
    while(fgets(input,BUFF_LEN,File_p) != NULL) {

/* allow for a comment starting with '#' at first or second position */
        if(input[0] == '#' || input[1] == '#')
            continue;

        if(!i)
            gomp_AtomTypes.AtomParams = gomp_AllocateVoidVector(sizeof(AtomData));
        else
            gomp_AtomTypes.AtomParams = gomp_ReallocateVoidVector(
                gomp_AtomTypes.AtomParams,
                ((i + 1) * sizeof(AtomData)));


        if(!sscanf(input,
                   "%d %f %f %f %c %f %f %f %c %s %f %d\n",
                   &gomp_AtomTypes.AtomParams[i].type,
                   &gomp_AtomTypes.AtomParams[i].bndrad,
                   &gomp_AtomTypes.AtomParams[i].vdwrad,
                   &gomp_AtomTypes.AtomParams[i].plurad,
                   &gomp_AtomTypes.AtomParams[i].global,
                   &gomp_AtomTypes.AtomParams[i].emin,
                   &gomp_AtomTypes.AtomParams[i].rmin,
                   &gomp_AtomTypes.AtomParams[i].patom,
                   &gomp_AtomTypes.AtomParams[i].hbond,
                   temp1,
                   &gomp_AtomTypes.AtomParams[i].mass,
                   &gomp_AtomTypes.AtomParams[i].ncharge))
            continue;

        gomp_CopyString(
            gomp_AtomTypes.AtomParams[i].atype,temp1,
            MAX_ATM_NAME_LEN );
        gomp_AtomTypes.AtomParams[i].cnct = 0;  /* not used now */
        ++i; 
    }

    gomp_AtomTypes.params = i;

#ifdef DEBUG   
    for(j = 0 ; j < i; j++) {

        printf("%d ",gomp_AtomTypes.AtomParams[j].type);
        printf("%f ",gomp_AtomTypes.AtomParams[j].bndrad);
        printf("%f ",gomp_AtomTypes.AtomParams[j].vdwrad);
        printf("%f ",gomp_AtomTypes.AtomParams[j].plurad);
        printf("%c ",gomp_AtomTypes.AtomParams[j].global);
        printf("%f ",gomp_AtomTypes.AtomParams[j].emin);
        printf("%f ",gomp_AtomTypes.AtomParams[j].rmin);
        printf("%f ",gomp_AtomTypes.AtomParams[j].mass);
        printf("%d ",gomp_AtomTypes.AtomParams[j].cnct);
        printf("%c ",gomp_AtomTypes.AtomParams[j].hbond);
        printf("%.4s\n",gomp_AtomTypes.AtomParams[j].atype);
    } 
#endif
    fclose(File_p);

#ifdef DEBUG
    gomp_PrintMessage(" done\n");
#endif 

    return (0);
}


/***********************************************************************/
int   gomp_GetNumberOfAtomParameters()
/***********************************************************************/
{
    return(gomp_AtomTypes.params);
}

/***********************************************************************/
int   gomp_GetAtomPointerForIndex(int Point)
/***********************************************************************/
{
    register int i;

    for(i = 0 ; i < gomp_GetNumberOfAtomParameters() ; i++) {

        if(Point == gomp_AtomTypes.AtomParams[i].type)
            return(i);
    }

    gomp_PrintERROR(" ** Can't find the given atom type in the atom type file");
    return(-1);

}
/***********************************************************************/
int   gomp_GetAtom_type(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].type);
}
/***********************************************************************/
float gomp_GetAtom_bndrad(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].bndrad);
}
/***********************************************************************/
float gomp_GetAtom_vdwrad(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].vdwrad);
}
/***********************************************************************/
float gomp_GetAtom_plurad(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].plurad);
}
/***********************************************************************/
char  gomp_GetAtom_global(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].global);
}
/***********************************************************************/
float gomp_GetAtom_emin(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].emin);
}
/***********************************************************************/
float gomp_GetAtom_rmin(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].rmin);
}
/***********************************************************************/
float gomp_GetAtom_patom(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].patom);
}
/***********************************************************************/
float gomp_GetAtom_mass(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].mass);
}
/***********************************************************************/
int   gomp_GetAtom_cnct(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].cnct);
}
/***********************************************************************/
char  gomp_GetAtom_hbond(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].hbond);
}
/***********************************************************************/
const char *gomp_GetAtom_atype(int Point)
/***********************************************************************/
{
    return(gomp_AtomTypes.AtomParams[Point].atype);
}
/***********************************************************************/
int gomp_GetAtomPointerForSymbol(const char *AtmSymbol)
/***********************************************************************/
{
    register int   i;
    static   char  Text[BUFF_LEN];

    for(i = 0 ; i < gomp_GetNumberOfAtomParameters() ; i++) {

        strncpy(Text,gomp_GetAtom_atype(i),MAX_ATM_NAME_LEN);

        if(!strncmp(AtmSymbol,Text,MAX_ATM_NAME_LEN) )
            return(i);
    }

    if( i < 1) {
        gomp_PrintWARNING("No parameters are available (read the first)");
        return(-1);
    }
    else {
        sprintf(Text,"** Can't find the given atom symbol '%s' in the atom type file",AtmSymbol);
        gomp_PrintERROR(Text);
        return(-1);
    }

}
