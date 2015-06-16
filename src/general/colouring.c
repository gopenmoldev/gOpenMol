/*  

Copyright (c) 1994 - 2005 by:
Leif Laaksonen , Center for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  

Enhancements 2003 - 2004 by:
Eero Häkkinen
*/

/*

Read colour table to be used by SCARECROW/gOpenMol

Leif Laaksonen 1995

RGB to grayscale:

luminosity = .299 red + .587 green + .114 blue

*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <tcl.h>

#include "colouring.h"
#include "gomenv.h"
#include "gommain.h"
#include "gomstring.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"

#include "stdafx.h"

#define COL_TBL_INC        50
#define MAX_COL_LEN       100   /* max num of characters for colour name */
#define COL_FILE_LINE_LEN 132   /* max number of characters on a line    */
#define COLOR_MATCH       1.0e-4

#define RABS(a)    ( ( a ) > 0.0 ? (a) : -(a))

static int col_tbl_max  = 0;  /* number of entries in the colour table */
static int col_tbl_high = 0;

typedef struct SCARE_COL_TABLE {
    float red;      /* read value */
    float green;    /* green value */
    float blue;     /* blue value */
    char name[MAX_COL_LEN];
} SCARE_col_table;

static SCARE_col_table *col_table;  /* colour table pointer */

static int PushColTable(int , int , int , const char *);

static void rainbow(double  , double  , double  ,
                    double  *, double  *, double  *);

static void dhsv2rgb(double  , double  , double  ,
                     double  *, double  *, double  *);

static int gomp_ColourDisplayType = 1; /* !=0 is colour, ==0 is grayscale */

/***************************************************************************/
int gomp_ReadColourTable(const char *tbl_name)
/*        *tbl_name;                  name of colour table file to be read */
/***************************************************************************/
{
    int i;
    int ColLoop;
    int Cr,Cg,Cb;
    int print_input;
    FILE *File_p;
    char chelp[BUFF_LEN];
    char inputl[COL_FILE_LINE_LEN];
    char OutText[BUFF_LEN];
     
/* look only into the data directory */
    if((File_p = fopen(tbl_name,"r")) == NULL) {

#if defined(WIN32)
        sprintf(chelp,"%s\\%s",gomp_ShowDataDir(),tbl_name);
#else
        sprintf(chelp,"%s/%s",gomp_ShowDataDir(),tbl_name);
#endif

        File_p = fopen(chelp,"r");
        if(File_p == NULL) {
            (void)gomp_PrintERROR(" ");
            sprintf(OutText,"?Can't open input file : '%s'",chelp);
            (void)gomp_PrintMessage(OutText);
            return(1); }
    }

    print_input=0;

/*   
     We are ready now to start reading.
     First line is a comment line. 
*/

#ifdef DEBUG
    printf("    Reading colour table file: %s ...\n",tbl_name);
#endif

/* calculate first number of entries in the file */
    ColLoop = 0;
    while(fgets(inputl,COL_FILE_LINE_LEN,File_p) != NULL)  ColLoop++;
    rewind(File_p);

    ColLoop--;

/* 1 line is a comment */
    if (fgets(inputl,COL_FILE_LINE_LEN,File_p) == NULL) {   /*  */
      gomp_PrintERROR("? Colour table file - ERROR");
      gomp_PrintMessage("$ File is empty");
      return(1);
    }
/*  read rest of the file (main loop)    */

    col_tbl_max = 0;

    for(i = 0 ; i < ColLoop ; i++) {

      if (fgets(inputl,COL_FILE_LINE_LEN,File_p) == NULL) {
	gomp_PrintERROR("? Colour table file - ERROR");
	gomp_PrintMessage("$ Colour entry missing");
	return(1);
      }

        *chelp = '\0';

        if(sscanf(inputl,"%d %d %d %s %s",&Cr,&Cg,&Cb,OutText,chelp) < 4) {
            gomp_PrintERROR(" ");
            gomp_PrintMessage("?In colour table set up ");
            return(1);
        }

        if(*chelp != '\0')
            sprintf(OutText,"%s %s",OutText,chelp);

        if(PushColTable(Cr , Cg , Cb , OutText)) {
            gomp_PrintERROR(" ");
            gomp_PrintMessage("?Can't push colour into colour table");
            return(1);
        }

    }


    return (0);
}
/***************************************************************************/
int  PushColTable(int Cr , int Cg , int Cb , const char *chelp)   
    /* push new colour into table */
/***************************************************************************/
{

    if(!col_tbl_high) { /* first time down here , grab space */
        col_table = malloc(COL_TBL_INC * sizeof(*col_table));
        col_tbl_high = COL_TBL_INC;

        if(col_table == NULL) {
            gomp_PrintMessage("?ERROR - can't allocate space for colour table");
            gomp_PrintEXIT("Can't continue");
        }
    }

    if(col_tbl_max == col_tbl_high) {
        col_tbl_high += COL_TBL_INC;
        col_table = (SCARE_col_table *) 
            realloc(col_table , 
                    col_tbl_high * sizeof(SCARE_col_table));

        if(col_table == NULL) {
            gomp_PrintMessage("?ERROR - can't reallocate space for colour table");
            gomp_PrintEXIT("Can't continue");
        }
    }
      
    col_table[col_tbl_max].red    = ((float)Cr)/255.0;
    col_table[col_tbl_max].green = ((float)Cg)/255.0;
    col_table[col_tbl_max].blue = ((float)Cb)/255.0;     

    strncpy(col_table[col_tbl_max].name,chelp,MAX_COL_LEN - 1);

    col_tbl_max++;

    return(0);
}

/*
  Check the colour name against names in the colour table
  If there is a hit return the Red, Green, Blue components.
*/

/*************************************************************************/
int gomp_ColourName2RGBSilent(const char *colour, float *r, float *g , float *b)      
    /* on return = 0 ok != 0 not in table */
/*************************************************************************/
{
    char chelp1[MAX_COL_LEN];
    char chelp2[MAX_COL_LEN];
    char ch;
    int  i;
    unsigned int  ir, ig, ib;
    static float sr = 1.0,sg = 1.0,sb = 1.0;

/* default colour is white */
    switch ( *colour ) {
    case '\0':
        *r = sr;
        *g = sg;
        *b = sb;
        return 0;
    case '#':
        {
            size_t clen = ( strlen(colour) - 1 ) / 3;
            if ( 1 <= clen && clen <= 4 ) {
                char format[] = "#%2x%2x%2x%c";
                sprintf(format, "#%%%zux%%%zux%%%zux%%c", clen, clen, clen);
                if ( sscanf(colour, format, &ir, &ig, &ib, &ch) == 3 ) {
                    /**
                     * Set max to 0xf, 0xff, 0xfff or to 0xffff
                     * according to a channel length.
                     */
                    unsigned short int max =
                        (((unsigned short int)1) << ( 4 * clen )) - 1;
                    *r = 1.0f * ir / max;
                    *g = 1.0f * ig / max;
                    *b = 1.0f * ib / max;
                    return 0;
                }
            }
        }
        break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
        if ( sscanf(colour, " %u %u %u %c", &ir, &ig, &ib, &ch) == 3 ) {
            if ( ir <= 255 && ig <= 255 && ib <= 255 ) {
                *r = 1.0f * ir / 255;
                *g = 1.0f * ig / 255;
                *b = 1.0f * ib / 255;
                return 0;
            }
            else
                return 1;
        }
        /* Fall through. */
    default:
        if ( sscanf(colour, " %f %f %f %c", r, g, b, &ch) == 3 ) {
            if ( 0.0 <= *r && *r <= 1.0 &&
                 0.0 <= *g && *g <= 1.0 &&
                 0.0 <= *b && *b <= 1.0 )
                return 0;
            else
                return 1;
        }
    }

/* check if the name is in the colour table */
    strncpy(chelp1,colour,MAX_COL_LEN);
    gomp_String2Lower(chelp1);

    for(i = 0 ; i < gomp_GetColourTableLength() ; i++) {
        strncpy(chelp2, col_table[i].name, MAX_COL_LEN);
        gomp_String2Lower(chelp2);
        if ( strncmp(chelp1, chelp2, strlen(chelp2)) == 0 ) {
            /* yes it is */
            *r = sr = col_table[i].red;
            *g = sg = col_table[i].green;
            *b = sb = col_table[i].blue;
            return 0;
        }
    } 

    return 1;
}

/*************************************************************************/
int gomp_ColourName2RGB(const char *colour, float *r, float *g , float *b)      
    /* on return = 0 ok != 0 not in table */
/*************************************************************************/
{
    if ( gomp_ColourName2RGBSilent(colour,r,g,b) == 0 )
        return(0);

    gomp_PrintWARNING("Problems assigning a colour");
    gomp_FormatMessage("?Unknown colour '%s' (will put it to white) ",colour);

    *r = *g = *b = 1.0f;

    return(1);
}
#if 0
/*************************************************************************/
int ColourName2RGBint(const char *colour, int *r, int *g , int *b)      
    /* on return = 0 ok != 0 not in table */
/*************************************************************************/
{
    char chelp1[MAX_COL_LEN];
    char chelp2[MAX_COL_LEN];
    char OutText[BUFF_LEN];
    float RedC,GreenC,BlueC;
    int   IRedC,IGreenC,IBlueC;
    int i;
    static int sr = 255,sg = 255,sb = 255;

/* default colour is white */
    if(colour[0] == '\0') {
        *r = sr;
        *g = sg;
        *b = sb;
        return(0);
    }

/* the colour name can also be a rgb coded name like "0.123 0.567 0.899"
   or "123 45 255". */

    if(gomp_IsStringAFloat(colour)) {
        sscanf(colour,"%f %f %f",&RedC,&GreenC,&BlueC);
        if(RedC   > 1.005 ||
           GreenC > 1.005 ||
           BlueC  > 1.005) { /* it's an integer list */
            sscanf(colour,"%d %d %d",&IRedC,&IGreenC,&IBlueC);
            *r = IRedC;
            *g = IGreenC;
            *b = IBlueC;
        }
        else {
            *r = (int)(RedC   * 255.0);
            *g = (int)(GreenC * 255.0);
            *b = (int)(BlueC  * 255.0);
        }
        return(0);
    }

    strncpy(chelp1,colour,MAX_COL_LEN);
    (void)gomp_String2Lower(chelp1);

/* check if the name is in the colour table */

    for(i = 0 ; i < gomp_GetColourTableLength() ; i++) {
        strncpy(chelp2,col_table[i].name,MAX_COL_LEN);
        (void)gomp_String2Lower(chelp2);
        if(gomp_Indexo(chelp1,chelp2) == 1) { /* yes it is */
            *r = sr = (int)(col_table[i].red   * 255.0);
            *g = sg = (int)(col_table[i].green * 255.0);
            *b = sb = (int)(col_table[i].blue  * 255.0);
            return(0);
        }
    } 

    gomp_PrintWARNING("Problems assigning a colour");
    sprintf(OutText,"?Unknown colour '%s' (will put it to white) ",colour);
    gomp_PrintMessage(OutText);

    *r = 255;
    *g = 255;
    *b = 255;
    return(1);
}
/*************************************************************************/
int gomp_ReturnColourFromIndex(int Which, char *Name , 
                               float *RedC , float *GreenC , float *BlueC)
/*************************************************************************/
{

    if(Which >= gomp_GetColourTableLength()) return(1);

    gomp_CopyString(Name,col_table[Which].name,BUFF_LEN);
    *RedC   = col_table[Which].red;
    *GreenC = col_table[Which].green;
    *BlueC  = col_table[Which].blue;

    return(0);
}
#endif
/*************************************************************************/
int gomp_GetColourTableLength()
/*************************************************************************/
{
    return(col_tbl_max);
}

/***********************************************************************/
void gomp_PreRainbow(double hh,float *rr,float *gg,float *bb)
/***********************************************************************/
{
    static double r,g,b,s=1.0,v=1.0;

/* truncate it to inside interval 0 ... 1 and reverse the colouring */

    if(hh > 1.0) hh = 1.0;
    if(hh < 0.0) hh = 0.0;

    hh = 1. - hh;   /* it's now from blue to red ... */

    rainbow(hh,s,v,&r,&g,&b);

    *rr = (float)r;
    *gg = (float)g;
    *bb = (float)b;
}     

/*   rainbow(h, s, v, r, g, b)
     double h, s, v, *r, *g, *b;
 
     This routine computes colors suitable for use in color level plots.
     Typically s=v=1 and h varies from 0 (red) to 1 (blue) in
     equally spaced steps.  (h=.5 gives green; 1<h<1.5 gives magenta.)
     To convert for frame buffer, use   R = floor(255.999*pow(*r,1/gamma))  etc.
     To get tables calibrated for other devices or to report complaints,
     contact  Eric Grosse   research!ehg    201-582-5828.
*/

static const double huettab[] = {
    0.0000, 0.0062, 0.0130, 0.0202, 0.0280, 0.0365, 0.0457, 0.0559, 0.0671, 0.0796,
    0.0936, 0.1095, 0.1275, 0.1482, 0.1806, 0.2113, 0.2393, 0.2652, 0.2892, 0.3119,
    0.3333, 0.3556, 0.3815, 0.4129, 0.4526, 0.5060, 0.5296, 0.5501, 0.5679, 0.5834,
    0.5970, 0.6088, 0.6191, 0.6281, 0.6361, 0.6430, 0.6490, 0.6544, 0.6590, 0.6631,
    0.6667, 0.6713, 0.6763, 0.6815, 0.6873, 0.6937, 0.7009, 0.7092, 0.7190, 0.7308,
    0.7452, 0.7631, 0.7856, 0.8142, 0.8621, 0.9029, 0.9344, 0.9580, 0.9755, 0.9889,
    1.0000
};
/* computed from the FMC-1 color difference formula */
/* Barco monitor, max(r,g,b)=1, n=61 magenta,  2 Jan 1986 */
 
/***********************************************************************/
void rainbow(double  h, double  s, double  v, 
             double  *r, double  *g, double  *b)
/***********************************************************************/
{
    static int i;
    static double  trash;
    h = 60*modf(h/1.5,&trash);
    i = (int)(floor(h)+0);
    h = huettab[i] + (huettab[i+1]-huettab[i])*(h-i);
    dhsv2rgb(h,s,v,r,g,b);
}
 
/***********************************************************************/
void dhsv2rgb(double  h, double  s, double  v, 
              double  *r, double  *g, double  *b)    /*...hexcone model...*/
    /* all variables in range [0,1] */
/***********************************************************************/
    /* here, h=.667 gives blue, h=0 or 1 gives red. */
{  /* see Alvy Ray Smith, Color Gamut Transform Pairs, SIGGRAPH '78 */
    static int i;
    static double f, m, n, k;
    static double  trash;
    h = 6*modf(h,&trash);
    i = (int)(floor(h)+0);
    f = h-i;
    m = (1-s);
    n = (1-s*f);
    k = (1-(s*(1-f)));
    switch(i){
    case 0: *r=1; *g=k; *b=m; break;
    case 1: *r=n; *g=1; *b=m; break;
    case 2: *r=m; *g=1; *b=k; break;
    case 3: *r=m; *g=n; *b=1; break;
    case 4: *r=k; *g=m; *b=1; break;
    case 5: *r=1; *g=m; *b=n; break;
    default: fprintf(stderr,"bad i: %f %d",h,i); gomp_Exit(1);
    }
    f = *r;
    if( f < *g ) f = *g;
    if( f < *b ) f = *b;
    f = v / f;
    *r *= f;
    *g *= f;
    *b *= f;
}


/************************************************************************/
int gomp_ColorByCharge(int Wstr , int ListEntries , const int *List,
                     double *Tmin, double *Tmax)
/************************************************************************/
{
    register int i;
    static   float rr,gg,bb;
    float    min_charge;
    float    max_charge;
    float    delta;
    double   step;
    char     OutText[BUFF_LEN];
    int      alt;
    const float   *atm_charge;

    if(!gomp_GetNumMolecStructs()) {
        gomp_PrintWARNING("no atom structures defined");
        return(1);
    }

    alt        = 0;
    atm_charge = gomp_GetAtomChargePointer(Wstr);

    min_charge =  1.e+30f;
    max_charge = -1.e+30f;

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
        if(atm_charge[i] < min_charge)
            min_charge = atm_charge[i];
        if(atm_charge[i] > max_charge)
            max_charge = atm_charge[i];
    }

    delta = max_charge - min_charge;

    if(RABS(delta) < 1.e-05) { /* charges most likely not assigned */
        gomp_PrintWARNING("charges not assigned (most likely)");
        return(1);
    }

    if(Tmin != NULL)
        min_charge = *Tmin;
    if(Tmax != NULL)
        max_charge = *Tmax;

    delta = max_charge - min_charge;

    sprintf(OutText,"Max charge: %7.3f (red)",max_charge);
    gomp_PrintMessage(OutText);
    sprintf(OutText,"Min charge: %7.3f (blue)",min_charge);
    gomp_PrintMessage(OutText);

    for(i = 0 ; i < ListEntries ; i++) {
        step = (atm_charge[List[i]] - min_charge) / delta;
        gomp_PreRainbow(step,&rr,&gg,&bb);
        gomp_PutAtomColour(Wstr, rr, gg, bb, List[i]);
    }

    return(0);
}
/************************************************************************/
int gomp_ColorByFourth(int StruL , int ListEntries , const int *List)
/************************************************************************/
{
    register int i,j;
    float        min_fourth;
    float        max_fourth;
    float        delta;
    double       step;
    static float rr,gg,bb;
    char         OutText[BUFF_LEN];
    int          Wstr;
    const float *bvalue;

    if(!gomp_GetNumMolecStructs()) {
        gomp_PrintWARNING("no atom structures defined");
        return(1);
    }

    min_fourth =  1.e+30f;
    max_fourth = -1.e+30f;

    Wstr       = StruL;
    bvalue     = gomp_GetAtomBValuePointer(Wstr);

    for(i = 0 ; i < gomp_GetNumAtomsInMolecStruct(Wstr) ; i++) {
        if(bvalue[i] < min_fourth)
            min_fourth = bvalue[i];
        if(bvalue[i] > max_fourth)
            max_fourth = bvalue[i];
    }

    if(RABS(max_fourth - min_fourth) < 1.e-10) {
        gomp_PrintMessage("?WARNING - fourth value is not set ");
        return(1);
    }

    sprintf(OutText,"Max fourth value: %7.3f (red)",max_fourth);
    gomp_PrintMessage(OutText);
    sprintf(OutText,"Min fourth value: %7.3f (blue)",min_fourth);
    gomp_PrintMessage(OutText);

    delta = max_fourth - min_fourth;

    for(i = 0 ; i < ListEntries ; i++) {
        j = List[i];
        step = (bvalue[j] - min_fourth) / delta;
        gomp_PreRainbow(step,&rr,&gg,&bb);
        gomp_PutAtomColour(Wstr, rr, gg, bb, j);
    }

    return(0);
}
/************************************************************************/
int gomp_ColorByResidueNumber(int StruL , int ListEntries , const int *List)
/************************************************************************/
{
    float delta;
    double step;
    int i,j;
    int maxres1,minres1;
    float rr,gg,bb;
    char OutText[BUFF_LEN];
    int        Wstr;
    const int *res1;

    if(!gomp_GetNumMolecStructs()) {
        gomp_PrintWARNING("no atom structures defined");
        return(1);
    }

    Wstr       = StruL;
    res1       = gomp_GetAtomResNum1Pointer(Wstr);
    maxres1    = gomp_GetMaxRes1Num();
    minres1    = gomp_GetMinRes1Num();

    sprintf(OutText,"Max res1: %d (red)",maxres1);
    gomp_PrintMessage(OutText);
    sprintf(OutText,"Min res1: %d (blue)",minres1);
    gomp_PrintMessage(OutText);

    delta = (float)(maxres1 - minres1);

    for(i = 0 ; i < ListEntries ; i++) {
        j = List[i];
        step = (double)(res1[j] - minres1) / delta;
        gomp_PreRainbow(step,&rr,&gg,&bb);
        gomp_PutAtomColour(Wstr, rr, gg, bb, j);
    }

    return(0);
}
/************************************************************************/
int gomp_SetDisplayColourType(int CType)
/************************************************************************/
{
    gomp_ColourDisplayType = CType;

    return(0);
}
/************************************************************************/
int gomp_GetDisplayColourType()
/************************************************************************/
{
    return(gomp_ColourDisplayType);
}
/************************************************************************/
int gomp_RGB2Grayscale(float *RedC, float *GreenC, float *BlueC)
/************************************************************************/
{
    float grayscale;

    grayscale = 0.299 * (*RedC) + 0.587 * (*GreenC) + 0.114 * (*BlueC);
    *RedC   = grayscale;
    *GreenC = grayscale;
    *BlueC  = grayscale;

    return(0);
}
