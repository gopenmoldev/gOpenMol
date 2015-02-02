
/*
                       Copyright (c) 2000 by:
        Leif Laaksonen , Centre for Scientific Computing , ESPOO, FINLAND
            Confidential unpublished property of Leif Laaksonen
                        All rights reserved


      This is a program to convert the Jaguar plt output to a plot format 
      recognized by gOpenMol.
x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-


!!!!!!!!!!!!!!!!!!O B S E R V E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
The cube data has to be given as an orthogonal x-, y-, z-coordinate system
in such a way that the x-axis comes first and the z-axis is given as the 
last one. This means that x is the slowest running coordinate and z is the
fastes running coordinate.
!!!!!!!!!!!!!!!!!!O B S E R V E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


*/

#include <stdio.h>
#include <ctype.h>
#include <sys/types.h>
#include <malloc.h>
#include <string.h>

#define BUFF_LEN        256
#define BOHR_RADIUS     0.52917715  /* conversion constant */
#define JAGUAR_TYPE   201
#define MAX_TITLE_LINES   5
#define Rabs(a)    ( ( a ) > 0.0 ? (a) : -(a))
#define SMALL 1.e-05

#define FWRITE(value_p , size)    { Items = \
                                 fwrite(value_p, size , 1 , Output_p);\
                                 if(Items < 1) {\
                     printf("?ERROR - in writing contour file (*)\n");\
                     return(1);}}

#define FWRITEN(value_p , num , size) { Items = \
                                fwrite(value_p, size , num , Output_p);\
                   if(Items < 1) {\
                     printf("?ERROR - in writing contour file (**)\n");\
                     return(1);}}

/* ................................................................... */
char *Usage = 
"This program converts the Jaguar output \n\
 to a format understandable to gopenmol.\n\
 Usage:\n\
 jaguar2plt -iinput.file -ooutput.file\n\
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\
 The input for Jaguar has to be defined as an orthogonal\n\
 coordinate system, where the x coordinate is the slowest running and\n\
 and z the fastest running coordinate\n";
/* ................................................................... */

   char *AtomSymbols={"\
Ac  Ag  Al  Am  As  Au  B   Ba  Be  Bi  Br  C   Ca  Cd  \
Ce  Cl  Co  Cr  Cs  Cu  D   Dy  Er  Eu  F   Fe  Ga  Gd  \
Ge  H   Hf  Hg  Ho  I   In  Ir  K   La  Li  Lu  Mg  Mn  \
Mo  N   Na  Nb  Nd  Ni  Np  O   Os  P   Pa  Pb  Pd  Pm  \
Po  Pr  Pt  Pu  Ra  Rb  Re  Rh  Ru  S   Sb  Sc  Se  Si  \
Sm  Sn  Sr  Ta  Tb  Tc  Te  Th  Ti  Tl  Tm  U   V   W   \
Y   Yb  Zn  Zr  He  Ne  Ar  Kr  Xe  Rn  "};


    int AtomSymbol_p[] =
       { 89 , 47 , 13 , 95 , 33 , 79 , 5  , 56 , 4  , 83 , 35 , 6  , 20 , 48 ,
         58 , 17 , 27 , 24 , 55 , 29 , 0  , 66 , 68 , 63 , 9  , 26 , 31 , 64 ,
         32 , 1  , 72 , 80 , 67 , 53 , 49 , 77 , 19 , 57 , 3  , 71 , 12 , 25 ,
         42 , 7  , 11 , 41 , 60 , 28 , 93 , 8  , 76 , 15 , 91 , 82 , 46 , 61 ,
         84 , 59 , 78 , 94 , 88 , 37 , 75 , 45 , 44 , 16 , 51 , 55 , 34 , 14 ,
         62 , 50 , 38 , 73 , 73 , 65 , 43 , 90 , 22 , 81 , 69 , 92 , 23 , 74 ,
         39 , 70 , 30 , 40 , 2  , 10 , 18 , 36 , 54 , 86};


/* input */
int TitleLines;                /* actual number of title lines */
int Natoms;                    /* number of atoms */
float Xorig,Yorig,Zorig;       /* x-, y- and z- origin */
int N1;                        /* # if incs in the slowest running direct */
float N1X1,N1Y1,N1Z1;          /* direction */
int N2;
float N2X1,N2Y1,N2Z1;          /* direction */
int N3;                        /* # if incs in the fastest running direct */
float N3X1,N3Y1,N3Z1;          /* direction */
float *Data;                   /* data */
int    ToSave;                    /* index into the orbital array */

/* output */
float Xmax,Ymax,Zmax;
int   TypeOfData = JAGUAR_TYPE;
int   ProduceCoordinateFile = 1;

char InputFile[BUFF_LEN];
char OutputFile[BUFF_LEN];
char CoordinateFile[BUFF_LEN];

/* functions */
int  ReadInputData(void);
int  WriteInputData(void);
int  WriteCoordinateFile(void);

/* externals */
extern char *Number2Name(int);

/**************************************************************************/
main(int args, char** argv)
/**************************************************************************/
{
    int i;

    printf("**********************************************************\n");
    printf("* Convert Jaguar output into the plot                    *\n");
    printf("* format known by gOpenMol (or SCARECROW).               *\n");
    printf("*                                                        *\n");
    printf("* Leif Laaksonen (CSC) 2000                              *\n");
    printf("* Email: Leif.Laaksonen@csc.fi                           *\n");
    printf("*                               Version: 10/08/2000      *\n");
    printf("**********************************************************\n\n");

    if(args < 2) {
    printf("%s",Usage);
    exit(0);}

    memset(OutputFile , 0 , BUFF_LEN);

/* go and hunt for the input switches */
    if(args > 1) {
       for(i = 1 ; i < args ; i++) {
                 /* check first that it is a '-' command */
                    if(argv[i][0] == '-' ) { /* yes it is */
/* '-o' output file name */
      if(argv[i][1] == 'o') {
      strncpy(OutputFile,&argv[i][2],BUFF_LEN);
      }
/* '-i' input file name */
      if(argv[i][1] == 'i') {
      strncpy(InputFile,&argv[i][2],BUFF_LEN);
      }
/* '-p' prevent coordinate file output */
      if(argv[i][1] == 'p') {
       ProduceCoordinateFile = 0;
      }
     }
    }
   }

    if(InputFile[0] == '\0') {
      printf("$ERROR - supply input file name\n");
      exit(1);
    }

    printf("File names:\n");
    printf("Input file:      '%s'\n",InputFile);
    
/*  Process the data ... */

    if(ReadInputData()) {
       printf("$ERROR - can't read input data\n");
       exit(1);
    }

    printf("Job done ...\n");
}
/**************************************************************************/
int ReadInputData()
/**************************************************************************/
{
    FILE *Input_p;
    char  Text[BUFF_LEN];
    int   i,j,k,l,ijk,ijkl;
    int   IsInteger;
    char  Temp1[BUFF_LEN];
    char  Temp2[BUFF_LEN];
    char  Temp3[BUFF_LEN];
    char  Temp4[BUFF_LEN];
    int   Hit;
    int   tatoms;
    float txorig;
    float tyorig;
    float tzorig;

    float textentx;
    float textenty;
    float textentz;

    int   tnptsx;
    int   tnptsy;
    int   tnptsz;

    int   tiplot;

    Input_p = fopen(InputFile , "r");
    if(Input_p == NULL) {
       printf("$ERROR - can't open input file '%s'\n",InputFile);
       return(1);
    }

    strncpy(Temp3 , OutputFile , BUFF_LEN);

    txorig = tyorig = tzorig = 0.0;

    while(fgets(Text,BUFF_LEN,Input_p) != NULL) { 

        if(strncmp(Text,"&plot",5)   == 0) {

          while(fgets(Text,BUFF_LEN,Input_p) != NULL) {

               if(strncmp(Text,"iplot=",6)   == 0) {
                   sscanf(Text,"%*s %d",&tiplot);
               }
               else if(strncmp(Text,"origin=",7)   == 0) {
                   sscanf(Text,"%*s %f %f %f",&txorig,&tyorig,&tzorig);
                   txorig *= BOHR_RADIUS;
                    tyorig *= BOHR_RADIUS;
                     tzorig *= BOHR_RADIUS;
               }
               else if(strncmp(Text,"extentx=",8)   == 0) {
                   sscanf(Text,"%*s %f %f %f",&N1X1,&N1Y1,&N1Z1);
                   N1X1 *= BOHR_RADIUS;
                    N1Y1 *= BOHR_RADIUS;
                     N1Z1 *= BOHR_RADIUS;               
               }
               else if(strncmp(Text,"extenty=",8)   == 0) {
                   sscanf(Text,"%*s %f %f %f",&N2X1,&N2Y1,&N2Z1);
                   N2X1 *= BOHR_RADIUS;
                    N2Y1 *= BOHR_RADIUS;
                     N2Z1 *= BOHR_RADIUS;               
               }
               else if(strncmp(Text,"extentz=",8)   == 0) {
                   sscanf(Text,"%*s %f %f %f",&N3X1,&N3Y1,&N3Z1);
                   N3X1 *= BOHR_RADIUS;
                    N3Y1 *= BOHR_RADIUS;
                     N3Z1 *= BOHR_RADIUS;               
               }
               else if(strncmp(Text,"npts=",5)   == 0) {
                   sscanf(Text,"%*s %d %d %d",&tnptsx,&tnptsy,&tnptsz);
               }
               else if(strncmp(Text,"&end",4)   == 0) {

/* check that this is a 'pure' x-coordinate */
                   if(Rabs(N1X1) < SMALL) {
                      printf("$ERROR - most likely your step in the x-direction (%f) is too small\n",N1X1);
                      exit(30);}
                   if(Rabs(N1Y1) > SMALL || Rabs(N1Z1) > SMALL) {
                      printf("$ERROR - first input has to be pure x-axis (y: %f , z: %f)\n",
                      N1Y1,N1Z1);
                      exit(31);}

/* check that this is a 'pure' y-coordinate */
                  if(Rabs(N2Y1) < SMALL) {
                     printf("$ERROR - most likely your step in the y-direction (%f) is too small\n",N2Y1);
                     exit(32);}
                  if(Rabs(N2X1) > SMALL || Rabs(N2Z1) > SMALL) {
                     printf("$ERROR - second input has to be pure y-axis (x: %f , z: %f)\n",
                     N2X1,N2Z1);
                     exit(33);}


/* check that this is a 'pure' z-coordinate */
                   if(Rabs(N3Z1) < SMALL) {
                      printf("$ERROR - most likely your step in the z-direction (%f) is too small\n",N3Z1);
                      exit(34);}
                   if(Rabs(N3X1) > SMALL || Rabs(N3Y1) > SMALL) {
                      printf("$ERROR - last input has to be pure z-axis (x: %f , y: %f)\n",
                      N3X1,N3Y1);
                      exit(35);}

/* look for a dot in the name ... */
                   {
                    char *pdest;
                    int result;

                        if(Temp3[0] == (char)NULL) {
                           strncpy(Temp1,InputFile,BUFF_LEN);
                        } else {
                           strncpy(Temp1,OutputFile,BUFF_LEN);
                        }

                        pdest = strrchr( Temp1, '.' );
                        result = pdest - Temp1;

                        strncpy(Temp2,Temp1,result);

                        if( pdest != NULL ) {
                            sprintf(&Temp2[result],"_gom_%d.plt\0",tiplot);
                        } else {
                            sprintf(Temp2,"%s_gom_%d.plt",Temp1,tiplot);
                        }
                   }

                   strncpy(OutputFile, Temp2 , BUFF_LEN);
                   printf("Output file:     '%s'\n",OutputFile);

                   Xorig = txorig ; 
                   Yorig = tyorig ; 
                   Zorig = tzorig ;

                   N1X1 = N1X1 / (tnptsx - 1);
                   N2Y1 = N2Y1 / (tnptsy - 1);
                   N3Z1 = N3Z1 / (tnptsz - 1);
                  
                   N1 = tnptsx;
                    N2 = tnptsy;
                     N3 = tnptsz;

                   Data = (float *)malloc(sizeof(float) * N1 * N2 * N3);
                   if(Data == NULL) exit(19);

                   for(i = 0   ; i < N1 ; i++)  { /* 1 */
                    for(j = 0  ; j < N2 ; j++)  { /* 2 */
                     for(k = 0 ; k < N3 ; k++)  { /* 3 */

                         ijk = i + N1 * j + N1 * N2 * k;

                         ijkl = fscanf(Input_p,"%f",&Data[ijk]);
                         if(ijkl != 1) {
                            printf("$ERROR - in reading the grid data\n");
                            printf("$ERROR - at ijk: %d %d %d\n",i,j,k);
                            exit(20);
                         }
                     } /* 3 */
                    }  /* 2 */
                   }   /* 1 */

                   if(WriteInputData()) {
                       printf("$ERROR - can't write the gOpenMol plt file\n");
                       exit(21);
                   }

                   strncpy(OutputFile, Temp3 , BUFF_LEN);
                   free(Data);

                   break;
               }
          }
        }
    }
/*
iplot=    1
iorb1a=    5
iorb2a=    6
iorb1b=    0
iorb2b=    0
origin=    0.000000    0.000000    0.000000
extentx=    4.000000    0.000000    0.000000
extenty=    0.000000    4.000000    0.000000
extentz=    0.000000    0.000000    4.000000
npts=   20   20   20
*/        
   fclose(Input_p);

   return(0);
}

/**************************************************************************/
int WriteInputData()
/**************************************************************************/
{
    FILE *Output_p;
    int   i,j,k,l,ijk,ijkl;
    int Items;
    float Help1;
    float Help2;
    float MinV;
	float MaxV;

    Output_p = fopen(OutputFile,"wb");
    if(Output_p == NULL) {
      printf("$ERROR - can't open output file: '%s' \n",OutputFile);
      return(1);
    }

	MinV =  1.0e+25;
	MaxV = -1.0e+25;

    i = 3;
    FWRITE(&i,sizeof(int));
    FWRITE(&TypeOfData , sizeof(int));

      FWRITE(&N3   , sizeof(int));
       FWRITE(&N2  , sizeof(int));
        FWRITE(&N1 , sizeof(int));

      Help1  = Zorig;
       Help2 = Zorig + (N3 - 1) * N3Z1;
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));
      Help1  = Yorig;
       Help2 = Yorig + (N2 - 1) * N2Y1;
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));
      Help1  = Xorig;
       Help2 = Xorig + (N1 - 1) * N1X1;
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));

   for(k = 0   ; k < N3 ; k++)  {
    for(j = 0  ; j < N2 ; j++)  { 
     for(i = 0 ; i < N1 ; i++)  {

     ijk = i + N1 * j + N1 * N2 * k;

     if(Data[ijk] < MinV) MinV = Data[ijk];
     if(Data[ijk] > MaxV) MaxV = Data[ijk];

     FWRITE(&Data[ijk] , sizeof(float));
     }
    }
 }

   fclose(Output_p);

   printf("Min value: %f\n",MinV);
   printf("Max value: %f\n",MaxV);

   return(0);
}
