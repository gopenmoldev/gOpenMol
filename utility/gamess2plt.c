
/*
                       Copyright (c) 2000 - 2004 by:
        Leif Laaksonen , Centre for Scientific Computing , ESPOO, FINLAND
            Confidential unpublished property of Leif Laaksonen
                        All rights reserved

      10.08.2001
      The program is extended to handle extended CUBE data.

      This is a program to convert the output using the "cube"
      command from the GAMESS program to a plot format recognized
      by gOpenMol.

      The PC Gamess PUNCH file can contain several "$cube" and
      $end" separated fields.

      The this program adds a running number (_#) after the dot
      or at the end of the file name(s) according to how many fields it finds
      in the "PUNCH" file.

      Example:

      Input file name:           "punch"
      Output plt file name(s):   punch_X_1.plt, punch_Y_2.plt, ...
      Output coord file name(s): punch_1.crd, punch_2.crd, ...

      Where the X and Y can be:
      * orb for an orbital density or orbital
      * den for electron density
      * grad_norm for gradient norm and
      * lap for Laplacian

      The running number corresponds to the coordinate file number.

      Three $cube groups:
      1) first for density
      2) second for gradient norm of the density,
      3) and finally the Laplacian of the density.

        They are easily distinguished because each $cube group contains an
        automatically-generated descriptor (it is always on the second
        line of the $cube group) like:
        UHF density
        Gradient norm of UHF density
        Laplacian of UHF density

        Thus you simply need to look at the keywords like
        "Gradient norm" and "Laplacian" to identify the
        type of the $cube contents.


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 This program converts the Gamess output using the 'cube'
 command to a form understandable to gopenmol.
 Usage:
 gamess2plt -iinput.cube -ooutput.plt
 Options:  -mXXX , where XXX is the molecular orbital number to be
                   placed in the plot file
           -p      prevent the output of the coordinate file
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      If you need any help please feel free to contact:

      Leif.Laaksonen@csc.fi



!!!!!!!!!!!!!!!!!!O B S E R V E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
The cube data has to be given as an orthogonal x-, y-, z-coordinate system
in such a way that the x-axis comes first and the z-axis is given as the
last one. This means that x is the slowest running coordinate and z is the
fastes running coordinate.

This program produces also a coordinate file in the CHARMM 'crd' format
which can be read by gOpenMol or SCARECROW.

!!!!!!!!!!!!!!!!!!O B S E R V E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


*/

#include <stdio.h>
#include <ctype.h>
#include <sys/types.h>
#include <malloc.h>
#include <string.h>

#define BUFF_LEN        256
#define BOHR_RADIUS     0.52917715  /* conversion constant */
#define GAMESS_TYPE   202
#define MAX_TITLE_LINES   5
#define Rabs(a)    ( ( a ) > 0.0 ? (a) : -(a))
#define SMALL 1.e-05
#define CUBE_HOOK_START1 " $cube"
#define CUBE_HOOK_END1   " $end"
#define CUBE_HOOK_START2 " $CUBE"
#define CUBE_HOOK_END2   " $END"

#define CUBE_TYPE_ORBITAL            "Orbitals"
#define CUBE_TYPE_ORBITAL_DENSITY    "Density of orbitals"
#define CUBE_TYPE_DENSITY            "density"
#define CUBE_TYPE_GRAD_NORM          "Gradient norm"
#define CUBE_TYPE_GRAD_FIELD         "Gradient field"
#define CUBE_TYPE_LAPLACIAN          "Laplacian of"
#define CUBE_TYPE_ESP                "ESP"

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
"This program converts the Gamess output using the 'cube'\n\
 command to a form understandable to gopenmol.\n\
 Usage:\n\
 gamess2plt -iinput.file -ooutput.file\n\
 Options:      -p      prevent the output of the coordinate file\n\
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\
 The 'Cube' input for GaussianXX has to be defined as an orthogonal\n\
 coordinate system, where the x coordinate is the slowest running and\n\
 and z the fastest running coordinate\n\
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\
 Please observe that a PCGAMESS CUBE file is NOT equivalent to a \n\
 GaussianXX CUBE file!!!!\n";
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
char TitleText[MAX_TITLE_LINES][BUFF_LEN];
int TitleLines;                /* actual number of title lines */
int Natoms;                    /* number of atoms */
float Xorig,Yorig,Zorig;       /* x-, y- and z- origin */
int N1;                        /* # if incs in the slowest running direct */
float N1X1,N1Y1,N1Z1;          /* direction */
int N2;
float N2X1,N2Y1,N2Z1;          /* direction */
int N3;                        /* # if incs in the fastest running direct */
float N3X1,N3Y1,N3Z1;          /* direction */
int *IA;                       /* atomic number arraypointer */
float *Chgn;                   /* charge array pointer */
float *XC,*YC,*ZC;             /* coordinate pointers  */
float *Data;                   /* data */
int    MolecularOrbitals = 0;  /* switch to handle molecular orbitals */
int   *MolecularOrbital_p;     /* index pointer */
int    MolecularOrbitalsDefined; /* orbitals defined in file */
int    MolecularOrbital2Plot;
int    ToSave;                    /* index into the orbital array */

/* output */
float Xmax,Ymax,Zmax;
int   TypeOfData = GAMESS_TYPE;
int   ProduceCoordinateFile = 1;
int   DataType = 0;

char InputFile[BUFF_LEN];
char OutputFile[BUFF_LEN];
char OutputFileTemp[BUFF_LEN];
char CoordinateFile[BUFF_LEN];

/* functions */
void MakeOutputFileName(char *);
int  ReadInputData(void);
int  WriteInputData(void);
int  WriteCoordinateFile(void);

int   FileLoop;

/* externals */
extern char *Number2Name(int);

/**************************************************************************/
main(int args, char** argv)
/**************************************************************************/
{
    int i;

    printf("******************************************************\n");
    printf("* Convert a 'cube' output from Gamess into the plot  *\n");
    printf("* format known by gOpenMol (or SCARECROW).           *\n");
    printf("*                                                    *\n");
    printf("* Leif Laaksonen (CSC) 2005                          *\n");
    printf("* Email: Leif.Laaksonen(at)csc.fi                    *\n");
    printf("*                               Version: 17/06/2005  *\n");
    printf("******************************************************\n\n");

    if(args < 2) {
    printf("%s",Usage);
    exit(0);}

/* go and hunt for the input switches */
    if(args > 1) {
       for(i = 1 ; i < args ; i++) {
                 /* check first that it is a '-' command */
                    if(argv[i][0] == '-' ) { /* yes it is */
/* '-m' orbital cube file                  */
      if(argv[i][1] == 'm') {

      printf("==> This option is not supported! <==\n    All available data will be used\n\n");          
/*
      sscanf(&argv[i][2],"%d",&MolecularOrbital2Plot);
      printf("Orbital number to plot: %d\n",MolecularOrbital2Plot);
*/
      }
/* '-o' output file name */
      if(argv[i][1] == 'o') {
      strncpy(OutputFileTemp,&argv[i][2],BUFF_LEN);
      }
/* '-i' input file name */
      if(argv[i][1] == 'i') {
      strncpy(InputFile,&argv[i][2],BUFF_LEN);
      }
/* '-p' prevent coordinate file output */
      if(argv[i][1] == 'p') {
       ProduceCoordinateFile = 0;
      }
/* '-t' type of file */
      if(argv[i][1] == 't') {
      sscanf(&argv[i][2],"%d",&DataType);
	  if (DataType == 1) {
         printf("** Orbital input data type\n");
	  } else if(DataType == 2) {
         printf("** Density input data type\n");
	  } else if(DataType == 3) {
         printf("** Gradient norm input data type\n");
	  } else if(DataType == 4) {
         printf("** Gradient field input data type\n");
	  } else if(DataType == 5) {
         printf("** Laplacian input data type\n");
	  } else {
	    printf("ERROR - unknown input data type '%d'\n",DataType);
	    exit(1);
	  }
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
void MakeOutputFileName(char *Filename)
/**************************************************************************/
{
     int i;
     char *pdest;
     int result;

     if(OutputFileTemp[0] == (char)NULL) {
/* look for a dot in the name ... */
       pdest = strrchr( Filename, '.' );
       result = pdest - Filename;
       if( pdest != NULL ) {
           strncpy(OutputFile,Filename,result);
           if(DataType == 1)
             sprintf(&OutputFile[result],"_orb_%%s_%d.plt\0",FileLoop);
           else if(DataType == 2)
             sprintf(&OutputFile[result],"_den_%%s_%d.plt\0",FileLoop);
           else if(DataType == 3)
             sprintf(&OutputFile[result],"_grad_norm_%%s_%d.plt\0",FileLoop);
           else if(DataType == 4)
             sprintf(&OutputFile[result],"_grad_field_%%s_%d.plt\0",FileLoop);
           else if(DataType == 5)
             sprintf(&OutputFile[result],"_lap_%%s_%d.plt\0",FileLoop);
           else if(DataType == 6)
             sprintf(&OutputFile[result],"_ESP_den_%%s_%d.plt\0",FileLoop);
           else
               exit(1);

           strncpy(CoordinateFile,Filename,result);
           sprintf(&CoordinateFile[result],"_%d.crd\0",FileLoop);
       }
       else {
           if(DataType == 1)
             sprintf(OutputFile,"%s_orb_%%s_%d.plt",Filename,FileLoop);
           else if(DataType == 2)
             sprintf(OutputFile,"%s_den_%%s_%d.plt",Filename,FileLoop);
           else if(DataType == 3)
             sprintf(OutputFile,"%s_grad_norm_%%s_%d.plt",Filename,FileLoop);
           else if(DataType == 4)
             sprintf(OutputFile,"%s_grad_field_%%s_%d.plt",Filename,FileLoop);
           else if(DataType == 5)
             sprintf(OutputFile,"%s_lap_%%s_%d.plt",Filename,FileLoop);
           else if(DataType == 6)
             sprintf(OutputFile,"%s_ESP_den_%%s_%d.plt",Filename,FileLoop);
           else
             exit(1);
           sprintf(CoordinateFile,"%s_%d.crd",Filename,FileLoop);
       }

     } else {

/* look for a dot in the name ... */
       pdest = strrchr( OutputFileTemp, '.' );
       result = pdest - OutputFileTemp;
       if( pdest != NULL ) {
           strncpy(OutputFile,OutputFileTemp,result);
           if(DataType == 1)
              sprintf(&OutputFile[result],"_orb_%%s_%d.plt\0",FileLoop);
           else if(DataType == 2)
              sprintf(&OutputFile[result],"_den_%%s_%d.plt\0",FileLoop);
           else if(DataType == 3)
              sprintf(&OutputFile[result],"_grad_norm_%%s_%d.plt\0",FileLoop);
           else if(DataType == 4)
              sprintf(&OutputFile[result],"_grad_field_%%s_%d.plt\0",FileLoop);
           else if(DataType == 5)
              sprintf(&OutputFile[result],"_lap_%%s_%d.plt\0",FileLoop);
           else if(DataType == 6)
              sprintf(&OutputFile[result],"_ESP_den_%%s_%d.plt\0",FileLoop);
           else
              exit(1);
           strncpy(CoordinateFile,OutputFileTemp,result);
           sprintf(&CoordinateFile[result],"_%d.crd\0",FileLoop);
       }
       else {
           if(DataType == 1)
             sprintf(OutputFile,"%s_orb_%%s_%d.plt",OutputFileTemp,FileLoop);
           else if(DataType == 2)
             sprintf(OutputFile,"%s_den_%%s_%d.plt",OutputFileTemp,FileLoop);
           else if(DataType == 3)
             sprintf(OutputFile,"%s_grad_norm_%%s_%d.plt",OutputFileTemp,FileLoop);
           else if(DataType == 4)
             sprintf(OutputFile,"%s_grad_field_%%s_%d.plt",OutputFileTemp,FileLoop);
           else if(DataType == 5)
             sprintf(OutputFile,"%s_lap_%%s_%d.plt",OutputFileTemp,FileLoop);
           else if(DataType == 6)
             sprintf(OutputFile,"%s_ESP_den_%%s_%d.plt",OutputFileTemp,FileLoop);
           else
             exit(1);
           sprintf(CoordinateFile,"%s_%d.crd",OutputFileTemp,FileLoop);
       }
     }

     return;
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
    char  InputFileTemp[BUFF_LEN];
    int   Hit;
    int   tatoms;
    float txorig;
    float tyorig;
    float tzorig;
	char *token;
    int   FilePosition;

    Input_p = fopen(InputFile , "r");
    if(Input_p == NULL) {
      printf("$ERROR - can't open input file '%s'\n",InputFile);
      return(1);
    }

    FileLoop = 1;

    while(!feof(Input_p)) {

/* save file postion, there is a need to backstep in analyzing multiline orbital lists */
       FilePosition = ftell(Input_p);

       fgets(Text,BUFF_LEN,Input_p);

       if(!strncmp(Text,CUBE_HOOK_START1,strlen(CUBE_HOOK_START1)) ||
          !strncmp(Text,CUBE_HOOK_START2,strlen(CUBE_HOOK_START2))) {

            printf("\nTitle in file (job title):\n");
/* first comes an unknow number of title lines (MAX lines is 5) */
            TitleLines = 0;

/* title loop # 1* == START */
            for(i = 0 ; i < MAX_TITLE_LINES ; i++) {
                fgets(Text,BUFF_LEN,Input_p);
                if(!strncmp(Text,CUBE_HOOK_END1,strlen(CUBE_HOOK_END1)) ||
                   !strncmp(Text,CUBE_HOOK_END2,strlen(CUBE_HOOK_END2))) {
                    break;
                }
/* try to figure out what the data is orbital, density, gradient norm, laplacian...*/
                if((i == 1) && (!DataType)) {
                    if(strstr(Text,CUBE_TYPE_ORBITAL)) {
                        DataType = 1;
                        printf("Found CUBE data of type orbital\n");

						strncpy(Temp1,Text,BUFF_LEN);
						MolecularOrbitalsDefined = 0;
                        token = strtok( Temp1, " " );
                        while( token != NULL )
                        {
                           /* While there are tokens in "string" */
                           /* Get next token: */
                           if(isdigit(*token) || (*token == '+') || (*token == '-'))
						        MolecularOrbitalsDefined++;
                           token = strtok( NULL, " " );

                           if(token == NULL) {
                             fgets(Temp1,BUFF_LEN,Input_p);
                               if(strchr(Temp1,'.') == NULL) {
                                  token = strtok(Temp1, " ");
                               }
                           }
                        }

                        if(fseek(Input_p , (long)FilePosition , SEEK_SET)) {
                           printf("ERROR - can't rewind data file at postion '1001'\n");
                           exit(1001);
                        }

                        MolecularOrbital_p = (int *)malloc(sizeof(int) * (MolecularOrbitalsDefined+1));
                        if(MolecularOrbital_p == NULL) exit(15);

                        fgets(Text,BUFF_LEN,Input_p);
                        strncpy(Temp1,Text,BUFF_LEN);
                        Hit = 0;
						token = strtok( Temp1, " " );
                        while( token != NULL )
                        {
							if(isdigit(*token) || (*token == '+') || (*token == '-')) {
                              sscanf(token,"%d", &MolecularOrbital_p[Hit]);
                              Hit++;
							}
                           /* Get next token: */
                           token = strtok( NULL, " " );

                           if(token == NULL) {
                             fgets(Temp1,BUFF_LEN,Input_p);
                               if(strchr(Temp1,'.') == NULL) {
                                  token = strtok(Temp1, " "); 
                               } else {
                                   strncpy(Text,Temp1,BUFF_LEN);
                               }
                           }
                        }

                        printf("Found '%d' orbital entries\n",MolecularOrbitalsDefined);

                    } else if(strstr(Text,CUBE_TYPE_ESP)) {
                        DataType = 6;
                        printf("Found CUBE data of type density\n");

						strncpy(Temp1,Text,BUFF_LEN);
						MolecularOrbitalsDefined = 0;
                        token = strtok( Temp1, " " );
                        while( token != NULL )
                        {
                           /* While there are tokens in "string" */
                           /* Get next token: */
						   if(isdigit(*token) || (*token == '+') || (*token == '-'))
						        MolecularOrbitalsDefined++;
                           token = strtok( NULL, " " );

                           if((MolecularOrbitalsDefined) && (token == NULL)) {
                             fgets(Temp1,BUFF_LEN,Input_p);
                               if(strchr(Temp1,'.') == NULL) {
                                  token = strtok(Temp1, " ");
                               }
                           }

                        }

						if(MolecularOrbitalsDefined) {
                          if(fseek(Input_p , (long)FilePosition , SEEK_SET)) {
                             printf("ERROR - can't rewind data file at postion '1006'\n");
                             exit(1006);
                          }

                          MolecularOrbital_p = (int *)malloc(sizeof(int) * (MolecularOrbitalsDefined+1));
                          if(MolecularOrbital_p == NULL) exit(15);

                          strncpy(Temp1,Text,BUFF_LEN);
                          Hit = 0;
						  token = strtok( Temp1, " " );
                          while( token != NULL )
                          {
							 if(isdigit(*token) || (*token == '+') || (*token == '-')) {
                                sscanf(token,"%d", &MolecularOrbital_p[Hit]);
                                Hit++;
							 }
                           /* Get next token: */
                             token = strtok( NULL, " " );

                           if(token == NULL) {
                             fgets(Temp1,BUFF_LEN,Input_p);
                               if(strchr(Temp1,'.') == NULL) {
                                  token = strtok(Temp1, " "); 
                               } else {
                                   strncpy(Text,Temp1,BUFF_LEN);
                               }
                           }
                          }

                           printf("Found '%d' orbital entries\n",MolecularOrbitalsDefined);
						}
					} else if(strstr(Text,CUBE_TYPE_LAPLACIAN)) {
                        DataType = 5;
                        printf("Found CUBE data of type Laplacian\n");

						strncpy(Temp1,Text,BUFF_LEN);
						MolecularOrbitalsDefined = 0;
                        token = strtok( Temp1, " " );
                        while( token != NULL )
                        {
                           /* While there are tokens in "string" */
                           /* Get next token: */
                           if(isdigit(*token) || (*token == '+') || (*token == '-'))
						        MolecularOrbitalsDefined++;
                           token = strtok( NULL, " " );

                           if(((MolecularOrbitalsDefined)) && (token == NULL)) {
                             fgets(Temp1,BUFF_LEN,Input_p);
                               if(strchr(Temp1,'.') == NULL) {
                                  token = strtok(Temp1, " ");
                               }
                           }
                        }

                        if(MolecularOrbitalsDefined) {

                          if(fseek(Input_p , (long)FilePosition , SEEK_SET)) {
                             printf("ERROR - can't rewind data file at postion '1005'\n");
                             exit(1005);
                          }

                          MolecularOrbital_p = (int *)malloc(sizeof(int) * (MolecularOrbitalsDefined+1));
                          if(MolecularOrbital_p == NULL) exit(15);

                          strncpy(Temp1,Text,BUFF_LEN);
                          Hit = 0;
		                  token = strtok( Temp1, " " );
                          while( token != NULL )
                          {
                              if(isdigit(*token) || (*token == '+') || (*token == '-')) {
                                sscanf(token,"%d", &MolecularOrbital_p[Hit]);
                                Hit++;
							  }
                           /* Get next token: */
                              token = strtok( NULL, " " );

                             if(token == NULL) {
                               fgets(Temp1,BUFF_LEN,Input_p);
                                 if(strchr(Temp1,'.') == NULL) {
                                    token = strtok(Temp1, " "); 
                                 } else {
                                     strncpy(Text,Temp1,BUFF_LEN);
                                 }
                             }
                           }

                           printf("Found '%d' orbital entries\n",MolecularOrbitalsDefined);
                        }
                    } else if(strstr(Text,CUBE_TYPE_GRAD_FIELD)) {
                        DataType = 4;
                        printf("Found CUBE data of type density gradient field\n");
                    } else if(strstr(Text,CUBE_TYPE_GRAD_NORM)) {
                        DataType = 3;
                        printf("Found CUBE data of type density gradient norm\n");
						strncpy(Temp1,Text,BUFF_LEN);
						MolecularOrbitalsDefined = 0;
                        token = strtok( Temp1, " " );
                        while( token != NULL )
                        {
                           /* While there are tokens in "string" */
                           /* Get next token: */
                           if(isdigit(*token) || (*token == '+') || (*token == '-'))
						        MolecularOrbitalsDefined++;
                           token = strtok( NULL, " " );

                           if((MolecularOrbitalsDefined) && (token == NULL)) {
                             fgets(Temp1,BUFF_LEN,Input_p);
                               if(strchr(Temp1,'.') == NULL) {
                                  token = strtok(Temp1, " ");
                               }
                           }
                        }

						if(MolecularOrbitalsDefined) {

                          if(fseek(Input_p , (long)FilePosition , SEEK_SET)) {
                             printf("ERROR - can't rewind data file at postion '1003'\n");
                             exit(1003);
                          }

                          MolecularOrbital_p = (int *)malloc(sizeof(int) * (MolecularOrbitalsDefined+1));
                          if(MolecularOrbital_p == NULL) exit(15);

                          strncpy(Temp1,Text,BUFF_LEN);
                          Hit = 0;
						  token = strtok( Temp1, " " );
                          while( token != NULL )
                          {
                             if(isdigit(*token) || (*token == '+') || (*token == '-')) {
                                sscanf(token,"%d", &MolecularOrbital_p[Hit]);
                                Hit++;
							 }
                           /* Get next token: */
                           token = strtok( NULL, " " );
                           if(token == NULL) {
                             fgets(Temp1,BUFF_LEN,Input_p);
                               if(strchr(Temp1,'.') == NULL) {
                                  token = strtok(Temp1, " "); 
                               } else {
                                   strncpy(Text,Temp1,BUFF_LEN);
                               }
                           }
                          }

                          printf("Found '%d' orbital entries\n",MolecularOrbitalsDefined);
						}
                    } else if(strstr(Text,CUBE_TYPE_DENSITY) || 
                              strstr(Text,CUBE_TYPE_ORBITAL_DENSITY)) {
                        DataType = 2;
                        printf("Found CUBE data of type density\n");
						strncpy(Temp1,Text,BUFF_LEN);
						MolecularOrbitalsDefined = 0;

                        token = strtok( Temp1, " " );
                        while( token != NULL )
                        {
                           /* While there are tokens in "string" */
                           /* Get next token: */
                           if(isdigit(*token) || (*token == '+') || (*token == '-'))
						        MolecularOrbitalsDefined++;
                           token = strtok( NULL, " " );

                           if((MolecularOrbitalsDefined) && (token == NULL)) {
                             fgets(Temp1,BUFF_LEN,Input_p);
                               if(strchr(Temp1,'.') == NULL) {
                                  token = strtok(Temp1, " ");
                               }
                           }
                        }

						if(MolecularOrbitalsDefined) {

                          if(fseek(Input_p , (long)FilePosition , SEEK_SET)) {
                             printf("ERROR - can't rewind data file at postion '1002'\n");
                             exit(1002);
                          }

                          MolecularOrbital_p = (int *)malloc(sizeof(int) * (MolecularOrbitalsDefined+1));
                          if(MolecularOrbital_p == NULL) exit(15);

                          strncpy(Temp1,Text,BUFF_LEN);
                          Hit = 0;
						  token = strtok( Temp1, " " );
                          while( token != NULL )
                          {
                             if(isdigit(*token) || (*token == '+') || (*token == '-')) {
                                sscanf(token,"%d", &MolecularOrbital_p[Hit]);
                                Hit++;
							 }
                           /* Get next token: */
                             token = strtok( NULL, " " );

                           if(token == NULL) {
                             fgets(Temp1,BUFF_LEN,Input_p);
                               if(strchr(Temp1,'.') == NULL) {
                                  token = strtok(Temp1, " "); 
                               } else {
                                   strncpy(Text,Temp1,BUFF_LEN);
                               }
                           }
                          }

                          printf("Found '%d' orbital entries\n",MolecularOrbitalsDefined);
                        } 
                    } else {
                        printf("%s",Text);
                    }
                }
                sscanf(Text,"%s %s %s %s",Temp1,Temp2,Temp3,Temp4);
/* if Temp1 is an integer and Temp2,temp3,temp4 are floats then the title
   lines are all read */
                IsInteger = 1;
                for(j = 0 ; j < strlen(Temp1) ; j++) {
                    if(isalpha(Temp1[j])) {
                       IsInteger = 0;
                       break;
                    }
                }
                if(!IsInteger) {
                    printf("%s",Text);
                    strncpy(TitleText[i],Text,BUFF_LEN);
                    TitleText[i][strlen(TitleText[i]) - 1] = (char)NULL;
                    TitleLines++;
                }
                else
                    break;
            }
/* title loop # 1* == STOP */

            if(!strncmp(Text,CUBE_HOOK_END1,strlen(CUBE_HOOK_END1))||
               !strncmp(Text,CUBE_HOOK_END2,strlen(CUBE_HOOK_END2))) {
                break;
            }

/* now starts the data ... */
            sscanf(Text,"%d %f %f %f",&Natoms,&Xorig,&Yorig,&Zorig);
            Xorig *= BOHR_RADIUS;
             Yorig *= BOHR_RADIUS;
              Zorig *= BOHR_RADIUS;
              printf("Number of atoms: %d, x-, y-, z-origin (in Anstrom): %f,%f,%f\n",
                      abs(Natoms),Xorig,Yorig,Zorig);

			 fscanf(Input_p,"%d %f %f %f",&N1,&N1X1,&N1Y1,&N1Z1);
             N1X1 *= BOHR_RADIUS;
              N1Y1 *= BOHR_RADIUS;
               N1Z1 *= BOHR_RADIUS;
               printf("Number of points: %d, in direction (x,y,z) %f %f %f\n",
                       N1,N1X1,N1Y1,N1Z1);

/* check that this is a 'pure' x-coordinate */
             if(Rabs(N1X1) < SMALL) {
                printf("$ERROR - most likely your step in the x-direction (%f) is too small\n",N1X1);
                exit(20);
             }
             if(Rabs(N1Y1) > SMALL || Rabs(N1Z1) > SMALL) {
                printf("$ERROR - first input has to be pure x-axis (y: %f , z: %f)\n",
                N1Y1,N1Z1);
                exit(21);
             }

             fscanf(Input_p,"%d %f %f %f",&N2,&N2X1,&N2Y1,&N2Z1);
             N2X1 *= BOHR_RADIUS;
              N2Y1 *= BOHR_RADIUS;
               N2Z1 *= BOHR_RADIUS;
               printf("Number of points: %d, in direction (x,y,z) %f %f %f\n",
                       N2,N2X1,N2Y1,N2Z1);

/* check that this is a 'pure' y-coordinate */
             if(Rabs(N2Y1) < SMALL) {
                printf("$ERROR - most likely your step in the y-direction (%f) is too small\n",N2Y1);
                exit(22);
             }
             if(Rabs(N2X1) > SMALL || Rabs(N2Z1) > SMALL) {
                printf("$ERROR - second input has to be pure y-axis (x: %f , z: %f)\n",
                N2X1,N2Z1);
                exit(23);
             }

             fscanf(Input_p,"%d %f %f %f",&N3,&N3X1,&N3Y1,&N3Z1);
             N3X1 *= BOHR_RADIUS;
              N3Y1 *= BOHR_RADIUS;
               N3Z1 *= BOHR_RADIUS;
               printf("Number of points: %d, in direction (x,y,z) %f %f %f\n",
                       N3,N3X1,N3Y1,N3Z1);

/* check that this is a 'pure' z-coordinate */
             if(Rabs(N3Z1) < SMALL) {
                printf("$ERROR - most likely your step in the z-direction (%f) is too small\n",N3Z1);
                exit(24);
             }
             if(Rabs(N3X1) > SMALL || Rabs(N3Y1) > SMALL) {
                printf("$ERROR - first input has to be pure z-axis (x: %f , y: %f)\n",
                N3X1,N3Y1);
                exit(25);
             }

             IA = (int *)malloc(sizeof(int) * Natoms);
             if(IA == NULL) exit(10);

             Chgn = (float *)malloc(sizeof(float) * Natoms);
             if(Chgn == NULL) exit(11);

             XC = (float *)malloc(sizeof(float) * Natoms);
             if(XC == NULL) exit(12);

             YC = (float *)malloc(sizeof(float) * Natoms);
             if(YC == NULL) exit(13);
             ZC = (float *)malloc(sizeof(float) * Natoms);
             if(ZC == NULL) exit(14);

/* atoms ... */

             printf("Atoms...\n");
             for(i = 0 ; i < Natoms ; i++)  {
                 fscanf(Input_p,"%d %f %f %f %f",&IA[i],&Chgn[i],&XC[i],&YC[i],&ZC[i]);
                 XC[i] *= BOHR_RADIUS;
                  YC[i] *= BOHR_RADIUS;
                   ZC[i] *= BOHR_RADIUS;
                   printf("Atomic number: %d, charge: %f, coord (x,y,z): %f %f %f\n",
                           IA[i],Chgn[i],XC[i],YC[i],ZC[i]);

             }

             if(DataType == 1) {    /* 1 */
/* molecular orbitals to handle ... */

                Data = (float *)malloc(sizeof(float) * N1 * N2 * N3 * MolecularOrbitalsDefined);
                if(Data == NULL) exit(17);

                for(i = 0   ; i < N1 ; i++)  { /* 2 */
                 for(j = 0  ; j < N2 ; j++)  { /* 3 */
                  for(k = 0 ; k < N3 ; k++)  { /* 4 */
                   for(l = 0; l < MolecularOrbitalsDefined ; l++) { /* 5 */

                       ijkl = i + N1 * j + N1 * N2 * k + N1 * N2 * N3 * l;

                       ijk = fscanf(Input_p,"%f",&Data[ijkl]);

                       if(ijk != 1) {
                          printf("$ERROR - in reading the grid data\n");
                          printf("$ERROR - at ijkl: %d %d %d %d\n",i,j,k,l);
                          exit(18);
                       }
                   }/* 5 */
                  } /* 4 */
                 }  /* 3 */
                }   /* 2 */
             } /* end *1* */
             else if((DataType == 2) || 
                     (DataType == 3) || 
                     (DataType == 5) ||
                     (DataType == 6)) {  /* start *1* */

				 if(MolecularOrbitalsDefined) {
                  Data = (float *)malloc(sizeof(float) * N1 * N2 * N3 * MolecularOrbitalsDefined);
                  if(Data == NULL) exit(17);

                    for(i = 0   ; i < N1 ; i++)  { /* 2 */
                     for(j = 0  ; j < N2 ; j++)  { /* 3 */
                      for(k = 0 ; k < N3 ; k++)  { /* 4 */
                       for(l = 0; l < MolecularOrbitalsDefined ; l++) { /* 5 */

                           ijkl = i + N1 * j + N1 * N2 * k + N1 * N2 * N3 * l;

                           ijk = fscanf(Input_p,"%f",&Data[ijkl]);

                           if(ijk != 1) {
                              printf("$ERROR - in reading the grid data\n");
                              printf("$ERROR - at ijkl: %d %d %d %d\n",i,j,k,l);
                              exit(18);
                           }
                       }/* 5 */
                      } /* 4 */
                     }  /* 3 */
                    }   /* 2 */
                 } else {
                   Data = (float *)malloc(sizeof(float) * N1 * N2 * N3);
                   if(Data == NULL) exit(19);

                   for(i = 0   ; i < N1 ; i++)  { /* 2 */
                    for(j = 0  ; j < N2 ; j++)  { /* 3 */
                     for(k = 0 ; k < N3 ; k++)  { /* 4 */

                         ijk = i + N1 * j + N1 * N2 * k;

                         ijkl = fscanf(Input_p,"%f",&Data[ijk]);
                         if(ijkl != 1) {
                            printf("$ERROR - in reading the grid data\n");
                            printf("$ERROR - at ijk: %d %d %d\n",i,j,k);
                            exit(20);
                         }
                     } /* 4 */
                    }  /* 3 */
                   }   /* 2 */
				 }
             } /* end *1* */
			 else {
				 printf("\n** ERROR ** this type '%d' not yet working\n",DataType);
                 DataType = 0;
                 continue;
			 }

            strncpy(InputFileTemp,InputFile,BUFF_LEN);
            MakeOutputFileName(InputFileTemp);
            if(ProduceCoordinateFile)
               printf("Coordinate file (in CHARMM 'crd' format): '%s'\n",CoordinateFile);
            else
               printf("No coordinate file will be written\n");
            if(ProduceCoordinateFile)
               (void)WriteCoordinateFile();
            (void)WriteInputData();

/* Reset DataType ! */
			DataType = 0;

            FileLoop++;
          }
    }

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
	int nn;
    float Help1;
    float Help2;
    float MinV;
	float MaxV;
	char  Temp[BUFF_LEN];
    char  OutFile[BUFF_LEN];

	if(DataType == 1) {

	   for(nn=0; nn <  MolecularOrbitalsDefined; nn++)
       {
          sprintf(Temp,"%d",MolecularOrbital_p[nn]);
		  sprintf(OutFile,OutputFile,Temp);
	      Output_p = fopen(OutFile,"wb");
	      if(Output_p == NULL) {
		   printf("$ERROR - can't open output file: '%s' \n",OutFile);
		   return(1);
	      }

		  printf("Saving file: '%s'\n",OutFile);

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

          for(k = 0   ; k < N3 ; k++)  {  /* 2 */
           for(j = 0  ; j < N2 ; j++)  {  /* 3 */
            for(i = 0 ; i < N1 ; i++)  {  /* 4 */

              ijk = i + N1 * j + N1 * N2 * k + N1 * N2 * N3 * nn;

              if(Data[ijk] < MinV) MinV = Data[ijk];
              if(Data[ijk] > MaxV) MaxV = Data[ijk];

              FWRITE(&Data[ijk] , sizeof(float));
            }                             /* 4 */
           }                              /* 3 */
          }                               /* 2 */
		   printf("Min value: %f\n",MinV);
           printf("Max value: %f\n",MaxV);
           fclose(Output_p);
        }
	} else if((DataType == 2) || 
              (DataType == 3) || 
              (DataType == 5) ||
              (DataType == 6)) {

	  if(MolecularOrbitalsDefined) {

	   for(nn=0; nn <  MolecularOrbitalsDefined; nn++)
       {
          sprintf(Temp,"%d",MolecularOrbital_p[nn]);
		  sprintf(OutFile,OutputFile,Temp);
	      Output_p = fopen(OutFile,"wb");
	      if(Output_p == NULL) {
		   printf("$ERROR - can't open output file: '%s' \n",OutFile);
		   return(1);
	      }

		  printf("Saving file: '%s'\n",OutFile);

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

          for(k = 0   ; k < N3 ; k++)  {  /* 2 */
           for(j = 0  ; j < N2 ; j++)  {  /* 3 */
            for(i = 0 ; i < N1 ; i++)  {  /* 4 */

              ijk = i + N1 * j + N1 * N2 * k + N1 * N2 * N3 * nn;

              if(Data[ijk] < MinV) MinV = Data[ijk];
              if(Data[ijk] > MaxV) MaxV = Data[ijk];

              FWRITE(&Data[ijk] , sizeof(float));
            }                             /* 4 */
           }                              /* 3 */
          }                               /* 2 */
		   printf("Min value: %f\n",MinV);
           printf("Max value: %f\n",MaxV);
           fclose(Output_p);
        }
	  } else {
          sprintf(Temp,"Total");
		  sprintf(OutFile,OutputFile,Temp);
	      Output_p = fopen(OutFile,"wb");
	      if(Output_p == NULL) {
		   printf("$ERROR - can't open output file: '%s' \n",OutFile);
		   return(1);
	      }

		  printf("Saving file: '%s'\n",OutFile);

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

          for(k = 0   ; k < N3 ; k++)  {  /* 2 */
           for(j = 0  ; j < N2 ; j++)  {  /* 3 */
            for(i = 0 ; i < N1 ; i++)  {  /* 4 */

              ijk = i + N1 * j + N1 * N2 * k;

              if(Data[ijk] < MinV) MinV = Data[ijk];
              if(Data[ijk] > MaxV) MaxV = Data[ijk];

              FWRITE(&Data[ijk] , sizeof(float));
            }                             /* 4 */
           }                              /* 3 */
          }                               /* 2 */
		   printf("Min value: %f\n",MinV);
           printf("Max value: %f\n",MaxV);
           fclose(Output_p);
	  }
	}
   return(0);
}

/**************************************************************************/
int WriteCoordinateFile()
/**************************************************************************/
{
    int   i,j;
    FILE *coord_p;
    char  AtomName[3];

    coord_p = fopen(CoordinateFile,"w");
    if(coord_p == NULL) {
      printf("$ERROR - can't open output file: '%s'",CoordinateFile);
      exit(1);
    }

    fprintf(coord_p,"* ++Automatic output generated from PC Gamess ++\n");
    fprintf(coord_p,"* ++using the '$cube' keyword                 ++\n");

    for(i = 0 ; i < TitleLines ; i++)
        fprintf(coord_p,"* %.132s\n",TitleText[i]);

    fprintf(coord_p,"*  \n");
    fprintf(coord_p,"%5d \n",Natoms);

    for(i = 0 ; i < Natoms ; i++) {

      for(j = 0 ; j < strlen(AtomSymbols)/4 ; j++) {
         if(IA[i] == AtomSymbol_p[j]) {
         strncpy(AtomName,&AtomSymbols[4 * j],2);
         AtomName[2] = '\0';
         break;}
      }
      fprintf(coord_p,
      "%5d%5d %-4.4s %-4.4s%10.5f%10.5f%10.5f %4.4s %-4d%10.5f \n",
      (i+1),1,"GAME",AtomName,XC[i],YC[i],ZC[i],"GAME",1,0.0);
      }

    fclose(coord_p);

    return(0);

}

