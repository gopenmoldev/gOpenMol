/*
From: Doug Moffatt <moffatt@ned1.sims.nrc.ca>


Doug Moffatt
Technical Officer
National Research Council of Canada
100 Sussex Driver 
Ottawa, Ontario
K1A 0R6

Tel.(613) 991-6376
Fax.(613) 954-8902


Hello Leif,

I am attaching the modified gcub2pl program; I call it g94cub2pl. It can
be as much as 10-15 times faster than gcub2pl if the user wishes to convert
many orbitals in the cube file. I do not have access to gopenmol on unix at the
moment so it has only been tested on nt. Maybe you could give it a quick test
on unix. I also include a small utility swbytes that
reverses the order of 4-byte units so that plt files generated on unix systems
can be viewed with the nt/95 vertion of gopenmol and vise-versa.
Do you think you will be distributing the source code soon. I just require the
algorithm to convert vertices into suitable triangles/quads to display the
orbital as a surface.

Regards,
Doug


                       Copyright (c) 1995 by:
        Leif Laaksonen , Centre for Scientific Computing , ESPOO, FINLAND
            Confidential unpublished property of Leif Laaksonen
                        All rights reserved


      This is a program to convert the output using the "cube"
      command from the Gaussian program to a plot format recognized
      by gOpenMol.

	Modified by Doug Moffatt, National Research Council of Canada
	This version is optimized to translate all orbitals in a gaussian cube file.


      Run this program on the output from GaussianXX using the
		'cube' keyword in GaussianXX.

      This program produces a coordinate file and a plot file, 
      that can be read into SCARECROW or gOpenMol.

      The running of this program is quite obvious (see later in this
      file).

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 This program converts the Gaussian output using the 'cube'
 command to a form understandable to gopenmol.
 Usage:
 g94cub2pl input output 
 Options:  -s  produce signed output (default)
           -a  produce absolute value output
	       -q  produce squared output	          

The program assumes an extension of .cub on input and appends .plt to the output name.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      If you need any help please feel free to contact:

      Leif.Laaksonen@csc.fi


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Taken from the;
                 Gaussian 94 User's Reference Manual
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

> Cube keyword


DESCRIPTION

The Cube properties keyword can be used to evaluate molecular orbitals, the
electrostatic potential, the electron density, density gradient, the norm of 
the density gradient, and Laplacian of the density over a 3 dimensional grid 
(cube) of points. By default, Cube evaluates the electron density (corresponding to
the Density option). Which density is used is controlled by the Density keyword; 
use Density=Current to evaluate the cube over the density from a correlated or
Cl-Singles wavefunction instead of the default Hartree-Fock density.

Note that only one of the available quantities can be evaluated within any one job 
step. Save the checkpoint file (using %Chk), and include Guess=(Reod,Only) 
Density=Checkpoint in the route section of a subsequent job (or job step) in 
order to evaluate a different quantity without repeating any of the other steps
 of the calculation.

Gaussian 94 provides reasonable defaults for grids, so Cube no longer requires
that the cube be specified by the user. However, the output filename must still
always be provided (see below).

Alternatively, Cube may be given a parameter specifying the number of points to
use per "side" (the default is 80). For example, Cube=100 specifies a grid of 
1,000,000 points ( 100 ), evenly distributed over the rectangular grid generated
by the program (which is not necessarily a cube). In addition, the input
format used by earlier versions of Gaussian is still supported; Cube=Cards
indicates that a grid will be input. It may be used to specify a grid of
arbitrary size and shape.

The files create by Cube can be manipulated using the cubman utility, described
in chapter 5.

Note that Pop=None will inhibit cube file creation.

INPUT FORMAT

When the user elects to provide it, the grid information is read from the input
stream. The first line-required for all Cube jobs-gives a file name for the cube
file. Subsequent lines, which are included only with Cube=Cards, must conform to
format (15,3F12.6), according to the following syntax:

Output-file-name              Required in all Cube jobs.
IFlag, X0, Y0, Z0             Output unit number and initial point.
N1, X1, Y1, Z1                Number of points and step-size in the X-direction.
N2, X2, Y2, Z2                Number of points and step-size in the Y-direction.
N3, X3, Y3, Z3                Number of points and step-size in the Z-direction.

If IFlag is positive, the output file is unformatted; if it is negative, 
the output file is formatted. If N1<O the input cube coordinates are assumed
to be in Bohr, otherwise, they are interpreted as Angstroms (|N1| is used as
the number of X-direction points in any case). Note that the three axes are 
used exactly as specified; they are not orthogonalized, so the grid need not 
be rectangular.

If the Orbitals option is selected, the cube filename (or cube filename and 
cube specification input) is immediately followed by a list of the orbitals 
to evaluate, in free-format, terminated by a blank line. In addition to 
numbers for the orbitals (with alpha orbitals numbered starting at N+l), the
following abbreviations can appear in the list:

HOMO             The highest occupied molecular orbital
LUMO             The lowest unoccupied molecular orbital
OCCA             All occupied (alpha) orbitals
OCCB             All beta occupied orbitals for UHF
ALL              All orbitals
VALENCE          All occupied non-core orbitals
VIRTUALS         All virtual orbitals


OUTPUT FILE FORMATS

Using the default input to Cube produces an unformatted output file (you can
use the cubman utility to convert it to a formatted version if you so desire;
see chapter 5). When the Cards option is specified, then the IFlag parameter's
sign determines the output file type. If IFlag>0, the output is unformatted. If
IFlag<0, the output is formatted. All values in the cube file are in atomic units,
regardless of the input units.

For density and potential grids, unformatted files have one row per record
(i.e., N1 * N2 records each of length N3). For formatted output, each row is
written out in format (6E13.5). In this case, if N3 is not a multiple of six,
then there may be blank space in some lines.

The norm of the density gradient is also a scalar (i.e., one value per point),
and is written out in the same manner. Density+gradient grids are similar, but
with two writes for each row (of lengths N3 and 3*N3). Density+gradient+Laplacian
grids have 3 writes per row (of lengths N3, 3*N3, and N3).

For example, for a density cube, the output file looks like this:

NAtoms, X-Origin, Y-Origin, Z-Origin
N1, X1, Y1, Z1              # of increments in the slowest running direction
N2, X2, Y2, Z2
N3, X3, Y3, Z3              # of increments in the fastest running direction
IA1, Chgl, X1, Y1, Z1       Atomic number, charge, and coordinates of thefirst atom
.....
IAn, Chgn, Xn, Yn, Zn       Atomic number, charge, and coordinates of the last atom

(N1 * N2) records, each of length N3  Values of the density at each point in the grid

Note that a separate write is used for each record.


For molecular orbital output, NAtoms will be less than zero, and an additional
record follows the data for the final atom (in format lOI5 if the file is formatted):

NMO,  ( MO ( I ), I = 1, NMO )                      Number of MOs and their numbers

If NMO orbitals were evaluated, then each record is NMo * N3 long and has the
values for all orbitals at each point together.

READING CUBE FILES WITH FORTRAN PROGRAMS

If one wishes to read the values of the density or potential back into an
array dimensioned X(N3,N2,N1) code like the following Fortran loop may be used:

		Do 10 I1 = 1,N1
		Do 10 I2 = 1,N2
			Read(n,'(6E13.5)') (X(I3,I2,I1),I3=1,N3)
10    Continue

where n is the unit number corresponding to the cube file.

If the origin is (X0,Y0,Z0), and the increments (X1,Y1,Z1), then point
(I1,I2,I3) has the coordinates:

X-co0rdinate: X0+(I1-1)*X1+(I2-1)*X2+(I3-1)* X3
Y-coordinate: Y0+(I1-1)*Y1+(I2-1)*Y2+(I3-1)* Y3
Z-coordinate: Z0+(I1-1)*Z1+(I2-1)*z2+(I3-1)* Z3

The output is similar if the gradient or gradient and Laplacian of the charge
density are also requested, except that in these cases there are two or three
records, respectively, written for each pair of I1, I2 values. Thus, if the
density, gradient, and Laplacian are to be read into arrays D(N3,N2,N1),
G(3,N3,N2,N1), RL(N3,N2,N1) from a formatted output file, a correct set of
Fortran loops would be:

		Do 10 I1 = 1, N1
		Do 10 I2 = 1, N2
		  Read(n,'(6F13.5)') (D(I3,I2,I1),I3=1,N3)
		  Read(n,'(6F13.5)') ((G(IXYZ,I3,I2,I1),IXYZ=1,3), I3=1,N3)
		  Read(n,'(6F13.5)') (RL(I3,I2,I1),I3=1,N3)
10    Continue

where again n is the unit number corresponding to the cube file.


OPTIONS

Density      Compute just the density values. This is the default.

Potential    Compute the electrostatic potential at each point.

Gradient     Compute the density and gradient.

Laplacian    Compute the density, gradient, and Laplacian ofthe density
				 ((nabla)2(rho)).

NormGradient Compute the norm of the density gradient at each point.

Orbitals     Compute the values of one or more molecular orbitals at each point.

FrozenCore   Remove the SCF core density. This is the default for the density,
             and is not allowed for the potential.

Full         Evaluate the density including all electrons.

Total        Use the total density. This is the default.

Alpha        Use only the alpha spin density.

Beta         Use only the beta spin density.

Spin         Use the spin density (difference between alpha and beta densities).

Cards        Read grid specification from the input stream (as described above).

Example 1: Produce the total electron density.

x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
#p rhf/6-31g* 5d test geom=modela  cube=(density,read) FormCheck=OptCart

Gaussian Test Job 257:
Density cube

0,1
o h f

input.cube
  -51        -2.0        -2.0        -1.0
	40         0.1         0.0         0.0
   40         0.0         0.1         0.0
   20         0.0         0.0         0.1
	   
x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-


Example 2: Produce the data to plot molecular orbitals

x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-
#p rhf/6-31g* 5d test geom=modela  cube=(orbitals) FormCheck=OptCart

Gaussian Test Job 257:
Density cube

0,1
o h f

input.cube
  -51        -2.0        -2.0        -1.0
   40         0.1         0.0         0.0
   40         0.0         0.1         0.0
   20         0.0         0.0         0.1
ALL

x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-


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
#define GAUSSIAN_TYPE   200
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
"This program converts all orbitals in a Gaussian cube file\n"
"to a form understandable to gopenmol.\n"
"Usage:\n"
"   g94cub2pl inputfile outputfile -[sqa]\n"
"program supplies input extension .cub output extension .plt\n"
"options: s - signed, q - squared, a - abs val(default)\n";
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
float *Data,*Datas;                   /* data */
int    MolecularOrbitals = 0;  /* switch to handle molecular orbitals */
int   *MolecularOrbital_p;     /* index pointer */
int    MolecularOrbitalsDefined; /* orbitals defined in file */

/* output */
float Xmax,Ymax,Zmax;
int   TypeOfData = GAUSSIAN_TYPE;

char InputFile[BUFF_LEN];
char OutputFile[BUFF_LEN];
char CoordinateFile[BUFF_LEN];

/* functions */
void MakeOutputFileName(char *);
int  ReadInputData(void);
int  WriteInputData(void);
int  WriteCoordinateFile(void);
int  outfmt=0;
char *outfmtstr[3] = {"abs. value","squared","signed"};
/* externals */
extern char *Number2Name(int);

/**************************************************************************/
main(int args, char** argv)
/**************************************************************************/
{
	 int i;

	 printf("**********************************************************\n");
	 printf("* Convert a 'cube' output from GaussianXX into the plot  *\n");
	 printf("* format known by gOpenMol (or SCARECROW).               *\n");
	 printf("*                                                        *\n");
	 printf("* Leif Laaksonen (CSC) 1995 - 2003                       *\n");
	 printf("* Email: Leif.Laaksonen@csc.fi                           *\n");
	 printf("*                               Version: 23/09/03        *\n");
	 printf("* modified by Doug Moffatt, National Research Council of *\n");
	 printf("* Canada, December, 1997                                 *\n");
	 printf("**********************************************************\n\n");

	 if(args < 3) {
	 printf("%s",Usage);
	 exit(0);}

		strncpy(OutputFile,argv[2],BUFF_LEN);
/*		sprintf(InputFile,"%s.cub",argv[1]); LUL1999-05-04*/
        sprintf(InputFile,"%s",argv[1]);

/* default signed */
    outfmt = 2;
	if(args > 3) {
     if(argv[3][0] == 'a') outfmt = 0;
	 if(argv[3][0] == 'q') outfmt = 1;
	 if(argv[3][0] == 's') outfmt = 2;
	}
	 printf("File names:\n");
	 printf("Input file:      '%s'\n",InputFile);

	 MakeOutputFileName(InputFile);

	 printf("Output file (plot file):  '%s'\n",OutputFile);
     printf("Coordinate file (in CHARMM 'crd' format): '%s'\n",CoordinateFile);

		/*  Process the data ... */

	 if(ReadInputData()) {
		 printf("$ERROR - can't read input data\n");
		 exit(1);
	 }
	 (void)WriteCoordinateFile();
	 (void)WriteInputData();

	 printf(" %s data stored\n",outfmtstr[outfmt]);

	 printf("Job done ...\n");
}
/**************************************************************************/
void MakeOutputFileName(char *Filename)
/**************************************************************************/
{
	  int i;
      int hit;

      if(OutputFile[0] != '\0') {
          hit = 0;
          for(i = 0 ; i < strlen(OutputFile) ; i++) {
			if(OutputFile[i] == '.') hit = i;
          }
          if(hit) {
			  sprintf(&OutputFile[hit],"\0");
			  strncpy(CoordinateFile,OutputFile,hit);
			  sprintf(&CoordinateFile[hit],".crd\0");
			  return;
          } else {
			  sprintf(CoordinateFile,"%s.crd\0",OutputFile);
              return;
          }
      }

/* look for a dot in the name ... */
      hit = 0;
	  for(i = 0 ; i < strlen(Filename) ; i++) {
			if(Filename[i] == '.') hit = i;
      }
      if(hit) {
         strncpy(OutputFile,Filename,hit);
         sprintf(&OutputFile[hit],"\0");
         strncpy(CoordinateFile,Filename,hit);
         sprintf(&CoordinateFile[hit],".crd\0");
         return;
      } else {
         sprintf(OutputFile,"%s\0",Filename);
         sprintf(CoordinateFile,"%s.crd\0",Filename);
      }
}

/**************************************************************************/
int ReadInputData()
/**************************************************************************/
{
	FILE *Input_p;
        char  Text[BUFF_LEN];
	int   i,j,k,ijk,ijkl,c0,c1,c2;
	register int l;
	int   IsInteger;
	char  Temp1[BUFF_LEN];
	char  Temp2[BUFF_LEN];
	char  Temp3[BUFF_LEN];
	char  Temp4[BUFF_LEN];
    int   Hit,kk;
    int   tatoms;
    float txorig;
    float tyorig;
    float tzorig;

    Input_p = fopen(InputFile , "r");
    if(Input_p == NULL) {
		printf("$ERROR - can't open input file '%s'\n",InputFile);
		return(1);
	 }

	printf("\nTitle in file (job title):\n");
/* first comes an unknow number of title lines (MAX lines is 5) */
	 TitleLines = 0;
	 for(i = 0 ; i < MAX_TITLE_LINES ; i++) {
		 fgets(Text,BUFF_LEN,Input_p);
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
		 TitleText[i][strlen(TitleText[i]) - 1] = '\0';
		 TitleLines++;}
	else
		 break;
	}

/* now starts the data ... */
	sscanf(Text,"%d %f %f %f",&Natoms,&Xorig,&Yorig,&Zorig);
	 Xorig *= BOHR_RADIUS;
	  Yorig *= BOHR_RADIUS;
		Zorig *= BOHR_RADIUS;
    printf("Number of atoms: %d, x-, y-, z-origin (in Anstrom): %f,%f,%f\n",
            abs(Natoms),Xorig,Yorig,Zorig);

/* check to correct for defect data from adf2cube program */ 
    if(!Natoms) { 
      printf("Your number of atoms is 0!\n");
      sscanf(Text,"%s",Temp1);
      if(Temp1[0] == '-') {
        MolecularOrbitals = 1;
      } else {
        MolecularOrbitals = 0;
      }
     } else {      

      MolecularOrbitals                = 0;
      if(Natoms < 0) MolecularOrbitals = 1;
    }

   Natoms = abs(Natoms);

	fscanf(Input_p,"%d %f %f %f",&N1,&N1X1,&N1Y1,&N1Z1);
	 N1X1 *= BOHR_RADIUS;
     N1Y1 *= BOHR_RADIUS;
      N1Z1 *= BOHR_RADIUS;
	 printf("Number of points: %d, in direction (x,y,z) %f %f %f\n",
            N1,N1X1,N1Y1,N1Z1);

/* check that this is a 'pure' x-coordinate */
   if(Rabs(N1X1) < SMALL) {
     printf("$ERROR - most likely your step in the x-direction (%f) is too small\n",N1X1);
     exit(20);}
	if(Rabs(N1Y1) > SMALL || Rabs(N1Z1) > SMALL) {
     printf("$ERROR - first input has to be pure x-axis (y: %f , z: %f)\n",
				N1Y1,N1Z1);
     exit(21);}

	fscanf(Input_p,"%d %f %f %f",&N2,&N2X1,&N2Y1,&N2Z1);
    N2X1 *= BOHR_RADIUS;
     N2Y1 *= BOHR_RADIUS;
		N2Z1 *= BOHR_RADIUS;
	 printf("Number of points: %d, in direction (x,y,z) %f %f %f\n",
            N2,N2X1,N2Y1,N2Z1);

/* check that this is a 'pure' y-coordinate */
	if(Rabs(N2Y1) < SMALL) {
     printf("$ERROR - most likely your step in the y-direction (%f) is too small\n",N2Y1);
     exit(22);}
	if(Rabs(N2X1) > SMALL || Rabs(N2Z1) > SMALL) {
     printf("$ERROR - second input has to be pure y-axis (x: %f , z: %f)\n",
				N2X1,N2Z1);
	  exit(23);}

	fscanf(Input_p,"%d %f %f %f",&N3,&N3X1,&N3Y1,&N3Z1);
    N3X1 *= BOHR_RADIUS;
	  N3Y1 *= BOHR_RADIUS;
      N3Z1 *= BOHR_RADIUS;
	 printf("Number of points: %d, in direction (x,y,z) %f %f %f\n",
            N3,N3X1,N3Y1,N3Z1);

/* check that this is a 'pure' z-coordinate */
   if(Rabs(N3Z1) < SMALL) {
     printf("$ERROR - most likely your step in the z-direction (%f) is too small\n",N3Z1);
     exit(24);}
	if(Rabs(N3X1) > SMALL || Rabs(N3Y1) > SMALL) {
     printf("$ERROR - first input has to be pure z-axis (x: %f , y: %f)\n",
				N3X1,N3Y1);
     exit(25);}

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
	if(MolecularOrbitals) {    /* 1 */
/* molecular orbitals to handle ... */
	fscanf(Input_p,"%d",&MolecularOrbitalsDefined);

	MolecularOrbital_p = (int *)malloc(sizeof(int) * MolecularOrbitalsDefined);
	if(MolecularOrbital_p == NULL) exit(15);

	for(i = 0 ; i < MolecularOrbitalsDefined ; i++) {
		 fscanf(Input_p,"%d",&MolecularOrbital_p[i]);
	}

	Data = (float *)malloc(sizeof(float) * N1 * N2 * N3 *
														MolecularOrbitalsDefined);
        if(Data == NULL) exit(17);
	c0 = N1*N2*N3;
	for(i = 0   ; i < N1 ; i++)  { /* 2 */
	 for(j = 0  ; j < N2 ; j++)  { /* 3 */
          c1 = i+N1*j;
	  for(k = 0 ; k < N3 ; k++)  { /* 4 */
	   c2 = N1*N2*k+c1;
	   Datas = &Data[c2];
	   for(l = 0; l < MolecularOrbitalsDefined ; l++) { /* 5 */

	  ijk = fscanf(Input_p,"%13e",Datas);
	  if(ijk != 1) {
		 printf("$ERROR - in reading the grid data\n");
		 printf("$ERROR - at ijkl: %d %d %d %d\n",i,j,k,l);
		 exit(18);
	  }
	  Datas += c0;
		 }/* 5 */
		} /* 4 */
	  }  /* 3 */
	if(!(i%10)) printf("%d of %d read\n",i+1,N1);

	 }   /* 2 */
	} /* end *1* */
	else {  /* start *1* */

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
	printf("%d of %d read\n",i+1,N1);
	}   /* 2 */
  } /* end *1* */
	fclose(Input_p);

	return(0);
}

/**************************************************************************/
int WriteInputData()
/**************************************************************************/
{
	 FILE *Output_p;
	 int   i,j,k,l,ijk,ijkl,nn;
	 int Items;
	 float Help1;
	 float Help2;
	 float MinV;
	float MaxV;
	char namo[256];

/* density */
    if(!MolecularOrbitalsDefined) {
        
	 sprintf(namo,"%s.plt",OutputFile);
	 Output_p = fopen(namo,"wb");
	 if(Output_p == NULL) {
		printf("$ERROR - can't open output file: '%s' \n",OutputFile);
		return(1);
	 }

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
	if(!outfmt) 
	for(k = 0   ; k < N3 ; k++)  {  /* 2 */
	 for(j = 0  ; j < N2 ; j++)  {  /* 3 */
	  for(i = 0 ; i < N1 ; i++)  {  /* 4 */
	  ijk = i + N1 * j + N1 * N2 * k;
	  if(Data[ijk] < 0) Data[ijk] = -Data[ijk];
	  FWRITE(&Data[ijk] , sizeof(float));
	  }                             /* 4 */
	 }                              /* 3 */
	}
											 /* 2 */
	if(outfmt == 1)  
	for(k = 0   ; k < N3 ; k++)  {  /* 2 */
	 for(j = 0  ; j < N2 ; j++)  {  /* 3 */
	  for(i = 0 ; i < N1 ; i++)  {  /* 4 */
	  ijk = i + N1 * j + N1 * N2 * k;
		Data[ijk] *= Data[ijk];
	  FWRITE(&Data[ijk] , sizeof(float));
	  }                             /* 4 */
	 }                              /* 3 */
	}
	if(outfmt == 2)   
	for(k = 0   ; k < N3 ; k++)  {  /* 2 */
	 for(j = 0  ; j < N2 ; j++)  {  /* 3 */
	  for(i = 0 ; i < N1 ; i++)  {  /* 4 */
	  ijk = i + N1 * j + N1 * N2 * k;
	  FWRITE(&Data[ijk] , sizeof(float));
	  }                             /* 4 */
	 }                              /* 3 */
	}


        fclose(Output_p);
        printf("%s written\n",namo);

      return(0);
   }

/* orbital/orbitals */
	for(nn=0; nn <  MolecularOrbitalsDefined; nn++)
     {
	 sprintf(namo,"%s%d.plt",OutputFile,MolecularOrbital_p[nn]);
	 Output_p = fopen(namo,"wb");
	 if(Output_p == NULL) {
		printf("$ERROR - can't open output file: '%s' \n",OutputFile);
		return(1);
	 }

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
	if(!outfmt) 
	for(k = 0   ; k < N3 ; k++)  {  /* 2 */
	 for(j = 0  ; j < N2 ; j++)  {  /* 3 */
	  for(i = 0 ; i < N1 ; i++)  {  /* 4 */
	  ijk = i + N1 * j + N1 * N2 * k + N1 * N2 * N3 * nn;
	  if(Data[ijk] < 0) Data[ijk] = -Data[ijk];
	  FWRITE(&Data[ijk] , sizeof(float));
	  }                             /* 4 */
	 }                              /* 3 */
	}
											 /* 2 */
	if(outfmt == 1)  
	for(k = 0   ; k < N3 ; k++)  {  /* 2 */
	 for(j = 0  ; j < N2 ; j++)  {  /* 3 */
	  for(i = 0 ; i < N1 ; i++)  {  /* 4 */
	  ijk = i + N1 * j + N1 * N2 * k + N1 * N2 * N3 * nn;
		Data[ijk] *= Data[ijk];
	  FWRITE(&Data[ijk] , sizeof(float));
	  }                             /* 4 */
	 }                              /* 3 */
	}
	if(outfmt == 2)   
	for(k = 0   ; k < N3 ; k++)  {  /* 2 */
	 for(j = 0  ; j < N2 ; j++)  {  /* 3 */
	  for(i = 0 ; i < N1 ; i++)  {  /* 4 */
	  ijk = i + N1 * j + N1 * N2 * k + N1 * N2 * N3 * nn;
	  FWRITE(&Data[ijk] , sizeof(float));
	  }                             /* 4 */
	 }                              /* 3 */
	}


	fclose(Output_p);
        printf("%s written\n",namo);
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

    if(!Natoms) {
      printf("No atom information available\n");
      return(1);
    }

	 coord_p = fopen(CoordinateFile,"w");
	 if(coord_p == NULL) {
		printf("$ERROR - can't open output file: '%s'",CoordinateFile);
		exit(1);
	 }

	 fprintf(coord_p,"* ++Automatic output generated from Gaussian++\n");
	 fprintf(coord_p,"* ++using the 'cube' keyword                ++\n");

	 for(i = 0 ; i < TitleLines ; i++)
		  fprintf(coord_p,"* %.70s\n",TitleText[i]);

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
      (i+1),1,"GAUS",AtomName,XC[i],YC[i],ZC[i],"GAUS",1,0.0);
      }

    fclose(coord_p);

    return(0);

}

#ifdef SWAP
/* swap bytes routine .... */
#include <stdio.h>
void swabi(void *,int);
void swordi(void *,int);
main (argc,argv)
int argc;
char **argv;
{
FILE *fpi,*fpo;
int i,j;
int buf[2048];
if(argc < 3) {printf("Reverses the order of 4-byte units.\n Usage: swbytes inputname outputname \n");
exit(0);}
fpi = fopen(argv[1],"rb");
fpo = fopen(argv[2],"wb");
i=0;
while(1)
{
j = fread(buf,4,1024,fpi);
if(j <= 0) break;
swabi(buf,j*4);
swordi(buf,j*2);
fwrite(buf,j,4,fpo);
i++;
if(!(i%100)) printf("%d\n",i);
}
fclose(fpi);
fclose(fpo);
}

/* swabi - - swap bytes in place */
void swabi(void *buf, int nbytes)
{
	long i,j;
	unsigned char *cbufx=buf,tmp;
	for(i=0; i < nbytes; i+=2)
	{tmp=cbufx[i]; cbufx[i]=cbufx[i+1]; cbufx[i+1]=tmp; }
}



/* swordi - - swap words in place */
void swordi(void *buf, int nwords)
{
	long i,j;
	unsigned short *cbufx=buf,tmp;
	for(i=0; i < nwords; i+=2)
	{tmp=cbufx[i]; cbufx[i]=cbufx[i+1]; cbufx[i+1]=tmp; }
}
#endif
