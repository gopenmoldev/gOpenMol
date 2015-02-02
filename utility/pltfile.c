
/*


      This is a small program to convert a plt file from
      unformatted to formatted and formatted to unformatted.

      Leif Laaksonen CSC 1996

      Compile:

      cc -O2 -o contman  contman.c

      Usage

      Unformatted to formatted:
      ./pltfile -uf -iinputfile1.plt -ooutputfile.txt

      Formatted to unformatted:
      ./pltfile -fu -iinputfile1.txt -ooutputfile.plt

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFF_LEN   256
#define BUFFER_LIM 100
#define STEP       1.e-05

#define Rabs(a)       ( ( a ) > 0.0 ? (a) : -(a))

#define MINUS      0
#define PLUS       1
#define NOTHING   -1


#define BUFF_LEN       256
#define MAX_CONT       100       /* maximum contouring levels               */
#define MAX_CONTOURS    50       /* max number of countours (files/entries) */
#define CONTOUR_TYPE_SOLID    0  /* display object as a solid surface       */
#define CONTOUR_TYPE_MESH     1  /* display object as a mesh surface        */

  FILE  *ISO_file;
      float ISO_xmin,ISO_xmax;
      float ISO_ymin,ISO_ymax;
      float ISO_zmin,ISO_zmax;
      float ISO_StepX;
      float ISO_StepY;
      float ISO_StepZ;
      int TypeOfSurface;        /* = 0 , unknown
                                   = 1 , VSS surface
                                   = 2 , orbital- or density surface
                                   = 3 , probe surface
                                */
char MY_NAME[80];
int VERBOSE = 0;

char   InputFile[BUFF_LEN];
FILE  *Input_p;
float  xmin,xmax;
float  ymin,ymax;
float  zmin,zmax;

char   OutputFile[BUFF_LEN];
FILE  *Output_p;

int Items;

int ContoursDefined;

int lulWriteContour(FILE *Output_p);

int      action = -1;
int      rank;
int      TypeOfSurface;
float    Atm2Ang  = 0.52917715;
/* contour structure */

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Remember that this structure is also defined as extern in 

   contdriver.c and contour.c in the graphics directory 

   Real bad programming practice 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

     int    ContoursDefined = 0;   /* number of contours defined */

     struct Contour_Struct {
     float *data;           /* contour data */
     float min;
     float max;
     int xdim;
     int ydim;
     int zdim;
     float Xtrans;
     float Ytrans;
     float Ztrans;
     float Xscale;
     float Yscale;
     float Zscale;
     float ColVal[MAX_CONT];
     int   NumVal;
     char  ColNam[MAX_CONT][BUFF_LEN];
     float RedC[MAX_CONT];
     float GreenC[MAX_CONT];
     float BlueC[MAX_CONT];
     int   DisplayType[MAX_CONT];
     char  ContFile[BUFF_LEN];
     float AlphaBlend[MAX_CONT];
     int   ContSmooth[MAX_CONT];
     char  Name[BUFF_LEN];
     int   Display[MAX_CONT];} ContourInfo[MAX_CONTOURS];

     int   gomContourInternalType;

main(int argc , char **argv)
{
    int i,j;
    char OutText[BUFF_LEN];

    for(i = 1 ; i < argc ; i++) {

      if(argv[i][0] == '-')     {

        switch(argv[i][1])       {

/* -i input contour files */
	case 'i':          /* read contour files */
                  if(sscanf(&argv[i][2],"%s",InputFile) 
                          != 1) {
                          printf("?: problems getting input file\n");
                          exit(1);
			}
                   printf("Input file:     '%s'\n",InputFile);
                   break;
/* -o output contour file */
	case 'o':          /* output contour     */
                   strncpy(OutputFile,&argv[i][2],BUFF_LEN);
                   printf("Output file:    '%s'\n",OutputFile);
                   break;
/* action */
        case 'u': 
                   if(argv[i][2] == 'f') {
                     action = 0;
		   }
                   else {
                   printf("$ERROR - unrecognized option\n");
                   exit(1);
		   }
                   break;
        case 'f': 
                   if(argv[i][2] == 'u') {
                     action = 1;
		   }
                   else {
                   printf("$ERROR - unrecognized option");
                   exit(1);
		   }
                   break;
 
/* -h print help ...        */
		     case 'h':
                       ;
                       exit(0);
		     }
           }
       }

      if(action < 0) {
         printf("$ERROR - please define what you want to do (u <=> f or f <=> u\n");
         exit(1);
      }

      ContoursDefined = 0;

      if(!action) {
        Input_p = fopen(InputFile,"rb");
           if(Input_p == NULL) {
           printf("?: can't open input file: '%s'\n",InputFile);
           exit(1);
         }

        if(ContourDriverUNFORM(InputFile, "1"))
           exit(1);

      }
      else {
        Input_p = fopen(InputFile,"r");
           if(Input_p == NULL) {
           printf("?: can't open input file: '%s'\n",InputFile);
           exit(1);
         }

        (void)lulReadContourFORM(Input_p);
      }

      xmax = ISO_xmax;
      xmin = ISO_xmin;
      ymax = ISO_ymax;
      ymin = ISO_ymin;
      zmax = ISO_zmax;
      zmin = ISO_zmin;

      printf("Xpoints: %d, YPoints: %d, Zpoints: %d\n",
              ContourInfo[0].xdim,
              ContourInfo[0].ydim,
              ContourInfo[0].zdim);
      printf("Xmin: %f, Xmax: %f, Ymin: %f, Ymax: %f, Zmin: %f, Zmax: %f\n",
              xmin,xmax,ymin,ymax,zmin,zmax);

      if(!action) {
        Output_p = fopen(OutputFile,"w");
           if(Output_p == NULL) {
           printf("?: can't open output file : '%s'\n",OutputFile);
           exit(1);
        }
        (void)lulWriteContourFORM(Output_p);
      }
      else {
        Output_p = fopen(OutputFile,"wb");
           if(Output_p == NULL) {
           printf("?: can't open output file : '%s'\n",OutputFile);
           exit(1);
        }
        (void)lulWriteContourUNFORM(Output_p);
      }

    return(0);
 }

#define FWRITE(value_p , size)    { Items = \
                                 fwrite(value_p, size , 1L , Output_p);\
                                 if(Items < 1) {\
                     printf("?ERROR - in writing contour file (*)\n");\
                     return(1);}}

#define FWRITEN(value_p , num , size) { Items = \
                                fwrite(value_p, size , num , Output_p);\
                   if(Items < 1) {\
                     printf("?ERROR - in writing contour file (**)\n");\
                     return(1);}}


/**************************************************************************/
int lulWriteContourUNFORM(FILE *Output_p)
/**************************************************************************/
{
      static int i,j,k,Loop,Items;
      static float xc,yc,zc,pot;
      static float Help1,Help2;
      static float Buffer[BUFFER_LIM];
      static float MaxV = -1.0e10;
      static float MinV =  1.0e10;

/* print out number of points in x- , y- and z-directions */

      FWRITE(&rank , sizeof(int));
      FWRITE(&TypeOfSurface , sizeof(int));
      FWRITE(&ContourInfo[0].zdim , sizeof(int));
       FWRITE(&ContourInfo[0].ydim , sizeof(int));
        FWRITE(&ContourInfo[0].xdim , sizeof(int));

      Help1 = ISO_zmin;
       Help2 = ISO_zmax;
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));
      Help1 = ISO_ymin;
       Help2 = ISO_ymax;
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));
      Help1 = ISO_xmin;
       Help2 = ISO_xmax;
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));

/*

      fprintf(Output_p," 3         ! Rank value\n");
      fprintf(Output_p," %d        ! number of points in the z-direction\n",npts[2]);
      fprintf(Output_p," %d        ! number of points in the y-direction\n",npts[1]);
      fprintf(Output_p," %d        ! number of points in the x-direction\n",npts[0]);
      fprintf(Output_p," %f  %f    ! zmin and zmax\n",xmin[2]*autang,xmax[2]*autang);
      fprintf(Output_p," %f  %f    ! ymin and ymax\n",xmin[1]*autang,xmax[1]*autang);
      fprintf(Output_p," %f  %f    ! xmin and xmax\n",xmin[0]*autang,xmax[0]*autang);

*/
      MaxV = -1.0e10;
      MinV =  1.0e10;

      Loop = 0;

      for(i = 0 ; 
          i < ContourInfo[0].xdim * ContourInfo[0].ydim * ContourInfo[0].zdim;
          i++) {

            Help1 = ContourInfo[0].data[i];

            if(Help1 < MinV) MinV = Help1;
            if(Help1 > MaxV) MaxV = Help1;

              Buffer[Loop] = Help1;
          
              FWRITE(&Help1 , sizeof(float));
      }      

      printf("Max value: %f \n",MaxV);
      printf("Min value: %f \n",MinV);
      printf("! Done.\n");

      fclose(Output_p);
}


/************************************************************************/
int  ContourDriverUNFORM(char *ContFile, char *ContName)
/************************************************************************/
{
     static float *Cdata;
     static int    Cxdim;
     static int    Cydim;
     static int    Czdim;
     static float  Cmax;
     static float  Cmin;
     static int i,j;
     static char OutText[BUFF_LEN];


     if(ContoursDefined == MAX_CONTOURS) {
        printf("?ERROR - Max number of contour files defined\n");
        return(-1);}

/* get the data from file ... */

  if (lulGetContourData(ContFile,&Cdata,
                       &Cxdim,
                       &Cydim,
                       &Czdim) == -1) {
     printf("Error from get_data\n");
     return(-1);
  }

/* save grabbed information */
     ContourInfo[ContoursDefined].xdim = Cxdim;
      ContourInfo[ContoursDefined].ydim = Cydim;
       ContourInfo[ContoursDefined].zdim = Czdim;
        ContourInfo[ContoursDefined].data = Cdata;

  strncpy(ContourInfo[ContoursDefined].ContFile,ContFile,BUFF_LEN);

/*
  get_max_min(ContourInfo[ContoursDefined].data,
                               ContourInfo[ContoursDefined].xdim,
                               ContourInfo[ContoursDefined].ydim,
                               ContourInfo[ContoursDefined].zdim,
                               &Cmax,
                               &Cmin);

     ContourInfo[ContoursDefined].max = Cmax;
      ContourInfo[ContoursDefined].min = Cmin;

     printf("Minimum value %f maximum value %f\n",Cmin,Cmax);
*/
/* ready with the file ...  */

/* look if there are some contours already defined */

       for(i = 0 ; i < ContourInfo[ContoursDefined].NumVal ; i++) {
       ContourInfo[ContoursDefined].Display[i]      = 0;
        ContourInfo[ContoursDefined].ContSmooth[i]  = 0;
         ContourInfo[ContoursDefined].AlphaBlend[i] = 1.;}

        ContourInfo[ContoursDefined].NumVal         = 0;

     if(ContName[0] == '\0') { /* no contour name defined */
        sprintf(OutText,"%d",ContoursDefined+1);
        strncpy(ContourInfo[ContoursDefined].Name,OutText,BUFF_LEN);}
     else {
        strncpy(ContourInfo[ContoursDefined].Name,ContName,BUFF_LEN);}

     ContoursDefined++;

     return(0);
}

/**************************** get_data ****************************/
/**************************** get_data ****************************/
/**************************** get_data ****************************/
/**************************** get_data ****************************/

          /* This subroutine will read in the data file. */
int lulGetContourData(char *filename,
                      float **data, int *xdim, int *ydim, int *zdim)
{

  int shape[3];
  int size;
  int ret;
  char OutText[BUFF_LEN];

  ret=SCAREgetdims(filename,&rank,shape,3);
  if (ret != 0) {
     printf("%s: error from SCAREgetdims %d for %s\n",
             MY_NAME,ret,filename);
     return -1;
  }

  if (rank != 3) {
     printf("%s: error from SCAREgetdims %d for %s\n",
             MY_NAME,ret,filename);
     return -1;
  }

  *xdim = shape[0];  *ydim = shape[1];  *zdim = shape[2];
  if (VERBOSE)
     printf("%s: data set size xdim=%d ydim=%d zdim=%d\n",
            MY_NAME,*xdim,*ydim,*zdim);

  size = (*xdim) * (*ydim) * (*zdim) * sizeof(float);

  if ((*data = (float *)malloc(size)) == NULL) {
    printf("%s: error, not enough memory for the data set\n",MY_NAME);
    printf("Need space for '%d * %d * %d = %d' entries\n",
                     *xdim,*ydim,*zdim,
                    (*xdim)*(*ydim)*(*zdim));
    fclose(ISO_file);
    return -1;
  }


  ret=SCAREgetdata(filename,rank,shape,*data);
  if (ret != 0) {
     printf("%s: error from SCAREgetdata %d file %s\n",
             MY_NAME,ret,filename);
     return -1;
  }

  return 0;
}

#define RECORD()   { Items = \
                   fread(&record, sizeof(int) , 1L , ISO_file);\
                   if(Items < 1) {\
                     printf("?ERROR - in reading contour file (*)\n");\
                     fclose(ISO_file);\
                     return(1);}}

#define FREAD(value_p , size)    { Items = \
                                 fread(value_p, size , 1L , ISO_file);\
                                 if(Items < 1) {\
                     printf("?ERROR - in reading contour file (**)\n");\
                     fclose(ISO_file);\
                     return(2);}}

#define FREADN(value_p , num , size) { Items = \
                                fread(value_p, size , num , ISO_file);\
                   if(Items != num) {\
                     printf("?ERROR - in reading contour file (***)");\
                     fclose(ISO_file);\
                     return(3);}}


/**************************** SCAREgetdims ****************************/
/**************************** SCAREgetdims ****************************/
/**************************** SCAREgetdims ****************************/
/**************************** SCAREgetdims ****************************/
int SCAREgetdims(filename,rank,shape,unknown)
      char *filename;
      int  *rank;
      int   shape[];
      int   unknown;
{
      char    OutText[BUFF_LEN];
      int     Items;
      int     record;
      
      ISO_file = fopen(filename,"rb");
      if(ISO_file == NULL) {
         printf("?ERROR - can't open file contour file '%s'\n",filename);
         return(1);}

/* first position is the rank (if this is a binary file) */
       FREAD(rank, sizeof(int));

/* pure binary file */
       if(*rank == 3) { /* this is a binary file */

       gomContourInternalType = 0;

/* type of surface            */
       FREAD(&TypeOfSurface, sizeof(int));

       if(VERBOSE) {
         switch(TypeOfSurface) {
         case 0: 
         printf("?WARNING - unknown surface type ");
         break;
          case 1:
          printf("=> Reading a VSS surface ");
          break;
           case 2:
           printf("=> Reading an orbital/density surface");
           break;
            case 3:
            printf("=> Reading a probe surface");
            break;}}

/* read shape ...             */
        FREAD(&shape[2] , sizeof(int));
        FREAD(&shape[1] , sizeof(int));
        FREAD(&shape[0] , sizeof(int));
       
       if(VERBOSE)
         printf("zdim: %d , ydim: %d , xdim: %d\n",shape[2],shape[1],shape[0]);
/* done with shape */

/* read min/max ... */
       FREAD(&ISO_zmin , sizeof(float));
       FREAD(&ISO_zmax , sizeof(float));
        FREAD(&ISO_ymin , sizeof(float));
        FREAD(&ISO_ymax , sizeof(float));
         FREAD(&ISO_xmin , sizeof(float));
         FREAD(&ISO_xmax , sizeof(float));

       }
/* this is a fortran unformatted file */
       else { /* this is a file produced with fortran */

        gomContourInternalType = 1;

        FREAD(rank, sizeof(int));
        FREAD(&TypeOfSurface, sizeof(int));

        RECORD();   /* control record */

        RECORD();   /* control record */
        FREAD(&shape[2] , sizeof(int));
        FREAD(&shape[1] , sizeof(int));
        FREAD(&shape[0] , sizeof(int));
        RECORD();   /* control record */

/* read min/max ... */
       RECORD();   /* control record */
       FREAD(&ISO_zmin , sizeof(float));
       FREAD(&ISO_zmax , sizeof(float));
        FREAD(&ISO_ymin , sizeof(float));
        FREAD(&ISO_ymax , sizeof(float));
         FREAD(&ISO_xmin , sizeof(float));
         FREAD(&ISO_xmax , sizeof(float));
       RECORD();   /* control record */
       }

       /* make some conversion (if necessary) */
       if(TypeOfSurface == 100) {  /* OpenMol data in atomic units */
         ISO_zmin *= Atm2Ang;
         ISO_zmax *= Atm2Ang;
         ISO_ymin *= Atm2Ang;
         ISO_ymax *= Atm2Ang;
         ISO_xmin *= Atm2Ang;
         ISO_xmax *= Atm2Ang;
       }

       if(VERBOSE)
         printf("x: %f %f \ny: %f %f \nz: %f %f\n",
         ISO_xmin,
         ISO_xmax,
         ISO_ymin,
         ISO_ymax,
         ISO_zmin,
         ISO_zmax);
          
/* done with min/max */

    return(0);
  }

/**************************** SCAREgetdata ****************************/
/**************************** SCAREgetdata ****************************/
/**************************** SCAREgetdata ****************************/
/**************************** SCAREgetdata ****************************/
int SCAREgetdata(char *filename,int rank,int *shape,float *data)
{   
      int Items;
      int il;
      int size;
      int record;
      float *Temp;

/* unformatted file */
      if(gomContourInternalType) {
       size = shape[1] * shape[0];
       Temp = data;
       for(il = 0 ; il < shape[2] ; il++) {
         RECORD();   /* control record */
         FREADN(Temp , size , sizeof(float));
         RECORD();   /* control record */
         Temp += size;
       }
      }
/* pure binary file */
      else {
        size = shape[2] * shape[1] * shape[0];

        FREADN(data , size , sizeof(float));
      }

      fclose(ISO_file);

      return(0);
}
void get_max_min(data, xdim, ydim, zdim, maxp, minp)
register float *data;
int xdim, ydim, zdim;
float *maxp, *minp;
/* This subroutine finds the maximum & minimum data values */
{
    float *enddata;
    float max, min;

    enddata = data + (xdim * ydim * zdim);
    max = min = *(data++);
    for (; data < enddata; data++) {
        if (*data > max) max = *data;
        if (*data < min) min = *data;
      }
    *maxp = max; *minp = min;
  }


/**************************************************************************/
int lulWriteContourFORM(FILE *Output_p)
/**************************************************************************/
{
      static int i,j,k,Loop,Loop1,Items;
      static float xc,yc,zc,pot;
      static float Help1,Help2;
      static float Buffer[BUFFER_LIM];
      static float MaxV = -1.0e10;
      static float MinV =  1.0e10;

/* print out number of points in x- , y- and z-directions */

      fprintf(Output_p,"%d %d\n",rank,TypeOfSurface);

      fprintf(Output_p,"%d %d %d\n",ContourInfo[0].zdim ,
                                    ContourInfo[0].ydim , 
                                    ContourInfo[0].xdim );

      fprintf(Output_p,"%e %e %e %e %e %e\n",zmin,zmax,ymin,ymax,xmin,xmax);


      MaxV = -1.0e10;
      MinV =  1.0e10;

      Loop = 0;

      for(i = 0 ; 
          i < ContourInfo[0].xdim * ContourInfo[0].ydim * ContourInfo[0].zdim;
          i++) {

            Help1 = ContourInfo[0].data[i];

            if(Help1 < MinV) MinV = Help1;
            if(Help1 > MaxV) MaxV = Help1;

            fprintf(Output_p,"%e ",Help1);

            Loop++;
            if(Loop == 6) {
              fprintf(Output_p,"\n");
              Loop = 0;}

      }      

      fprintf(Output_p,"\n");

      printf("Max value: %f \n",MaxV);
      printf("Min value: %f \n",MinV);
      printf("! Done.\n");

      fclose(Output_p);

      return(0);
}


/**************************************************************************/
int lulReadContourFORM(FILE *Input_p)
/**************************************************************************/
{
      static int i,j,k,Items,Entries;
      static float xc,yc,zc,pot;
      static float Help1,Help2;
      static float Buffer[BUFFER_LIM];

/* print out number of points in x- , y- and z-directions */
      fscanf(Input_p,"%d %d",&rank,&TypeOfSurface);

      fscanf(Input_p,"%d %d %d",&ContourInfo[0].zdim ,
                                &ContourInfo[0].ydim , 
                                &ContourInfo[0].xdim );

      fscanf(Input_p,"%e %e %e %e %e %e",&zmin,&zmax,&ymin,&ymax,&xmin,&xmax);

       ISO_zmin = zmin;
       ISO_zmax = zmax;
       ISO_ymin = ymin;
       ISO_ymax = ymax;
       ISO_xmin = xmin;
       ISO_xmax = xmax;

      Entries = ContourInfo[0].xdim * ContourInfo[0].ydim * ContourInfo[0].zdim;

      ContourInfo[0].data = (float *)malloc(Entries * sizeof(float));
        if(ContourInfo[0].data == NULL) {
           printf("$ERROR - can't allocate memory\n");
           exit(1);
	}

      for(i = 0 ; i < Entries; i++) {

            fscanf(Input_p,"%e",&Help1);

            ContourInfo[0].data[i] = Help1;

      }      

      fclose(Input_p);

      return(0);
}

