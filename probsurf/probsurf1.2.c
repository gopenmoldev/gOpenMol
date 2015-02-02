/*
  A probe surface program based on the ideas by R. Voorintholt,
  M.T. Koster, G. Vegter, G. Vriend and W.G.J. Hol in J. Mol. Graphics
  7 (1989) 243.

  Copyright  Leif Laaksonen 1991,1992

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <malloc.h>

#define LINE_LEN 80
#define BUFF_LEN 500
#define BUFFER_LIM 500
#define CONTOUR_FILE "probesurf.plt"
#define SEP_STRING ","
#define LOOP_TRIGGER 1.0e-4

#define FWRITE(value_p , size)    { Items = \
                                 fwrite(value_p, size , 1L , PROBE_p);\
                                 if(Items < 1) {\
                     printf("?ERROR - in writing contour file (*)\n");\
                     return(1);}}

#define FWRITEN(value_p , num , size) { Items = \
                                fwrite(value_p, size , num , PROBE_p);\
                   if(Items < 1) {\
                     printf("?ERROR - in writing contour file (**)\n");\
                     return(1);}}

#define CVALUE   100.
#define RABS(a)    ( ( a ) > 0.0 ? (a) : -(a))

      float *ar3pot; /* potential array */
      float *x,*y,*z; /* x,y and z-coordinate pointers */
      float *VdwRad;
      int    nat;      /* total number of atoms */
      int npts[3];     /* number of points in x,y and z */
      int BuffMult;    /* number of inner loops in buffer */
      float xmin[3];  /* min value */
      float xmax[3];  /* max value */
      float rinc[3];  /* step length */
      float autang = 0.52917715;   /*             */
      float ProbeRad; /* probe vdW radius         */
      float ProbeRad2; /* probe vdW radius**2     */
      char PROBE_file[BUFF_LEN] = CONTOUR_FILE;
      char PROBE_limits[BUFF_LEN];

           char *announce =
"      **********************************************************\n\
      *                                                        *\n\
      *    Molecular Probe Surface program                     *\n\
      *    PROBSURF                                            *\n\
      *                                                        *\n\
      *    Leif Laaksonen    Version 1.2                       *\n\
      *                      1998                              *\n\
      *                                                        *\n\
      **********************************************************\n\n";

           char *usage =
"Usage:\n\
 probesurf -ofile.name -llimits < input.file > output.file\n\
                                                          \n\
           -h this text\n\
           -o defines the name of the contour file (Default 'probesurf.plt')\n\
           -l defines the limits as:\n\
              numx,numy,numz,xmin,xmax,ymin,ymax,zmin,zmax\n\
              where numx,numy,numz are the number of points in\n\
              the x,y and z directions and the rest are the coordinate\n\
              limits of the surface\n\n";


float DoubleCheck(float , float , float);
int   uotm(float , float , float , float *);                                             
int   slpmap(void);
int   ParsePROBElimits(int * , float * , float *);

main(argc,argv)
    int argc;
    char **argv;
{
      float angtau;
      int iu,nh,i;
      char ctitle[LINE_LEN];
      char SegmentName[BUFF_LEN];
      char ResidueName[BUFF_LEN];
      char AtomName[BUFF_LEN];

      printf("%s",announce);

/* parse input line */
       if(argc > 1) {
               for(i = 1 ; i < argc ; i++) {
                  if(argv[i][0] == '-') {
                     switch(argv[i][1]) {
		         case 'h':    /* help text        */
                                  printf("%s",usage);
                                  exit(0);
                         case 'o':    /* output file name */
                                  strncpy(PROBE_file,&argv[i][2],BUFF_LEN);
                                  break;
			 case 'l':    /* points and coordinates limits
                         numx,numy,numz,xmin,xmax,ymin,ymax,zmin,zmax */
                                  strncpy(PROBE_limits,&argv[i][2],BUFF_LEN);
                                  break;
                         default:
                                  printf("?ERROR - unknown parameter '-%s'\n",
                                  argv[i][1]);
                                  break;
                      } /* switch */
                   } /* if */
                } /* for */
        } /* if */ 

/* read title. */
      gets(ctitle);
      printf("Title: %s \n",ctitle);                                  
      gets(ctitle);
      printf("Title: %s \n",ctitle);                                  
/* read total number of atoms, number of hydrogens and units flag. */
      gets(ctitle);
      sscanf(ctitle,"%d",&nat);

      printf("Number of atoms    : %d \n",nat);

/* read atom coordinates and number of valence electrons */
      x      = (float *) malloc(nat * sizeof(float));
      y      = (float *) malloc(nat * sizeof(float));
      z      = (float *) malloc(nat * sizeof(float));
      VdwRad = (float *) malloc(nat * sizeof(float));

      if( x == NULL ||
          y == NULL ||
          z == NULL ||
          VdwRad == NULL) {
            printf("?ERROR - can't allocate memory\n");
            exit(1);}

      for(i = 0; i < nat ; i++ ) {
      gets(ctitle);
      sscanf(ctitle,"%*d %s %s %s %f %f %f %f",
              SegmentName,ResidueName,AtomName,
             &x[i],&y[i],&z[i],&VdwRad[i]);
      
      printf("'%d' (%s:%s:%s) x: %f y: %f z: %f vdWrad: %f \n",i+1,
                                    SegmentName,ResidueName,AtomName,
                                    x[i],y[i],z[i],VdwRad[i]);  
      }

/* read coordinate bounds and number of intervals along the axes. */
      gets(ctitle);
      sscanf(ctitle,"%f %f %f %f %f %f",
      &xmin[0],&xmax[0],&xmin[1],&xmax[1],&xmin[2],&xmax[2]);
      gets(ctitle);
      sscanf(ctitle,"%d %d %d",
      &npts[0],&npts[1],&npts[2]);
      gets(ctitle);
      sscanf(ctitle,"%f",&ProbeRad);

      if(PROBE_limits[0] != '\0') ParsePROBElimits(npts,xmin,xmax);

      printf("Xmin: %f , Xmax: %f\n",xmin[0],xmax[0]);
      printf("Ymin: %f , Ymax: %f\n",xmin[1],xmax[1]);
      printf("Zmin: %f , Zmax: %f\n",xmin[2],xmax[2]);
      printf("Xpts: %d , Ypts: %d , Zpts: %d\n",npts[0],npts[1],npts[2]);
      printf("Probe radius: %f\n",ProbeRad);

/*    check inner loop (x) */
      BuffMult = BUFFER_LIM / npts[0];
      if(BuffMult < 1) {
       printf("?ERROR - internal buffer (BUFFER_LIM) has to be deeper than number of points in x direction\n");
       exit(1);}

      ProbeRad2 = ProbeRad * ProbeRad;

/* calculate step-lengths  */
        for(i = 0 ; i < 3 ; i++) {        
          if(npts[i] > 1) 
            rinc[i] = (xmax[i]-xmin[i])/(npts[i]-1);
          else
            rinc[i]=0.0;
	}
/* generate the binary map file and listing file */
      slpmap(); /* go and get it now ....*/ 

      exit(0);
      }


int slpmap()
{
      static int i,j,k,Loop,Loop1,Items;
      static float xc,yc,zc,pot;
      static float Help1,Help2;
      FILE *PROBE_p;
      static float Buffer[BUFFER_LIM];
      static float MaxV = -1.0e10;
      static float MinV =  1.0e10;
      static int TypeOfSurface;

/* print out number of points in x- , y- and z-directions */

      PROBE_p = fopen(PROBE_file,"wb");
      if(PROBE_p == NULL) {
        printf("?ERROR - can't open contour file '%s' \n",PROBE_file);
        return(1);}

      printf("\n==> Writing contour to file '%s' \n",PROBE_file);


      i = 3;
      FWRITE(&i , sizeof(int));
      TypeOfSurface = 3;
      FWRITE(&TypeOfSurface , sizeof(int));
      FWRITE(&npts[2] , sizeof(int));
       FWRITE(&npts[1] , sizeof(int));
        FWRITE(&npts[0] , sizeof(int));

      Help1 = xmin[2];
       Help2 = xmax[2];
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));
      Help1 = xmin[1];
       Help2 = xmax[1];
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));
      Help1 = xmin[0];
       Help2 = xmax[0];
        FWRITE(&Help1 , sizeof(float));
        FWRITE(&Help2 , sizeof(float));

/*

      fprintf(PROBE_p," 3         ! Rank value\n");
      fprintf(PROBE_p," %d        ! number of points in the z-direction\n",npts[2]);
      fprintf(PROBE_p," %d        ! number of points in the y-direction\n",npts[1]);
      fprintf(PROBE_p," %d        ! number of points in the x-direction\n",npts[0]);
      fprintf(PROBE_p," %f  %f    ! zmin and zmax\n",xmin[2]*autang,xmax[2]*autang);
      fprintf(PROBE_p," %f  %f    ! ymin and ymax\n",xmin[1]*autang,xmax[1]*autang);
      fprintf(PROBE_p," %f  %f    ! xmin and xmax\n",xmin[0]*autang,xmax[0]*autang);

*/
      MaxV = -1.0e10;
      MinV =  1.0e10;

      Loop = 0;
       Loop1 = 0;

      for(i = 0 ; i < npts[2] ; i++) { /* loop over z */
            zc = xmin[2] + rinc[2] * i;
         for(j = 0 ; j < npts[1] ; j++) { /* loop over y */
            yc = xmin[1] + rinc[1] * j;
            for(k = 0; k < npts[0] ; k++) { /* loop over x */
            xc = xmin[0] + rinc[0] * k;

            uotm(xc,yc,zc,&pot);

            Help1 = pot;

            if(Help1 < MinV) MinV = Help1;
            if(Help1 > MaxV) MaxV = Help1;

              Buffer[Loop] = Help1;
              Loop++;
          
	   }      /* x */

            Loop1++;

            if(Loop1 == BuffMult){
              FWRITEN(Buffer , Loop , sizeof(float));
                Loop = 0;
                 Loop1 = 0;}

	  }       /* y */
      }       /* z */

      if(Loop)
            FWRITEN(Buffer , Loop , sizeof(float));

      printf("Max value: %f \n",MaxV);
      printf("Min value: %f \n",MinV);
      printf("! Done.\n");

      fclose(PROBE_p);

      return(0);
}


int uotm(float xh, float yh, float zh, float *pot)                                             
{
      static float dist,dist2;
      static float xi,yi,zi;
      static float xf,yf,zf;
      static float rado,rado2,radi,radi2;
      static float value,temp;
      register int    i;
      register float *X_v;
      register float *Y_v;
      register float *Z_v;

      value = 0.0;
      X_v   = x;
       Y_v  = y;
        Z_v = z;

      for(i = 0 ; i < nat ; i++) {

/* if we reached a value close enough to CVALE stop looping */
       if(fabs(CVALUE - value) < LOOP_TRIGGER) {
          *pot = value;
           return(0);
       }

       xi= *(X_v + i);
       xf = xi - xh;

       dist2 = xf * xf;

        yi= *(Y_v + i); 
        yf = yi - yh;

        dist2 = yf * yf + dist2;

         zi= *(Z_v + i);
         zf = zi - zh;

         dist2 = zf * zf + dist2;

         rado  = ProbeRad + VdwRad[i];
         rado2 = rado * rado;

        if(dist2 < rado2) {
          
          radi  = VdwRad[i];
          dist  = sqrt(dist2);

          if(dist >= radi) {

            radi2 = radi * radi;

            temp = CVALUE * (rado2 - dist2)/(rado2 - radi2);
          } else {
            temp = CVALUE;
          }

          if(temp > value) value = temp;
        }

     }

     *pot = value;

     return(0);
}


int ParsePROBElimits(int *GetNpts , float *GetMin , float *GetMax)
{
     char *find;
     int loop;

     printf(".... > Redefining input limits < ....\n");

     find = strtok(PROBE_limits,SEP_STRING);

     if(find ==NULL) return(1);

     GetNpts[0] = atoi(find);

     for(loop = 1; loop < 9 ; loop++) {

     find = strtok(NULL,SEP_STRING);

       if(find == NULL) return(0);

       switch(loop) {
       case 1: /* ypts */
              GetNpts[1] = atoi(find);
              break;
       case 2: /* zpts */
              GetNpts[2] = atoi(find);
              break;
       case 3: /* xmin */
              GetMin[0] = atof(find);
              break;
       case 4: /* xmax */
              GetMin[1] = atof(find);
              break;
       case 5: /* ymin */
              GetMin[2] = atof(find);
              break;
       case 6: /* ymax */
              GetMin[3] = atof(find);
              break;
       case 7: /* zmin */
              GetMin[4] = atof(find);
              break;
       case 8: /* zmax */
              GetMin[5] = atof(find);
              break;}
   }

   return(0);
}
float DoubleCheck(float xh , float yh , float zh)
{
      static int CloseIndex;
      static float Dmin,dist,Dmin2,dist2;
      static float xi,yi,zi;
      static float xf,yf,zf;
      static float value;
      static float Dtrigger;
      register int    i;
      register float *X_v;
      register float *Y_v;
      register float *Z_v;

      X_v = x;
       Y_v = y;
        Z_v = z;

      value    = 0.0;
      Dtrigger = 1.0e+20;

      for(i = 0 ; i < nat ; i++) {

       Dmin  = ProbeRad + VdwRad[i];
       Dmin2 = Dmin * Dmin;

       xi= *(X_v + i);
       xf = xi - xh;

       dist2 = xf * xf;
       if(dist2 > Dmin2) continue;

        yi= *(Y_v + i); 
        yf = yi - yh;

        dist2 = yf * yf + dist2;
        if(dist2 > Dmin2) continue;

         zi= *(Z_v + i);
         zf = zi - zh;

         dist2 = zf * zf + dist2;

        if((dist2 < Dmin2) && (dist2 < Dtrigger) ) {

          dist = sqrt(dist2);

          if(dist > VdwRad[i] && dist <= Dmin) {
            value = CVALUE * (Dmin2 - dist2)/(Dmin2 - VdwRad[i] * VdwRad[i]);
          } else {
            value = CVALUE;
          }

          Dtrigger = dist2;
        }

     }

     return(value);
}






