
/*
     This is a program for generating orbital and density plots
     of molecules. As the program is now it only accepts the
     wave function from an Extended Huckel calculation 

     Leif Laaksonen 1990

*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <malloc.h>

#define VERBOSE    0
#define BUFF_LEN 500
#define BUFFER_LEN 256
#define SEP_STRING ","

/*
#define RECORD()   { Items = \
                   fread(&record, sizeof(int) , 1L , ISO_file);\
                   if(Items < 1) {\
                     PrintMessage("?ERROR - in reading wave function file");\
                     return(1);}}

#define FREAD(value_p , size)    { Items = \
                                 fread(value_p, size , 1L , ISO_file);\
                                 if(Items < 1) {\
                     PrintMessage("?ERROR - in reading wave function file");\
                     return(1);}}

#define FREADN(value_p , num , size) { Items = \
                                fread(value_p, size , num , point_f);\
                   if(Items < 1) {\
                     PrintMessage("?ERROR - in reading trajectory file");\
                     return(1);}}
*/

#define FWRITE(value_p , size)    { Items = \
                                 fwrite(value_p, size , 1L , point_f);\
                                 if(Items < 1) {\
                     printf("?ERROR - in writing contour file\n");\
                     return(1);}}

#define FWRITEN(value_p , num , size) { Items = \
                                fwrite(value_p, size , num , point_f);\
                   if(Items < 1) {\
                     printf("?ERROR - in writing contour file\n");\
                     return(1);}}


#define RABS(a)      ((a) > 0.0 ? (a) : -(a) )
#define MAXEXP          20
#define PI              3.14159265
#define CONTOUR_FILE "orbital.plt"

#define SUM_Tools(i)  {        /* 2s-orbital */\
                               sum += dist * eigv[orb][j] * \
                               rnorms[i] * exp(-exps[i] * dist);\
                               j++;\
                               /* 2p-orbital */\
                               phelp = exp(-expp[i] * dist);\
                               sum += (deltz * eigv[orb][j+2] +\
                                       delty * eigv[orb][j+1] +\
                                       deltx * eigv[orb][j])  *\
                                       rnormp[i] * phelp;\
                                       j += 3;}

#define SUM_Tools1(i) {        /* 3s-orbital */\
                               sum += dist*dist * eigv[orb][j] * \
                               rnorms[i] * exp(-exps[i] * dist);\
                               j++;\
                               /* 3p-orbital */\
                               phelp = exp(-expp[i] * dist);\
                               sum += (deltz * dist * eigv[orb][j+2] +\
                                       delty * dist * eigv[orb][j+1] +\
                                       deltx * dist * eigv[orb][j])  *\
                                       rnormp[i] * phelp;\
                                       j += 3;\
                               /* 3d-orbital */\
                               phelp = exp(-expp[i] * dist);\
                               sum += ((deltx * deltx - delty * delty) * \
                                                        eigv[orb][j] +\
                                       (3.0 * deltz * deltz - dist *dist) *\
                                                        eigv[orb][j+1] +\
                                       deltx * delty *  eigv[orb][j+2]   +\
                                       deltx * deltz *  eigv[orb][j+3]   +\
                                       delty * deltz *  eigv[orb][j+4])  *\
                                        rnormd[i] * phelp;\
                                       j += 5;}

#define SUM_Tools2(i)  {       /* 4s-orbital */\
                               dist2 = dist * dist;\
                               sum += dist * dist2 * eigv[orb][j] * \
                               rnorms[i] * exp(-exps[i] * dist);\
                               j++;\
                               /* 4p-orbital */\
                               phelp = exp(-expp[i] * dist);\
                               sum += (deltz * eigv[orb][j+2] +\
                                       delty * eigv[orb][j+1] +\
                                       deltx * eigv[orb][j])  *\
                                       rnormp[i] * dist2 * phelp;\
                                       j += 3;}

           char *announce =
"      **********************************************************\n\
      *                                                        *\n\
      *    Molecular orbital and density program               *\n\
      *    PSI/DENSITY                                         *\n\
      *                                                        *\n\
      *    Leif Laaksonen    Version 1.0                       *\n\
      *                      1991                              *\n\
      *                                                        *\n\
      **********************************************************\n\n";

           char *symbol = 
" C N O FSi P SClLiBe BNaMgAl HCACBCGCDCECHCZNANBNGNHNDNENZOAOBOCOGODOEOZOH\
SASBSGSDHHHGHZHEHAHBHDHNCa";

           int ihelpv[60] = {
                   1, 2, 3, 4, 5, 6, 7,
                   8, 9,10,11,12,13,14,15, 
                   1, 1, 1, 1, 1, 1, 1,
                   2, 2, 2, 2, 2, 2, 2,
                   3, 3, 3, 3, 3, 3, 3, 3,
                   7, 7, 7, 7,
                  15,15,15,15,15,15,15,15,16};
     
           int velec[MAXEXP] = {4,5,6,7,4,5,6,7,1,2,3,5,2,3,1,2};

           float atmu = 0.529167;  

#ifdef NEW

/* s-orbitals ... */
           int ns[MAXEXP] = { 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 2 , 2 , 2 ,
                          2 , 3 , 3 , 1 , 4};

           float exps[MAXEXP] = {1.608,1.924,2.246,2.564,1.634,
           1.881,2.122,2.356,0.65,0.975,1.3,0.733,0.95,1.167,1.3,1.2};

/* p-orbitals ... */
           int np[MAXEXP] = { 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 
                          2 , 3 , 3 , 0 , 4};

           float expp[MAXEXP] = {1.568,1.917,2.227,2.55,1.428,
           1.628,1.827,2.039,0.65,0.975,1.3,0.733,0.95,1.167,0.0,1.2};

/* d-orbitals ... */
           int nd[MAXEXP] = { 0 , 0 , 0 , 0 , 3 , 3 , 3 , 3};

           float expd[MAXEXP] = {0.0,0.0,0.0,0.0,1.383,1.4,1.5,2.033};

#else

/* s-orbitals ... */
           int ns[MAXEXP] = { 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 2 , 2 , 2 ,
                          2 , 3 , 3 , 1 , 4};

           float exps[MAXEXP] = {1.625,1.950,2.275,2.425,1.383,
           1.6,1.817,2.033,0.65,0.975,1.3,0.733,0.95,1.167,1.3,1.2};

/* p-orbitals ... */
           int np[MAXEXP] = { 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 
                          2 , 3 , 3 , 0 , 4};

           float expp[MAXEXP] = {1.625,1.950,2.275,2.425,1.383,
           1.6,1.817,2.033,0.65,0.975,1.3,0.733,0.95,1.167,1.2};

/* d-orbitals ... */
           int nd[MAXEXP] = { 0 , 0 , 0 , 0 , 3 , 3 , 3 , 3};

           float expd[MAXEXP] = {0.0,0.0,0.0,0.0,1.383,1.4,1.5,2.033};

#endif

/* detailed section ...                   */

      int numat;           /* number of atoms */
      int numorb;          /* number of orbitals */
      int numelec;         /* number of electrons */
      int CalcType = 0;    /* Calclulation type
                              = 0 , orbital
                              = 1 , orbital density
                              = 2 , total electron density */

 char eigv_file[BUFFER_LEN];     /* name of eigen vector matrix file */
 char coord_file[BUFFER_LEN];    /* name of coordinate file (EHT input file) */
 char DENSITY_limits[BUFF_LEN];

      float rnorms[MAXEXP];  /* normalization constatant */
      float rnormp[MAXEXP];
      float rnormd[MAXEXP];

      float dx,dy,dz; /* step length */
      int plot_orb;   /* orbital to be plotted */

      struct {
      float minx;
      float maxx;
      float miny;
      float maxy;
      float minz;
      float maxz;} limits;

      struct {
      int x;
      int y;
      int z;} steps; 

      float *x;  /* pointers to x,y and z coordinates */
      float *y;
      float *z;

      float *orbocc; /* pointer to orbital occupation number */

      float **eigv; /* pointer to eigenvectors */

      char *atnam; /* pointer to atom names (only two characters long each) */

      char *atmtype; /* pointer to "atom type" */

      FILE *point_f;  /* file pointer to the point file */
      char point_file[BUFFER_LEN] = CONTOUR_FILE;

      struct {
      int ReadType;        /* = 0 , EHT file name supplied
                              = 1 , coordinates supplied
                           */
              } ControlBlock;

/****************************************************************************/
main(int argc , char *argv[])
/****************************************************************************/
{
       int i;
/* parse input line */
       if(argc > 1) {
               for(i = 1 ; i < argc ; i++) {
                  if(argv[i][0] == '-') {
                     switch(argv[i][1]) {
                         case 'o':    /* output file name */
                                  strncpy(point_file,&argv[i][2],BUFFER_LEN);
                                  break;
			 case 'l':    /* points and coordinates limits
                         numx,numy,numz,xmin,xmax,ymin,ymax,zmin,zmax */
                                  strncpy(DENSITY_limits,&argv[i][2],BUFF_LEN);
                                  break;
                         default:
                                  printf("?ERROR - unknown parameter '-%s'\n",
                                  argv[i][1]);
                                  break;
                      } /* switch */
                   } /* if */
                } /* for */
        } /* if */ 

       printf("%s",announce);

       read_input();

       if(!ControlBlock.ReadType)
                        read_eht_input(); /* get atom coordinates */

       count_num_elec(); /* count numbers of electrons and assign atom types */
   
       read_eig_vec();   /* read eigenvectors */

       convert(); /* convert units to A */

       switch(CalcType) {

       case 0:
       psi(plot_orb);
       break;

       case 1:
       psi2(plot_orb);
       break;

       case 2:
       psit();
       break;}
}
/* convert from Angstrom units to atomic units */
/****************************************************************************/
convert()
/****************************************************************************/
{

      static int i;
      static float srt3,srt4,srt5;

/* convert coordinates */
      for(i = 0 ; i < numat ; i++) {
        x[i] /= atmu;
         y[i] /= atmu;
          z[i] /= atmu;}

/* convert limits */
      limits.maxx /= atmu;
      limits.minx /= atmu;
      limits.maxy /= atmu;
      limits.miny /= atmu;
      limits.maxz /= atmu;
      limits.minz /= atmu;

      srt3 = sqrt(3.);
      srt4 = sqrt(2. / 15.);
      srt5 = srt4/srt3;

/* calculate normalization constants */
/* s-orbitals */
/* p-orbitals */
/* d-orbitals */
/* C */
       rnormd[0] = 0.0;
       rnormp[0] = sqrt(pow(exps[0],5.) / PI);
       rnorms[0] = rnormp[0] / srt3;
/* N */
       rnormd[1] = 0.0;
       rnormp[1] = sqrt(pow(exps[1],5.) / PI);
       rnorms[1] = rnormp[1] / srt3;
/* O */
       rnormd[2] = 0.0;
       rnormp[2] = sqrt(pow(exps[2],5.) / PI);
       rnorms[2] = rnormp[2] / srt3;
/* F */
       rnormd[3] = 0.0;
       rnormp[3] = sqrt(pow(exps[3],5.) / PI);
       rnorms[3] = rnormp[3] / srt3;
/* Si */
       rnormd[4] = sqrt(pow(exps[4],7.) / PI);
       rnormp[4] = rnormd[4] * srt4;
       rnorms[4] = rnormd[4] * srt5;
/* P */
       rnormd[5] = sqrt(pow(exps[5],7.) / PI);
       rnormp[5] = rnormd[5] * srt4;
       rnorms[5] = rnormd[5] * srt5;
/* S */
       rnormd[6] = sqrt(pow(exps[6],7.) / PI);
       rnormp[6] = rnormd[6] * srt4;
       rnorms[6] = rnormd[6] * srt5;
/* Cl */
       rnormd[7] = sqrt(pow(exps[7],7.) / PI);
       rnormp[7] = rnormd[7] * srt4;
       rnorms[7] = rnormd[7] * srt5;
/* Li */
       rnormd[8] = 0.0;
       rnormp[8] = sqrt(pow(exps[8],5.) / PI);
       rnorms[8] = rnormp[8] / srt3;
/* Be */
       rnormd[9] = 0.0;
       rnormp[9] = sqrt(pow(exps[9],5.) / PI);
       rnorms[9] = rnormp[9] / srt3;
/* B */
       rnormd[10] = 0.0;
       rnormp[10] = sqrt(pow(exps[10],5.) / PI);
       rnorms[10] = rnormp[10] / srt3;
/* Na*/ 
       rnormd[11] = sqrt(pow(exps[11],7.) / PI);
       rnormp[11] = rnormd[11] * srt4;
       rnorms[11] = rnormd[11] * srt5;
/* Mg */
       rnormd[12] = sqrt(pow(exps[12],7.) / PI);
       rnormp[12] = rnormd[12] * srt4;
       rnorms[12] = rnormd[12] * srt5;
/* Al */
       rnormd[13] = sqrt(pow(exps[13],7.) / PI);
       rnormp[13] = rnormd[13] * srt4;
       rnorms[13] = rnormd[13] * srt5;
/* H */
       rnorms[14] = sqrt(pow(exps[14],3.) / PI);
/* Ca */
       rnormd[15] = 0.0;
       rnormp[15] = sqrt(pow(expp[15],9.) / (105. * PI));
       rnorms[15] = sqrt(pow(exps[15],9.) / (315. * PI));

}
/* read eht input file for coordinates and atom names */
/****************************************************************************/
read_eht_input()
/****************************************************************************/
{
    FILE *eht_p;

    int i,j,k;
    char inputs[BUFFER_LEN];
    char NumAtm[3];
    char TempC;

    eht_p = fopen(coord_file,"r");

/* file error ... */
    if(eht_p == NULL) {
    printf("?ERROR - can't open EHT input file '%s' \n",coord_file);
    exit(1);}

/* first title card */
  printf("*** EHT input file title: \n");
  loop:;
    fgets(inputs,120,eht_p);
    printf("%s\n",inputs);
    if(strlen(inputs) > 71) 
       if(inputs[71] == '*') goto loop;

    fgets(inputs,120,eht_p);
     strncpy(NumAtm,inputs+3,3);
      sscanf(NumAtm,"%d",&numat);
       printf("Number of atoms: %d \n",numat);

/* go and get space for the coordinate arrays and the atom name array */

    x = (float *) malloc(numat * sizeof(float));
    y = (float *) malloc(numat * sizeof(float));
    z = (float *) malloc(numat * sizeof(float));

    i = (numat / 40 + 1) * 80;
    atnam = (char *) malloc( i );

/* check that everything is ok */

    if(x == NULL || y == NULL || z == NULL || atnam == NULL) {
      printf("?ERROR - memory allocation problems \n");
      exit(1);}

/* start reading coordinates */

    for(i = 0 ; i < numat ; i++) {
       fgets(inputs,120,eht_p);
       sscanf(inputs,"%f %f %f",&x[i],&y[i],&z[i]);

#ifdef DEBUG
printf("x %f y %f z %f \n",x[i],y[i],z[i]);
#endif

     }

/* atom labels */

       i = 2 * numat;
       j = 0;
       while(1){ 
        TempC = getc(eht_p);
         if(i == j) break;
         if(TempC == 13 || TempC == 10) continue;
          atnam[j] = TempC;
          j++;}

#ifdef DEBUG
    for(i = 0 ; i < numat ; i++) printf(" => %d :%.2s:\n",i+1,atnam+2*i);
#endif

/* ok close file and return */

    fclose(eht_p);      

}

/****************************************************************************/
float **matrix(dim)
        int dim;
/****************************************************************************/
{
        static int i;
        float **m;

        m = (float **) malloc(dim * sizeof(float *));
        if(m == NULL) {
        printf("?ERROR - allocation failure in matrix \n");
        exit(1);}

        for(i = 0 ; i < dim ; i++) {
           m[i] = (float *) malloc(dim * sizeof(float));
           if(m[i] == NULL) {
             printf("?ERROR - allocation failure in matrix \n");
             exit(1);}
        }

        return(m);
}
/****************************************************************************/
count_num_elec()
/****************************************************************************/
{
     static int i,j,k;
     static int hit;
     static int collect;

     k = strlen(symbol)/2;

/* grabb space for 'atmtype' array */
     atmtype = (char *) malloc(numat);
     if(atmtype == NULL) {
       printf("?ERROR - allocation failure in 'count_num_elec' \n");
       exit(1);}
/* done .. */

     numorb = 0;
     collect = 0;
     hit = 0;

     for(i = 0; i < numat ; i++) {
        hit = 0;
        for(j = 0 ; j < k ; j++) {
           if(strncmp(atnam+2*i,symbol+2*j,2) == 0) {
             atmtype[i] = ihelpv[j]-1;
             collect += velec[atmtype[i]];
             if(ns[atmtype[i]] > 0) /* s-type contribution */
                numorb += 1;
             if(np[atmtype[i]] > 0) /* p-type contribution */
                numorb += 3;
             if(nd[atmtype[i]] > 0) /* d-type contribution */
                numorb += 5;
             hit = 1;
             break;}
        }
      if(!hit) {
       printf("?ERROR - unresolved atom name '%.2s' at: %d\n",atnam+2*i,(i+1));
        exit(1);}        
      }
        numelec = collect;
        printf("Number of orbitals: %d number of electrons; %d \n",numorb,numelec);

/* orbital occupation array */
     orbocc = (float *) malloc(numorb * sizeof(float));
     if(orbocc == NULL) {
       printf("?ERROR - allocation failure in 'count_num_elec' \n");
       exit(1);}

       k = numelec / 2;
       for(i = 0 ; i < numorb ; i++) orbocc[i] = 0.0;

       for(i = 0; i < k ; i++) orbocc[i] = 2.;

       if((numelec- k*2) > 0) orbocc[k] = 1.0;

}
/****************************************************************************/
read_eig_vec()
/****************************************************************************/
{
     FILE *eigv_p;
     static int i,j,k,l;
     static char inputs[BUFFER_LEN];
     static float sum;

/* grabb first the space .. */
     eigv = matrix(numorb); 

/* then read the eigenvector matrix */
     eigv_p = fopen(eigv_file,"r");
     if(eigv_p == NULL) {
     printf("?ERROR - can't open eigen vector matrix file '%s' \n",eigv_file);
     exit(1);}

/* molecular orbitals in rows */

     fgets(inputs,BUFFER_LEN,eigv_p);
     printf("%s\n",inputs);

     for(i = 0 ; i < numorb ; i++) {
         k = numorb - i - 1;
        for(j = 0 ; j < numorb ; j++) {
        if(!fscanf(eigv_p,"%f",&eigv[k][j])) {
          printf("?ERROR - can't read wave function at (%d:%d)\n",i,j);
          exit(1);}
          }
        }

     fclose(eigv_p); /* close file */

/* check normalisation */
     for(i = 0 ; i < numorb ; i++) {

        sum = 0;
        for(j = 0 ; j < numorb ; j++) 
           sum += eigv[i][j] * eigv[i][j];

        if(RABS(1.-sqrt(sum)) > 1.e-5) { /* not normalized */
        printf("Orbital nr: %d not normalized (%f) , will do it now \n",
                (i+1),sqrt(sum));

        sum = 1. /sqrt(sum);
        for(j = 0 ; j < numorb ; j++) {
           eigv[i][j] *=sum;
#ifdef DEBUG
        printf("> %d %d %f \n",i,j,eigv[i][j]);
#endif
        }
      }
      }
}
/* calculate orbital nr orb */
/****************************************************************************/
psi(orb)
   int orb;
/****************************************************************************/
{
     static int i,j,k,Loop,Items;
     static int ix,iy,iz;
     static float sum,dist,dist2;
     static float xc,yc,zc;
     static float deltx,delty,deltz;
     static float phelp,dhelp;
     static int shape[3];
     static float Buffer[BUFF_LEN];
     static float MaxV = -1.0e10;
     static float MinV =  1.0e10;

     dx = (limits.maxx - limits.minx) / (float)(steps.x - 1);
     dy = (limits.maxy - limits.miny) / (float)(steps.y - 1);
     dz = (limits.maxz - limits.minz) / (float)(steps.z - 1);
  
/* main loops */

/* open output file */
   printf("   Open output file '%s' for the points\n",point_file);

   shape[0] = steps.z;
    shape[1] = steps.y;
     shape[2] = steps.x;

   printf("\n==> Writing contour information to file '%s'\n",point_file);

/* print out number of points in x- , y- and z-directions */

   i = SCAREwritedims(point_file , 3 , shape , 0);
       if(i) {
           printf("?ERROR - in writing to contour file \n");
           exit(1);}

      MaxV = -1.0e10;
      MinV =  1.0e10;

   Loop = 0; 

   for(iz = 0 ; iz < steps.z ; iz++) {
      zc = dz * (float)iz + limits.minz;
      for(iy = 0 ; iy < steps.y ; iy++) {
      yc = dy * (float)iy + limits.miny;
         for(ix = 0 ; ix < steps.x ; ix++) {
         xc = dx * (float)ix + limits.minx;

/* now loop the atoms */
        j = 0;
        sum = 0.0;
        for(i = 0 ; i < numat ; i++) {

        deltx = xc - x[i];
        delty = yc - y[i];
        deltz = zc - z[i];

        dist = sqrt(deltx*deltx+delty*delty+deltz*deltz);

        switch(atmtype[i]) {

           case 0: /* C */

               SUM_Tools(0);

               break;

	     case 14: /* H */
                /* s-orbital */
               sum += eigv[orb][j] * rnorms[14] * exp(-exps[14] * dist);
                    j++;
               break;

            case 1: /* N */

               SUM_Tools(1);

               break;

	     case 2: /* O */

               SUM_Tools(2);

               break;

	     case 3: /* F */

               SUM_Tools(3);

               break;

             case 4: /* Si */

               SUM_Tools1(4);

               break;

             case 5: /* P */

               SUM_Tools1(5);

               break;

             case 6: /* S */

               SUM_Tools1(6);

               break;

             case 7: /* Cl */

               SUM_Tools1(7);

               break;

             case 8: /* Li */
 
                /* s-orbital */
               sum += dist * eigv[orb][j] * rnorms[8] * exp(-exps[8] * dist);
                    j++;

               break;

             case 10: /* B */

               SUM_Tools(10);

               break;

             case 11: /* Na */

               sum += dist*dist * eigv[orb][j] * 
               rnorms[11] * exp(-exps[11] * dist);
               j++;

               break;

             case 13: /* Al */

               SUM_Tools1(13);

               break;

             case 15: /* Ca */

               SUM_Tools2(15);

               break;

             default:
             printf("Atom type %d not yet implemented \n",atmtype[i]);
	     }   /* end of switch */
             }

            if(sum < MinV) MinV = sum;
            if(sum > MaxV) MaxV = sum;

            if(Loop == BUFF_LEN){
              FWRITEN(Buffer , BUFF_LEN , sizeof(float));
                Loop = 0;}

              Buffer[Loop] = sum;
              Loop++;
/*
             fprintf(point_f," %f \n",sum);
*/
           } /* x */
          }  /* y */
         }   /* z */


      if(Loop)
            FWRITEN(Buffer , Loop , sizeof(float));

  fclose(point_f);

  printf("Max value: %f \n",MaxV);
  printf("Min value: %f \n",MinV);
  printf("! Done.\n");
} 

/****************************************************************************/
read_input()
/****************************************************************************/
{

       int i;
       char inputs[BUFFER_LEN];

       ControlBlock.ReadType = 0;

/* first card is title */
       gets(inputs);
       printf("   Title: %s\n",inputs);
/* coordinate file name */
       gets(inputs);
/* if name starts with '@' it is a file otherwise the coordinates follow */
       if( inputs[0] == '@') {
         strncpy(coord_file,inputs+1,BUFFER_LEN);
         printf("   Coordinate file name: %s\n",coord_file);}
       else { /* coordinates follow */
            ControlBlock.ReadType = 1;
            GetCoordinates();}

/* eigenvector file name */
       gets(eigv_file);
       printf("   Eigenvector file name: %s\n",eigv_file);
/* minx , maxx , miny , maxy , minz , maxz */ 
       gets(inputs);
       sscanf(inputs,"%f %f %f %f %f %f",&limits.minx,&limits.maxx,
                                         &limits.miny,&limits.maxy,
                                         &limits.minz,&limits.maxz);
/* number of points in x,y and z direction */
       gets(inputs);
       sscanf(inputs,"%d %d %d",&steps.x,&steps.y,&steps.z);

       if(DENSITY_limits[0] != '\0') ParseDensitylimits();

       printf("   xmin: %f xmax: %f\n\
   ymin: %f ymax: %f \n\
   zmin: %f zmax: %f \n",
                                 limits.minx,limits.maxx,
                                 limits.miny,limits.maxy,
                                 limits.minz,limits.maxz);
       printf("Number of points in x: %d y: %d z: %d directions \n",
               steps.x,steps.y,steps.z);
/* orbital to plot                          */
       gets(inputs);
       sscanf(inputs,"%d %d",&plot_orb,&CalcType);
       printf("Orbital nr: %d will be plotted \n",plot_orb);
       printf("Calc type %d \n",CalcType);
       plot_orb--;

}

/****************************************************************************/
psi2(orb)
   int orb;
/****************************************************************************/
{
     static int i,j,k,Loop,Items;
     static int ix,iy,iz;
     static float sum,dist,dist2;
     static float xc,yc,zc;
     static float deltx,delty,deltz;
     static float phelp,dhelp;
     static int shape[3];
     static float Buffer[BUFF_LEN];
     static float MaxV = -1.0e10;
     static float MinV =  1.0e10;

     dx = (limits.maxx - limits.minx) / (float)(steps.x - 1);
     dy = (limits.maxy - limits.miny) / (float)(steps.y - 1);
     dz = (limits.maxz - limits.minz) / (float)(steps.z - 1);
  
/* main loops */


/* open output file */
   printf("   Open output file '%s' for the points\n",point_file);

   shape[0] = steps.z;
    shape[1] = steps.y;
     shape[2] = steps.x;

   printf("\n==> Writing contour information to file '%s'\n",point_file);

   i = SCAREwritedims(point_file , 3 , shape , 0);
       if(i) {
           printf("?ERROR - in writing to contour file \n");
           exit(1);}

      MaxV = -1.0e10;
      MinV =  1.0e10;

   Loop = 0;

   for(iz = 0 ; iz < steps.z ; iz++) {
      zc = dz * iz + limits.minz;
      for(iy = 0 ; iy < steps.y ; iy++) {
      yc = dy * iy + limits.miny;
         for(ix = 0 ; ix < steps.x ; ix++) {
         xc = dx * ix + limits.minx;

/* now loop the atoms */
        j = 0;
        sum = 0.0;
        for(i = 0 ; i < numat ; i++) {

        deltx = xc - x[i];
        delty = yc - y[i];
        deltz = zc - z[i];

        dist = sqrt(deltx*deltx+delty*delty+deltz*deltz);

        switch(atmtype[i]) {

           case 0: /* C */

               SUM_Tools(0);

               break;

	     case 14: /* H */
                /* s-orbital */
               sum += eigv[orb][j] * rnorms[14] * exp(-exps[14] * dist);
                    j++;
               break;

            case 1: /* N */

               SUM_Tools(1);

               break;

	     case 2: /* O */

               SUM_Tools(2);

               break;

	     case 3: /* F */

               SUM_Tools(3);

               break;

             case 4: /* Si */

               SUM_Tools1(4);

               break;

             case 5: /* P */

               SUM_Tools1(5);

               break;

             case 6: /* S */

               SUM_Tools1(6);

               break;

             case 7: /* Cl */

               SUM_Tools1(7);

               break;

             case 8: /* Li */
 
                /* s-orbital */
               sum += dist * eigv[orb][j] * rnorms[8] * exp(-exps[8] * dist);
                    j++;

               break;

             case 10: /* B */

               SUM_Tools(10);

               break;

             case 11: /* Na */

               sum += dist*dist * eigv[orb][j] * 
               rnorms[11] * exp(-exps[11] * dist);
               j++;

               break;

             case 13: /* Al */

               SUM_Tools1(13);

               break;

             case 15: /* Ca */

               SUM_Tools2(15);

               break;

             default:
             printf("Atom type %d not yet implemented \n",atmtype[i]);
	     }   /* end of switch */
             }

            if(sum < MinV) MinV = sum;
            if(sum > MaxV) MaxV = sum;

            if(Loop == BUFF_LEN){
              FWRITEN(Buffer , BUFF_LEN , sizeof(float));
                Loop = 0;}

              Buffer[Loop] = sum*sum*orbocc[orb];
              Loop++;
/*
             fprintf(point_f," %f \n",sum*sum*orbocc[orb]);
*/
           } /* x */
          }  /* y */
         }   /* z */


      if(Loop)
            FWRITEN(Buffer , Loop , sizeof(float));

  fclose(point_f);

  printf("Max value: %f \n",MaxV);
  printf("Min value: %f \n",MinV);
  printf("! Done.\n");
} 

/****************************************************************************/
psit()
/****************************************************************************/
{
     static int i,j,k,orb,Loop,Items;
     static int ix,iy,iz;
     static float sum,dist,sumt,dist2;
     static float xc,yc,zc;
     static float deltx,delty,deltz;
     static float phelp,dhelp;
     static int shape[3];
     static float Buffer[BUFF_LEN];
     static float MaxV = -1.0e10;
     static float MinV =  1.0e10;

     dx = (limits.maxx - limits.minx) / (float)(steps.x - 1);
     dy = (limits.maxy - limits.miny) / (float)(steps.y - 1);
     dz = (limits.maxz - limits.minz) / (float)(steps.z - 1);
  
/* main loops */

/* do the contour level plot for both +- level                  */

/* open output file */
   printf("   Open output file '%s' for the points\n",point_file);

   shape[0] = steps.z;
    shape[1] = steps.y;
     shape[2] = steps.x;

   printf("\n==> Writing contour information to file '%s'\n",point_file);

   i = SCAREwritedims(point_file , 3 , shape , 0);
       if(i) {
           printf("?ERROR - in writing to contour file \n");
           exit(1);}

#ifdef JUNK
/* print out number of points in x- , y- and z-directions */

     fprintf(point_f," 3         ! Rank value\n");
     fprintf(point_f," %d        ! number of points in the z-direction\n",steps.z);
     fprintf(point_f," %d        ! number of points in the y-direction\n",steps.y);
     fprintf(point_f," %d        ! number of points in the x-direction\n",steps.x);
     fprintf(point_f," %f  %f    ! zmin and zmax\n",limits.minz*atmu,limits.maxz*atmu);
     fprintf(point_f," %f  %f    ! ymin and ymax\n",limits.miny*atmu,limits.maxy*atmu);
     fprintf(point_f," %f  %f    ! xmin and xmax\n",limits.minx*atmu,limits.maxx*atmu);

#endif

      MaxV = -1.0e10;
      MinV =  1.0e10;
 
   Loop = 0;

   for(iz = 0 ; iz < steps.z ; iz++) {
      zc = dz * iz + limits.minz;
      for(iy = 0 ; iy < steps.y ; iy++) {
      yc = dy * iy + limits.miny;
         for(ix = 0 ; ix < steps.x ; ix++) {
         xc = dx * ix + limits.minx;

       sumt = 0;

       for(orb = 0 ; orb < numorb ; orb++) {

           if(orbocc[orb] < 1.0) continue;

/* now loop the atoms */
        j = 0;
        sum = 0.0;
        for(i = 0 ; i < numat ; i++) {

        deltx = xc - x[i];
        delty = yc - y[i];
        deltz = zc - z[i];

        dist = sqrt(deltx*deltx+delty*delty+deltz*deltz);

        switch(atmtype[i]) {

           case 0: /* C */

               SUM_Tools(0);

               break;

	     case 14: /* H */
                /* s-orbital */
               sum += eigv[orb][j] * rnorms[14] * exp(-exps[14] * dist);
                   j++;
               break;

            case 1: /* N */

               SUM_Tools(1);

               break;

	     case 2: /* O */

               SUM_Tools(2);

               break;

	     case 3: /* F */

               SUM_Tools(3);

               break;

             case 4: /* Si */

               SUM_Tools1(4);

               break;

             case 5: /* P */

               SUM_Tools1(5);

               break;

             case 6: /* S */

               SUM_Tools1(6);

               break;

             case 7: /* Cl */

               SUM_Tools1(7);

               break;

             case 8: /* Li */
 
                /* s-orbital */
               sum += dist * eigv[orb][j] * rnorms[8] * exp(-exps[8] * dist);
                    j++;

               break;

             case 10: /* B */

               SUM_Tools(10);

               break;

             case 11: /* Na */

               sum += dist*dist * eigv[orb][j] * 
               rnorms[11] * exp(-exps[11] * dist);
               j++;

               break;

             case 13: /* Al */

               SUM_Tools1(13);

               break;

             case 15: /* Ca */

               SUM_Tools2(15);

               break;

             default:
             printf("Atom type %d not yet implemented \n",atmtype[i]);
	     }   /* end of switch */
             }
             sumt += sum*sum*orbocc[orb]; 

	    }

            if(sumt < MinV) MinV = sumt;
            if(sumt > MaxV) MaxV = sumt;

            if(Loop == BUFF_LEN){
              FWRITEN(Buffer , BUFF_LEN , sizeof(float));
                Loop = 0;}

              Buffer[Loop] = sumt;
              Loop++;

           } /* x */
          }  /* y */
         }   /* z */



      if(Loop)
            FWRITEN(Buffer , Loop , sizeof(float));

  fclose(point_f);

  printf("Max value: %f \n",MaxV);
  printf("Min value: %f \n",MinV);
  printf("! Done.\n");
} 


int SCAREwritedims(filename,rank,shape,unknown)
      char *filename;
      int   rank;
      int   shape[];
      int   unknown;
{
      char OutText[BUFF_LEN];
      int Items;
      float Help;
      int TypeOfSurface;
      
      point_f = fopen(filename,"w");
      if(point_f == NULL) {
         printf("?ERROR - can't open contour file '%s'\n",filename);
         return(1);}

/* first position is the rank */
         FWRITE(&rank, sizeof(int));
/* type of surface            */
         TypeOfSurface = 2;
         FWRITE(&TypeOfSurface, sizeof(int));

/* write shape ...             */
        FWRITE(&shape[0] , sizeof(int));
        FWRITE(&shape[1] , sizeof(int));
        FWRITE(&shape[2] , sizeof(int));
       
       if(VERBOSE)
         printf("zdim: %d , ydim: %d , xdim: %d\n",shape[0],shape[1],shape[2]);
/* done with shape */

/* write min/max ... */
       Help = limits.minz*atmu;
       FWRITE(&Help , sizeof(float));
        Help = limits.maxz*atmu;
        FWRITE(&Help , sizeof(float));
         Help = limits.miny*atmu;
         FWRITE(&Help , sizeof(float));
          Help = limits.maxy*atmu;
          FWRITE(&Help , sizeof(float));
           Help = limits.minx*atmu;
           FWRITE(&Help , sizeof(float));
            Help = limits.maxx*atmu;
            FWRITE(&Help , sizeof(float));

/* done with min/max */

    return(0);
  }

GetCoordinates()
{
     char inputs[BUFFER_LEN];
     char Temp[2];
     int  i;

     gets(inputs);
     sscanf(inputs,"%d",&numat);

/* go and get space for the coordinate arrays and the atom name array */

    x = (float *) malloc(numat * sizeof(float));
    y = (float *) malloc(numat * sizeof(float));
    z = (float *) malloc(numat * sizeof(float));

       printf("Number of atoms: %d \n",numat);

    atnam = (char *) malloc( 2 * numat );

/* check that everything is ok */

    if(x == NULL || y == NULL || z == NULL || atnam == NULL) {
      printf("?ERROR - memory allocation problems \n");
      exit(1);}

/* start reading coordinates */

    for(i = 0 ; i < numat ; i++) {
       gets(inputs);
       sscanf(inputs,"%f %f %f %s",&x[i],&y[i],&z[i],Temp);

       if(Temp[1] == '\0') {
         atnam[2*i] = ' ';
          atnam[2*i + 1] = Temp[0];}
       else
         strncpy(atnam+2*i , Temp , 2);

       printf(":%d: x %f y %f z %f '%.2s'\n",(i+1),x[i],y[i],z[i],atnam+2*i);

     }

}

ParseDensitylimits()
{
     char *find;
     int loop;

     printf(".... > Redefining input limits < ....\n");

     find = strtok(DENSITY_limits,SEP_STRING);

     if(find ==NULL) return;

     steps.x = atoi(find);     

     for(loop = 1; loop < 9 ; loop++) {

        find = strtok(NULL,SEP_STRING);

        if(find == NULL) return;

        switch(loop) {
        case 1: /* ypts */
               steps.y = atoi(find);
               break;      
        case 2: /* zpts */
               steps.z = atoi(find);
               break;
        case 3: /* xmin */
               limits.minx = atof(find);
               break;      
	case 4: /* xmax */
               limits.maxx = atof(find);
               break;
        case 5: /* ymin */
               limits.miny = atof(find);
               break;      
	case 6: /* ymax */
               limits.maxy = atof(find);
               break;
        case 7: /* zmin */
               limits.minz = atof(find);
               break;      
	case 8: /* zmax */
               limits.maxz = atof(find);
               break;}
      }
}
