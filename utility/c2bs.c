
/* filter for converting a CHARMM file to Ball&Stick format

   Leif Laaksonen 1991
   
*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#define BUFF_LEN 256

  char *Tell = 
"*********************************************************\n\
 * CHARMM to Ball and Stick file converter               *\n\
 * Leif Laaksonen 1991                                   *\n\
 *********************************************************\n";
 
main(argc , argv)
    int   argc;
    char *argv[];
{

    FILE *Input_FP;   /* input file pointer */
    static char  InputFileName[BUFF_LEN];
    FILE *Output_FP;  /* output file pointer */
    static char  OutputFileName[BUFF_LEN] = "BallStick.input";
    int  NumAtoms,i,j;
    char InputLine[BUFF_LEN];
    char AtmName[BUFF_LEN];
    float Xc,Yc,Zc;
    
    printf("%s",Tell);
    
    printf("Give input file name : ");
    gets(InputFileName);
    
    Input_FP = fopen(InputFileName,"r");
    if(Input_FP == NULL) {
       printf("?ERROR - can't open input file '%s' \n",InputFileName);
       exit(1);}
       
    Output_FP = fopen(OutputFileName,"w");
    if(Output_FP == NULL) {
       printf("?ERROR - can't open output file '%s' \n",OutputFileName);
       exit(1);}

    printf("Output file name '%s' \n",OutputFileName);

again:;
    
    fgets(InputLine , BUFF_LEN , Input_FP);

    if(InputLine[0] == '*') {
       printf("%s\n",InputLine);
        goto again;}
    
    NumAtoms = atoi(InputLine);
    
    printf("Converter found %d atoms\n",NumAtoms);
    
    fprintf(Output_FP,"TITLE (Just a dummy title)\n");
    
    for(i = 0 ; i< NumAtoms ; i++) {
        fgets(InputLine , BUFF_LEN , Input_FP);
        sscanf(InputLine,"%*d %*d %*s %s %f %f %f",AtmName,&Xc,&Yc,&Zc);
        
        printf("%.4s \t %f \t %f \t %f \n",AtmName,Xc,Yc,Zc);
        
        fprintf(Output_FP,"%.4s \t %f \t %f \t %f\n",AtmName,Xc,Yc,Zc);}
        
    fprintf(Output_FP,"END\n");      
    printf("Done ...\n");
}
