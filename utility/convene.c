#include <stdio.h>
#include <math.h>

#define BUFF_LEN 120

/*

    This program reads in (from standard input) an energy file (*.ene).
    The CHARMM22 energy file contains the information on 3 separate
    lines (not on one as before). This program writes out the 3 lines
    on one line on standard output.

    Leif Laaksonen CSC 1992

*/

/************************************************************************/
main() /* Convert energy file from CHARMM22 to contain all info on one line */
/************************************************************************/

{


   int et_count;
   char InBuff[BUFF_LEN];
   int i;


   et_count = 0;

/* ready to continue ... */

    
    while(gets(InBuff) != NULL) {
     printf("%s ",InBuff);
      et_count++;
      if(!(et_count%3)) printf("\n");}
}
