/*Copyright 2002-3 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
/*Modified 2003 by Eero Häkkinen, CSC
  Replaced conio functions by stdio funkcions.*/
/*Really only for debugging purposes.*/
#include "keys.h"

namespace Plugin {
namespace Fchk {
    
void PauseForAction()
{
/* Waits for user to press key.  Diagnostic tool. */
    int g;

    g=getchar();
}

void Assert(void *flag,char* text,int Pause)
{
    int* iflag=(int*)flag;
    if(!(*iflag))
        printf("\n%s\n",text);
    if(Pause)
        PauseForAction();
}

} // namespace Fchk
} // namespace Plugin
