/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
#include <stdlib.h>
#include <stdio.h>

#ifndef parser_included
namespace Plugin {
namespace Fchk {
    
int ReadToEOL(FILE* ftr,char* rs);
int SplitString(char* IString,char* OString,int Start);
int cmp(const char* c1,const char* c2);

} // namespace Fchk
} // namespace Plugin
#define parser_included 1
#endif
