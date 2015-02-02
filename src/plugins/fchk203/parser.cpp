/*Copyright 2003 Kevin J. Boyd, University of New Orleans
Permission granted to distribute and modify this code for personal, educational, and
research use. */
/*Contains text-handling routines used in several locations*/
#include "parser.h"

namespace Plugin {
namespace Fchk {
    
int ReadToEOL(FILE* ftr,char* rs)
{
    char c1='\000';
    char EOLn='\012';
    int nc=0;

    while((c1!=EOLn)&&(!feof(ftr)))
    {
        fscanf(ftr,"%c",&c1);
        rs[nc++]=c1;
    }
    rs[nc]='\000';
/*  if(c1==EOLn)
    {
        fscanf(ftr,"%c",&c1);
        nc++;
    }*/
    if(!feof(ftr))
        return ++nc;
    else
        return -(++nc);
}


char StringDelim[11]=" ,\n\";[]@()";

int SplitString(char* IString,char* OString,int Start)
{
    char c1;
    int PosInString=0;
    int flag=0;
    int ctr;
    
    while(!flag)
    {
        c1=IString[PosInString+Start];
        for(ctr=0;ctr<11;ctr++)
        {
            if(c1==StringDelim[ctr])
            {
                flag=1;
            }
        }
        if(!flag)
        {
            OString[PosInString++]=c1;
        }
    }
    OString[PosInString++]='\000';
    return PosInString+Start;
}

int cmp(const char* c1,const char* c2)
{
    int i;

    for(i=0;;i++)
    {
        if(c1[i]!=c2[i])
                return i+1;
        if(c1[i]=='\000')
                return 0;
    }
}

} // namespace Fchk
} // namespace Plugin
