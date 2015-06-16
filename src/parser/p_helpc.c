/*

Copyright (c) 1996 - 2004 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include <ctype.h>
#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>
#include <sys/types.h>

#include "gomenv.h"
#include "gomhelp.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "printmsg.h"
#include "tclutils.h"

#include "stdafx.h"

#define MAIN_HELP        "main_help.html"
#define ATOM_HELP        "atom_help.html"
#define CALCULATE_HELP   "calculate_help.html"
#define CENTER_HELP      "center_help.html"
#define CONTOUR_HELP     "contour_help.html"
#define DEFINE_HELP      "define_help.html"
#define DISPLAY_HELP     "display_help.html"
#define EDIT_HELP        "edit_help.html"
#define EXPORT_HELP      "export_help.html"
#define FILL_HELP        "fill_help.html"
#define FIND_HELP        "find_help.html"
#define HARDCOPY_HELP    "hardcopy_help.html"
#define HELP_HELP        "main_help.html"
#define IMPORT_HELP      "import_help.html"
#define MANIPULATE_HELP  "manipulate_help.html"
#define MONITOR_HELP     "monitor_help.html"
#define MOPEN_HELP       "mopen_help.html"
#define MSAVE_HELP       "msave_help.html"
#define PLOT_HELP        "plot_help.html"
#define PLUMBER_HELP     "plumber_help.html"
#define RESET_HELP       "reset_help.html"
#define ROTATE_HELP      "rotate_help.html"
#define RUN_HELP         "run_help.html"
#define SCALE_HELP       "scale_help.html"
#define SELECT_HELP      "select_help.html"
#define SHOW_HELP        "show_help.html"
#define TRACE_HELP       "trace_help.html"
#define TRAJECTORY_HELP  "trajectory_help.html"
#define TRANSLATE_HELP   "translate_help.html"
#define WINDOW_HELP      "window_help.html"

static struct _NetHelp {
    char HURL[BUFF_LEN];
    char WebBrowser[BUFF_LEN];
} NetHelp;

static int         PrintHelpInfo(const char *,const char *);

/*********************************************************************/
int gomp_HelpCommand(ClientData clientdata, Tcl_Interp *interp,
                   int argc, const char **argv)
/*********************************************************************/
{
    static char Text1[BUFF_LEN];
    static char Text2[BUFF_LEN];

/* #1   help  ... */
    strncpy(Text1,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text1 , "help")) {

        gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
        gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
/* node */
        if(gomp_StringMatch(Text1 , "?") ||
           (strlen(Text1) == 0)     ||
           gomp_StringMatch(Text1 , "help")) {
            if(!PrintHelpInfo(MAIN_HELP,""))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
        else if(gomp_StringMatch(Text1 , "atom")) {
            if(!PrintHelpInfo(ATOM_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "calc$ulate")) {
            if(!PrintHelpInfo(CALCULATE_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "cent$er")) {
            if(!PrintHelpInfo(CENTER_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "conto$ur")) {
            if(!PrintHelpInfo(CONTOUR_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "defi$ne")) {
            if(!PrintHelpInfo(DEFINE_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "disp$lay")) {
            if(!PrintHelpInfo(DISPLAY_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "edit")) {
            if(!PrintHelpInfo(EDIT_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "expo$rt")) {
            if(!PrintHelpInfo(EXPORT_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "fill")) {
            if(!PrintHelpInfo(FILL_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "find")) {
            if(!PrintHelpInfo(FIND_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "hard$copy")) {
            if(!PrintHelpInfo(HARDCOPY_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "help")) {
            if(!PrintHelpInfo(HELP_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "impo$rt")) {
            if(!PrintHelpInfo(IMPORT_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "mani$pulate")) {
            if(!PrintHelpInfo(MANIPULATE_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "moni$tor")) {
            if(!PrintHelpInfo(MONITOR_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "mope$n")) {
            if(!PrintHelpInfo(MOPEN_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "msav$e")) {
            if(!PrintHelpInfo(MSAVE_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "plot")) {
            if(!PrintHelpInfo(PLOT_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "plum$ber")) {
            if(!PrintHelpInfo(PLUMBER_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "rese$t")) {
            if(!PrintHelpInfo(RESET_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "rota$te")) {
            if(!PrintHelpInfo(ROTATE_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "run")) {
            if(!PrintHelpInfo(RUN_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "scal$e")) {
            if(!PrintHelpInfo(SCALE_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "sele$ct")) {
            if(!PrintHelpInfo(SELECT_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "show")) {
            if(!PrintHelpInfo(SHOW_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "trac$e")) {
            if(!PrintHelpInfo(TRACE_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "traj$ectory")) {
            if(!PrintHelpInfo(TRAJECTORY_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "tran$slate")) {
            if(!PrintHelpInfo(TRANSLATE_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
        else if(gomp_StringMatch(Text1 , "wind$ow")) {
            if(!PrintHelpInfo(WINDOW_HELP,Text2))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }   
/* use a web browser to get the help */
        else if(gomp_StringMatch(Text1 , "webb$rowser")) {
	  if(strlen(gomp_GetWebBrowser()) == 0) {
                gomp_PrintERROR("web browser is not defined");
                return(TCL_ERROR);
            }
	  if(strlen(gomp_GetHURL()) == 0) {
                gomp_PrintERROR("URL to the help files is not defined");
                return(TCL_ERROR);
            }
            sprintf(Text1,"%s %s",gomp_GetWebBrowser(),gomp_GetHURL());
            return(system(Text1));
        }   
    }
    else  {
        gomp_PrintERROR("'help' command not recognized");
        return(TCL_ERROR);
    }

/*  E R R O R command not recognized         */
    gomp_PrintERROR("'help' command not recognized");
    return(TCL_ERROR);

}

/*********************************************************************/
int PrintHelpInfo(const char *FileName,const char *anchor)
/*********************************************************************/
{
    static FILE *File_p;
    static char  Path[3*BUFF_LEN];
    char *line,*tag=NULL,*p,*q,*r;
    size_t length;
    int    PrintAnchors;

#if defined(WIN32)
    sprintf(Path,"%s\\%s",gomp_ShowHelpDir(),FileName);
#else
    sprintf(Path,"%s/%s",gomp_ShowHelpDir(),FileName);
#endif

    line = NULL;
    File_p = fopen(Path,"r");
    if(File_p == NULL) {
        gomp_PrintERROR("can't open help file:");
        gomp_PrintERROR(FileName);
        return(1);
    }

    /* Find <PRE> tag. */
    while(fgets(Path,BUFF_LEN,File_p) != NULL) {
        for( tag = Path - 1; (tag = strchr(tag+1,'<')); ) {
            if( toupper((unsigned char)tag[1]) == 'P' &&
                toupper((unsigned char)tag[2]) == 'R' &&
                toupper((unsigned char)tag[3]) == 'E' &&
                tag[4] == '>')
                break;
        }
        if( tag ) {
            /* We found the <PRE> tag! */
            if( tag[5] && tag[5]!='\n' )
                line = tag + 5;
            break;
        }
    }

    if( anchor && *anchor ) {
        if( gomp_StringMatch(anchor, "anch$ors") ||
            gomp_StringMatch(anchor, "help") ) {
            PrintAnchors = 1;
            gomp_PrintMessage(
                "**************************************************************************");
        }
        else
            PrintAnchors = 0;

        while(fgets(Path,BUFF_LEN,File_p) != NULL) {
            for( tag = Path - 1; (tag = strchr(tag+1,'<')); ) {
                if( toupper((unsigned char)tag[1]) == 'A' &&
                    isspace((unsigned char)tag[2])) {
                    int i;
                    for( i=3; isspace((unsigned char)tag[i]); i++) 
                        ;
                    if( toupper((unsigned char)tag[i+0]) == 'N' &&
                        toupper((unsigned char)tag[i+1]) == 'A' &&
                        toupper((unsigned char)tag[i+2]) == 'M' &&
                        toupper((unsigned char)tag[i+3]) == 'E' &&
                        tag[i+4]=='=' &&
                        tag[i+5]=='\"') {
                        if( strncmp( tag+i+6 , anchor , strlen(anchor) ) == 0 &&
                            tag[i+6+strlen(anchor)]=='\"') {
                            /* We found the <A name="$$"> tag! */
                            line = strchr(tag,'>');
                            if( line )
                                line++;
                            break;
                        }
                        else if( PrintAnchors ) {
                            /* We found the <A name="$$"> tag, but */
                            /* a non-matching one.                 */
                            line = strchr(tag+i+6,'\"');
                            if(line) {
                                *line = '\0';
                                gomp_PrintMessage(tag+i+6);
                                *line = '\"';
                            }
                        }
                    }
                }   
            }
            if( tag )
                /* We found the <A name="$$"> tag! */
                break;
        }
        if( !tag ) {
            if( PrintAnchors ) {
                gomp_PrintMessage(
                    "**************************************************************************");
                fclose(File_p);
                return(0);
            }
            else {
                sprintf(Path,"Anchor not found in help file: %s#%s",FileName,anchor);
                gomp_PrintERROR(Path);
                fclose(File_p);
                return(1);
            }
        }
        gomp_PrintMessage(
            "**************************************************************************");
    }

    while(1) {
        if( !line ) {
            if( fgets(Path,BUFF_LEN,File_p) == NULL )
                break;
            line = Path;
        }
        q = strchr(line,'<');
        r = strchr(line,'&');
        if( q  || r) {
            if( r && ( !q || r<q) )
                q = r;
            p = q;

            while( *q ) {
                if( *q == '<' ) {
                    /* Omit tags. */
                    if( tag[1] == '/' &&
                        toupper((unsigned char)tag[2]) == 'P' &&
                        toupper((unsigned char)tag[3]) == 'R' &&
                        toupper((unsigned char)tag[4]) == 'E' &&
                        tag[5] == '>') {

                        fclose(File_p);
                        return(0);
                    }
                    if( anchor && *anchor &&
                        toupper((unsigned char)tag[1]) == 'A' &&
                        isspace((unsigned char)tag[2]) ) {

                        fclose(File_p);
                        gomp_PrintMessage(
                            "**************************************************************************");
                        return(0);
                    }
                    q = strchr(q,'>') + 1;
                    if( !(q-1) || !*q )
                        break;
                }
                else if( *q == '&' ) {
                    /* Handle some entities. */
                    if( strncmp( q , "&lt;" , 4 ) ) {
                        *p++ = '<';
                        q   += 4;
                    }
                    else if( strncmp( q , "&gt;" , 4 ) ) {
                        *p++ = '>';
                        q   += 4;
                    }
                    else
                        *p++ = *q++;
                }
                else
                    *p++ = *q++;
            }
            *p = '\0';
        }

        length = strlen(line);
        if( length && line[length-1] == '\n' )
            line[length-1] = '\0';
        gomp_PrintMessage(line);
        line = NULL;
    }

    fclose(File_p);

    return(0);
}

/*********************************************************************/
int   gomp_SetHURL(const char *HURL)
/*********************************************************************/
{
    gomp_CopyString(NetHelp.HURL,HURL,BUFF_LEN);

    return(0);
}
/*********************************************************************/
int   gomp_SetWebBrowser(const char *WebBrowser)
/*********************************************************************/
{
    gomp_CopyString(NetHelp.WebBrowser,WebBrowser,BUFF_LEN);

    return(0);
}
/*********************************************************************/
const char *gomp_GetHURL()
/*********************************************************************/
{
    return(NetHelp.HURL);
}

/*********************************************************************/
const char *gomp_GetWebBrowser()
/*********************************************************************/
{
    return(NetHelp.WebBrowser);
}

