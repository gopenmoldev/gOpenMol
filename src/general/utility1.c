/*  

Copyright (c) 1994 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing, ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved
  


Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include <ctype.h>
#include <stdarg.h>
#include "gomstdio.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <tcl.h>

#if defined(WIN32)
#include <windows.h>
#include <winsock.h>
#include <time.h>
#include <process.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#if defined(IRIX)
#include <sys/procfs.h>
#endif

#include "colouring.h"
#include "coord_file.h"
#include "g_status.h"
#include "gomcast.h"
#include "gomenv.h"
#include "gomfile.h"
#include "gomlog.h"
#include "gommain.h"
#include "gomproc.h"
#include "gomtext.h"
#include "gomversion.h"
#include "gomstring.h"
#include "memalloc.h"
#include "model_file.h"
#include "molecule.h"
#include "molecstruct.h"
#include "printmsg.h"
#include "projview.h"
#include "text_stack.h"
#include "tclutils.h"

#include "stdafx.h"

#define MAX_FLOAT  1.e+35    /* This max float is just a guess */
#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))
#define Rabs(a)    ( ( a ) > 0   ? (a) : -(a))
#define Fabs(a)    ( ( a ) > 0.0 ? (a) : -(a))
#define SCALE(a,b,c)     scaleO(a,b,c) /* use own scale function */

/* functions */
static int Get_wall_secs(float *);
static int Get_env(void); /* go and hunt for environment parameters */
static const char *GetModuleName(void);
static const char *AtomConversionFile(void);
static const char *AtomParametersFile(void);
static const char *ElementParametersFile(void);
static const char *ColourTableFile(void); 
static int InformIObytes(void);
static int JoinPath(char *,const char *, const char *);

/* environment related stuff */
static const char *GetGomTempDir(void);
static const char *GetGomDirectory(const char *);
static void SetGomDirectory(const char *, const char *, int);
static int  SetValuesFromEnvironment(void);
static struct {                     /* environment file structure */
    char File[BUFF_LEN];
} Environment;

static const char *GetEnvironmentFile(void);

/* structures */

static struct {                     /* log file structure */
    char file[BUFF_LEN];
    FILE *log_p;
    int ok;
    int lines;
} LogFile = { "", NULL , 1 };

static struct {
    char ElementParameters[BUFF_LEN];
    char AtomConversion[BUFF_LEN];
    char AtomParameters[BUFF_LEN];
    char ColourTable[BUFF_LEN];
} FileNames = 
    {ELEM_PARAM , ATOM_EQ , ATOM_PARAM , COLOUR_TABLE};
#if 0
static struct {
    char TERM[BUFF_LEN];             /* terminal */
    char DISPLAY[BUFF_LEN];          /* display  */
} Terminal;
#endif

/* coordinate and trajectory default file types */
static struct {
    int  Length;
    int  AMBERdefault;
    char AMBERfileType[BUFF_LEN];
    int  CHARMMdefault;
    char CHARMMfileType[BUFF_LEN];
    int  GAMESSdefault;
    char GAMESSfileType[BUFF_LEN];
    int  GAUSSIANdefault;
    char GAUSSIANfileType[BUFF_LEN];
    int  GROMACSdefault;
    char GROMACSfileType[BUFF_LEN];
    int  HYPERCHEMdefault;
    char HYPERCHEMfileType[BUFF_LEN];
    int  INSIGHTdefault;
    char INSIGHTfileType[BUFF_LEN];
    int  JAGUARdefault;
    char JAGUARfileType[BUFF_LEN];
    int  MOL2default;
    char MOL2fileType[BUFF_LEN];
    int  MOPACgraphdefault;
    char MOPACgraphfileType[BUFF_LEN];
    int  MUMODdefault;
    char MUMODfileType[BUFF_LEN];
    int  OPENMOLdefault;
    char OPENMOLfileType[BUFF_LEN];
    int  PDBdefault;
    char PDBfileType[BUFF_LEN];
    int  TINKERdefault;
    char TINKERfileType[BUFF_LEN];
    int  UHBDdefault;
    char UHBDfileType[BUFF_LEN];
    int  XMOLdefault;
    char XMOLfileType[BUFF_LEN];
    int  XYZdefault;
    char XYZfileType[BUFF_LEN];
    int  YASPdefault;
    char YASPfileType[BUFF_LEN];
    int  OPENMOL_center_default;
    char OPENMOL_center_fileType[BUFF_LEN];
} CoordFileType = { 18 ,
                    0 ,
                    "" ,     /* amber */
                    1 ,
                    "crd" ,  /* charmm */
                    0 ,
                    "dat" ,  /* gamess */
                    0 ,
                    "" ,     /* gaussian */
                    0 ,
                    "gro" ,  /* gromac */
                    0 ,
                    "" ,     /* hyperchem */
                    0 ,
                    "car" ,  /* insight */
                    0 ,
                    "" ,     /* jaguar */
                    0 ,
                    "mol2" , /* mol2 */
                    0 ,
                    "gpt" ,  /* mopac graphics file */
                    0 ,
                    "",      /* mumod */
                    0 ,
                    "" ,     /* openmol */
                    0 ,
                    "pdb" ,  /* pdb */
                    0 ,
                    "xyz" ,  /* TINKER */
                    0 ,
                    "qcd" ,  /* UHBD */
                    0 ,
                    "xmol",  /* xmol */
                    0 ,
                    "xyz",   /* xyz */
                    0 ,
                    "",      /* yasp */
                    0 ,
                    "cbi__binary_input.data"};

static char ModuleName[BUFF_LEN];

static int  SuppressErrors = 0;

/* atm, res and seg name stacks */ 
static char *atnam_stack;    /* stack to contain the different atnam, resnam */
static int   atnam_stack_deep=0;
static int   atnam_stack_max=100;
static int *atnam_stack_num;
static char *resnam_stack;   /* and segment names */
static int   resnam_stack_deep=0; 
static int   resnam_stack_max=50;
static int *resnam_stack_num;
static char *segment_stack;
static int   segment_stack_deep=0;
static int   segment_stack_max=30;
static int *segment_stack_num;

/* B A N N E R                              */

static const char *gOpenMolBanner[] = {
"***************************************************************",
"",
"* * * * *   g O p e n M o l   * * * * *",
"",
"",
GOPENMOL_COPYRIGHT,
"",
"\t", /* Version */
"",
"***************************************************************",
NULL};

static const char *gOpenMolUsage[] = {
"  *******************************************************************",
"  *                                                                 *",
"  *  Help page to get you off the ground.                           *",
"  *                                                                 *",
"  *  Following options allowed on the command line:                 *",
"  *                                                                 *",
"  *  gopenmol -cICoord.Ext -h -t -ln -mModel.gom -pFile.Name File.n *",
"  *                                                                 *",
"  *  -mModel.gom start GopenMol with the model file Model.gom       *",
"  *  -cICoord.Ext start gOpenMol with coordinate file               *",
"  *  -pFile.Name start gOpenMol with parameter file 'File.Name'     *",
"  *  -t start gOpenMol without any graphics                         *",
"  *  -ln do not write a log file  (GOPENMOL.LOG)                    *",
"  *  -ly force the write of the log file (GOPENMOL.LOG)             *",
"  *  -s hardware stereo (if available)                              *",
"  *  -h print this help                                             *",
"  *  File.n start gOpenMol with the coordinate file File.n          *",
"  *                                                                 *",
"  *******************************************************************",
NULL};

/* end of declarations */

/***********************************************************************/
const char *gomp_ShowRootDir()
/***********************************************************************/
{
    return GetGomDirectory("GOM_ROOT");
}
/***********************************************************************/
const char *gomp_ShowHomeDir()
/***********************************************************************/
{
    return GetGomDirectory("GOM_HOME");
}
/***********************************************************************/
const char *gomp_ShowBinDir()
/***********************************************************************/
{
    return GetGomDirectory("GOM_BIN");
}
/***********************************************************************/
const char *gomp_ShowDataDir()
/***********************************************************************/
{
    return GetGomDirectory("GOM_DATA");
}
/***********************************************************************/
const char *gomp_ShowHelpDir()
/***********************************************************************/
{
    return GetGomDirectory("GOM_HELP");
}
/***********************************************************************/
const char *gomp_ShowTempDir()
/***********************************************************************/
{
    return GetGomDirectory("GOM_TEMP");
}
/***********************************************************************/

/*
  Get the time (in seconds) since last called this function.

  Returns first time when it is called 0.0 (resets itself)
  and after that the delta time since first time called.

  Leif Laaksonen (1994)
*/
/***********************************************************************/
int Get_wall_secs(float *TimeSecs)
/***********************************************************************/
{
    static int Switch = -1;

    static time_t TimeAtStart = 0;
    static time_t GetTime = 1;

    (void)time(&GetTime);

    *TimeSecs = 0.0;

    if(Switch < 0) {
        TimeAtStart = GetTime;
        Switch = 1;
    }
    else {
        *TimeSecs = (float)(GetTime - TimeAtStart);
    }

    return(0);
}
#if 0
/***********************************************************************/
int gomp_GetWallSeconds()
/***********************************************************************/
{      
    return((int)time(0));
}
#endif
/***********************************************************************/
int gomp_StartUpStuff(int argc , const char *argv[])
/***********************************************************************/
{
    float TimeNow = 0.0;

/* setup the interrupts         */
    gomp_DefineSignals();

/* set timer and go ... */
    Get_wall_secs(&TimeNow);

/* The execution starts here ... */

/* get the environment parameters */
    if(Get_env())
        gomp_PrintEXIT("can't get environment parameters");
/*
  if(gomp_ReadCommandStack())
  gomp_PrintEXIT("$Can't read old command stack");

  if(gomp_ParseStartUpLine(argc , argv))
  return(1);
*/
/* Read atom parameter file */
    if(gomp_ReadElementParams(ElementParametersFile()))
        gomp_PrintEXIT("$Can't read element parameter file");

/* Read atom parameter file */
    if(gomp_ReadAtomParams(AtomParametersFile()))
        gomp_PrintEXIT("$Can't read atom parameter file");

/* Read atom conversion table */
    if(gomp_ReadAtomConversionTable(AtomConversionFile()))
        gomp_PrintEXIT("$Can't read atom conversion file");

/* Read colour table */
    if(gomp_ReadColourTable(ColourTableFile()))
        gomp_PrintEXIT("$Can't read colour table file");

    return(0);
}


/**************************************************************************/
int gomp_OpenLogFile()
/**************************************************************************/
{
#if defined(WIN32)
    char Text[BUFF_LEN];
    const char *Temp;

    if(LogFile.log_p) {
        gomp_PrintERROR("log file is already open");
        return(1);
    }

    Temp = GetGomTempDir();
    if( *Temp == '\0' ) {
        gomp_PrintERROR("Can't open log file");
        LogFile.ok = 1;
        return(1);
    }

    strncpy(LogFile.file,LOGFILE,strlen(LOGFILE));
    sprintf(Text,"%s\\%s",Temp,LogFile.file);
    printf("    Log file is : '%s'\n",Text);
    

    LogFile.lines = 0;

    if(LogFile.ok) return(0);  /* don't use log file */

    LogFile.log_p = fopen(Text,"w");
    
    if(LogFile.log_p == (FILE *)NULL) {
        gomp_PrintERROR("?Can't open log file ");
        LogFile.ok = 1;
        return(1);
    }

#else
    strncpy(LogFile.file,LOGFILE,strlen(LOGFILE));

    LogFile.lines = 0;

    if(LogFile.ok) return(0);  /* don't use log file */

    LogFile.log_p = fopen(LogFile.file,"w");

    if(LogFile.log_p == NULL) {
        gomp_PrintMessage("?Can't open log file ");
        LogFile.ok = 1;
        return(1);
    }

#endif

    LogFile.ok = 0;

    return(0);


}
/**************************************************************************/
int gomp_WriteToLogFile(const char *text)
/**************************************************************************/
{
    if(LogFile.ok) return(0);  /* don't write to log file */

    if(*text == (char)NULL) return(0);

    if(fprintf(LogFile.log_p,"%s\n",text) < 0) {
        printf("?Can't write to log file\n");
        return(1);
    }

    fflush(LogFile.log_p);

    return(0);
}
/**************************************************************************/
int gomp_CloseLogFile()
/**************************************************************************/
{

    if(LogFile.ok) return(0);  /* no log file open ... */

    if(fclose(LogFile.log_p) == EOF) {
        gomp_PrintERROR("can't close log file");
        return(1);
    }

    return(0);
}
/***********************************************************************/
int gomp_PrintMessage(const char *Text)
/***********************************************************************/
{
    if ( printf("%s\n",Text) == EOF )
        return(1);

    if ( gomp_PushText2Stack("",Text) )
        return(1);

    if ( gomp_WriteToLogFile(Text) )
        return(1);

    return(0);
}
/***********************************************************************/
int gomp_PrintERROR(const char *Text)
/***********************************************************************/
{
    static char SavedText[BUFF_LEN];
    int         Code;

    if ( Text )
        gomp_SendTclReturn(Text);

/* originally implemented as:     
    if ( ! SuppressErrors ) {

    However, I can't get it to work preperly because gOpenMol comes
    here with SuppressErrors == 1 and prevents the error output!
    
    Only if executed from inside eval {}.
    
    Anyway, errors MUST be shown then errors are not to be suppressed.
*/

    if ( ! Text )
        Text = Tcl_GetStringFromObj(Tcl_GetObjResult(gomp_GetTclInterp()),NULL);

    if ( ! SuppressErrors ) {
        if ( printf("!gOpenMol - ERROR: %s\n",Text) == EOF )
            return(1);

        if ( gomp_GetTermType() == GRAPHICS_AVAILABLE &&
             strncmp(Text,SavedText,BUFF_LEN) ) {
            /* Save the real result. */
            Tcl_SavedResult resultState;
            Tcl_Obj *cmd;
            Tcl_SaveResult(gomp_GetTclInterp(),&resultState);
            /* Generate command as Tcl list. */
            cmd  = gomp_CreateTclList("%s %s","lulErrorDialog",Text);

            Tcl_IncrRefCount(cmd);
            Code = Tcl_GlobalEvalObj(gomp_GetTclInterp(),cmd);
            Tcl_DecrRefCount(cmd);

            strncpy(SavedText,Text,BUFF_LEN);
            if(Code != TCL_OK) {
                gomp_PrintMessage("can't put this error message to widget:");
                gomp_PrintMessage("internal error:");
                gomp_PrintMessage(Tcl_GetStringResult(gomp_GetTclInterp()));
                gomp_PrintMessage("error message:");
                gomp_PrintMessage(Text);
                /* Restore the real result. */
                Tcl_RestoreResult(gomp_GetTclInterp(),&resultState);
                return(1);
            }
            /* Restore the real result. */
            Tcl_RestoreResult(gomp_GetTclInterp(),&resultState);
        }
    }

    if ( gomp_PushText2Stack("!gOpenMol - ERROR: ",Text) ||
         gomp_WriteToLogFile(Text) )
        return(1);

    return(0);
}

/***********************************************************************/
int gomp_PrintWARNING(const char *Text)
/***********************************************************************/
{
    if ( printf("!gOpenMol - WARNING: %s\n",Text) == EOF )
        return(1);

    if ( gomp_PushText2Stack("!gOpenMol - WARNING: ",Text) )
        return(1);

    if ( gomp_WriteToLogFile(Text) )
        return(1);

    return(0);
}

/***********************************************************************/
int gomp_PrintEXIT(const char *Text)
/***********************************************************************/
{
    if(printf("$gOpenMol - EXIT: %s\n",Text) == EOF)
        exit(1);

    if(gomp_WriteToLogFile(Text))
        exit(1);

    exit(0);
}
/************************************************************************/
int gomp_SuppressErrors(int suppress)
/************************************************************************/
{
    int oldSuppress = SuppressErrors;
    SuppressErrors  = suppress;
    return oldSuppress;
}

/************************************************************************/
int gomp_Get_date(int alt)        /* get and print date and time */
/************************************************************************/
{
    time_t t;
    const char *date_time;
    float secs,secs1;

    time(&t);
    if(t == -1) {
        gomp_PrintERROR("Can't get the time ");
        return(1);
    }

    date_time = ctime(&t);

#if !defined(WIN32)
    if(gomp_RunStatistics() == -1) {
        gomp_PrintERROR("?Can't get information about the resource utilization");
    }
#endif

    switch(alt) {

    case 1:  /* print date and time at start */
        gomp_FormatMessage("    ****** Starting on host '%s' on %s ",gomp_LongHostName(),date_time);
        break;

    case 2: /* print date and time at the end */
        (void)gomp_Get_cpu_secs(&secs , &secs1);
        gomp_FormatMessage("    ****** Exit on host '%s' on %s ",gomp_LongHostName(),date_time);
        gomp_FormatMessage("   Wasted cpu secs: %.2f (Children: %.2f secs)",secs,secs1);
        (void)Get_wall_secs(&secs1);
        gomp_FormatMessage("   In %.3f hours of wall clock time (%.3f%c)\n",
                secs1/3600.,(secs/secs1)*100.,'%');
        (void) InformIObytes();

        break;    
 
    case 3: /* print general date and time */
        gomp_FormatMessage("    ****** Host  '%s' on %s",gomp_LongHostName(),date_time);
        (void)gomp_Get_cpu_secs(&secs , &secs1);
        gomp_FormatMessage("   Wasted cpu secs: %.2f (Children: %.2f secs)",secs,secs1);
        (void)Get_wall_secs(&secs1);
        gomp_FormatMessage("   In %.3f hours of wall clock time (%.3f%c)",
                secs1/3600.,(secs/secs1)*100.,'%');
        (void) InformIObytes();
        break; 

    }

    return(0);
}

/***********************************************************************/
int InformIObytes()
/***********************************************************************/
{
#if HAVE_GETRUSAGE

    gomp_PrintMessage("     ****  I / O  I N F O  ****");
    gomp_FormatMessage(
        "     Write: %ld Gbytes %ld Kbytes",
        (long int)gomp_RUsage.ru_oublock,
        (long int)gomp_RUsage.ru_oublock / 1000);
    gomp_FormatMessage(
        "     Read : %ld Gbytes %ld Kbytes",
        (long int)gomp_RUsage.ru_inblock,
        (long int)gomp_RUsage.ru_inblock / 1000);

#elif defined(IRIX)

    gomp_PrintMessage("     ****  I / O  I N F O  ****");
    gomp_FormatMessage(     
        "     Write: %ld Gbytes %ld Kbytes",
        (long int)gomp_ProcessInfo.pu_gbwrit,
        (long int)gomp_ProcessInfo.pu_bwrit / 1000);
    gomp_FormatMessage(
        "     Read : %ld Gbytes %ld Kbytes",
        (long int)gomp_ProcessInfo.pu_gbread,
        (long int)gomp_ProcessInfo.pu_bread / 1000);

#endif

    return(0);
}

/***********************************************************************/
int gomp_SaveModuleName(const char *MName)
/***********************************************************************/
{
    strncpy(ModuleName , MName , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
const char *GetModuleName()
/***********************************************************************/
{
    return(ModuleName);
}
/************************************************************************/
int Get_env()   /* go and hunt for environment parameters */
/************************************************************************/
{
    const char *dir, *root;
    char   Temp[BUFF_LEN];

/* GOM_ROOT */
    root = getenv("GOM_ROOT");
    if ( ! root || ! *root ) {
        /* sorry root directory not defined */
        gomp_PrintEXIT(
            "Path to root directory not defined. Program has to halt!");
    }
    SetGomDirectory("GOM_ROOT",root,1);

    if ( ! *GetEnvironmentFile() ) {
        JoinPath(Temp, root, "environment.txt");
        gomp_SetEnvironmentFile(Temp);
    }

    SetValuesFromEnvironment();

/* check GOM_HOME */
    if ( ! *GetGomDirectory("GOM_HOME") ) {
        dir = getenv("HOME");
        if ( ! dir || !*dir  )
            dir = root;
        SetGomDirectory("GOM_HOME",dir,0);
    }
    
/* check GOM_TEMP */
    if ( ! *GetGomDirectory("GOM_TEMP") )
        SetGomDirectory("GOM_TEMP",GetGomTempDir(),1);

#if defined(JUNK)
/* check TCL_LIBRARY  */
    CTemp = GetValueFromEnvironment(
        "TCL_LIBRARY value is missing!",
        "TCL_LIBRARY",NULL);

    sprintf(Temp,"TCL_LIBRARY=%s",CTemp);
    if(putenv(Temp))
        gomp_PrintEXIT("Can't get memory to exand the environment");

/* check TK_LIBRARY  */
    CTemp = GetValueFromEnvironment(
        "TK_LIBRARY value is missing!",
        "TK_LIBRARY",NULL);

#endif
/* ............... */
#if 0
/* get "term"    */
#if defined(WIN32)
    (void)gomp_SetTerminal("WIN32");
#else
    envp = getenv("TERM");
    if(envp != NULL) {
/*
  sprintf(OutText,"    Terminal:  '%s'",envp);
  gomp_PrintMessage(OutText);
*/
        (void)gomp_SetTerminal(envp);
    }
    else
        gomp_PrintMessage("   ** Unknown terminal type **");
#endif
/* get "display" */
    envp = getenv("DISPLAY");
    if(envp != NULL) {
/*
  sprintf(OutText,"    Display:   '%s'",envp);
  gomp_PrintMessage(OutText);
*/
        (void)gomp_SetDisplay(envp);
    }
    else 
        (void)gomp_SetDisplay("WIN32 Display");
#endif

    return(0);
}
/**************************************************************************/
const char *GetGomTempDir()
/**************************************************************************/
{
    const char *dir;
    
    dir = GetGomDirectory("GOM_TEMP");
    
    if ( ! dir || ! *dir ) {
#ifdef WIN32
        /* GetTempPath honours TEMP and TMP. */
        static char TempDir[MAX_PATH];
        DWORD length = GetTempPath(MAX_PATH, TempDir);
        if ( 0 < length && length < MAX_PATH )
            dir = TempDir;
#else
        dir = getenv("TMPDIR");
        if ( ! dir || ! *dir )
            dir = getenv("TMP");
        if ( ! dir || ! *dir ) {
            if ( access( "/tmp/.", R_OK | W_OK | X_OK | F_OK ) == 0 )
                dir = "/tmp";
        }
#endif

        if ( ! dir || ! *dir )
            gomp_PrintEXIT(
                "Path to temp directory not defined. Program has to halt");
    }
    
    return dir;
}
/***********************************************************************/
const char *GetGomDirectory(const char *name)
/***********************************************************************/
{
    const char *value;
    value = Tcl_GetVar2(gomp_GetTclInterp(),"gomEnv",name,TCL_GLOBAL_ONLY);
    return value ? value : "";
}
/***********************************************************************/
void SetGomDirectory(const char *name, const char *dir, int check)
/***********************************************************************/
{
    char Temp[BUFF_LEN];
    if ( Tcl_GetPathType( dir ) != TCL_PATH_ABSOLUTE ) {
        if ( ! *dir )
            dir = gomp_ShowRootDir();
        else {
            JoinPath(Temp, gomp_ShowRootDir(), dir);
            dir = Temp;
        }
    }
    if ( check && strpbrk(dir, " \t") )
        gomp_FormatEXIT(
            "%s\n%s",
            "gOpenMol is installed in a directory with a name containing "
            "space(s)\n"
            "Please reinstall in a directory with a name without space(s)",
            dir);
    if ( ! Tcl_SetVar2(
        gomp_GetTclInterp(),"gomEnv",name,dir,TCL_GLOBAL_ONLY) )
        gomp_FormatEXIT("can't set 'gomEnv(%s)' environment value",name);
}
/***********************************************************************/
int SetValuesFromEnvironment()
/***********************************************************************/
{
    int argc, code;
    const char *value;
    const char **argv;
    FILE *EnvFile;
    char InputLine[BUFF_LEN];

    if ( ! *Environment.File )
        return(0);

    EnvFile = fopen(Environment.File,"r");
    
    if ( ! EnvFile )
        gomp_FormatEXIT("Can't open environment file: %s", Environment.File);

    while ( fgets(InputLine,BUFF_LEN,EnvFile) != NULL) { 
   
        code = Tcl_SplitList(gomp_GetTclInterp(), InputLine, &argc, &argv);

        if ( code != TCL_OK )
            gomp_FormatEXIT("Can't parse environment file.\n"
                            "File: %s\n"
                            "Line: %s",
                            Environment.File,
                            InputLine);

        if ( argc >= 2 && argv[0][0] != '#' ) {
            /**
             * Environment values take precedence over the values from the
             * environment file.
             */
            value = getenv(argv[0]);
            SetGomDirectory(argv[0], value ? value : argv[1], 1);
        }

        if ( argc )
            Tcl_Free((char *)CONST_CAST(char **, argv));
    }

    fclose(EnvFile);

    return(0);
}
/**************************************************************************/
int gomp_Get_proc_info()       /* get process memory */
/**************************************************************************/
{

    unsigned long pagesize;    /* system page size */
    int pid_num;     /* process pid number */
    int i,str_long;
    const char *Working;
    int t_size;  /* total size of process */
    int r_size;  /* total resident size of process */

#if defined(WIN32)
    SYSTEM_INFO SystemInfo;
#endif

    Working = GetModuleName();

    if(Working[0] == (char)NULL) {
        gomp_PrintWARNING("?Lost the program name");
        return(1);
    }

    str_long = strlen(Working);
    for(i = str_long - 1 ; i >= 0 ; i--) {
        if(Working[i] == '/') break;
    }

#if defined(WIN32)
    GetSystemInfo(&SystemInfo);
    pagesize = SystemInfo.dwPageSize;
    pid_num  = _getpid();
#else
#ifdef _SC_PAGESIZE
    pagesize = sysconf(_SC_PAGESIZE);
#else
    pagesize = getpagesize();
#endif
    pid_num  = getpid();
#endif
    t_size = 0;
    r_size = 0;

#if !defined(WIN32)
    if(gomp_RunStatistics() == -1) {
        gomp_PrintERROR(
            "Can't get information about the resource utilization");
        return(1);
    }
#endif
    gomp_FormatMessage(
        "    Program name:                    %s ",Working+i+1);
    gomp_FormatMessage(
        "    Process pid:                     %d ",pid_num);
    gomp_FormatMessage(
        "    System pagesize:                 %ld \t(bytes)",pagesize);
    {
        int major;
        int minor; 
        int patchLevel; 
/*       Tcl_ReleaseType type;*/ 
        Tcl_GetVersion(&major, &minor, &patchLevel, NULL);
        gomp_FormatMessage(
            "    Tcl/Tk library version:          %d.%d.%d ",
            major,minor,patchLevel);
        
    }

#if HAVE_GETRUSAGE
    t_size = gomp_RUsage.ru_maxrss;
    gomp_FormatMessage(
        "    Total size of process:           %d \t(Kbytes)",t_size);
#elif defined(IRIX)
    t_size = gomp_ProcessInfo.pu_size    * pagesize / 1000;
    r_size = gomp_ProcessInfo.pu_rss     * pagesize / 1000;
    gomp_FormatMessage(
        "    Total size of process:           %d \t(Kbytes)",t_size);
    gomp_FormatMessage(
        "    Total resident size of process:  %d \t(Kbytes)",r_size);
#endif

    return(0);
}
/**************************************************************************/
const char *AtomConversionFile()
/**************************************************************************/
{
    return(FileNames.AtomConversion);
}
/**************************************************************************/
const char *AtomParametersFile()
/**************************************************************************/
{
    return(FileNames.AtomParameters);
}
/**************************************************************************/
const char *ElementParametersFile()
/**************************************************************************/
{
    return(FileNames.ElementParameters);
}
/**************************************************************************/
int gomp_SetAtomParametersFile(const char *FileName)
/**************************************************************************/
{
    gomp_CopyString(FileNames.AtomParameters,FileName,BUFF_LEN);

    return(0);
}

/**************************************************************************/
const char *ColourTableFile()
/**************************************************************************/
{
    return(FileNames.ColourTable);
}
#if 0
/**************************************************************************/
int gomp_SetTerminal(const char *Term)
/**************************************************************************/
{
    gomp_CopyString(Terminal.TERM,Term,BUFF_LEN);

    return(0);
}
/**************************************************************************/
int gomp_SetDisplay(const char *Display)
/**************************************************************************/
{
    gomp_CopyString(Terminal.DISPLAY,Display,BUFF_LEN);

    return(0);
}
/**************************************************************************/
const char *gomp_GetTerminal()
/**************************************************************************/
{
    return(Terminal.TERM);
}
/**************************************************************************/
const char *gomp_GetDisplay()
/**************************************************************************/
{
    return(Terminal.DISPLAY);
}
#endif
/**************************************************************************/
static void PrintBannerLine(const char *line, int width)
/**************************************************************************/
{
    static const char* Spaces =
        "                                                                 ";
    size_t length;

    length = strlen(line);

    /* Pad the line with spaces. */
    /* Print the paddes line. */
    gomp_FormatMessage(
        "    *%.*s%s%.*s*",
        (int)( ( width - length ) / 2 ),
        Spaces,
        line,
        (int)( width - length - ( width - length ) / 2 ),
        Spaces);
}
/**************************************************************************/
int gomp_PrintBANNER()
/**************************************************************************/
{
    static char VersionLine[BUFF_LEN];
    int width = strlen(gOpenMolBanner[0]);
    int i;
    
    for ( i = 0 ; gOpenMolBanner[i] ; i++ ) {
        switch (gOpenMolBanner[i][0]) {
        case '\t':
            /* Print version line. */
            sprintf(VersionLine,"Version %d.%02d  %s  (r. %s)",
                    GOPENMOL_VERSION / 100,GOPENMOL_VERSION % 100,
                    OS_NAME, GOPENMOL_RELEASE);
            PrintBannerLine(VersionLine,width);
            break;
        default:
            PrintBannerLine(gOpenMolBanner[i],width);
        }
    }

    return(0);
}
/**************************************************************************/
int gomp_PrintUSAGE()
/**************************************************************************/
{
    int i;
    
    for( i = 0 ; gOpenMolUsage[i] ; i++ )
        printf("%s\n",gOpenMolUsage[i]);

    return(0);
}
/**************************************************************************/
const char *gomp_LongHostName()
/**************************************************************************/
{
    static const char *this_host;
    static char    This_Host[BUFF_LEN];

    if(!(this_host = Tcl_GetHostName())) {
        strcpy(This_Host,"MyHost");
        return(This_Host);
    }

    return(this_host);
}
#if 0
/***********************************************************************/
int gomp_PrintEnvironment()
/***********************************************************************/
{
    char OutText[BUFF_LEN];

    sprintf(OutText,"    Data path: %s ",gomp_ShowDataDir());
    gomp_PrintMessage(OutText);

    sprintf(OutText,"    Bin path:  %s ",gomp_ShowBinDir());
    gomp_PrintMessage(OutText);


    sprintf(OutText,"    Home path: %s ",gomp_ShowHomeDir());
    gomp_PrintMessage(OutText);

    sprintf(OutText,"    Terminal:  '%s'",gomp_GetTerminal());
    gomp_PrintMessage(OutText);

    sprintf(OutText,"    Display:   '%s'",gomp_GetDisplay());
    gomp_PrintMessage(OutText);

    return(0);
}
#endif
/***********************************************************************/
int gomp_GetFileSize(const char *Name)
/***********************************************************************/
{

#if defined(WIN32)
    static int FileSize;
    static int Temp;
    static FILE *File_p;

    File_p = fopen(Name , "r");
    if(File_p == NULL) {
        gomp_PrintWARNING("Can't get the file size! It's not a file?");
        return(0);
    }
    rewind(File_p);
    Temp = ftell(File_p);
    fseek(File_p,0L,SEEK_END);
    FileSize = ftell(File_p) - Temp;
    fclose(File_p);
    return(FileSize);
#else
    static struct stat FileInfo;

    if(stat(Name, &FileInfo)) return(-1);

    return((int)FileInfo.st_size);
#endif
}
#if 0
/***********************************************************************/
int gomp_PutText2StatusLine1(const char *Text)
/***********************************************************************/
{
    if(gomp_GetTermType() == GRAPHICS_AVAILABLE)
        return(gomp_PutText2Info1(Text));
    else
        return(0);
}
#endif
/***********************************************************************/
int gomp_PutText2StatusLine2(int Value)
/***********************************************************************/
{
    if(gomp_GetTermType() == GRAPHICS_AVAILABLE)
        return(gomp_PutChars2Info2(Value));
    else
        return(0);
}
/***********************************************************************/
char *gomp_Fgets (char *s, int n, FILE *stream)
/***********************************************************************/
{
    char *retv;
    int   i;

    retv = fgets ( s, n, stream);

    i = strlen(s) - 1;

    if(s[i] == '\n') s[i] = '\0' ;

    return(retv);
}

/***********************************************************************/
const char *gomp_GetAMBERcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.AMBERfileType);
}

/***********************************************************************/
int   gomp_PutAMBERcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.AMBERfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetAMBERdefault()
/***********************************************************************/
{
    return(CoordFileType.AMBERdefault);
}
/***********************************************************************/
int   gomp_PutAMBERdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.AMBERdefault = 1;

    return(0);
}
/***********************************************************************/
const char *gomp_GetCHARMMcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.CHARMMfileType);
}

/***********************************************************************/
int   gomp_PutCHARMMcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.CHARMMfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetCHARMMdefault()
/***********************************************************************/
{
    return(CoordFileType.CHARMMdefault);
}
/***********************************************************************/
int   gomp_PutCHARMMdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.CHARMMdefault = 1;

    return(0);
}
/***********************************************************************/
const char *gomp_GetGAMESScoordFileType()
/***********************************************************************/
{
    return(CoordFileType.GAMESSfileType);
}

/***********************************************************************/
int   gomp_PutGAMESScoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.GAMESSfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetGAMESSdefault()
/***********************************************************************/
{
    return(CoordFileType.GAMESSdefault);
}
/***********************************************************************/
int   gomp_PutGAMESSdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.GAMESSdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetGAUSSIANcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.GAUSSIANfileType);
}

/***********************************************************************/
int   gomp_PutGAUSSIANcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.GAUSSIANfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetGAUSSIANdefault()
/***********************************************************************/
{
    return(CoordFileType.GAUSSIANdefault);
}
/***********************************************************************/
int   gomp_PutGAUSSIANdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.GAUSSIANdefault = 1;

    return(0);
}
/***********************************************************************/
const char *gomp_GetGROMACScoordFileType()
/***********************************************************************/
{
    return(CoordFileType.GROMACSfileType);
}

/***********************************************************************/
int   gomp_PutGROMACScoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.GROMACSfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetGROMACSdefault()
/***********************************************************************/
{
    return(CoordFileType.GROMACSdefault);
}
/***********************************************************************/
int   gomp_PutGROMACSdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.GROMACSdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetHYPERCHEMcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.HYPERCHEMfileType);
}

/***********************************************************************/
int   gomp_PutHYPERCHEMcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.HYPERCHEMfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetHYPERCHEMdefault()
/***********************************************************************/
{
    return(CoordFileType.HYPERCHEMdefault);
}
/***********************************************************************/
int   gomp_PutHYPERCHEMdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.HYPERCHEMdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetINSIGHTcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.INSIGHTfileType);
}

/***********************************************************************/
int   gomp_PutINSIGHTcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.INSIGHTfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetINSIGHTdefault()
/***********************************************************************/
{
    return(CoordFileType.INSIGHTdefault);
}
/***********************************************************************/
int   gomp_PutINSIGHTdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.INSIGHTdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetJAGUARcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.JAGUARfileType);
}

/***********************************************************************/
int   gomp_PutJAGUARcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.JAGUARfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetJAGUARdefault()
/***********************************************************************/
{
    return(CoordFileType.JAGUARdefault);
}
/***********************************************************************/
int   gomp_PutJAGUARdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.JAGUARdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetMOL2coordFileType()
/***********************************************************************/
{
    return(CoordFileType.MOL2fileType);
}

/***********************************************************************/
int   gomp_PutMOL2coordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.MOL2fileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetMOL2default()
/***********************************************************************/
{
    return(CoordFileType.MOL2default);
}
/***********************************************************************/
int   gomp_PutMOL2default()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.MOL2default = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetMOPACgraphcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.MOPACgraphfileType);
}

/***********************************************************************/
int   gomp_PutMOPACgraphcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.MOPACgraphfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetMOPACgraphdefault()
/***********************************************************************/
{
    return(CoordFileType.MOPACgraphdefault);
}
/***********************************************************************/
int   gomp_PutMOPACgraphdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.MOPACgraphdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetMUMODcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.MUMODfileType);
}

/***********************************************************************/
int   gomp_PutMUMODcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.MUMODfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetMUMODdefault()
/***********************************************************************/
{
    return(CoordFileType.MUMODdefault);
}
/***********************************************************************/
int   gomp_PutMUMODdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.MUMODdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetOPENMOLcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.OPENMOLfileType);
}

/***********************************************************************/
int   gomp_PutOPENMOLcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.OPENMOLfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetOPENMOLdefault()
/***********************************************************************/
{
    return(CoordFileType.OPENMOLdefault);
}
/***********************************************************************/
int   gomp_PutOPENMOLdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.OPENMOLdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetPDBcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.PDBfileType);
}

/***********************************************************************/
int   gomp_PutPDBcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.PDBfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetPDBdefault()
/***********************************************************************/
{
    return(CoordFileType.PDBdefault);
}
/***********************************************************************/
int   gomp_PutPDBdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.PDBdefault = 1;

    return(0);
}
/***********************************************************************/
const char *gomp_GetTINKERcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.TINKERfileType);
}

/***********************************************************************/
int   gomp_PutTINKERcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.TINKERfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetTINKERdefault()
/***********************************************************************/
{
    return(CoordFileType.TINKERdefault);
}
/***********************************************************************/
int   gomp_PutTINKERdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.TINKERdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetUHBDcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.UHBDfileType);
}

/***********************************************************************/
int   gomp_PutUHBDcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.UHBDfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetUHBDdefault()
/***********************************************************************/
{
    return(CoordFileType.UHBDdefault);
}
/***********************************************************************/
int   gomp_PutUHBDdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.UHBDdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetXMOLcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.XMOLfileType);
}

/***********************************************************************/
int   gomp_PutXMOLcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.XMOLfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetXMOLdefault()
/***********************************************************************/
{
    return(CoordFileType.XMOLdefault);
}
/***********************************************************************/
int   gomp_PutXMOLdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.XMOLdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetXYZcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.XYZfileType);
}

/***********************************************************************/
int   gomp_PutXYZcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.XYZfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetXYZdefault()
/***********************************************************************/
{
    return(CoordFileType.XYZdefault);
}
/***********************************************************************/
int   gomp_PutXYZdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.XYZdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetYASPcoordFileType()
/***********************************************************************/
{
    return(CoordFileType.YASPfileType);
}

/***********************************************************************/
int   gomp_PutYASPcoordFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.YASPfileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetYASPdefault()
/***********************************************************************/
{
    return(CoordFileType.YASPdefault);
}
/***********************************************************************/
int   gomp_PutYASPdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.YASPdefault = 1;

    return(0);
}

/***********************************************************************/
const char *gomp_GetOPENMOLcenterFileType()
/***********************************************************************/
{
    return(CoordFileType.OPENMOL_center_fileType);
}

/***********************************************************************/
int   gomp_PutOPENMOLcenterFileType(const char *FileType)
/***********************************************************************/
{
    strncpy(CoordFileType.OPENMOL_center_fileType , FileType , BUFF_LEN-1);

    return(0);
}
/***********************************************************************/
int   gomp_GetOPENMOLcenterdefault()
/***********************************************************************/
{
    return(CoordFileType.OPENMOL_center_default);
}
/***********************************************************************/
int   gomp_PutOPENMOLcenterdefault()
/***********************************************************************/
{
    (void)gomp_ResetAllDefaultFileTypes();

    CoordFileType.OPENMOL_center_default = 1;

    return(0);
}

/***********************************************************************/
int   gomp_ResetAllDefaultFileTypes()
/***********************************************************************/
{
    CoordFileType.AMBERdefault           = 0;
    CoordFileType.CHARMMdefault          = 0;
    CoordFileType.GAUSSIANdefault        = 0;
    CoordFileType.HYPERCHEMdefault       = 0;
    CoordFileType.INSIGHTdefault         = 0;
    CoordFileType.MOL2default            = 0;
    CoordFileType.MOPACgraphdefault      = 0;
    CoordFileType.MUMODdefault           = 0;
    CoordFileType.OPENMOLdefault         = 0;
    CoordFileType.PDBdefault             = 0;
    CoordFileType.XMOLdefault            = 0;
    CoordFileType.XYZdefault             = 0;
    CoordFileType.YASPdefault            = 0;
    CoordFileType.OPENMOL_center_default = 0;

    return(0);
}
/***********************************************************************/
int gomp_SetLogFileStatus(int StatusValue)
/***********************************************************************/
{
    LogFile.ok = StatusValue;

    return(0);
}
/***********************************************************************/
int gomp_GetLogFileStatus()
/***********************************************************************/
{
    return(LogFile.ok);
}

/***********************************************************************/
void   gomp_Exit(int StatusValue)
/***********************************************************************/
{
    gomp_CloseLogFile();

    Tcl_Exit(StatusValue);
}
/***********************************************************************/
int gomp_SetEnvironmentFile(const char *FileName)
/***********************************************************************/
{
    gomp_CopyString(Environment.File,FileName,BUFF_LEN);

    return(0);
}
/***********************************************************************/
const char *GetEnvironmentFile()
/***********************************************************************/
{
    return(Environment.File);
}
/***********************************************************************/
int JoinPath(char *dst, const char *prefix, const char *file)
/***********************************************************************/
{
    gomp_CopyString(dst, prefix, BUFF_LEN);

    if ( *file ) {
        const char *const end = dst + BUFF_LEN - 1;
        char *p;
        p = dst + strlen(dst);
        
        /* Copy the file name and convert all directory separators. */
        if ( p != end ) {
#ifdef WIN32
            *p++ = '\\';
#else
            *p++ = '/';
#endif
        }
        while ( *file && p != end ) {
            switch ( *file ) {
#ifdef WIN32
            case '/': *p++ = '\\'; break;
#else
            case '\\': *p++ = '/'; break;
#endif
            default: *p++ = *file;
            }
            file++;
        }
        
        *p = '\0';
    }

    return(0);
}
/***********************************************************************/
int gomp_SplitFile(char *dst,const char *file)
/***********************************************************************/
{
    int argc;
    const char **argv;

    Tcl_SplitPath(file, &argc, &argv);

    if(!argc)
        strcpy(dst , "Unknown");
    else
        gomp_CopyString(dst, argv[argc - 1], BUFF_LEN);

    Tcl_Free((char *)CONST_CAST(char **, argv));

    return(0);
}
/*
  alt #1: atom name
  #2: residue name
  #3: segment name
*/
/*************************************************************************/
int gomp_UpdateNameStack(int alt, const char *temp)
/*************************************************************************/
{

    int i;

    switch(alt) {

    case 1:
        if(atnam_stack_deep) {
            for(i = 0 ; i < atnam_stack_deep ; i++ ) {
                if(!strncmp(atnam_stack+MAX_ATM_NAME_LEN*i,temp,MAX_ATM_NAME_LEN)) {
                    atnam_stack_num[i]++; /* increment counter */
                    return 0; /* yes it is in the stack */
                }
            }
        } else {
            atnam_stack     = gomp_AllocateCharVector(atnam_stack_max * MAX_ATM_NAME_LEN);
            atnam_stack_num = gomp_AllocateIntVector(atnam_stack_max);
        }
        /* no it is not so put it there */
        if(atnam_stack_deep == atnam_stack_max) {
            atnam_stack_max = atnam_stack_max + 100;
            atnam_stack =     realloc(atnam_stack, 
                                               atnam_stack_max * MAX_ATM_NAME_LEN);
            atnam_stack_num = realloc(atnam_stack_num, 
                                               atnam_stack_max * sizeof(*atnam_stack_num));
        }
        strncpy(atnam_stack+MAX_ATM_NAME_LEN*atnam_stack_deep,temp,
                MAX_ATM_NAME_LEN);
        atnam_stack_num[atnam_stack_deep] = 1;
        atnam_stack_deep++;
        break;

    case 2:
        if(resnam_stack_deep) {
            for(i = 0 ; i < resnam_stack_deep ; i++ ) {
                if(!strncmp(resnam_stack+MAX_RES_NAME_LEN*i,temp,MAX_RES_NAME_LEN)) {
                    resnam_stack_num[i]++;
                    return 0; /* yes it is in the stack */
                }
            }
        } else {
            resnam_stack     = gomp_AllocateCharVector(resnam_stack_max * MAX_RES_NAME_LEN);
            resnam_stack_num = gomp_AllocateIntVector(resnam_stack_max);
        }
        /* no it is not so put it there */
        if(resnam_stack_deep == resnam_stack_max) {
            resnam_stack_max = resnam_stack_max + 50;
            resnam_stack =     realloc(resnam_stack, 
                                                resnam_stack_max * MAX_RES_NAME_LEN);
            resnam_stack_num = realloc(resnam_stack_num, 
                                                resnam_stack_max * sizeof(*resnam_stack_num));
        }
        strncpy(resnam_stack+MAX_RES_NAME_LEN*resnam_stack_deep,temp,
                MAX_RES_NAME_LEN);
        resnam_stack_num[resnam_stack_deep] = 1;
        resnam_stack_deep++;
        break;

    case 3:
        if(segment_stack_deep) {
            for(i = 0 ; i < segment_stack_deep ; i++ ) {
                if(!strncmp(segment_stack+MAX_SEG_NAME_LEN*i,temp,MAX_SEG_NAME_LEN)) {
                    segment_stack_num[i]++;
                    return 0; /* yes it is in the stack */
                }
            }
        } else {
            segment_stack     = gomp_AllocateCharVector(segment_stack_max * MAX_RES_NAME_LEN);
            segment_stack_num = gomp_AllocateIntVector(segment_stack_max);
        }
        /* no it is not so put it there */
        if(segment_stack_deep == segment_stack_max) {
            segment_stack_max = segment_stack_max + 30;
            segment_stack =     realloc(segment_stack, 
                                                 segment_stack_max * MAX_SEG_NAME_LEN);
            segment_stack_num = realloc(segment_stack_num, 
                                                segment_stack_max * sizeof(*segment_stack_num));
        }
        strncpy(segment_stack+MAX_SEG_NAME_LEN*segment_stack_deep,temp,
                MAX_SEG_NAME_LEN);
        segment_stack_num[segment_stack_deep] = 1;
        segment_stack_deep++;
        break;
    }

    return 1;
}
/*************************************************************************/
int gomp_ResetNameStack()
/*************************************************************************/
{

    if(atnam_stack_deep) {
        free(atnam_stack);    /* stack to contain the different atnam, resnam */
        free(atnam_stack_num);
        atnam_stack_deep = 0;
    }
    if(resnam_stack_deep) {
        free(resnam_stack);
        free(resnam_stack_num);
        resnam_stack_deep = 0;
    }
    if(segment_stack_deep) {
        free(segment_stack);
        free(segment_stack_num);
        segment_stack_deep = 0;
    }

    return 0;
}
/*************************************************************************/
const char *gomp_ShowAtomNameStack()
/*************************************************************************/
{
    static char *AtomNameList = (char *)NULL;
    char Text[BUFF_LEN];
    int  i;

    if(AtomNameList) {
        free(AtomNameList);
    }

    if(!atnam_stack_deep) {

        sprintf(Text,"{ } { }");
        AtomNameList = gomp_AllocateCharVector(strlen(Text) + 1);
        strncpy(AtomNameList,Text,strlen(Text));
        AtomNameList[strlen(Text)] = (char)NULL;

    } else {
        for(i = 0 ; i < atnam_stack_deep ; i++) {
            sprintf(Text,"{%.4s %d} ",atnam_stack+i*MAX_ATM_NAME_LEN,atnam_stack_num[i]);
            if(i) {
                AtomNameList = gomp_ReallocateCharVector(AtomNameList , strlen(AtomNameList) + strlen(Text) + 1);
                strncat(AtomNameList,Text,strlen(Text));
                AtomNameList[strlen(AtomNameList)] = (char) NULL;
            } else {
                AtomNameList = gomp_AllocateCharVector(strlen(Text) + 1);
                strncpy(AtomNameList,Text,strlen(Text));
                AtomNameList[strlen(Text)] = (char)NULL;
            }
        }
    }

    return(AtomNameList);
}
/*
  Show the numer of atoms 
  1) with the atom labels
  2) under the different residue name tags
  3) under the different segment name tags
*/
/*************************************************************************/
const char *gomp_ShowResidueNameStack()
/*************************************************************************/
{
    static char *ResidueNameList = (char *)NULL;
    char Text[BUFF_LEN];
    int  i;

    if(ResidueNameList) {
        free(ResidueNameList);
    }

    if(!resnam_stack_deep) {

        sprintf(Text,"{ } { }");
        ResidueNameList = gomp_AllocateCharVector(strlen(Text) + 1);
        strncpy(ResidueNameList,Text,strlen(Text));
        ResidueNameList[strlen(Text)] = (char)NULL;

    } else {
        for(i = 0 ; i < resnam_stack_deep ; i++) {
            sprintf(Text,"{%.4s %i} ",resnam_stack+i*MAX_RES_NAME_LEN,resnam_stack_num[i]);
            if(i) {
                ResidueNameList = gomp_ReallocateCharVector(ResidueNameList , strlen(ResidueNameList) + strlen(Text) + 1);
                strncat(ResidueNameList,Text,strlen(Text));
                ResidueNameList[strlen(ResidueNameList)] = (char) NULL;
            } else {
                ResidueNameList = gomp_AllocateCharVector(strlen(Text) + 1);
                strncpy(ResidueNameList,Text,strlen(Text));
                ResidueNameList[strlen(Text)] = (char)NULL;
            }
        }
    }

    return(ResidueNameList);
}
/*************************************************************************/
const char *gomp_ShowSegmentNameStack()
/*************************************************************************/
{
    static char *SegmentNameList = (char *)NULL;
    char Text[BUFF_LEN];
    int  i;

    if(SegmentNameList) {
        free(SegmentNameList);
    }

    if(!segment_stack_deep) {

        sprintf(Text,"{ } { }");
        SegmentNameList = gomp_AllocateCharVector(strlen(Text) + 1);
        strncpy(SegmentNameList,Text,strlen(Text));
        SegmentNameList[strlen(Text)] = (char)NULL;

    } else {
        for(i = 0 ; i < segment_stack_deep ; i++) {
            sprintf(Text,"{%.4s %i} ",segment_stack+i*MAX_SEG_NAME_LEN,segment_stack_num[i]);
            if(i) {
                SegmentNameList = gomp_ReallocateCharVector(SegmentNameList , strlen(SegmentNameList) + strlen(Text) + 1);
                strncat(SegmentNameList,Text,strlen(Text));
                SegmentNameList[strlen(SegmentNameList)] = (char) NULL;
            } else {
                SegmentNameList = gomp_AllocateCharVector(strlen(Text) + 1);
                strncpy(SegmentNameList,Text,strlen(Text));
                SegmentNameList[strlen(Text)] = (char)NULL;
            }
        }
    }

    return(SegmentNameList);
}

static FILE *chkdsk;

/*************************************************************************/
int gomp_OpenPipe2Program(const char *Program , const char *Mode)
/*************************************************************************/
{
    if( chkdsk )
        pclose(chkdsk);

    chkdsk = popen( Program, Mode );
    if ( ! chkdsk )
        return( 1 );

    return(0);

}
/*************************************************************************/
int gomp_ClosePipe2Program()
/*************************************************************************/
{
    /* Close pipe and print return value of CHKDSK. */
    if(chkdsk == (FILE *)NULL) return(0);
#if defined(WIN32)
    _pclose( chkdsk );
#else
    pclose( chkdsk );
#endif
    chkdsk = (FILE *)NULL;

    return(0);
}
/*************************************************************************/
int gomp_SendPipe2Program(const char *Command)
/*************************************************************************/
{

    if(chkdsk == (FILE *)NULL) {
        gomp_PrintERROR("no pipe open to write into");
        return(1);
    }

    if(fprintf(chkdsk , "%s\n", Command) < 0) return(1);

    return(0);
}
/*   Utility function to check last characters in a line for a text file
     UNIX end a text file line with     \n
     Windows end a text file line with  \r\n
*/
/***********************************************************************/
int gomp_CheckTextFileLineEnding(const char *Filename)
/***********************************************************************/
{
    FILE *fp;

    fp = fopen(Filename,"rb");
    if(fp == NULL) {
        gomp_FormatERROR("Can't open input file: %s", Filename);
        return(1);
    }

    for ( ;; ) {
        switch ( getc(fp) ) {
        case '\n':
            /* UNIX line ending. */
            fclose(fp);
#ifdef WIN32
            return 1;
#else
            return 0;
#endif  
        case '\r':
#ifdef WIN32
            if ( getc(fp) == '\n' ) {
                /* Windows line ending. */
                fclose(fp);
                return 0;
            }
#endif
            /* Windows or Macintosh line ending. */
            fclose(fp);
            return 1;
        case EOF:
            fclose(fp);
            gomp_PrintERROR("Could not decide on the text file correctness");
            return 1;
        }
    }
}
