/*

Copyright (c) 2003 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Coded by: Eero Häkkinen
*/

#include "maindefs.h"

#include <limits.h>
#include <string.h>
#include <tcl.h>

#include "colouring.h"
#include "gomtcl.h"
#include "memalloc.h"
#include "printmsg.h"
#include "parser.h"
#include "tclutils.h"

#include "stdafx.h"

/*********************************************************************/
static Tcl_Obj *CreateTclList(const char **format,va_list *args)
/*********************************************************************/
{
    Tcl_Obj    *list = Tcl_NewListObj(0,NULL);
    const char *text;
    char        ch;

    for ( ;; ) {
        switch ( **format ) {
        case '\0':
            /* end of list */
            return list;
        case ')':
            /* end of sub list */
            ++*format;
            return list;
        case '(':
            /* sub list */
            ++*format;
            Tcl_ListObjAppendElement(NULL,list,
                                     CreateTclList(format,args));
            continue;
        case ' ':
        case '\t':
            /* Visual groupping. */
            break;
        case '%':
            ++*format;
            switch ( **format ) {
            case 'p':
                /* object */
                Tcl_ListObjAppendElement(
                    NULL,list,va_arg(*args,Tcl_Obj*));
                break;
            case 's':
                /* string */
                text = va_arg(*args,const char *);
                Tcl_ListObjAppendElement(NULL,list,Tcl_NewStringObj(text,-1));
                break;
            case 'c':
                /* character as int */
                ch = (char)va_arg(*args,int);
                Tcl_ListObjAppendElement(
                    NULL,list,Tcl_NewStringObj(&ch,1));
                break;
            case 'd':
            case 'i':
                /* int */
                if ( (*format)[1] == 'b' ) {
                    ++*format;
                    /* boolean */
                    Tcl_ListObjAppendElement(
                        NULL,list,Tcl_NewBooleanObj(va_arg(*args,int)));
                }
                else
                    Tcl_ListObjAppendElement(
                        NULL,list,Tcl_NewIntObj(va_arg(*args,int)));
                break;
            case 'l':
                ++*format;
                switch ( **format ) {
                case 'd':
                case 'i':
                    /* long int */
                    Tcl_ListObjAppendElement(
                        NULL,list,Tcl_NewLongObj(va_arg(*args,long int)));
                    break;
                default:
                    gomp_FormatEXIT(
                        "invalid format specification '%10s' "
                        "for gomp_CreateTclList",*format-2);
                }
                break;
            case 'f':
                /* double  */
                Tcl_ListObjAppendElement(
                    NULL,list,Tcl_NewDoubleObj(va_arg(*args,double)));
                break;
            default:
                gomp_FormatEXIT(
                    "invalid format specification '%10s' "
                    "for gomp_CreateTclList",*format-1);
            }
            break;
        default:
            gomp_FormatEXIT(
                "invalid format specification '%10s' "
                "for gomp_CreateTclList",*format);
        }
        ++*format;
    }           
}
/*********************************************************************/
Tcl_Obj *gomp_CreateTclList(const char *format, ...)
/*********************************************************************/
{
    Tcl_Obj *list;
    va_list args;
    va_start(args,format);
    list = CreateTclList(&format,&args);
    va_end(args);
    return list;
}
/*********************************************************************/
Tcl_Obj *gomp_ParserCreateList(const char *format, ...)
/*********************************************************************/
{
    Tcl_Obj *list;
    va_list args;
    va_start(args,format);
    list = CreateTclList(&format,&args);
    va_end(args);
    return list;
}

int gomp_ParserAddArgDescriptionToErrorMsg(Tcl_Interp *interp,
                                           const char *desc)
{
    if ( interp ) {
        /* Change the error message to be */
        /* <desc>: <old error message>    */
        Tcl_Obj *err = Tcl_NewStringObj(desc,-1);
        Tcl_AppendToObj(err,": ",2);
        Tcl_AppendObjToObj(err,Tcl_GetObjResult(interp));
        Tcl_SetObjResult(interp,err);
    }
    return TCL_ERROR;
}

int gomp_ParserParseColour(Tcl_Interp *interp,Tcl_Obj *obj,
                           gom_FloatColour *fc)
{
    const char *colour = Tcl_GetStringFromObj(obj,NULL);
    if ( gomp_ColourName2RGBSilent(colour,&fc->red,&fc->green,&fc->blue) == 0 )
        return 1;
    else {
        gomp_ParserFormatReturn(
            interp,"unknown colour \"%s\"",colour);
        return 0;
    }
}

int gomp_ParserCheckDoubleRange(Tcl_Interp *interp,double value,
                                double min,double max)
{
    if ( min <= value && value <= max )
        return 1;
    else {
        gomp_ParserFormatReturn(
            interp,"value %f out of allowed range: %f ... %f",
            value,min,max);
        return 0;
    }
}

int gomp_ParserCheckIntRange(Tcl_Interp *interp,int value,
                             int min,int max)
{
    if ( min <= value && value <= max )
        return 1;
    else {
        gomp_ParserFormatReturn(
            interp,"value %d out of allowed range: %d ... %d",
            value,min,max);
        return 0;
    }
}

int gomp_ParserParseIntEnum(Tcl_Interp *interp,Tcl_Obj *obj,
                            const gom_ParserIntEnum entries[],
                            int *valuePtr)
{
    int index;
    if ( Tcl_GetIndexFromObjStruct(
             interp,obj,entries,sizeof(*entries),"option",0,&index) != TCL_OK )
        return 0;
    *valuePtr = entries[index].value;
    return 1;
}

static void ParseSyntaxFromArgumentLists(
    Tcl_DString *syntaxes,
    Tcl_Interp *interp,int objc,Tcl_Obj *const objv[],
    const char *msg,
    const gom_ParserArgumentList *list)
{
    char        syntax[BUFF_LEN];
    int         FirstRound = 1;
    const char *command;

    if ( interp )
        Tcl_DStringAppend(syntaxes,"wrong # args: should be one of\n",-1);

    for ( ; list->action.list || list->action.fnctn ; list++ ) {

        command = list->command;
        /* Handle all command aliases. */
        do {
            if ( ! interp && objc > 0 ) {
                /* Filter commands. */
                const char *name = Tcl_GetStringFromObj(objv[0],NULL);
                if ( ! command )
                    /* There is no command, so command doesn't match. */
                    break;
                if ( strncmp(command,name,strlen(name)) != 0 ) {
                    /* Continue to the next command alias. */
                    command += strlen(command) + 1;
                    continue;
                }
            }
            if ( command ||
                 ( list->args.msg && *list->args.msg) ) {
                /* We need at least one argument more.            */
                /* All previous arguments are therefor mandatory. */
                /* Remove question marks.                         */
                const char *m = msg;
                char       *s = syntax;
                for ( ; *m ; m++ ) {
                    if ( *m != '?' )
                        *s++ = *m;
                }
                *s = '\0';
            }
            else
                strcpy(syntax,msg);

            if ( command ) {
                /* There is a command name. */
                if ( *syntax )
                    strcat(syntax," ");
                strcat(syntax,command);
                /* Go to the next alias. */
                command += strlen(command) + 1;
            }
            if ( list->args.msg && *list->args.msg ) {
                /* There are arguments. */
                if ( *syntax )
                    strcat(syntax," ");
                strcat(syntax,list->args.msg);
            }
            if ( interp ) {
                const char *r;
                if ( list->action.list ) {
                    /* This isn't the final entry. */
                    if ( *syntax )
                        strcat(syntax," ...");
                    else
                        strcat(syntax,"...");
                }
                Tcl_ResetResult(interp);
                Tcl_WrongNumArgs(interp,objc,objv,*syntax ? syntax : NULL);
                r = Tcl_GetStringResult(interp);
                syntax[0] = '\t';
                while ( *r != '\"' )
                    r++;
                strcpy(syntax + 1, r);
                if ( ! FirstRound )
                    Tcl_DStringAppend(syntaxes,"\n",1);
                Tcl_DStringAppend(syntaxes,syntax,-1);
            }
            else {
                if ( list->action.list ) {
                    /* This isn't the final entry. */
                    if ( ! FirstRound )
                        Tcl_DStringAppend(syntaxes,"\n",1);
                    ParseSyntaxFromArgumentLists(
                        syntaxes,interp,
                        objc - (interp ? 0 : 1),
                        objv + (interp ? 0 : 1),
                        syntax,list->action.list);
                }
                else if ( objc <= 1 ) {
                    /* This is the final entry. */
                    if ( ! FirstRound )
                        Tcl_DStringAppend(syntaxes,"\n",1);
                    Tcl_DStringAppend(syntaxes,syntax,-1);
                }
                else
                    /* Command doesn't match. No sub command available. */
                    continue;
            }
            FirstRound = 0;
        } while ( command && *command );
    }
}

static void ReportCorrectSyntax(
    Tcl_Interp *interp,int objc,Tcl_Obj *const objv[],
    const char *msg,
    const gom_ParserArgumentList *list)
{
    Tcl_DString syntaxes;
    Tcl_DStringInit(&syntaxes);
    /* Get syntax. */
    ParseSyntaxFromArgumentLists(&syntaxes,interp,objc,objv,msg,list);
    /* Set error info. */
    Tcl_ResetResult(interp);
    Tcl_DStringResult(interp,&syntaxes);
    Tcl_DStringFree(&syntaxes);
}

int gomp_ParserCommand(ClientData clientData,Tcl_Interp *interp,
                       int objc,Tcl_Obj *const objv[])
{
    const gom_ParserArgumentList *list =
        (const gom_ParserArgumentList *)clientData;
    const gom_ParserArgumentList *list_start;
    int       arg_entry_count;
    int       new_objc  = 0;
    Tcl_Obj **new_objv  = NULL;
    int       obj_index = 1; /* Skip the main command name. */
    intptr_t  Flags     = 0;
    int       i;

    list_start      = list;
    arg_entry_count = 0;

    while ( list->action.list || list->action.fnctn ) {
        if ( list->command ) {
            /* Find a matching command name. */
            int index;
            const char *command;
            const char *alias;

            if ( obj_index >= objc )
                /* Command name is missing. */
                Tcl_WrongNumArgs(interp,obj_index,objv,"cmd ?arg? ...");
            else if ( Tcl_GetIndexFromObjStruct(
                     interp,objv[obj_index],list,sizeof(*list),
                     "command",0,&index) != TCL_OK ) {
                Tcl_Obj *result = Tcl_GetObjResult(interp);
                /* Check if there is a command alias. */
                command = Tcl_GetStringFromObj(objv[obj_index],NULL);

                while ( list->command ) {
                    alias = list->command + strlen(list->command) + 1;
                    while ( *alias ) {
                        if ( strncmp(command, alias, strlen(alias) ) == 0 ) {
                            /* We found a command alias entry. */
                            Tcl_ResetResult(interp);
                            index = 0;
                            goto cmd_entry_found;
                        }
                        if ( ! arg_entry_count )
                            /* Append command alias to the error message.
                             * Error message will be regenerated if there are
                             * argument entries.
                             */
                            Tcl_AppendStringsToObj(
                                result,", or ",alias,(char *)NULL);
                        alias += strlen(alias) + 1;
                    }
                    list++;
                }
            }
            else
                /* We found the entry. */
                goto cmd_entry_found;

            /* Check if we can continue to the plain argument entries. */
            while ( list->command )
                list++;
            if ( list->action.list || list->action.fnctn ) {
                /* We can continue. Clear error message. */
                Tcl_ResetResult(interp);
                continue;
            }
            else {
                /* There are no more entries left to be checked.     */
                /* Make list->command to be non-NULL.                */
                /* That tells we already have an error message.      */
                list--;
                break;
            }

          cmd_entry_found:
            /* We found the entry. */
            obj_index++;
            /* Select matching command entry. */
            list = &list[index];
        }
        else
            arg_entry_count++;

        /* Check argument count. */
        if ( ( objc - obj_index < list->args.count.min ) ||
             ( objc - obj_index > list->args.count.max &&
               list->action.fnctn ) ) {
            /* We don't have enough arguments or                    */
            /* we have too many arguments.                          */
            /* Argument count may exceed the maximum count if this  */
            /* is not the final entry (list->action.fnctn is NULL). */
            if ( list->command ) {
                /* Command name did match so argument count should match. */
                /* We can't continue.                                     */
                if ( list->action.list )
                    ReportCorrectSyntax(interp,obj_index,objv,
                        list->args.msg,list->action.list);
                else
                    Tcl_WrongNumArgs(interp,obj_index,objv,list->args.msg);
                break;
            }
            /* Continue to the next entry to   */
            /* use different set of arguments. */
            list++;
            continue;
        }

        Flags &= list->action.Flags.inherit;
        Flags ^= list->action.Flags.toggle;

        if ( ! new_objv ) {
            if ( list->action.fnctn )
                /* This is the last argument list part. We don't have
                 * previous arguments and all the rest should be passed
                 * to the function. We don't have to reorganize the list.
                 */
                return list->action.fnctn(
                    interp,objc-obj_index,&objv[obj_index],Flags);
            if ( list->args.count.max > 0 && obj_index < objc ) {
                /* We have to reorganize the arguments (well it is
                 * possible that the next entries are ARGS entries but
                 * we don't know that yet).
                 */
                new_objv = gomp_AllocateVoidVector(
                    ( objc - obj_index ) * sizeof(*new_objv) );
                if ( ! new_objv ) {
                    gomp_PrintERROR(NULL);
                    return TCL_ERROR;
                }
            }   
        }
        
        /* Append mandatory arguments to the array. */
        for ( i = 0 ; i < list->args.count.min ; i++ )
            new_objv[new_objc++] = objv[obj_index++];
        /* Append optional arguments to the array. */
        for ( ; i < list->args.count.max && obj_index < objc ; i++ )
            new_objv[new_objc++] = objv[obj_index++];

        if ( list->action.fnctn ) {
            /* This is the last argument list part. */
            /* All arguments are processed.         */
            int code;
            code = list->action.fnctn(interp,new_objc,new_objv,Flags);
            gomp_FreeVector(new_objv);
            return code;
        }
        else {
            /* Proceed to the next part. */
            list_start = list = list->action.list;
            arg_entry_count   = 0;
        }
    }

    gomp_FreeVector(new_objv);

    if ( ! ( list->command && arg_entry_count == 0 ) )
        /* Create complate syntax message. */
        ReportCorrectSyntax(interp,obj_index,objv,"",list_start);
    /* else */
    /*     Use previous error message. */

    gomp_PrintERROR(NULL);
    return TCL_ERROR;
}

gom_ParserIntEnum gomp_ParserBooleanEnumTable[] = {
    {"false",0},
    {"true" ,1},
    {"no"   ,0},
    {"yes"  ,1},
    {"off"  ,0},
    {"on"   ,1},
    GOM_PARSER_END_INT_ENUM_LIST
};

static int GetCommandSyntax(GOM_PARSER_ARGLIST,ptrdiff_t Flags)
{
    Tcl_DString syntaxes;
    Tcl_CmdInfo info;
    const char *Cmd;

    GOM_PARSER_START_ARGS
        GOM_PARSER_ARG(Cmd)
        GOM_PARSER_ARG(SubCmds)
    GOM_PARSER_END_ARGS;
    GOM_PARSER_VERIFY_GOM_PARSER(
        "command name",
        GOM_PARSER_RETRIEVE_STRING(Cmd) );
    
    
    if ( ! Tcl_GetCommandInfo(interp,Cmd,&info) ) {
        Tcl_DStringFree(&syntaxes);
        gomp_FormatERROR("not a command name: '%s'",Cmd);
        return TCL_ERROR;
    }
    if ( ! info.isNativeObjectProc ||
         info.objProc != gomp_ParserCommand ) {
        Tcl_DStringFree(&syntaxes);
        gomp_FormatERROR("not a gOpenMol command: '%s'",Cmd);
        return TCL_ERROR;
    }

    Tcl_DStringInit(&syntaxes);
    ParseSyntaxFromArgumentLists(
        &syntaxes, NULL, objc - 1, objv + 1, Cmd,
        (const gom_ParserArgumentList *)info.objClientData);
    Tcl_DStringResult(interp,&syntaxes);
    Tcl_DStringFree(&syntaxes);

    GOM_PARSER_SUCCEEDED;
}

const gom_ParserArgumentList gomp_SyntaxCommand[] = {
    GOM_PARSER_FINAL_ARGS(GOM_PARSER_NEED_ARGS_RANGE(
        "Cmd ?SubCmd ...?",1,INT_MAX),
                          GetCommandSyntax,
                          GOM_PARSER_UNUSED_VALUE),
    GOM_PARSER_END_ARGUMENT_LIST
};
