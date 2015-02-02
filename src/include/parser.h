#ifndef INC_GOPENMOL_GOM_PARSER
#define INC_GOPENMOL_GOM_PARSER

#include <tcl.h>
#ifdef STDC_HEADERS
#  include <stddef.h>
#endif
#ifdef HAVE_STDINT_H
#  include <stdint.h>
#endif
#ifdef HAVE_INTTYPES_H
#  include <inttypes.h>
#endif

#include "gomcast.h"
#include "gomcext.h"
#include "parser_types.h"
#include "printmsg.h"

/****************************************************************/
/* Use GOM_PARSER_PART and GOM_PARSER_FINAL macros to init      */
/* elements of arrays of gom_ParserArgumentList type.           */
/****************************************************************/
typedef struct gom_ParserArgumentList_struct {
    const char *command;
    struct {
        const char *msg;
        struct {
            int min;
            int max;
        } count;
    } args;
    struct {
        const struct gom_ParserArgumentList_struct *list;
        int (*fnctn)(Tcl_Interp *, int, Tcl_Obj * const[], intptr_t);
        struct {
            intptr_t toggle;
            intptr_t inherit;
        } Flags;
    } action;
} gom_ParserArgumentList;

/****************************************************************/
/* If CMD macro is used:                                        */
/* - The next argument is treated as a command name.            */
/* - If there is an entry for whom command name equal to cmd    */
/*   field  or is an unique abbreviation of the cmd field, that */
/*   entry is selected. GOM_PARSER_CMD_SEPARATOR item may be    */
/*   separate commands into different sets (ie. "color" and     */
/*   "colour" shouldn't be in the same set).                    */
/* - Argument count is checked. If check fails, an error is     */
/*   reported.                                                  */
/* - If the argument is missing or no match occurs, all         */
/*   immediately following CMD entries are omitted. If          */
/*   there is no more entries, an error is reported. Otherwise  */
/*   ARGS entry will be processed.                              */
/* if ARGS macro is used:                                       */
/* - Entries are went through until argument count check        */
/*   succeeds or entry is not an ARGS entry.                    */
/* - If check succeeds, the entry is selected.                  */
/* - If a CMD entry is reached, CMD entry will be processed.    */
/* - If the end of the array is reached, an error is reported.  */
/*                                                              */
/* In the both cases:                                           */
/* - GOM_PARSER_END_ARGUMENT_LIST marks the end of the list.    */
/* - If PART macro is used, the actual argument count may       */
/*   exceed the maximum argument count. Arguments over the      */
/*   maximum count are left to the next part.                   */
/* - Arguments are stored to a Tcl_Obj* array.                  */
/* - If PART macro is used, the next part (array) is processed. */
/* - If FINAL macro  is used, the function is called  and the   */
/*   return value of the function is returned.                  */
/* - Function should be declared as:                            */
/*   int function_name(Tcl_Interp *interp,                      */
/*                     int objc,Tcl_Obj * const objv[],         */
/*                     int userValue)                           */
/* - Stored arguments are in objv array which have objc         */
/*   elements.                                                  */
/* - Some Parser macros require variables interp, objc and      */
/*   objv to exist, so it is better to use these names.         */
/* - Value of userValue is controlled by flags fields.          */
/****************************************************************/
#define GOM_PARSER_CMD_PART(cmd,args,array,flags) \
    {cmd "\0",args,{array,NULL,flags}}
#define GOM_PARSER_FINAL_CMD(cmd,args,fnctn,flags) \
    {cmd "\0",args,{NULL,fnctn,flags}}
#define GOM_PARSER_ARGS_PART(args,array,flags) \
    {NULL,args,{array,NULL,flags}}
#define GOM_PARSER_FINAL_ARGS(args,fnctn,flags) \
    {NULL,args,{NULL,fnctn,flags}}

#define GOM_PARSER_END_ARGUMENT_LIST \
    {NULL,{NULL,{0,0}},{NULL,NULL,{0,0}}}


/****************************************************************/
/* May be used to give aliases for command names.               */
/*                                                              */
/* Adding similar gom_ParserArgumentList entries is more        */
/* effient way to implement aliases. Use this only in the case  */
/* the aliases have a common prefix like the commands "colour"  */
/* and "color".                                                 */
/****************************************************************/
#define GOM_PARSER_CMD_ALIAS                    "\0"

/****************************************************************/
/* May be used as args field for GOM_PARSER_PART and            */
/* GOM_PARSER_FINAL macros.                                     */
/****************************************************************/
#define GOM_PARSER_NO_MORE_ARGS                 {NULL,{0,0}}
#define GOM_PARSER_NEED_ARGS(msg,count)         {msg,{count,count}}
#define GOM_PARSER_NEED_ARGS_RANGE(msg,min,max) {msg,{min,max}}

/****************************************************************/
/* May be used as flags field for GOM_PARSER_PART and           */
/* GOM_PARSER_FINAL macros.                                     */
/****************************************************************/
#define GOM_PARSER_UNUSED_VALUE               {0,0}
#define GOM_PARSER_INHERIT_VALUE              {0,~(intptr_t)0}
#define GOM_PARSER_SET_VALUE(value)           {value,0}
#define GOM_PARSER_SET_POINTER_VALUE(pointer) \
    {(intptr_t)(const char *)&*(pointer),0}
#define GOM_PARSER_ADD_FLAGS(flags)           {flags,~(intptr_t)flags}
#define GOM_PARSER_TOGGLE_FLAGS(flags)        {flags,~(intptr_t)0}
#define GOM_PARSER_REMOVE_FLAGS(flags)        {0,~(intptr_t)flags}

#define GOM_PARSER_GET_POINTER_VALUE(value)  ((void *)value)

typedef struct {
    const char *key;
    int         value;
} gom_ParserIntEnum;

#define GOM_PARSER_END_INT_ENUM_LIST       {NULL,0}

extern int gomp_ParserAddArgDescriptionToErrorMsg(
    Tcl_Interp *,const char *);
extern int gomp_ParserCheckDoubleRange(Tcl_Interp *,double,double,double);
extern int gomp_ParserCheckIntRange(Tcl_Interp *,int,int,int);
extern int gomp_ParserParseIntEnum(Tcl_Interp *,Tcl_Obj *,
                                   const gom_ParserIntEnum [], int *);
extern int gomp_ParserParseColour(Tcl_Interp *interp,Tcl_Obj *obj,
                                  gom_FloatColour *fc);
extern int gomp_ParserFormatReturn(Tcl_Interp *, const char *, ...)
                CHECK_FORMAT_PRINTF_1_EXTRA;

/****************************************************************/
/* Format: Only following sequences are allowed:                */
/*     space,tab    Will be omitted.                            */
/*     %s           string      (const char *)                  */
/*     %c           character   (int)                           */
/*     %f           real number (double)                        */
/*     %d, %i       integer     (int)                           */
/*     %ld,%li,     long        (long int)                      */
/*     %db,%ib      boolean     (int)                           */
/*     %p           object      (Tcl_Obj*)                      */
/*     (...)        sublist                                     */
/****************************************************************/
extern Tcl_Obj* gomp_ParserCreateList(const char * /* format */, ...)
                CHECK_FORMAT_PRINTF;

extern gom_ParserIntEnum gomp_ParserBooleanEnumTable[];

/****************************************************************/
/* Parse arguments of a Tcl command.                            */
/* Client data must be address to a constant                    */
/* TclParserArgumentList array whose last element is            */
/* GOM_PARSER_END_INT_ENUM_LIST.                                */
/****************************************************************/
extern int gomp_ParserCommand(ClientData,Tcl_Interp *,int,Tcl_Obj *const []);

/****************************************************************/
/* Declare command arguments. Some of the arguments may be      */
/* optional.                                                    */
/* GOM_PARSER_START_ARGS                                        */
/*  GOM_PARSER_ARG(Seg) GOM_PARSER_ARG(Res) GOM_PARSER_ARG(Atm) */
/*  GOM_PARSER_ARG(Fradius)                                     */
/* GOM_PARSER_END_ARGS;                                         */
/****************************************************************/
#define GOM_PARSER_START_ARGS enum {
#define GOM_PARSER_ARG(arg)       GomParser_I_##arg,
#define GOM_PARSER_SELECTIONLIST_ARG(num) \
                                  GomParser_I_Seg##num, \
                                  GomParser_I_Res##num, \
                                  GomParser_I_Atm##num,
#define GOM_PARSER_END_ARGS       GomParser_I_TclParser_Not_Used \
                              }

/****************************************************************/
/* Verify the result of GOM_PARSER macro expression.            */
/* If expression is false, adds argument description to the Tcl */
/* error message, generates a gOpenMol error and returns        */
/* TCL_ERROR.                                                   */
/****************************************************************/
#define GOM_PARSER_VERIFY_GOM_PARSER(description,expr) \
    if ( ! (expr) ) { \
        gomp_ParserAddArgDescriptionToErrorMsg( \
            interp,description); \
        gomp_PrintERROR(NULL); \
        return TCL_ERROR; \
    }

/****************************************************************/
/* Verify the result of a parser function.                      */
/* This can be used to share common part of parser functions.   */
/* Return TCL_ERROR if child function fails.                    */
/****************************************************************/
#define GOM_PARSER_VERIFY_CHILD(expr) \
    if ( (expr) != TCL_OK ) \
        return TCL_ERROR;

/****************************************************************/
/* Verify the result of an expression (should be true).         */
/* Expression must generate a gOpenMol error in case of         */
/* failure.                                                     */
/****************************************************************/
#define GOM_PARSER_VERIFY(expr) \
    if ( ! (expr) ) \
        return TCL_ERROR;

/****************************************************************/
/* All macros descripted below evaluate to nonzero on success.  */
/* On failure they set the Tcl error message and evaluate to    */
/* zero.                                                        */
/* Some macros have interp argument. Passing NULL as that       */
/* argument will leave Tcl error message unchanged.             */
/* That can be useful then wanting to keep the previous Tcl     */
/* error message.                                               */
/****************************************************************/

/****************************************************************/
/* Retrieve boolean value from the object declared by           */
/* GOM_PARSER_ARG(intArg).                                      */
/* The result is stored to a variable intArg which must be      */
/* declared as int.                                             */
/* Nonzero numbers and values "true", "yes" and "on" are        */
/* considered to be true.                                       */
/* Number zero and values "false", "no" and "off" are           */
/* considered to be false.                                      */
/* The object must exist (see GOM_PARSER_DEFAULT).              */
/****************************************************************/
#define GOM_PARSER_RETRIEVE_BOOLEAN(intArg) \
    (gomp_ParserParseIntEnum(\
         interp,GOM_PARSER_GET_OBJECT(intArg), \
         gomp_ParserBooleanEnumTable,&intArg) || \
     Tcl_GetBooleanFromObj( \
         NULL,GOM_PARSER_GET_OBJECT(intArg),&intArg)  == TCL_OK)
#define GOM_PARSER_RETRIEVE_BOOLEAN_SILENT(intArg) \
    (gomp_ParserParseIntEnum(\
         NULL,GOM_PARSER_GET_OBJECT(intArg), \
         gomp_ParserBooleanEnumTable,&intArg) || \
     Tcl_GetBooleanFromObj( \
         NULL,GOM_PARSER_GET_OBJECT(intArg),&intArg)  == TCL_OK)

/****************************************************************/
/* Retrieve double value from the object declared by            */
/* GOM_PARSER_ARG(doubleArg).                                   */
/* The result is stored to a variable doubleArg which must be   */
/* declared as double.                                          */
/* The object must exist (see GOM_PARSER_DEFAULT).              */
/****************************************************************/
#define GOM_PARSER_RETRIEVE_DOUBLE(doubleArg) \
    (Tcl_GetDoubleFromObj(\
        interp,GOM_PARSER_GET_OBJECT(doubleArg),&doubleArg) == TCL_OK)

/****************************************************************/
/* Retrieve integer value from the object declared by           */
/* GOM_PARSER_ARG(intArg).                                      */
/* The result is stored to a variable intArg which must be      */
/* declared as int.                                             */
/* The object must exist (see GOM_PARSER_DEFAULT).              */
/****************************************************************/
#define GOM_PARSER_RETRIEVE_INT(intArg) \
    (Tcl_GetIntFromObj(\
        interp,GOM_PARSER_GET_OBJECT(intArg),&intArg) == TCL_OK)
#define GOM_PARSER_RETRIEVE_INT_SILENT(intArg) \
    (Tcl_GetIntFromObj(\
        NULL,GOM_PARSER_GET_OBJECT(intArg),&intArg) == TCL_OK)

/****************************************************************/
/* Retrieve string value from the object declared by            */
/* GOM_PARSER_ARG(stringArg).                                   */
/* The result is stored to a variable stringArg which must be   */
/* declared as const char *.                                    */
/* The object must exist (see GOM_PARSER_DEFAULT).              */
/****************************************************************/
#define GOM_PARSER_RETRIEVE_STRING(stringArg) \
    (stringArg = Tcl_GetStringFromObj( \
        GOM_PARSER_GET_OBJECT(stringArg),NULL))

/****************************************************************/
/* Retrieve string value from the object declared by            */
/* GOM_PARSER_SELECTIONLIST_ARG(num).                           */
/* The result is stored to a variable listVar which must be     */
/* declared as SelectionList.                                   */
/* The objects do not have to exist.                            */
/****************************************************************/
#define GOM_PARSER_RETRIEVE_SELECTIONLIST(num,listVar) \
    ( (listVar).Segment = \
      GOM_PARSER_HAS_ARG(Seg##num) ? \
      Tcl_GetStringFromObj(GOM_PARSER_GET_OBJECT(Seg##num),NULL) : "*", \
      (listVar).Residue = GOM_PARSER_HAS_ARG(Res##num) ? \
      Tcl_GetStringFromObj(GOM_PARSER_GET_OBJECT(Res##num),NULL) : "*", \
      (listVar).Atom    = GOM_PARSER_HAS_ARG(Atm##num) ? \
      Tcl_GetStringFromObj(GOM_PARSER_GET_OBJECT(Atm##num),NULL) : "*", \
      1 )

/****************************************************************/
/* Retrieve colour value from the object declared by            */
/* GOM_PARSER_ARG(colourArg).                                   */
/* Argument may be colour name or a RGB triplet which items are */
/* intergers between 0 and 255 or floats between 0.0 and 1.0.   */
/* The result is stored to a variable colourArg which must be   */
/* declared as FloatColour.                                     */
/* The object must exist (see GOM_PARSER_DEFAULT).              */
/****************************************************************/
#define GOM_PARSER_RETRIEVE_COLOUR(colourArg) \
    gomp_ParserParseColour(interp, \
        GOM_PARSER_GET_OBJECT(colourArg),&colourArg)

/****************************************************************/
/* Return Tcl object declared by GOM_PARSER_ARG(arg).           */
/* The object must exist (see GOM_PARSER_DEFAULT and            */
/* GOM_PARSER_HAS_ARG).                                         */
/****************************************************************/
#define GOM_PARSER_GET_OBJECT(arg) objv[GomParser_I_##arg]

/****************************************************************/
/* Return array of Tcl objects starting from the object         */
/* declared by GOM_PARSER_ARG(arg).                             */
/****************************************************************/
#define GOM_PARSER_GET_OBJECT_LIST(arg) &objv[GomParser_I_##arg]

/****************************************************************/
/* Return true (nonzero) value if argument arg is supplied.     */
/****************************************************************/
#define GOM_PARSER_HAS_ARG(arg)    (GomParser_I_##arg<objc)

/****************************************************************/
/* Assign default value to the variable if corresponding object */
/* doesn't exist.                                               */
/* Return nonzero if value is assigned and zero otherwise.      */
/* So it is possible to use expressions like:                   */
/*     int Iatom;                                               */
/*     GOM_PARSER_START_ARGS                                    */
/*         GOM_PARSER_ARG(Iatom)                                */
/*     GOM_PARSER_END_ARGS;                                     */
/*     GOM_PARSER_VERIFY_GOM_PARSER("atom index",               */
/*         GOM_PARSER_DEFAULT(Iatom,1) ||                       */
/*         GOM_PARSER_RETRIEVE_INT(Iatom) );                    */
/****************************************************************/
#define GOM_PARSER_DEFAULT(arg,defval) \
    (GOM_PARSER_HAS_ARG(arg) ? 0 : (arg = (defval), 1))

/****************************************************************/
/* If argument <arg> does exist, sets pointer p<arg> to the     */
/* address of <arg> and returns zero.                           */
/* If argument <arg> does not exist, sets pointer p<arg> to     */
/* NULL and returns nonzero value.                              */
/* So it is possible to use expressions like:                   */
/*     int Iatom,*pIatom;                                       */
/*     GOM_PARSER_START_ARGS                                    */
/*         GOM_PARSER_ARG(Iatom)                                */
/*     GOM_PARSER_END_ARGS;                                     */
/*     GOM_PARSER_VERIFY_GOM_PARSER("atom index",               */
/*         GOM_PARSER_NULL_POINTER(Iatom) ||                    */
/*         GOM_PARSER_RETRIEVE_INT(Iatom) );                    */
/****************************************************************/
#define GOM_PARSER_NULL_POINTER(arg) \
    (GOM_PARSER_HAS_ARG(arg) ? (p##arg = &arg, 0) : (p##arg = NULL, 1))

/****************************************************************/
/* Retrieve enumerated integer value from the object declared   */
/* by GOM_PARSER_ARG(intArg).                                   */
/* The result is stored to a variable intArg which must be      */
/* declared as int.                                             */
/* Entries is an array of type                                  */
/* const GOM_PARSER_END_INT_ENUM_LIST[] whose last element is   */
/* GOM_PARSER_END_INT_ENUM_LIST.                                */
/* The object must exist (see GOM_PARSER_DEFAULT).              */
/****************************************************************/
#define GOM_PARSER_PARSE_INT_ENUM(intArg,entries) \
    gomp_ParserParseIntEnum( \
        interp,GOM_PARSER_GET_OBJECT(intArg),entries,&intArg)
#define GOM_PARSER_PARSE_INT_ENUM_NO_ERROR(intArg,entries) \
    gomp_ParserParseIntEnum( \
        NULL,GOM_PARSER_GET_OBJECT(intArg),entries,&intArg)

/****************************************************************/
/* Check the value of the variable to be in the range.          */
/* Value of the variable must be set.                           */
/****************************************************************/
#define GOM_PARSER_CHECK_DOUBLE_RANGE(doubleValue,min,max) \
    gomp_ParserCheckDoubleRange(interp,doubleValue,min,max)

/****************************************************************/
/* Check the value of the variable to be in the range.          */
/* Value of the variable must be set.                           */
/****************************************************************/
#define GOM_PARSER_CHECK_INT_RANGE(intValue,min,max) \
    gomp_ParserCheckIntRange(interp,intValue,min,max)

/****************************************************************/
/* These are macros compine RETRIEVE and CHECK_RANGE macros.    */
/****************************************************************/
#define GOM_PARSER_RETRIEVE_DOUBLE_CHECK_RANGE(doubleArg,min,max) \
    (GOM_PARSER_RETRIEVE_DOUBLE(doubleArg) && \
     GOM_PARSER_CHECK_DOUBLE_RANGE(doubleArg,min,max))
#define GOM_PARSER_RETRIEVE_INT_CHECK_RANGE(intArg,min,max) \
    (GOM_PARSER_RETRIEVE_INT(intArg) && \
     GOM_PARSER_CHECK_INT_RANGE(intArg,min,max))

#define GOM_PARSER_ARGLIST \
    Tcl_Interp *interp, int objc, Tcl_Obj * const objv[]
#define GOM_PARSER_PASS_ARGS interp,objc,objv

#define GOM_PARSER_RETURN_INT(value) \
    Tcl_SetObjResult(interp, Tcl_NewIntObj(value))
#define GOM_PARSER_RETURN_BOOLEAN(value) \
    Tcl_SetObjResult(interp, Tcl_NewBooleanObj((value)!=0))
#define GOM_PARSER_RETURN_DOUBLE(value) \
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value))
#define GOM_PARSER_RETURN_STRING(value) \
    Tcl_SetResult(interp,CONST_CAST(char*,value),TCL_VOLATILE)
#define GOM_PARSER_RETURN_LIST(list) \
    Tcl_SetObjResult(interp,gomp_ParserCreateList list)

/****************************************************************/
/* List operations.                                             */
/* Declaration:                                                 */
/*     gom_ParserList list;                                     */
/* GOM_PARSER_LIST_INIT                                         */
/*     Init the list.                                           */
/*     Returns a true value on success.                         */
/* GOM_PARSER_LIST_APPEND_*                                     */
/*     Append an element to the list.                           */
/*     Returns a true value on success.                         */
/* GOM_PARSER_LIST_APPEND_LIST                                  */
/*     Append a sub list to the list. Do not change or free     */
/*     the sub list if appending succees.                       */
/*     Returns a true value on success.                         */
/* GOM_PARSER_LIST_RETURN                                       */
/*     Return the list. Do not free the list after this.        */
/* GOM_PARSER_LIST_FREE                                         */
/*     Free the list. May be called even if                     */
/*     GOM_PARSER_LIST_INIT has failed.                         */
/****************************************************************/
typedef Tcl_Obj * gom_ParserList;
#define GOM_PARSER_LIST_INIT(list) ((list) = Tcl_NewListObj(0,NULL))
#define GOM_PARSER_LIST_APPEND_INT(list,value) \
    (Tcl_ListObjAppendElement(interp, list, Tcl_NewIntObj(value))==TCL_OK)
#define GOM_PARSER_LIST_APPEND_DOUBLE(list,value) \
    (Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(value))==TCL_OK)
#define GOM_PARSER_LIST_APPEND_STRING(list,value) \
    (Tcl_ListObjAppendElement(interp, list,\
                              Tcl_NewStringObj(value,-1))==TCL_OK)
#define GOM_PARSER_LIST_APPEND_LIST(list,sublist) \
    (Tcl_ListObjAppendElement(interp, list, sublist)==TCL_OK)
#define GOM_PARSER_LIST_RETURN(list) \
    Tcl_SetObjResult(interp,list)
#define GOM_PARSER_LIST_FREE(list) \
    if ( list ) { \
        Tcl_IncrRefCount(list); \
        Tcl_DecrRefCount(list); \
    }

#define GOM_PARSER_SUCCEEDED  return TCL_OK
#define GOM_PARSER_FAILED     return TCL_ERROR

typedef struct {
    float (*Func)(void);
} FloatFunc;

typedef struct {
    int (*Func)(void);
} IntFunc;

typedef struct {
    int (*Func)(void);
    struct {
        const char *false_value;
        const char *true_value;
    } strings;
} BooleanEnumFunc;

typedef struct {
    const char *(*Func)(void);
} StringFunc;

extern int gomp_ParserReturnFloatValueFromFunc(GOM_PARSER_ARGLIST,intptr_t);
extern int gomp_ParserReturnIntValueFromFunc(GOM_PARSER_ARGLIST,intptr_t);
extern int gomp_ParserReturnBooleanValueFromFunc(GOM_PARSER_ARGLIST,intptr_t);
extern int gomp_ParserReturnBooleanEnumValueFromFunc(GOM_PARSER_ARGLIST,intptr_t);
extern int gomp_ParserReturnStringValueFromFunc(GOM_PARSER_ARGLIST,intptr_t);
extern int gomp_ParserVerifyFuncResult(GOM_PARSER_ARGLIST,intptr_t);

#define GOM_PARSER_FINAL_FLOAT_CMD(cmd,pointer) \
GOM_PARSER_FINAL_CMD(cmd,\
                     GOM_PARSER_NO_MORE_ARGS,\
                     gomp_ParserReturnFloatValueFromFunc,\
                     GOM_PARSER_SET_POINTER_VALUE(pointer))

#define GOM_PARSER_FINAL_INT_CMD(cmd,pointer) \
GOM_PARSER_FINAL_CMD(cmd,\
                     GOM_PARSER_NO_MORE_ARGS,\
                     gomp_ParserReturnIntValueFromFunc,\
                     GOM_PARSER_SET_POINTER_VALUE(pointer))

#define GOM_PARSER_FINAL_BOOLEAN_CMD(cmd,pointer) \
GOM_PARSER_FINAL_CMD(cmd,\
                     GOM_PARSER_NO_MORE_ARGS,\
                     gomp_ParserReturnBooleanValueFromFunc,\
                     GOM_PARSER_SET_POINTER_VALUE(pointer))

#define GOM_PARSER_FINAL_BOOLEAN_ENUM_VALUE_CMD(cmd,pointer) \
GOM_PARSER_FINAL_CMD(cmd,\
                     GOM_PARSER_NO_MORE_ARGS,\
                     gomp_ParserReturnBooleanEnumValueFromFunc,\
                     GOM_PARSER_SET_POINTER_VALUE(pointer))

#define GOM_PARSER_FINAL_STRING_CMD(cmd,pointer) \
GOM_PARSER_FINAL_CMD(cmd,\
                     GOM_PARSER_NO_MORE_ARGS,\
                     gomp_ParserReturnStringValueFromFunc,\
                     GOM_PARSER_SET_POINTER_VALUE(pointer))

#define GOM_PARSER_FINAL_VERIFY_FUNC_CMD(cmd,pointer) \
GOM_PARSER_FINAL_CMD(cmd,\
                     GOM_PARSER_NO_MORE_ARGS,\
                     gomp_ParserVerifyFuncResult,\
                     GOM_PARSER_SET_POINTER_VALUE(pointer))

#endif /* INC_GOPENMOL_GOM_PARSER */
