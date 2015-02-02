/*

Copyright (c) 1996 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

in collaboration with

OpenMol Molecular Astrophysics Group
Max-Planck-Institut  fuer Astrophysik
Garching, GERMANY

Enhancements 2002 - 2003 by:
Eero Häkkinen
*/

#include "maindefs.h"

#define YES 1
#define NO  0
#define RABS(a)    ( ( a ) > 0.0 ? (a) : -(a))

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>

#if defined(GLUT)
#include <GL/glut.h>
#else
#include <GL/glaux.h>
#endif
#endif /* ENABLE_GRAPHICS */

#include <tcl.h>

#ifdef ENABLE_GRAPHICS
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>
#endif /* ENABLE_GRAPHICS */

#include "colouring.h"
#include "gomfile.h"
#include "gomstring.h"
#include "gomtext.h"
#include "memalloc.h"
#include "printmsg.h"
#include "projview.h"
#include "tclutils.h"

#include "stdafx.h"

#define MAX_FONTS  1000

/* text or symbol handling */

static int     txt_deep;          /* elements in text stack */

typedef struct {
    int     txt_dim;           /* text to be placed in 2D or 3D space */
    char   *txt_msg;           /* text stack                          */
    char   *txt_fnt;           /* font for text in stack              */
    float   txt_col[3];        /* text colour                         */
    float   txt_size;          /* font size                           */
    float   txt_xc;            /* x coordinate stack                  */
    float   txt_yc;            /* y coordinate stack                  */
    float   txt_zc;            /* z coordinate stack                  */
    float   txt_Xstep;         /* x step for moving text              */
    float   txt_Ystep;         /* y step for moving text              */
    float   txt_Zstep;         /* Z step for moving text              */
    float   txt_Xlim;          /* place in x for the text to stop     */
    float   txt_Ylim;          /* place in y for the text to stop     */
    float   txt_Zlim;          /* place in z for the text to stop     */
} TextManager_t;
   
static TextManager_t *TextManager;

static const char *GetTextEntryFromStack(int);
static const char *GetTextFontFromStack(int);
static float GetTextXCoordFromStack(int);
static float GetTextYCoordFromStack(int);
static float GetTextZCoordFromStack(int);
static float GetTextRedColourFromStack(int);
static float GetTextGreenColourFromStack(int);
static float GetTextBlueColourFromStack(int);
static int   GetTextPlotDimension(int);

/* ....................... */

#ifndef GLUT
static int  DispListBase;
static int  DispListBaseLength;
#endif

/**************************************************************************/
void gomp_PrintString(const char *s, const char *Font)
/**************************************************************************/
{
#ifdef ENABLE_GRAPHICS

#if defined(GLUT)
    int i;
    int il = strlen(s);
    const char *Value;

    for(i = 0 ; i < il ; i++) {

        if(Font[0] == '*')
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12 , s[i]);
        else
            if(gomp_StringMatch(Font ,        "BITMAP_8_BY_13"))  {
                glutBitmapCharacter(GLUT_BITMAP_8_BY_13        , s[i]);
            } else if(gomp_StringMatch(Font , "BITMAP_9_BY_15"))  {
                glutBitmapCharacter(GLUT_BITMAP_9_BY_15        , s[i]);
            } else if(gomp_StringMatch(Font , "BITMAP_TIMES_ROMAN_10"))  {
                glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10 , s[i]);
            } else if(gomp_StringMatch(Font , "BITMAP_TIMES_ROMAN_24"))  {
                glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24 , s[i]);
            } else if(gomp_StringMatch(Font , "BITMAP_HELVETICA_10"))  {
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10   , s[i]);
            } else if(gomp_StringMatch(Font , "BITMAP_HELVETICA_12"))  {
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12   , s[i]);
            } else if(gomp_StringMatch(Font , "BITMAP_HELVETICA_18"))  {
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18   , s[i]);
            } else {
                gomp_PrintERROR("unknown font name.\nResetting font to 'BITMAP_HELVETICA_12'");
                glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12 , s[i]);

                Value = Tcl_SetVar(gomp_GetTclInterp(),"gomTextFont","BITMAP_HELVETICA_12",TCL_GLOBAL_ONLY);
                if(!Value) {
                    gomp_PrintERROR("can't set tcl variable 'gomTextFont'");
                }
            }

    }

#else

#if defined(WIN32)

    glPushAttrib(GL_LIST_BIT);

/* display a string: 
   indicate start of glyph display lists */ 
    glListBase (DispListBase); 
/* now draw the characters in a string */
    glCallLists (strlen(s), GL_UNSIGNED_BYTE, (unsigned const char *)s); 

    glPopAttrib();
#else

    glPushAttrib(GL_LIST_BIT);
    glListBase(DispListBase);
    glCallLists(strlen(s) , GL_UNSIGNED_BYTE , (unsigned const char *)s);
    glPopAttrib();
#endif
#endif
#endif /* ENABLE_GRAPHICS */
}
/**************************************************************************/
int gomp_PushText2AnnotateStack(int   Index , const char *Text , const char *Font , 
                              float Red   , float Green, float Blue ,
                              float Xc    , float Yc   , float Zc)
/**************************************************************************/
{
    const float *sumXYZ;


    sumXYZ = gomp_GetTranslateArray();

    if(!txt_deep) {

        TextManager = gomp_AllocateVoidVector(sizeof(*TextManager));

        TextManager[0].txt_msg = gomp_AllocateCharVector(strlen(Text) + 1);
        strncpy(TextManager[0].txt_msg , Text , strlen(Text));
        TextManager[0].txt_msg[strlen(Text)] = (char)NULL;
        TextManager[0].txt_fnt = gomp_AllocateCharVector(strlen(Font) + 1);
        TextManager[0].txt_fnt[strlen(Font)] = (char)NULL;

        strncpy(TextManager[0].txt_fnt , Font , strlen(Font));

        TextManager[0].txt_col[0]    = Red;
        TextManager[0].txt_col[1]    = Green;
        TextManager[0].txt_col[2]    = Blue;

        if(Index == 3) {
            TextManager[0].txt_xc     = Xc - sumXYZ[0];
            TextManager[0].txt_yc     = Yc - sumXYZ[1];
            TextManager[0].txt_zc     = Zc - sumXYZ[2];
        }
        else {
            TextManager[0].txt_xc     = Xc;
            TextManager[0].txt_yc     = Yc;
            TextManager[0].txt_zc     = 0.0;
        }

        TextManager[0].txt_dim       = Index;

        txt_deep = 1;
    }
    else {
        TextManager = 
            gomp_ReallocateVoidVector(TextManager , 
                          (txt_deep + 1) * sizeof(*TextManager));

        TextManager[txt_deep].txt_msg = gomp_AllocateCharVector(strlen(Text) + 1);
        strncpy(TextManager[txt_deep].txt_msg , Text , strlen(Text));
        TextManager[txt_deep].txt_msg[strlen(Text)] = (char)NULL;

        TextManager[txt_deep].txt_fnt = gomp_AllocateCharVector(strlen(Font) + 1);
        strncpy(TextManager[txt_deep].txt_fnt , Font , strlen(Font));
        TextManager[txt_deep].txt_fnt[strlen(Font)] = (char)NULL;

        TextManager[txt_deep].txt_col[0]    = Red;
        TextManager[txt_deep].txt_col[1]    = Green;
        TextManager[txt_deep].txt_col[2]    = Blue;

        if(Index == 3) {
            TextManager[txt_deep].txt_xc     = Xc - sumXYZ[0];
            TextManager[txt_deep].txt_yc     = Yc - sumXYZ[1];
            TextManager[txt_deep].txt_zc     = Zc - sumXYZ[2];
        }
        else {
            TextManager[txt_deep].txt_xc     = Xc;
            TextManager[txt_deep].txt_yc     = Yc;
            TextManager[txt_deep].txt_zc     = 0.0;
        }

        TextManager[txt_deep].txt_dim       = Index;

        txt_deep++;
    }

    return(0);
}

/**************************************************************************/
int gomp_DeleteTextStack()
/**************************************************************************/
{
    int i;

    if(!txt_deep) {
/*       gomp_PrintWARNING("no text stack available to be deleted");*/
        return(0);
    }

    for(i = 0 ; i < txt_deep ; i++) {
        free(TextManager[i].txt_msg);
        free(TextManager[i].txt_fnt);
    }
   
    free(TextManager);

    txt_deep = 0;

    return(0);
}

/**************************************************************************/
int   gomp_GetEntriesInTextStack()
/**************************************************************************/
{
    return(txt_deep);
}
/**************************************************************************/
const char *GetTextEntryFromStack(int Index)
/**************************************************************************/
{
    if(Index < 0 || Index >= txt_deep) {
        gomp_PrintERROR("index to text stack is out of range");
        return("");
    }

    return(TextManager[Index].txt_msg);
}
/**************************************************************************/
const char *GetTextFontFromStack(int Index)
/**************************************************************************/
{
    if(Index < 0 || Index >= txt_deep) {
        gomp_PrintERROR("index to text stack is out of range");
        return("");
    }

    return(TextManager[Index].txt_fnt);
}

/**************************************************************************/
float GetTextXCoordFromStack(int Index)
/**************************************************************************/
{
    if(Index < 0 || Index >= txt_deep) {
        gomp_PrintERROR("index to text stack is out of range");
        return((float)-1.e-20);
    }

    return(TextManager[Index].txt_xc);
}

/**************************************************************************/
float GetTextYCoordFromStack(int Index)
/**************************************************************************/
{
    if(Index < 0 || Index >= txt_deep) {
        gomp_PrintERROR("index to text stack is out of range");
        return((float)-1.e-20);
    }

    return(TextManager[Index].txt_yc);
}

/**************************************************************************/
float GetTextZCoordFromStack(int Index)
/**************************************************************************/
{
    if(Index < 0 || Index >= txt_deep) {
        gomp_PrintERROR("index to text stack is out of range");
        return((float)-1.e-20);
    }

    return(TextManager[Index].txt_zc);
}

/**************************************************************************/
float GetTextRedColourFromStack(int Index)
/**************************************************************************/
{
    if(Index < 0 || Index >= txt_deep) {
        gomp_PrintERROR("index to text stack is out of range");
        return((float)0.0);
    }

    return(TextManager[Index].txt_col[0]);
}

/**************************************************************************/
float GetTextGreenColourFromStack(int Index)
/**************************************************************************/
{
    if(Index < 0 || Index >= txt_deep) {
        gomp_PrintERROR("index to text stack is out of range");
        return((float)0.0);
    }

    return(TextManager[Index].txt_col[1]);
}
/**************************************************************************/
float GetTextBlueColourFromStack(int Index)
/**************************************************************************/
{
    if(Index < 0 || Index >= txt_deep) {
        gomp_PrintERROR("index to text stack is out of range");
        return((float)0.0);
    }

    return(TextManager[Index].txt_col[2]);
}

/**************************************************************************/
int   gomp_PlotTextStack()
/**************************************************************************/
{
#ifdef ENABLE_GRAPHICS
    static int    i;
/*  static int   mm;*/
    static float xc;
    static float yc;
    static float zc;
    static float col[3];
    static const char *Value;
    static char  font[BUFF_LEN];
    static const float *RotMP;

    glDisable(GL_LIGHTING);

    Value = Tcl_GetVar(gomp_GetTclInterp(),"gomTextFont",TCL_GLOBAL_ONLY);

    if(Value) {
        gomp_CopyString(font,Value,BUFF_LEN);
    } else {
        gomp_CopyString(font,"BITMAP_HELVETICA_12",BUFF_LEN);
    }

    for (i = 0 ; i < txt_deep ; i++) {

        switch(GetTextPlotDimension(i)) {

        case 2: /* 2D plot */
            glMatrixMode(GL_PROJECTION);
            glPushMatrix();
            glLoadIdentity();
            glOrtho(0.0, 1.0 , 0.0, 1.0 , 0.0 , 1.0);
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
/*           glMatrixMode(mm);*/

            xc = GetTextXCoordFromStack(i);
            yc = GetTextYCoordFromStack(i);
            zc = GetTextZCoordFromStack(i);

            col[0] = GetTextRedColourFromStack(i);
            col[1] = GetTextGreenColourFromStack(i);
            col[2] = GetTextBlueColourFromStack(i);

            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&col[0] , &col[1] , &col[2]);

            glColor3fv(col);
            glRasterPos3f( xc , yc , 0.0);
            gomp_PrintString(GetTextEntryFromStack(i) , font);

            glPopMatrix();
            glMatrixMode(GL_PROJECTION);
            glPopMatrix();
            glMatrixMode(GL_MODELVIEW);

            break;
        case 3: /* 3D plot */


            RotMP = gomp_GetSavedModelViewMatrixMT(0);
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();

            if(gomp_GetProjectionTransformation() == ORTHOGRAPHIC_VIEW) {
            } else if(gomp_GetProjectionTransformation() == PERSPECTIVE_VIEW) {
                gluLookAt(0.0 , 0.0 , 
                          (gomp_GetPerspectiveNear() + 0.5 * gomp_GetPerspectiveWindow()) ,
                          0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0);
            } else {
                gomp_PrintERROR("undefined Projection Transformation defined");
                return(1);
            }

            glMultMatrixf(RotMP);

            xc = GetTextXCoordFromStack(i);
            yc = GetTextYCoordFromStack(i);
            zc = GetTextZCoordFromStack(i);

            col[0] = GetTextRedColourFromStack(i);
            col[1] = GetTextGreenColourFromStack(i);
            col[2] = GetTextBlueColourFromStack(i);

            if(!gomp_GetDisplayColourType())
                (void)gomp_RGB2Grayscale(&col[0] , &col[1] , &col[2]);

            glColor3fv(col);
            glRasterPos3f( xc , yc , zc);
            gomp_PrintString(GetTextEntryFromStack(i) , font);

            glPopMatrix();

            break;
        }
    }

    glEnable(GL_LIGHTING);
#endif /* ENABLE_GRAPHICS */

    return(0);
}
/**************************************************************************/
int   GetTextPlotDimension(int Index)
/**************************************************************************/
{

    if(Index < 0 || Index >= txt_deep) {
        gomp_PrintERROR("index to text stack is out of range");
        return(0);
    }

    return(TextManager[Index].txt_dim);

}
/**************************************************************************/
int    gomp_WriteText2ModelFile(FILE *Model_f)
/**************************************************************************/
{
    int    i;
    const float *sumxyz;

    if(gomp_GetEntriesInTextStack() < 1) {
        return(0);
    }

    sumxyz = gomp_GetTranslateArray();

/* start the writing ... */

    fprintf(Model_f , "[Text]\n");

    fprintf(Model_f , "%d\n",gomp_GetEntriesInTextStack());

    for(i = 0 ; i < gomp_GetEntriesInTextStack() ; i++) {

        if(GetTextPlotDimension(i) == 2) {
            fprintf(Model_f , "%f %f %f %f %f %f %d\n",
                    GetTextXCoordFromStack(i),
                    GetTextYCoordFromStack(i),
                    GetTextZCoordFromStack(i),
                    GetTextRedColourFromStack(i),
                    GetTextGreenColourFromStack(i),
                    GetTextBlueColourFromStack(i),
                    GetTextPlotDimension(i));
        } else {
            fprintf(Model_f , "%f %f %f %f %f %f %d\n",
                    GetTextXCoordFromStack(i)+sumxyz[0],
                    GetTextYCoordFromStack(i)+sumxyz[1],
                    GetTextZCoordFromStack(i)+sumxyz[2],
                    GetTextRedColourFromStack(i),
                    GetTextGreenColourFromStack(i),
                    GetTextBlueColourFromStack(i),
                    GetTextPlotDimension(i));
        }

        fprintf(Model_f , "%s\n",
                GetTextFontFromStack(i));
        fprintf(Model_f , "%s\n",
                GetTextEntryFromStack(i));
    } 
       
    return(0);
}
/**************************************************************************/
int    gomp_ReadText2ModelFile(FILE *Model_f)
/**************************************************************************/
{
    int  TextEntries;
    int  dime;
    int  i;
    char Text[BUFF_LEN];
    char Text1[BUFF_LEN];
    float xc;
    float yc;
    float zc;
    float redc;
    float greenc;
    float bluec;
       
/* start the process ... */

    (void)gomp_DeleteTextStack();

    gomp_Fgets(Text,BUFF_LEN,Model_f);

    sscanf(Text ,"%d", &TextEntries);

    for(i = 0 ; i < TextEntries ; i++) {
        gomp_Fgets(Text,BUFF_LEN,Model_f);
        sscanf(Text ,"%f %f %f %f %f %f %d",
               &xc,&yc,&zc,&redc,&greenc,&bluec,&dime);
        gomp_Fgets(Text, BUFF_LEN,Model_f);
        gomp_Fgets(Text1,BUFF_LEN,Model_f);
        (void)gomp_PushText2AnnotateStack(
            dime,Text1,Text,redc,greenc,bluec,xc,yc,zc); 
    }

    return(0);
}
