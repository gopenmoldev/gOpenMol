/*

Copyright (c) 1994 - 2005 by:
Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND
Confidential unpublished property of 
Leif Laaksonen
All rights reserved

Enhancements 2002 - 2005 by:
Eero Häkkinen
*/

#include "maindefs.h"

#include "gomstdio.h"
#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#if !defined(WIN32)
#include <sys/types.h>
#endif

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef ENABLE_GRAPHICS
#include <GL/gl.h>
#include <GL/glu.h>
#endif /* ENABLE_GRAPHICS */

#include "bond.h"
#include "cell.h"
#include "colouring.h"
#include "contour.h"
#include "coord_file.h"
#include "gomhelp.h"
#include "gomstring.h"
#include "gomtcl.h"
#include "gomwindow.h"
#include "light_model.h"
#include "molecoord.h"
#include "molecule.h"
#include "molecstruct.h"
#include "picking.h"
#include "plot_molec.h"
#include "plot.h"
#include "printmsg.h"
#include "projview.h"
#include "selection.h"
#include "stereo.h"
#include "tclutils.h"
#include "text_stack.h"
#include "trajectory.h"

#include "stdafx.h"

/*********************************************************************/
int gomp_DefineCommand(ClientData clientdata, Tcl_Interp *interp,
                     int argc, const char **argv)
/*********************************************************************/
{
    static int   Type;
    static char  Text[BUFF_LEN];
    static char  Text1[BUFF_LEN];
    static char  Text2[BUFF_LEN];
    static char  Text3[BUFF_LEN];
    static char  Text4[BUFF_LEN];
    static char  Text5[BUFF_LEN];
    static int   ITemp;
    static int   ITemp1;
    static int   i;
    static float FTemp1;
    static float FTemp2;
    static float FTemp3;
    static int   Wstr;
    static float RedC;
    static float GreenC;
    static float BlueC;
    static const char *Value;
    static int   gomp_InputView;

/* #1 define  ... */
    strncpy(Text,gomp_GetNextFromParserStack(argc , argv),BUFF_LEN-1);
    if(gomp_StringMatch(Text , "defi$ne")) {

        Wstr = 0;

        gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
/* node */
        if(gomp_StringMatch(Text , "node")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("node name missing");
                return(TCL_ERROR);
            }
            if(!gomp_PushNode2Stack(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* help function used to move text between c-code and tcl */
        else if(gomp_StringMatch(Text , "gtext")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(!gomp_SetGlobalTextString(Text1)) 
                return(TCL_OK);
            else
                return(TCL_ERROR);
       
        }
/* colour type */
        else if(gomp_StringMatch(Text , "colourt$ype") || 
                gomp_StringMatch(Text , "colort$ype")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text , "gray$scale") ||
               gomp_StringMatch(Text , "grey$scale")) {
                (void)gomp_SetDisplayColourType(OFF);
                (void)gomp_ResetPrepare1DTexture(); 
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text , "colo$ur") ||
                    gomp_StringMatch(Text , "colo$r")) {
                (void)gomp_SetDisplayColourType(ON);
                (void)gomp_ResetPrepare1DTexture(); 
                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("wrong colour type defined (grayscale/colour)");
                return(TCL_ERROR);
            }
        }
/* projection */
        else if(gomp_StringMatch(Text , "proj$ection")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text , "orth$ographic")) {
                (void)gomp_SetProjectionTransformation(ORTHOGRAPHIC_VIEW);
                (void)gomp_ResetView();
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text , "pers$pective")) {
                (void)gomp_SetProjectionTransformation(PERSPECTIVE_VIEW);
                (void)gomp_ResetView();
                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("wrong type of Projection Type defined");
                return(TCL_ERROR);
            }
        }
#ifdef ENABLE_GRAPHICS
#if defined(GLUT)
/* windowing */
        else if(gomp_StringMatch(Text , "wind$owing")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text , "sing$le")) {
                (void)gomp_SetWindowingStyle(SINGLE_WINDOWING);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text , "mult$i")) {
                (void)gomp_SetWindowingStyle(MULTI_WINDOWING);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text , "upda$te") ||
                    gomp_StringMatch(Text , "redr$aw")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_StringMatch(Text1 , "auto$matic")) {
                    if(gomp_SetUpdateDisplayMode(ON)) {
                        return(TCL_ERROR);
                    } else {
                        Value = Tcl_SetVar(gomp_GetTclInterp(),"gomInstantDisplay","1",0);
                        if(!Value) {
                            gomp_PrintERROR("can't set instant display trigger into variable 'gomInstantDisplay'");
                            return(TCL_ERROR);
                        }
                        return(TCL_OK);
                    }
                } else if(gomp_StringMatch(Text1 , "manu$al")) {
                    if(gomp_SetUpdateDisplayMode(OFF)) {
                        return(TCL_ERROR);
                    } else {
                        Value = Tcl_SetVar(gomp_GetTclInterp(),"gomInstantDisplay","0",0);
                        if(!Value) {
                            gomp_PrintERROR("can't set instant display trigger into variable 'gomInstantDisplay'");
                            return(TCL_ERROR);
                        }
                        return(TCL_OK);
                    }
                } else {
                    gomp_PrintERROR("wrong type of windowing update (automatic/manual)");
                    return(TCL_ERROR);
                }
            }
            else {
                gomp_PrintERROR("wrong type of windowing defined (single/multi/redraw)");
                return(TCL_ERROR);
            }
        }
#endif
#endif /* ENABLE_GRAPHICS */
/* default coordinate extension */
        else if(gomp_StringMatch(Text , "defa$ult")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text , "coor$dinate")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_StringMatch(Text1 , "type")) {
                    gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                    if(gomp_StringMatch(Text2 , "ambe$r")) {
                        (void)gomp_PutAMBERdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "char$mm")) {
                        (void)gomp_PutCHARMMdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "gaus$sian")) {
                        (void)gomp_PutGAUSSIANdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "hype$rchem")) {
                        (void)gomp_PutHYPERCHEMdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "insi$ght")) {
                        (void)gomp_PutINSIGHTdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "mol2")) {
                        (void)gomp_PutMOL2default();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "mumo$d")) {
                        (void)gomp_PutMUMODdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "open$mol")) {
                        (void)gomp_PutOPENMOLcenterdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "pdb")) {
                        (void)gomp_PutPDBdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "xmol")) {
                        (void)gomp_PutXMOLdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "xyz")) {
                        (void)gomp_PutXYZdefault();
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "yasp")) {
                        (void)gomp_PutYASPdefault();
                        return(TCL_OK);
                    }
                    else {
                        gomp_PrintERROR("unknown file type");
                        return(TCL_ERROR);
                    }
                }
                else if(gomp_StringMatch(Text1 , "exte$nsion")) {
                    gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                    gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                    if(gomp_StringMatch(Text2 , "ambe$r")) {
                        (void)gomp_PutAMBERcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "char$mm")) {
                        (void)gomp_PutCHARMMcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "gaus$sian")) {
                        (void)gomp_PutGAUSSIANcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "hype$rchem")) {
                        (void)gomp_PutHYPERCHEMcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "insi$ght")) {
                        (void)gomp_PutINSIGHTcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "mol2")) {
                        (void)gomp_PutMOL2coordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "mumo$d")) {
                        (void)gomp_PutMUMODcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "open$mol")) {
                        (void)gomp_PutOPENMOLcenterFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "pdb")) {
                        (void)gomp_PutPDBcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "xmol")) {
                        (void)gomp_PutXMOLcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "xyz")) {
                        (void)gomp_PutXYZcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else if(gomp_StringMatch(Text2 , "yasp")) {
                        (void)gomp_PutYASPcoordFileType(Text1);
                        return(TCL_OK);
                    }
                    else {
                        gomp_PrintERROR("unknown file type");
                        return(TCL_ERROR);
                    }
                }
            }
/* ................................................ */
/* default trajectory extension */
            else if(gomp_StringMatch(Text , "traj$ectory")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(gomp_StringMatch(Text1 , "type")) {
                    int Type = gomp_ParseTrajType(
                        gomp_GetNextFromParserStack(argc,NULL));
                    if ( Type > 0 ) {
                        gomp_PutTrajectoryTypeDefault(Type);
                        return(TCL_OK);
                    }
                    else
                        return(TCL_ERROR);
                }
                else if(gomp_StringMatch(Text1 , "exte$nsion")) {
                    int Type = gomp_ParseTrajType(
                        gomp_GetNextFromParserStack(argc,NULL));
                    gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                    if ( Type > 0 ) {
                        gomp_PutTrajectoryTypeFileExtension(Type, Text1);
                        return(TCL_OK);
                    }
                    else
                        return(TCL_ERROR);
                }
            }
            gomp_PrintERROR("unrecognized 'define trajectory' command");
            return(TCL_ERROR);
        }

/* draw buffer */
        else if(gomp_StringMatch(Text , "drawb$uffer")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("buffer value missing (front/back)");
                return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text , "fron$t")) {
                (void)gomp_SetDrawBuffer(1);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text , "back")) {
                (void)gomp_SetDrawBuffer(0);
                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("undefined buffer name");
                return(TCL_ERROR);
            }
        }
/* licorice sphere */
        else if(gomp_StringMatch(Text , "licos$phere")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("licorice sphere value missing");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                ITemp1 = 1;
            }
            else if(gomp_StringMatch(Text1, "all") || gomp_StringMatch(Text1, "*")) {
                if(!gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("no structures available");
                    return(TCL_ERROR);
                } 
                for ( i = 0 ; i < gomp_GetNumMolecStructs(); i++) {
                    if(gomp_PutAtomLicoRadS(i , atof(Text))) {
                        gomp_PrintERROR("can't set licorice sphere radius");
                        return(TCL_ERROR);
                    }
                }
                return(TCL_OK);
            }
            else
                ITemp1 = atoi(Text1);
            if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("structure number out of range");
                return(TCL_ERROR);
            } 
            ITemp1--;

            if(!gomp_PutAtomLicoRadS(ITemp1 , atof(Text)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* licorice cylinder */
        else if(gomp_StringMatch(Text , "licoc$ylinder")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("licorice cylinder value missing");
                return(TCL_ERROR);
            }
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                ITemp1 = 1;
            }
            else if(gomp_StringMatch(Text1, "all") || gomp_StringMatch(Text1, "*")) {
                if(!gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("no structures available");
                    return(TCL_ERROR);
                } 
                for ( i = 0 ; i < gomp_GetNumMolecStructs(); i++) {
                    if(gomp_PutAtomLicoRadC(i , atof(Text))) {
                        gomp_PrintERROR("can't set licorice cylinder radius");
                        return(TCL_ERROR);
                    }
                }
                return(TCL_OK);
            }
            else
                ITemp1 = atoi(Text1);
            if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                gomp_PrintERROR("structure number out of range");
                return(TCL_ERROR);
            } 
            ITemp1--;

            if(!gomp_PutAtomLicoRadC(ITemp1 , atof(Text)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* display list */
        else if(gomp_StringMatch(Text , "disp$laylists")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("parameter value missing");
                return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text , "on")) {
                (void)gomp_SetDisplayListState(ON);
                gomp_PrintMessage("Display lists are used");
            } else if(gomp_StringMatch(Text , "off")) {
                (void)gomp_SetDisplayListState(OFF);
                gomp_PrintMessage("Display lists are not used");
            } else {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("parameter value missing");
                    return(TCL_ERROR);
                }
                if(gomp_StringMatch(Text1 , "on"))
                    Type = ON;
                else if(gomp_StringMatch(Text1 , "off"))
                    Type = OFF;
                else {
                    gomp_PrintERROR("value has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
                if(gomp_StringMatch(Text , "defa$ult"))
                    (void)gomp_SetDefaultDisplayListState(Type);
                else if(gomp_ParseObjectDisplayListState(Text,Type)) {
                    gomp_PrintERROR("unrecognized type for 'define displaylists' command");
                    return(TCL_ERROR);
                }
            }

            return(TCL_OK);
        }

/* molecule line width */
        else if(gomp_StringMatch(Text , "mlin$ewidth")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("molecule line value missing");
                return(TCL_ERROR);
            }
            if(!gomp_SetMoleculeLineWidth(atoi(Text)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* contour line width */
        else if(gomp_StringMatch(Text , "clin$ewidth")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("contour line value missing");
                return(TCL_ERROR);
            }
            if(!gomp_SetContourLineWidth(atoi(Text)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* sphere quality */
        else if(gomp_StringMatch(Text , "sphe$requality")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("sphere quality value missing");
                return(TCL_ERROR);
            }
            if(!gomp_SetSphereQuality(atoi(Text)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }

/* cylinder quality */
        else if(gomp_StringMatch(Text , "cyli$nderquality")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("cylinder quality value missing");
                return(TCL_ERROR);
            }
            if(!gomp_SetCylinderQuality(atoi(Text)))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }

/* help url */
        else if(gomp_StringMatch(Text , "hurl")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("URL to the help files is missing");
                return(TCL_ERROR);
            }
/* update the widget ... */
/*      (void)gomp_SetHelpWebURL(Text); */
            if(!gomp_SetHURL(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* web browser */
        else if(gomp_StringMatch(Text , "webb$rowser")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("web browser name missing");
                return(TCL_ERROR);
            }
/* update the widget ... */
/*      (void)gomp_SetHelpWebBrowser(Text);*/
            if(!gomp_SetWebBrowser(Text))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* background colour */
        else if(gomp_StringMatch(Text , "bgco$lour") ||
                gomp_StringMatch(Text , "bgco$lor")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text) == 0) {
                gomp_PrintERROR("colour name/value is missing");
                return(TCL_ERROR);
            }

            if(gomp_ColourName2RGB(Text , &RedC , &GreenC , &BlueC)) {
                gomp_PrintERROR("can't convert colour name/value");
                return(TCL_ERROR);
            }
           
            if(!gomp_SetBGColor(RedC , GreenC , BlueC))
                return(TCL_OK);
            else
                return(TCL_ERROR);
        }
/* near clipping plane */
        else if(gomp_StringMatch(Text , "near$plane")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            if(gomp_StringMatch(Text , "step")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("near clipping plane step value is missing");
                    return(TCL_ERROR);
                }
                (void)gomp_SetPerspectiveNearStep(atof(Text1));
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text , "valu$e")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("near clipping plane value is missing");
                    return(TCL_ERROR);
                }
                {
                    float Twindow;
                    float Tnear;
                    float TnearNew;
                    float TnearOld;
                    float Tfar;

                    TnearOld = gomp_GetPerspectiveNear();
                    TnearNew = atof(Text1);
                    Twindow  = gomp_GetPerspectiveWindow() + 2.0 * (TnearOld - TnearNew);
                    Tnear    = TnearNew;
/* This is done too fool the far plane */
/*         Tfar     = Tnear + Twindow;*/
                    Tfar     = gomp_GetPerspectiveFar();
/*
  Tnear   = atof(Text1);
  Tfar    = gomp_GetPerspectiveFar();
  Twindow = Tfar - Tnear; 
*/
                    if(Tfar < Tnear) {
                        gomp_PrintERROR("Far clipping plane < Near clipping plane");
                        return(TCL_ERROR);
                    }
                    (void)gomp_SetPerspectiveNear(Tnear);
                    (void)gomp_SetPerspectiveFar(Tfar);
                    (void)gomp_SetPerspectiveWindow(Twindow);
#ifdef ENABLE_GRAPHICS
                    (void)gomp_ResetProjection();
#endif /* ENABLE_GRAPHICS */
                }
                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("unrecognized 'define nearplane' command");
                return(TCL_ERROR);
            }
        }
/* far clipping plane */
        else if(gomp_StringMatch(Text , "farp$lane")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            if(gomp_StringMatch(Text , "step")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("near clipping plane step value is missing");
                    return(TCL_ERROR);
                }
                (void)gomp_SetPerspectiveFarStep(atof(Text1));
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text , "valu$e")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("far clipping plane value is missing");
                    return(TCL_ERROR);
                }
                {
                    float Tnear;
                    float Tfar;
#if defined(CLIPPING_JUNK)
                    float Twindow;
                    float Tnear;
                    float Tfar;
                    float TfarNew;
                    float TfarOld;
                    float Diff;

                    TfarOld = gomp_GetPerspectiveFar();
                    TfarNew = atof(Text1);
                    Diff    = TfarNew - TfarOld;
                    Twindow = gomp_GetPerspectiveWindow() + 2.0 * Diff;
                    Tnear   = gomp_GetPerspectiveNear() - Diff;
                    Tfar    = TfarNew;

                    if(Tfar < Tnear) {
                        gomp_PrintERROR("Far clipping plane < Near clipping plane");
                        return(TCL_ERROR);
                    }

/*
  Twindow = Tfar - Tnear;
*/
                    (void)gomp_SetPerspectiveNear(Tnear - Diff);
                    (void)gomp_SetPerspectiveFar(Tfar);
                    (void)gomp_SetPerspectiveWindow(Twindow);
#ifdef ENABLE_GRAPHICS
                    (void)gomp_ResetProjection();
#endif /* ENABLE_GRAPHICS */
#endif
                    Tfar    = atof(Text1);
                    Tnear   = gomp_GetPerspectiveNear();
                    if(Tfar < Tnear) {
                        gomp_PrintERROR("Far clipping plane < Near clipping plane");
                        return(TCL_ERROR);
                    }
                    (void)gomp_SetPerspectiveFar(Tfar);
                }

                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("unrecognized 'define farplane' command");
                return(TCL_ERROR);
            }
        }
/* viewing angle */
        else if(gomp_StringMatch(Text , "viewa$ngle")) {
            gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            (void)gomp_SetPerspectiveAngle(atof(Text));
#ifdef ENABLE_GRAPHICS
            if(gomp_PrepareDisplay(0)) {
                return(TCL_ERROR);
            }
#endif /* ENABLE_GRAPHICS */
            return(TCL_OK);
        }
/* material */
        else if(gomp_StringMatch(Text , "mate$rial")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            if(gomp_StringMatch(Text1 , "spec$ular")) {

                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
         
/* red */
                if(gomp_StringMatch(Text2 , "red")) {

                    strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                            BUFF_LEN-1);
     
                    if(strlen(Text3) == 0) {
                        gomp_PrintERROR("red material specular value missing");
                        return(TCL_ERROR);
                    }
                    if(!gomp_SetMaterialSpecularRed(atof(Text3)))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
/* green */
                else if(gomp_StringMatch(Text2 , "gree$n")) {

                    strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                            BUFF_LEN-1);
     
                    if(strlen(Text3) == 0) {
                        gomp_PrintERROR("green material specular value missing");
                        return(TCL_ERROR);
                    }
                    if(!gomp_SetMaterialSpecularGreen(atof(Text3)))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);

                }
/* blue */
                else if(gomp_StringMatch(Text2 , "blue")) {

                    strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                            BUFF_LEN-1);
     
                    if(strlen(Text3) == 0) {
                        gomp_PrintERROR("blue material specular value missing");
                        return(TCL_ERROR);
                    }
                    if(!gomp_SetMaterialSpecularBlue(atof(Text3)))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);

                }
                else {
                    gomp_PrintERROR("unrecognized colour, should be red, green or blue"); 
                    return(TCL_ERROR);
                }
            }
            else if(gomp_StringMatch(Text1 , "shin$iness")) {

                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("material shininess value missing");
                    return(TCL_ERROR);
                }

                if(!gomp_SetMaterialShininess(atof(Text2)))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("unrecognized 'define material' command");
                return(TCL_ERROR);
            }
        }
/* node ... */
        else if(gomp_StringMatch(Text , "-nod$e")) {
            gomp_PrintERROR("not yet implemented");
            return(TCL_ERROR);
        }
/* very elaborate stuff */
        else if(gomp_StringMatch(Text , "stru$cture")) {
/* structure name */
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("structure name is missing");
                return(TCL_ERROR);
            }

/* atoms to be in the structure */
            gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text2) == 0) {
                gomp_PrintERROR("number of atoms is missing");
                return(TCL_ERROR);
            }
            ITemp = atoi(Text2);
            if(ITemp <= 0) {
                gomp_PrintERROR("number of atoms has to be > 0");
                return(TCL_ERROR);
            }
/* action parameter */
            gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text3) == 0 || gomp_StringMatch(Text3 , "new"))
                return(gomp_CreateMolecStruct(Text1 , ITemp , NEW) >= 0 ? TCL_OK : TCL_ERROR);
            else if(gomp_StringMatch(Text3 , "appe$nd"))
                return(gomp_CreateMolecStruct(Text1 , ITemp , APPEND) >= 0 ? TCL_OK : TCL_ERROR);
            else {
                gomp_PrintERROR("option has to be 'new' or 'append'");
                return(TCL_ERROR);
            }
        }
/* atom */
        else if(gomp_StringMatch(Text , "atom")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "labe$l") ||
               gomp_StringMatch(Text1 , "symb$ol")) {
                gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text) == 0) {
                    gomp_PrintERROR("atom symbol is missing");
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("atom index is missing");
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text2) == 0) {
                    ITemp1 = 1;
                }
                else
                    ITemp1 = atoi(Text2);
                if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure number out of range");
                    return(TCL_ERROR);
                }
                ITemp1--;
                ITemp = atoi(Text1);
                if(ITemp < 1 || ITemp > gomp_GetNumAtomsInMolecStruct(ITemp1)) {
                    gomp_PrintERROR("atom number is out of range");
                    return(TCL_ERROR);
                }
                ITemp--;

                if(gomp_GetAtomSegName(ITemp1 , ITemp) == NULL)
                    (void)gomp_PutAtomSegName(ITemp1 , "gom" , ITemp);
                if(gomp_GetAtomResName(ITemp1 , ITemp) == NULL)
                    (void)gomp_PutAtomResName(ITemp1 , "gom" , ITemp);
                (void)gomp_PutAtomAtmName(ITemp1 , Text ,  ITemp);
                return(TCL_OK);
            }
/* define cross size (length) */
            else if(gomp_StringMatch(Text1 , "cros$s")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("cross length value missing");
                    return(TCL_ERROR);
                }
                if(!gomp_SetCrossLen(atof(Text2))) 
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
       
            }
/* atom coordinates */
            else if(gomp_StringMatch(Text1 , "coor$dinates")) {
                const float *sumxyz = gomp_GetTranslateArray();
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text5,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text4) == 0) {
                    gomp_PrintERROR("atom index is missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text5) == 0) {
                    ITemp1 = 1;
                }
                else
                    ITemp1 = atoi(Text5);
                if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure number out of range");
                    return(TCL_ERROR);
                }
                ITemp1--;
                ITemp = atoi(Text4);
                if(ITemp < 1 || ITemp > gomp_GetNumAtomsInMolecStruct(ITemp1)) {
                    gomp_PrintERROR("atom number is out of range");
                    return(TCL_ERROR);
                }
                ITemp--;
                (void)gomp_PutAtomXCoord(ITemp1 , atof(Text1) - sumxyz[0] , ITemp);
                (void)gomp_PutAtomYCoord(ITemp1 , atof(Text2) - sumxyz[1] , ITemp);
                (void)gomp_PutAtomZCoord(ITemp1 , atof(Text3) - sumxyz[2] , ITemp);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "sele$ction")) {

                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
          
                if(gomp_StringMatch(Text2 , "on")) {
                    if(!gomp_SetSelectionModeStatus(ON)) {
                        Value = Tcl_SetVar(gomp_GetTclInterp(),"gomSelectionState","1",TCL_GLOBAL_ONLY);
                        if(!Value) {
                            gomp_PrintERROR("can't set tcl variable 'gomSelectionState'");
                            return(TCL_ERROR);
                        }
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else if(gomp_StringMatch(Text2 , "off")) {
                    if(!gomp_SetSelectionModeStatus(OFF)) {
                        Value = Tcl_SetVar(gomp_GetTclInterp(),"gomSelectionState","0",TCL_GLOBAL_ONLY);
                        if(!Value) {
                            gomp_PrintERROR("can't set tcl variable 'gomSelectionState'");
                            return(TCL_ERROR);
                        }
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else {
                    gomp_PrintERROR("option has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
            }
/* atom selection */
            else if(gomp_StringMatch(Text1 , "iden$tify")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(gomp_StringMatch(Text2 , "on")) {
                    if(!gomp_SetIdentifyAtomActive(ON))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else if(gomp_StringMatch(Text2 , "off")) {
                    if(!gomp_SetIdentifyAtomActive(OFF))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else {
                    gomp_PrintERROR("option has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
            }
/* atom connectivity */
            else if(gomp_StringMatch(Text1 , "reco$nnectivity")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(gomp_StringMatch(Text2 , "on")) {
                    if(!gomp_SetBondReconnectivityState(ON))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else if(gomp_StringMatch(Text2 , "off")) {
                    if(!gomp_SetBondReconnectivityState(OFF))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else {
                    gomp_PrintERROR("option has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
            }
/* atom hydrogen bond connectivity */
            else if(gomp_StringMatch(Text1 , "hbre$connectivity")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(gomp_StringMatch(Text2 , "on")) {
                    if(!gomp_SetHBondReconnectivityState(ON))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else if(gomp_StringMatch(Text2 , "off")) {
                    if(!gomp_SetHBondReconnectivityState(OFF))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else {
                    gomp_PrintERROR("option has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
            }
/* max atom connectivity */
            else if(gomp_StringMatch(Text1 , "maxc$onnectivity")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("atom connectivity max value missing");
                    return(TCL_ERROR);
                }
                ITemp = atoi(Text2);

                if(!gomp_SetMaxAtomConnections(ITemp))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
/* covar radius */
            else if(gomp_StringMatch(Text1 , "cova$r")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(gomp_ParseSetAtomCovar(Text2,Text3))
                    return(TCL_ERROR);
                else
                    return(TCL_OK);
            }
/* atom window */
            else if(gomp_StringMatch(Text1 , "windo$w")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("search window value missing");
                    return(TCL_ERROR);
                }
                if(gomp_StringMatch(Text2 , "max")) {
                    { int ii;
                      ITemp = 0;
                      for(ii = 0 ; ii < gomp_GetNumMolecStructs() ; ii++) {
                         if(gomp_GetNumAtomsInMolecStruct(ii) > ITemp)
                             ITemp = gomp_GetNumAtomsInMolecStruct(ii);
                      }
/* if no system defined can't set the value from max number of atoms */
                      if(!ITemp) return(TCL_ERROR);
                    }
                } else {
                  ITemp = atoi(Text2);
                }
                if(!gomp_SetSearchWindow(ITemp))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
/* charge */
            else if(gomp_StringMatch(Text1 , "charg$e")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("atom index is missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text3) == 0) {
                    ITemp1 = 1;
                }
                else
                    ITemp1 = atoi(Text3);
                if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure number out of range");
                    return(TCL_ERROR);
                }
                ITemp1--;
                ITemp = atoi(Text2);
                if(ITemp < 1 || ITemp > gomp_GetNumAtomsInMolecStruct(ITemp1)) {
                    gomp_PrintERROR("atom number is out of range");
                    return(TCL_ERROR);
                }
                ITemp--;
                (void)gomp_PutAtomCharge(ITemp1 , atof(Text1) , ITemp);
                return(TCL_OK);
            }
/* VDW */
            else if(gomp_StringMatch(Text1 , "vdw")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("atom index is missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text3) == 0) {
                    ITemp1 = 1;
                }
                else
                    ITemp1 = atoi(Text3);
                if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure number out of range");
                    return(TCL_ERROR);
                }
                ITemp1--;
                ITemp = atoi(Text2);
                if(ITemp < 1 || ITemp > gomp_GetNumAtomsInMolecStruct(ITemp1)) {
                    gomp_PrintERROR("atom number is out of range");
                    return(TCL_ERROR);
                }
                ITemp--;
                (void)gomp_PutAtomVdwRad(ITemp1 , atof(Text1) , ITemp);
                return(TCL_OK);
            }
/* atom colour */
            else if(gomp_StringMatch(Text1 , "colo$r") ||
                    gomp_StringMatch(Text1 , "colo$ur")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("atom index is missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text3) == 0) {
                    ITemp1 = 1;
                }
                else
                    ITemp1 = atoi(Text3);
                if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure number out of range");
                    return(TCL_ERROR);
                }
                ITemp1--;
                ITemp = atoi(Text2);
                if(ITemp < 1 || ITemp > gomp_GetNumAtomsInMolecStruct(ITemp1)) {
                    gomp_PrintERROR("atom number is out of range");
                    return(TCL_ERROR);
                }
                ITemp--;
                if(gomp_ColourName2RGB(Text1 , &RedC , &GreenC , &BlueC)) {
                    gomp_PrintERROR("can't convert colour name/value");
                    return(TCL_ERROR);
                }
                gomp_PutAtomColour(ITemp1, RedC, GreenC, BlueC, ITemp);
                return(TCL_OK);
            }
/* residue number */
            else if(gomp_StringMatch(Text1 , "resn$umber")) {
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("atom index is missing");
                    return(TCL_ERROR);
                }
                if(strlen(Text3) == 0) {
                    ITemp1 = 1;
                }
                else
                    ITemp1 = atoi(Text3);
                if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure number out of range");
                    return(TCL_ERROR);
                }
                ITemp1--;
                ITemp = atoi(Text2);
                if(ITemp < 1 || ITemp > gomp_GetNumAtomsInMolecStruct(ITemp1)) {
                    gomp_PrintERROR("atom number is out of range");
                    return(TCL_ERROR);
                }
                ITemp--;
                (void)gomp_PutAtomResNum1(ITemp1 , atoi(Text1) , ITemp);
                (void)gomp_PutAtomResNum2(ITemp1 , atoi(Text1) , ITemp);
                return(TCL_OK);
            }
        }
/* atom ...... */
/* residue */
        else if(gomp_StringMatch(Text , "resi$due")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "labe$l") ||
               gomp_StringMatch(Text1 , "symb$ol")) {
                gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text) == 0) {
                    gomp_PrintERROR("residue name is missing");
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("residue index is missing");
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text2) == 0) {
                    ITemp1 = 1;
                }
                else
                    ITemp1 = atoi(Text2);
                if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure number out of range");
                    return(TCL_ERROR);
                }
                ITemp1--;
                ITemp = atoi(Text1);
                if(ITemp < 1 || ITemp > gomp_GetNumAtomsInMolecStruct(ITemp1)) {
                    gomp_PrintERROR("residue number is out of range");
                    return(TCL_ERROR);
                }
                ITemp--;
                (void)gomp_PutAtomResName(ITemp1 , Text , ITemp);
                return(TCL_OK);
            }
        }
/* residue ...... */
/* segment */
        else if(gomp_StringMatch(Text , "segm$ent")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "labe$l") ||
               gomp_StringMatch(Text1 , "symb$ol")) {
                gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text) == 0) {
                    gomp_PrintERROR("segment name is missing");
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text1) == 0) {
                    gomp_PrintERROR("segment index is missing");
                    return(TCL_ERROR);
                }
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text2) == 0) {
                    ITemp1 = 1;
                }
                else
                    ITemp1 = atoi(Text2);
                if(ITemp1 < 1 || ITemp1 > gomp_GetNumMolecStructs()) {
                    gomp_PrintERROR("structure number out of range");
                    return(TCL_ERROR);
                }
                ITemp1--;
                ITemp = atoi(Text1);
                if(ITemp < 1 || ITemp > gomp_GetNumAtomsInMolecStruct(ITemp1)) {
                    gomp_PrintERROR("atom number is out of range");
                    return(TCL_ERROR);
                }
                ITemp--;
                (void)gomp_PutAtomSegName(ITemp1 , Text , ITemp);
                return(TCL_OK);
            }
        }
/* stereo pair ...... */
        else if(gomp_StringMatch(Text , "spai$r")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

            if(gomp_StringMatch(Text1 , "on")) {
                (void)gomp_SetStereoPlotState(ON);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "off")) {
                (void)gomp_SetStereoPlotState(OFF);
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "dist$ance")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text2) != 0) 
                    (void)gomp_SetStereoPlotTranslate(atof(Text2));
                else {
                    gomp_PrintERROR("translate value missing");
                    return(TCL_ERROR);
                }
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "angl$e")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text2) != 0) 
                    (void)gomp_SetStereoPlotAngle(atof(Text2));
                else {
                    gomp_PrintERROR("angle value missing");
                    return(TCL_ERROR);
                }
                return(TCL_OK);
            }
        }
/* QUAD stereo ...... */
        else if(gomp_StringMatch(Text , "quad$stereo")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
      
            if(gomp_StringMatch(Text1 , "on")) {
                gomp_QuadStereoOn();
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "off")) {
                gomp_QuadStereoOff();
                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "ang$le")) {
                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text2) != 0) 
                    gomp_SetQuadStereoHalfAngle(0.5*atof(Text2));
                else {
                    gomp_PrintERROR("angle value missing");
                    return(TCL_ERROR);
                }
                return(TCL_OK);
            }
        }
        else if(gomp_StringMatch(Text , "cell")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "dime$nsions")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text4,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                gomp_SetCellA(atof(Text2));
                gomp_SetCellB(atof(Text3));
                gomp_SetCellC(atof(Text4));

                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "angl$es")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text4,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                gomp_SetCellAlpha(atof(Text2));
                gomp_SetCellBeta(atof(Text3));
                gomp_SetCellGamma(atof(Text4));

                return(TCL_OK);
            }
            else if(gomp_StringMatch(Text1 , "line$width")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("line width value missing");
                    return(TCL_ERROR);
                }
                if(!gomp_SetCellLinewidth(atoi(Text2)))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "colo$urs") ||
                    gomp_StringMatch(Text1 , "colo$rs")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("colour name or code missing");
                    return(TCL_ERROR);
                }
                {
                    float Red;
                    float Green;
                    float Blue;

                    if(gomp_ColourName2RGB(Text2, &Red , &Green , &Blue)) {
                        return(TCL_ERROR);
                    }
                    if(!gomp_SetCellColour(Red, Green,Blue))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            }
            else if(gomp_StringMatch(Text1 , "tran$slation") ||
                    gomp_StringMatch(Text1 , "plac$e")) {
                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                strncpy(Text4,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                (void)gomp_SetCellXtrans(atof(Text2));
                (void)gomp_SetCellYtrans(atof(Text3));
                (void)gomp_SetCellZtrans(atof(Text4));

                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("error in 'define cell' command");
                return(TCL_ERROR);
            }
        }
/* rotation */
        else if(gomp_StringMatch(Text , "rota$tion")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            if(gomp_StringMatch(Text1 , "stat$e")) {

                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(gomp_StringMatch(Text2 , "on")) {
                    if(!gomp_SetRotationState(ON)) {
                        Value = Tcl_SetVar(gomp_GetTclInterp(),"gomManipulationState","0",TCL_GLOBAL_ONLY);
                        if(!Value) {
                            gomp_PrintERROR("can't set tcl variable 'gomManipulationState'");
                            return(TCL_ERROR);
                        }
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else if(gomp_StringMatch(Text2 , "off")) {
                    if(!gomp_SetRotationState(OFF)) {
                        Value = Tcl_SetVar(gomp_GetTclInterp(),"gomManipulationState","1",TCL_GLOBAL_ONLY);
                        if(!Value) {
                            gomp_PrintERROR("can't set tcl variable 'gomManipulationState'");
                            return(TCL_ERROR);
                        }
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else {
                    gomp_PrintERROR("option has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
            }
            else {
                gomp_PrintERROR("unrecognized 'define rotation' command");
                return(TCL_ERROR);
            }
        }
/* translation */
        else if(gomp_StringMatch(Text , "tran$slation")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            if(gomp_StringMatch(Text1 , "stat$e")) {

                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(gomp_StringMatch(Text2 , "on")) {
                    if(!gomp_SetTranslationState(ON)) {
                        Value = Tcl_SetVar(gomp_GetTclInterp(),"gomManipulationState","1",TCL_GLOBAL_ONLY);
                        if(!Value) {
                            gomp_PrintERROR("can't set tcl variable 'gomManipulationState'");
                            return(TCL_ERROR);
                        }
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else if(gomp_StringMatch(Text2 , "off")) {
                    if(!gomp_SetTranslationState(OFF)) {
                        Value = Tcl_SetVar(gomp_GetTclInterp(),"gomManipulationState","0",TCL_GLOBAL_ONLY);
                        if(!Value) {
                            gomp_PrintERROR("can't set tcl variable 'gomManipulationState'");
                            return(TCL_ERROR);
                        }
                        return(TCL_OK);
                    }
                    else {
                        return(TCL_ERROR);
                    }
                }
                else {
                    gomp_PrintERROR("option has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
            }
            else {
                gomp_PrintERROR("unrecognized 'define translation' command");
                return(TCL_ERROR);
            }
        }
/* gromos96 */
        else if(gomp_StringMatch(Text , "gromos96")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
            if(gomp_StringMatch(Text1 , "coorda$mplifier")) {

                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("coordinate amplifier value missing");
                    return(TCL_ERROR);
                }

                (void)gomp_SetGROMOS96CoordAmplifier(atof(Text2));
                return(TCL_OK);
            }
            else {
                gomp_PrintERROR("unrecognized 'define gromos96' command");
                return(TCL_ERROR);
            }
        }
        else if(gomp_StringMatch(Text , "bondst$yle")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);

            if(gomp_StringMatch(Text1 , "smoo$th")) {
                if(!gomp_SetBondDisplayStyle(0))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "half") ||
                    gomp_StringMatch(Text1 , "defa$ult")) {
                if(!gomp_SetBondDisplayStyle(1))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("option has to be 'smooth' or 'half'");
                return(TCL_ERROR);
            }
        }
/* system */
        else if(gomp_StringMatch(Text , "syst$em")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(gomp_StringMatch(Text1 , "tran$slation")) {
                gomp_CopyString(Text,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                if(strlen(Text) == 0) {
                    gomp_PrintERROR("translation state missing (on/off)");
                    return(TCL_ERROR);
                }
                if(gomp_StringMatch(Text , "on")) {
                    if(!gomp_SetSystemTranslateState(1))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else if(gomp_StringMatch(Text , "off")) {
                    if(!gomp_SetSystemTranslateState(0))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
                else {
                    gomp_PrintERROR("option has to be 'on' or 'off'");
                    return(TCL_ERROR);
                }
            }
        }
/* light mapping type */
        else if(gomp_StringMatch(Text , "colorm$apping") ||
                gomp_StringMatch(Text , "colourm$apping")) {
            strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
          
            if(gomp_StringMatch(Text2 , "text$ure")) {
                if(!gomp_SetColorMappingType(COLOR_MAPPING_TYPE_TEXTURE))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text2 , "rain$bow")) {
                if(!gomp_SetColorMappingType(COLOR_MAPPING_TYPE_RAINBOW))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("option has to be 'texture' or 'rainbow'");
                return(TCL_ERROR);
            }
        }
/* light position */
        else if(gomp_StringMatch(Text , "ligh$t")) {

            strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);
          
            if(gomp_StringMatch(Text2 , "posi$tion")) {

                strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(gomp_StringMatch(Text3 , "c")) {
                    if(!gomp_SetLightPosition(0))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                } else if(gomp_StringMatch(Text3 , "nw")) {
                    if(!gomp_SetLightPosition(8))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                } else if(gomp_StringMatch(Text3 , "ne")) {
                    if(!gomp_SetLightPosition(2))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                } else if(gomp_StringMatch(Text3 , "n")) {
                    if(!gomp_SetLightPosition(1))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                } else if(gomp_StringMatch(Text3 , "w")) {
                    if(!gomp_SetLightPosition(7))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                } else if(gomp_StringMatch(Text3 , "e")) {
                    if(!gomp_SetLightPosition(3))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                } else if(gomp_StringMatch(Text3 , "sw")) {
                    if(!gomp_SetLightPosition(6))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                } else if(gomp_StringMatch(Text3 , "se")) {
                    if(!gomp_SetLightPosition(4))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                } else if(gomp_StringMatch(Text3 , "s")) {
                    if(!gomp_SetLightPosition(5))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
            } else if(gomp_StringMatch(Text2 , "diff$use")) {

                gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
         
/* red */
                if(gomp_StringMatch(Text2 , "red")) {
  
                    strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                            BUFF_LEN-1);
     
                    if(strlen(Text3) == 0) {
                        gomp_PrintERROR("red light diffuse value missing");
                        return(TCL_ERROR);
                    }
                    if(!gomp_SetLightDiffuseRed(atof(Text3)))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);
                }
/* green */
                else if(gomp_StringMatch(Text2 , "gree$n")) {

                    strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                            BUFF_LEN-1);
     
                    if(strlen(Text3) == 0) {
                        gomp_PrintERROR("green light diffuse value missing");
                        return(TCL_ERROR);
                    }
                    if(!gomp_SetLightDiffuseGreen(atof(Text3)))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);

                }
/* blue */
                else if(gomp_StringMatch(Text2 , "blue")) {

                    strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                            BUFF_LEN-1);
     
                    if(strlen(Text3) == 0) {
                        gomp_PrintERROR("blue light diffuse value missing");
                        return(TCL_ERROR);
                    }
                    if(!gomp_SetLightDiffuseBlue(atof(Text3)))
                        return(TCL_OK);
                    else
                        return(TCL_ERROR);

                } else {
                    gomp_PrintERROR("option has to be 'red, green or blue'");
                    return(TCL_ERROR);
                }
            }
            else {
                gomp_PrintERROR("option has to be 'position or diffuse'");
                return(TCL_ERROR);
            }
        }
/* define trajectory */
        else if(gomp_StringMatch(Text , "traj$ectory")) {
            strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);

            if(gomp_StringMatch(Text2 , "fid")) {
                strncpy(Text3,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);

                if(gomp_StringMatch(Text3 , "on")) {
                    (void)gomp_SetDisplayRunningFrameNumberState(ON);
                    return(TCL_OK);
                } else if(gomp_StringMatch(Text3 , "off")) {
                    (void)gomp_SetDisplayRunningFrameNumberState(OFF);
                    return(TCL_OK);
                } else {
                    gomp_PrintERROR("command 'define trajectory fid' has to include 'on' or 'off'");
                    return(TCL_ERROR);
                }
            }
        }
        else if(gomp_StringMatch(Text , "cutp$lane")) {
            strncpy(Text1,gomp_GetNextFromParserStack(argc,NULL),
                    BUFF_LEN-1);

            if(gomp_StringMatch(Text1 , "damp$ing")) {

                strncpy(Text2,gomp_GetNextFromParserStack(argc,NULL),
                        BUFF_LEN-1);
                if(strlen(Text2) == 0) {
                    gomp_PrintERROR("damping factor is missing");
                    return(TCL_ERROR);
                }
                (void)gomp_SetCutPlaneDamping(atof(Text2));
                return(TCL_OK);
            } else {
                gomp_PrintERROR("undefined option in 'define cutplane'");
                return(TCL_ERROR);
            }        
        }
/* define redisplay mode fast/slow */
        else if(gomp_StringMatch(Text , "redi$splay")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("redisplay mode has to be (fast/slow)");
                return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text1 , "fast")) {
                if(!gomp_SetSystemRedisplayMode(0))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "slow")) {
                if(!gomp_SetSystemRedisplayMode(1))
                    return(TCL_OK);
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("option has to be 'fast' or 'slow'");
                return(TCL_ERROR);
            }
        }
/* define global/local scaling */
        else if(gomp_StringMatch(Text , "scal$ing")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("scaling parameter has to be local/global");
                return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text1 , "glob$al")) {
                if(!gomp_SetAllowIndividualScaling(0)) {
                    return(TCL_OK);
                }
                else {
                    return(TCL_ERROR);
                }
            }
            else if(gomp_StringMatch(Text1 , "loca$l")) {
                if(!gomp_SetAllowIndividualScaling(1)) {
                    return(TCL_OK);
                }
                else {
                    return(TCL_ERROR);
                }
            }
            else {
                gomp_PrintERROR("option has to be 'global' or 'local'");
                return(TCL_ERROR);
            }
        }
/* define global/local manipulation */
        else if(gomp_StringMatch(Text , "transf$ormation")) {
            gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
            if(strlen(Text1) == 0) {
                gomp_PrintERROR("manipulation parameter has to be local/global");
                return(TCL_ERROR);
            }
            if(gomp_StringMatch(Text1 , "glob$al")) {
/* its already in global mode */
                if(!gomp_GetObjectCenterType()) return(TCL_OK);
                if(!gomp_SetObjectCenterType(0)) {

                    Value = gomp_SetGlobalTransformationState();

                    if(!Value) {
                        gomp_InputView = 0;
                    } else {
                        if(atoi(Value) > 0) {
                            gomp_InputView = 1;
                        } else {
                            gomp_InputView = 0;
                        }
                    }

                    if(gomp_InputView) {
                        if(gomp_DoViewingTransformationOverStructures( -1 )) 
                            return(TCL_ERROR);
                    }

                    if(gomp_ResetView()) {
                        gomp_PrintERROR("cant't reset view (global), the result will be unreliable");
                    }
                    return(TCL_OK);
                }
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "loca$l")) {
/* its already in local mode */
                if(gomp_GetObjectCenterType()) return(TCL_OK);
                if(!gomp_SetObjectCenterType(1)) {

                    Value = gomp_SetLocalTransformationState();

                    if(!Value) {
                        gomp_InputView = 0;
                    } else {
                        if(atoi(Value) > 0) {
                            gomp_InputView = 1;
                        } else {
                            gomp_InputView = 0;
                        }
                    }

                    if(gomp_InputView) {
                        if(gomp_DoViewingTransformationOverStructures( -1 )) 
                            return(TCL_ERROR);
                    }

                    if(gomp_ResetView()) {
                        gomp_PrintERROR("cant't reset view (local), the result will be unreliable");
                    }
                    return(TCL_OK);
                }
                else
                    return(TCL_ERROR);
            }
            else if(gomp_StringMatch(Text1 , "user")) {
/* its already in local mode */
                if(gomp_GetObjectCenterType()) return(TCL_OK);
                if(!gomp_SetObjectCenterType(1)) {
                    Value = Tcl_SetVar(gomp_GetTclInterp(),"lulTraceTransformationValue","2",TCL_GLOBAL_ONLY);
                    if(!Value) {
                        gomp_PrintERROR("can't set tcl variable 'lulTraceTransformationValue'");
                        return(TCL_ERROR);
                    }

                    if(gomp_ResetView()) {
                        gomp_PrintERROR("cant't reset view (local), the result will be unreliable");
                    }

                    gomp_CopyString(Text1,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                    gomp_CopyString(Text2,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                    gomp_CopyString(Text3,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);
                    gomp_CopyString(Text4,gomp_GetNextFromParserStack(argc,NULL),BUFF_LEN);

                    ITemp  = atoi(Text1);
                    FTemp1 = atof(Text2);
                    FTemp2 = atof(Text3);
                    FTemp3 = atof(Text4);

                    if(ITemp < 1 || ITemp > gomp_GetNumMolecStructs()) {

                        gomp_PrintERROR("structure index out of allowe range");
                        return(TCL_ERROR);
                    }

                    sprintf(Text,"Will apply a MT translation (x,y,z) for structure #%d: %f %f %f",
                            ITemp,FTemp1,FTemp2,FTemp3);
                    gomp_PrintMessage(Text);

                    (void)gomp_SaveTranslateArrayMT(   ITemp - 1 , FTemp1 , FTemp2 , FTemp3); 

                    return(TCL_OK);
                }
                else
                    return(TCL_ERROR);
            }
            else {
                gomp_PrintERROR("option has to be 'global', 'local' or 'user'");
                return(TCL_ERROR);
            }
        }
        else {
            gomp_PrintERROR("'define' command not recognized");
            return(TCL_ERROR);
        }
    }
/*  E R R O R command not recognized         */
    gomp_PrintERROR("'define' command not recognized");

    return(TCL_ERROR);

}

/*********************************************************************/
const char *gomp_SetGlobalTransformationState()
/*********************************************************************/
{
    static const char *Value;

    Value = Tcl_SetVar(gomp_GetTclInterp(),"lulTraceTransformationValue","0",TCL_GLOBAL_ONLY);
    if(!Value) {
        gomp_PrintERROR("can't set tcl variable 'lulTraceTransformationValue'");
        return(NULL);
    }
    Value  = Tcl_GetVar(gomp_GetTclInterp() , 
                                "lulDoViewTransGlobalLocal", 
                                TCL_GLOBAL_ONLY);

    return(Value);
}
/*********************************************************************/
const char *gomp_SetLocalTransformationState()
/*********************************************************************/
{
    static const char *Value;

    Value = Tcl_SetVar(gomp_GetTclInterp(),"lulTraceTransformationValue","1",TCL_GLOBAL_ONLY);
    if(!Value) {
        gomp_PrintERROR("can't set tcl variable 'lulTraceTransformationValue'");
        return(NULL);
    }

    Value  = Tcl_GetVar(gomp_GetTclInterp() , 
                                "lulDoViewTransGlobalLocal", 
                                TCL_GLOBAL_ONLY);

    return(Value);
}
