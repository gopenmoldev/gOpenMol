#########################################################################
# ****************** graphical user interface ***************************
#                    g O p e n M o l 
# This is the Tcl/Tk graphical user interface part to gOpenMol
#
# Copyright Leif Laaksonen  2004 version 2.31
#
# Date:  2004-03-24
#
#########################################################################
#
#
#####################################
# Command stack
# Variable is gomCommandStack

set gomCommandStackDeepMax 30
set gomCommandStackPeek     0
set gomCommandStackFill     0
set gomCommandStackNext     0
set gomCommandStackPrev     0
set gomCommandStackCounter  0

#####################################
#
# Function keys
#
#####################################
set F1  "display;# display scene"
set F2  "window 1 fullscreen;# use fullscreen mode"
set F3  "window 1 resize 500 500;window 1 move 10 20;# return from fullscreen mode to a 500x500 window"
set F4  ""
set F5  ""
set F6  ""
set F7  ""
set F8  ""
set F9  ""
set F10 "htmlShowHelp gOpenMol_begin.html#COMMANDS;# help on commands"
set F11 "htmlShowHelp gopenmol.html;# main help"
set F12 "exit;# exit gOpenMol"

# announce it
gomPrint "Function keys 1...12"
gomPrint "F1:  $F1"
gomPrint "F2:  $F2"
gomPrint "F3:  $F3"
gomPrint "F4:  $F4"
gomPrint "F5:  $F5"
gomPrint "F6:  $F6"
gomPrint "F7:  $F7"
gomPrint "F8:  $F8"
gomPrint "F9:  $F9"
gomPrint "F10: $F10"
gomPrint "F11: $F11"
gomPrint "F12: $F12"

#####################################
# setup default parameters
# set default gOpenMol file name
#####################################
set gomFileName             ""
set gomNumStructures        0
set gomAtomColour           #ffffff
set gomControlFont8         {Arial 8 bold}
set gomControlFont          {Arial 10 bold}
set gomControlFont12        {Arial 12 bold}
set gomControlFontHuge      {Arial 20 bold}
set gomControlFontFixed     {Courier 10 bold}
font create gomfont8    -family arial -size 8  -weight bold
font create gomfont     -family arial -size 10 -weight bold
font create gomfont12   -family arial -size 12 -weight bold
font create gomfonthuge -family arial -size 20 -weight bold
# default rotation step in degrees
set gomDefaultRotaStep      5.0

# default translation step in Ångströms
set gomDefaultTranStep      1.0

# cut plane animation state
set gomCutplaneLoopState    0

# default values for licorice sphere and cylinder
set   gomLicoSphere   0.3
set   gomLicoCylinder 0.3

#cylinder arrow control (hat/cylinder radius ratio and
#                        length/hat length ratio 
#set gomCylinderArrowControl "1.3 0.7"

#
# Socket server state
#
set gomServerState 0
set gomSocketPort  2309

#   variable for fchk plugin (controls popup of contour widget)
set lulContourControlStartup 1
#
# Most grahics cards are able to give good OpenGL performance
# so the default redisplay type is show all graphics objects
# while rotating/translating molecule
# 2003-06-09/LUL
#
define redisplay slow
#####################################
# hydrogen bond defaults
#####################################
set gomHbondColour [lulColourFloat2Hex 0.0 1.0 1.0]

#####################################
# gradient display default color
#
#set gomDefaultGradienColor {1.0 0.0 0.0}
#set gomVectorDefaultType  0        ;# line display mode
#set gomVectorDefaultType "1 0.1"   ;# licorice display mode with a radius of 0.1
set gomCurrentVectorDisplayStatus 0 ;# no current display

##############################
# GLUT bitmap fonts
#
# BITMAP_8_BY_13
# BITMAP_9_BY_15
# BITMAP_TIMES_ROMAN_10
# BITMAP_TIMES_ROMAN_24
# BITMAP_HELVETICA_10
# BITMAP_HELVETICA_12
# BITMAP_HELVETICA_18
#
##############################
  
# atom label font
set gomAtomLabelFont        BITMAP_HELVETICA_12
#
# i is structure number and j is running atom index number
#
set gomAtomLabelStringDef      {_[show atom atomname $j $i]}
set gomAtomLabelStringExt      {_[show atom segment  $j $i]:[show atom resnumber $j $i]:[show atom atomname $j $i]($j)}
set gomAtomLabelNumberLimit    200

# text font
set gomTextFont BITMAP_HELVETICA_12

# monitor distance/angle/torsion precision
set gomMonitorDistancePrecision "%.2f"
set gomMonitorAnglePrecision    "%.2f"
set gomMonitorTorsionPrecision  "%.2f"

# no viewing transformation when switching between Global/Local transformation
set lulDoViewTransGlobalLocal 0

# viewing transformation is applied to "export input" coordinates
set lulInputView 0
 
# default export program input type
set gomExportInputType probesurf

set lulGlobalEnteringWidgetValue ""

set gomDefaultRMSDmagnifier 1.0

# switch to control the calculation of connectivity when a new file is readinto gOpenMol
# options are yes/no
set lulCalculateAtomConnectivity yes

# the help files are here ...
set gomHelpDir [show helpdirectory]
# value here has to match with the next line (manual == 0 , automatic == 1)
set gomInstantDisplay 1
# manual update of window
define window update automatic
# used as a global to exchange an URL to gOpenMol
set gomURLFileName ""
# put system translation on as default
define system translation on

######################################################################################

# default trajectory file type is defined here!
set gomTrajType            "charmm"
set gomTrajFileName ""
set gomTrajTypes "amber cerius2 charmm discover famber frame gromos gromos96a gromos96b hyperchem mumod xmol xplor yasp"
set gomTrajExt(amber)      "*"
set gomTrajExt(cerius2)    "*"
set gomTrajExt(charmm)     ".dcd .DCD"
set gomTrajExt(discover)   "*"
set gomTrajExt(dl_poly)    "*"
set gomTrajExt(famber)     "*"
set gomTrajExt(Gromacs)    ".trj .trr .xtc .TRJ .TRR .XTC"
set gomTrajExt(Gromos)     "*"
set gomTrajExt(Gromos96A)  ".dat"
set gomTrajExt(Gromos96B)  ".dat"
set gomTrajExt(hyperchem)  ".snp .SNP"
set gomTrajExt(mumod)      "*"
set gomTrajExt(tinker)     ".arc .ARC"
#set gomTrajExt(openmol)
set gomTrajExt(xmol)       ".xmol .XMOL"
set gomTrajExt(xplor)      ".dcd .DCD"
set gomTrajExt(yasp)       "*"
# default value (change this to change the coordinates and colour)
#set gomCurrentFrameProperty {0.75 0.9 red BITMAP_HELVETICA_12}
# display loop sleep value
#set gomDisplayLoopSleep 0
set gomHardcopyType            ""
set gomHardcopyExt(bmp)        ".bmp .BMP"
set gomHardcopyExt(jpeg)       ".jpg .JPG"
set gomHardcopyExt(postscript) ".ps .PS"
set gomHardcopyExt(rgb)        ".rgb .RGB"
set gomHardcopyExt(tga)        ".tga .TGA"
set gomHardcopyExt(xwd)        ".xwd .XWD"

# default JPEG quality
#set gomJPEGquality 75


set gomCellColour          #ff0000

# address to gOpenMol Home Page
set WebHome [file join $gomEnv(GOM_HELP) gopenmol.html]

# working scalars...
set gomAppend               0
set gomFileSize             ""
set gomNumAtoms             ""
set gomSwitch               0
set gomPickSwitch           0
set gomStructure            0
set gomCutX                 0
set gomCutY                 0
set gomCutZ                 0
set gomSwitchDistMon        0
set gomSwitchAngMon         0
set gomSwitchTorsMon        0
set gomAppendDistList       0
set gomAppendAngList        0
set gomAppendTorsList       0

##################################################################
# help files

set gomHelpFile(atomtype)            "atom_type_widget.html"
set gomHelpFile(atomtree)            "atom_tree_widget.html"
set gomHelpFile(atomdisplaymask)     "atom_display_mask_widget.html"
set gomHelpFile(atomlabeldisplay)    "atom_label_selection_widget.html"
set gomHelpFile(calccorrelation)     "calc_correlation_widget.html"
set gomHelpFile(calcatomconn)        "calc_atom_conn.html"
set gomHelpFile(centersystem)        "center_widget.html"
set gomHelpFile(colourselection)     "colour_selection_widget.html"
set gomHelpFile(contmapping)         "contour_mapping_widget.html"
set gomHelpFile(contglueing)         "contour_glueing_widget.html"
set gomHelpFile(contour_details)     "contour_details_widget.html"
set gomHelpFile(contourplane)        "contour_plane_widget.html"
set gomHelpFile(contourxyzplane)     "contour_xyzplane_widget.html"
set gomHelpFile(contourclipplane)    "contour_clip_plane_widget.html"
set gomHelpFile(contourxyzclipplane) "contour_xyzclip_plane_widget.html"
set gomHelpFile(contplane_animate)   "cont_plane_anim_widget.html"
set gomHelpFile(copytimeseries)      "copy_timeseries_widget.html"
set gomHelpFile(copytimeseries2disk) "copy_timeseries2disk_widget.html"
set gomHelpFile(cutplane_plottype)   "cutplane_type_widget.html"
set gomHelpFile(editmol)             "edit_mol_widget.html"
set gomHelpFile(definecell)          "define_cell_widget.html"
set gomHelpFile(displayprops)        "display_props_widget.html"
set gomHelpFile(editdisplay)         "edit_disp_attr_widget.html"
set gomHelpFile(editlightprop)       "edit_light_widget.html"
set gomHelpFile(editmatprop)         "edit_material_widget.html"
set gomHelpFile(elementinfo)         "element_info_widget.html"
set gomHelpFile(exportcoordinates)   "export_coordinates_widget.html"
set gomHelpFile(exportinput)         "export_input_widget.html"
set gomHelpFile(hardcopy)            "hardcopy_widget.html"
set gomHelpFile(hbondcriteria)       "hbond_criteria_widget.html"
set gomHelpFile(hydrogenbonds)       "hydrogen_bonds_widget.html"
set gomHelpFile(identifyatom)        "identify_atom_widget.html"
set gomHelpFile(joingamessirc)       "join_gamess_irc_widget.html"
set gomHelpFile(importcoordinates)   "import_coordinates_widget.html"
set gomHelpFile(importdictionary)    "import_dictionary_widget.html"
set gomHelpFile(importgbasis)        "import_gbasis_widget.html"
set gomHelpFile(importmodel)         "import_model_widget.html"
set gomHelpFile(importvector)        "import_vector_widget.html"
set gomHelpFile(maintrajectory)      "main_trajectory_widget.html"
set gomHelpFile(makemovie)           "make_animation_widget.html"
set gomHelpFile(manipulate)          "manipulate_widget.html"
set gomHelpFile(mantimeseries)       "mantimeseries_widget.html"
set gomHelpFile(meansqrdisp)         "msd_widget.html"
set gomHelpFile(measuregeom)         "measuregeom_widget.html"
set gomHelpFile(monitordistance)     "mondist_widget.html"
set gomHelpFile(monitorangle)        "monang_widget.html"
set gomHelpFile(monitortorsion)      "montors_widget.html"
set gomHelpFile(plotldp)             "ldp_widget.html"
set gomHelpFile(periodictable)       "periodic_table_widget.html"
set gomHelpFile(plotaxis)            "plot_axis_widget.html"
set gomHelpFile(plotcolourscale)     "plot_cscale_widget.html"
set gomHelpFile(plotrmsd)            "plot_rmsd_widget.html"
set gomHelpFile(plotstereopair)      "plot_stereo_pair_widget.html"
set gomHelpFile(plotvectorfile)      "plot_vector_file_widget.html"
set gomHelpFile(quadstereo)          "plot_quadstereo_widget.html"
set gomHelpFile(plumber)             "plumber_widget.html"
set gomHelpFile(plumberlist)         "plumber_list_widget.html"
set gomHelpFile(plumberseco)         "plumber_seco_widget.html"
set gomHelpFile(rottranscontrol)     "rotate_widget.html"
set gomHelpFile(runcontman)          "runcontman_widget.html"
set gomHelpFile(rungcube2plt)        "rungcube2plt_widget.html"
set gomHelpFile(runjaguar2plt)       "runjaguar2plt_widget.html"
set gomHelpFile(runpltfile)          "runpltfile_widget.html"
set gomHelpFile(runprobesurf)        "runprobesurf_widget.html"
set gomHelpFile(runtclscript)        "runtclscript_widget.html"
set gomHelpFile(runkont2plt)         "runkont2plt_widget.html"
set gomHelpFile(rungamess2plt)       "rungamess2plt_widget.html"
set gomHelpFile(runautodock2plt)     "runautodock2plt_widget.html"
set gomHelpFile(runuhbd2plt)         "runuhbd2plt_widget.html"
set gomHelpFile(runturbomole2plt)    "runturbomole2plt_widget.html"
set gomHelpFile(runxvibs)            "runxvibs_widget.html"
set gomHelpFile(selectmask)          "selmask_widget.html"
set gomHelpFile(selectatom)          "select_widget.html"
set gomHelpFile(selclusteratom)      "calculate_cluster_widget.html"
set gomHelpFile(selrdfatom)          "calculate_rdf_widget.html"
set gomHelpFile(simpose)             "simpose_widget.html"
set gomHelpFile(traceatoms)          "trace_driver_widget.html"
set gomHelpFile(writecluster)        "write_cluster_widget.html"
set gomHelpFile(writeinput)          "write_input_widget.html"
set gomHelpFile(writemodel)          "write_model_widget.html"
set gomHelpFile(writermsd)           "write_rmsd_widget.html"
set gomHelpFile(writetimes)          "write_times_widget.html"


##################################################################

# window manager........................
set platform  $tcl_platform(platform)
if {$platform == "unix"} {
wm geometry . 520x500+520+5
} else {
wm geometry . 510x500+510+5
}
wm title . "gOpenMol Tcl/Tk interface (version [format "%.2f" [expr [string trim [show version]]/100.]])"
# if this window is closed close down the whole program!
wm protocol . WM_DELETE_WINDOW  {catch {quit}}
#.......................................
# Put the appname to be used with the send command
tk appname gopenmol
#
# Menu bar
frame .frame -borderwidth 2 -relief raised -bd 2
pack  .frame -side top   -fill x

frame .dummy -width 15c -height 20c
pack  .dummy -side top   -fill both
frame .dummy.left -borderwidth 2 -relief flat -width 8c -height 20c
pack  .dummy.left   -side left  -fill y
# 
label   .dummy.left.gomlabel  -text "gOpenMol version $gomVersion" \
        -font {Helvetica 14 bold}
pack .dummy.left.gomlabel  -side top
#
frame .dummy.left.structures -borderwidth 2 -relief ridge -bd 2
pack  .dummy.left.structures -side top -anchor w
text  .dummy.left.structures.text \
       -yscrollcommand ".dummy.left.structures.scrollbar set"  \
       -xscrollcommand ".dummy.left.structures.scrollbarx set" \
       -width 20 -height 5 -bg gray
#
scrollbar .dummy.left.structures.scrollbar \
           -command ".dummy.left.structures.text yview"
scrollbar .dummy.left.structures.scrollbarx \
           -orient horizontal \
           -command ".dummy.left.structures.text xview"
pack  .dummy.left.structures.scrollbar  -side right  -fill y
pack  .dummy.left.structures.scrollbarx -side bottom -fill x
pack  .dummy.left.structures.text -side top -anchor w
#
frame .dummy.left.structures.text.entry -borderwidth 2 -relief ridge
pack  .dummy.left.structures.text.entry -side top -anchor w
#
label .dummy.left.structures.text.entry.nothing \
       -text "*No structure(s) defined*" -bg red
pack  .dummy.left.structures.text.entry.nothing -side top -anchor w
.dummy.left.structures.text window create 1.0 \
  -window .dummy.left.structures.text.entry
#
# put in radio button for global/local transformations
#

frame   .dummy.left.buttons -borderwidth 2 -relief raised
pack    .dummy.left.buttons -side top -anchor w 
frame   .dummy.left.buttons.top 
frame   .dummy.left.buttons.bottom
pack    .dummy.left.buttons.top .dummy.left.buttons.bottom -side top -anchor w 

label   .dummy.left.buttons.top.text -text "Transformation:"
pack    .dummy.left.buttons.top.text -side top -anchor w

radiobutton .dummy.left.buttons.bottom.global    -text "Global"   -value 0   \
             -variable lulTraceTransformationValue                           \
             -command  {define transformation global}   
radiobutton .dummy.left.buttons.bottom.individual -text "Local" -value 1 \
             -variable lulTraceTransformationValue                            \
             -command  {define transformation local}
pack    .dummy.left.buttons.bottom.global .dummy.left.buttons.bottom.individual -side left -anchor w

  if {[show transformation type] == "global"} {
      .dummy.left.buttons.bottom.global     select
  } else {
      .dummy.left.buttons.bottom.individual select
  }

#trace variable lulTraceTransformationValue w {lulModifyTransformationRadiobutton}

#
frame .dummy.right -width 15c -height 20c
pack  .dummy.right  -side right -anchor n -fill y

# File menu ...
menubutton .frame.file -text "File" -menu .frame.file.menu -underline 0
pack .frame.file -side left
set m .frame.file.menu
menu .frame.file.menu
#$m add command -label "New..."     -command {}
$m add command -label "Open..."    -command {lulImportModelFile .dummy.left;lulInvalidateDisplay}
#$m add command -label "Close"      -command {}
$m add command -label "Save"       -command {lulSaveModelFile .dummy.left}
$m add command -label "Save As..." -command {lulExportModelFile .dummy.left}
$m add separator
$m add cascade -label "Import"      -menu $m.import
$m add cascade -label "Export"      -menu $m.export
$m add separator
$m add cascade -label "Reset..."    -menu $m.reset
$m add separator
$m add command -label "Hardcopy..." -command {lulMakeHardcopy}
$m add separator
$m add command -label "Exit" -command quit

menu  $m.import
set w $m.import
$w add command -label "Cluster..."      -command "lulImportClusterFile $w"
$w add command -label "Coords..."       -command "lulFile::ImportCoordFile"
$w add command -label "Dict..."         -command "lulImportDictionaryFile $w"
$w add command -label "GBasis..."       -command "lulImportGBasisFile"

$w add cascade -label "Vector Data..."  -menu $w.vectors 
menu  $w.vectors
set u $w.vectors
$u add command -label "CHARMm..."          -command "lulImportVectorFile 1 $w"
$u add command -label "Vector File..."     -command "lulImportVectorFile 2 $w"

$w add command -label "Tcl script..."   -command "lulDisplayEditRunTCLScript"
$w add separator
$w add cascade -label "Atom charges..." -menu $w.charges

menu  $w.charges
set u $w.charges
$u add cascade -label "GAUSSIAN..."     -menu $u.gaussian 
$u add command -label "ICON8..."        -command "lulImportAtomChargesICON8  $w"
$u add command -label "MOPAC6..."       -command "lulImportAtomChargesMOPAC6 $w"
$u add command -label "USER..."         -command "lulImportAtomChargesUSER   $w"

menu  $u.gaussian
set w $u.gaussian
$w add command -label "Natural Population" -command "lulImportAtomChargesGAUSSIAN  $w 1"
$w add command -label "Merz-Kollmann"      -command "lulImportAtomChargesGAUSSIAN  $w 2"
$w add command -label "Chelpg"             -command "lulImportAtomChargesGAUSSIAN  $w 3"

menu  $m.export
set w $m.export
$w add command -label "Cluster..."      -command "lulExportClusterFile $w"
$w add command -label "Coords..."       -command "lulFile::ExportCoordFile"
$w add command -label "Correlation..."  -command "lulExportCorrelationFile $w"
$w add cascade -label "Input"           -menu    $w.probesurf
$w add command -label "RDF..."          -command "lulExportRDFfile $w"
$w add command -label "Timeseries..."   -command "lulExportTimeSeriesFile"
#$w add command -label "RMSD..."        -command {puts RMSD}

menu   $m.export.probesurf
set  w $m.export.probesurf
$w add command -label "GAMESS..."          -command "set gomExportInputType gamess;lulCreateExportProgramInputFile"
$w add command -label "ICON8"              -command "set gomExportInputType icon8;lulCreateExportProgramInputFile"
$w add command -label "Mopac6..."          -command "set gomExportInputType mopac;lulCreateExportProgramInputFile"
$w add command -label "Probesurf..."       -command "set gomExportInputType probesurf;lulCreateExportProgramInputFile"
$w add command -label "User"               -command "set gomExportInputType user;lulCreateExportProgramInputFile"

menu  $m.reset
set w $m.reset
$w add command -label "Atom colours"     -command {reset atomcolour}
$w add command -label "gOpenmol"         -command {reset gopenmol}
$w add command -label "View"             -command {reset view}

# Edit menu ...
menubutton .frame.edit -text "Edit" -menu .frame.edit.menu -underline 0
pack .frame.edit -side left
set m .frame.edit.menu
menu .frame.edit.menu

# this works only for Windows so far ...
if {$tcl_platform(platform) == "windows"} {
 $m add cascade -label "Copy" -menu $m.copy
 menu $m.copy 
 #$m.copy add command -label "Atoms"         -command {}
 $m.copy add command -label "Bitmap"        -command {copy bitmap}
 $m.copy add command -label "Correlation"   -command {copy correlation}
 $m.copy add command -label "MSD"           -command {copy msdisplacement}
 $m.copy add command -label "RDF"           -command {copy rdf}
 $m.copy add command -label "Timeseries..." -command {lulCopyTimeSeries}
}

#$m add command -label "Annotate..."       -command {error "this is just a demo: no action has been defined for the \"Open ...\" entry"}
$m add command -label  "Cell..."                -command {lulDefineCell}
$m add command -label  "Center..."              -command {lulCenterSystem}
$m add command -label  "Identify atom..."       -command {lulIdentifyAtom}
#$m add command -label "Manipulate..."     -command {error "this is just a demo: no action has been defined for the \"Print ...\" entry"}
$m add command -label  "Light properties..."    -command {lulEditLightProperties}
$m add command -label  "Material properties..." -command {lulEditMaterialProperties}
$m add command -label  "Merge structures"       -command {atom structure merge}
$m add command -label  "Molecule..."            -command {lulEditMolecule}
$m add command -label  "Rotate/Translate..."    -command {lulRotaTransControl}
$m add command -label  "Select..."              -command {lulSelectAtoms}
$m add cascade -label  "Stereo..."              -menu $m.stereo 
$m add separator
#$m add command -label "Display attribs..."    -command "lulEditDisplay"
$m add command -label "Display properties..."  -command "lulChangeDisplayProperties"
$m add command -label "Display list properties..." -command "lulChangeDisplayListProperties"
$m add separator
$m add command -label "3D =>2D"            -command "lulMol3Dto2D"
#$m add command -label "Setup..."          -command exit

menu  $m.stereo
set w $m.stereo
$w add command -label "Pair"               -command "lulDisplayStereoPairControl"
if {[show quadstereo available]} {
    $w add command -label "Quad"           -command "lulQuadStereoControl"
}
# View menu ...
menubutton .frame.view -text "View" -menu .frame.view.menu -underline 0
pack .frame.view -side left
set m .frame.view.menu
menu .frame.view.menu

$m add command -label "Atom colour..."     -command {lulDefineAtomColour}
$m add command -label "Atom label(s)..."   -command {lulDisplayAtomLabels} 
$m add command -label "Atom tree..."       -command {lulAtomTree::create}
$m add command -label "Atom mask..."       -command {lulDefineDisplayAtoms} 
$m add command -label "Atom type..."       -command {lulDefineADType}
$m add cascade -label "Background colour"  -menu $m.background

menu  $m.background
set w $m.background
$w add command -label "Colour..." -command DefineBGColour
$w add command -label "RGB..."    -command DefineBGColourRGB

$m add command -label "Measure..."         -command {lulMeasureGeometry}
$m add command -label "Periodic table..."  -command {lulDefinePeriodicTable}
$m add cascade -label "Stereo..."          -menu $m.stereo 
$m add separator
$m add command -label "Text output..."     -command {show output}

menu  $m.stereo
set w $m.stereo
$w add command -label "Pair"               -command "lulDisplayStereoPairControl"
if {[show quadstereo available]} {
    $w add command -label "Quad"           -command "lulQuadStereoControl"
}

# Calculate menu ...
menubutton .frame.calculate -text "Calculate" -menu .frame.calculate.menu -underline 4
pack .frame.calculate -side left
set m .frame.calculate.menu
menu .frame.calculate.menu
$m add command -label "Average structure..." -command "calculate avstruct;lulInvalidateDisplay"
$m add command -label "Cluster..."           -command "lulSelectClusterAtoms"
$m add command -label "Connectivity..."      -command "lulCalculateConnectivity"
$m add command -label "Hydrogen bonds..."    -command "lulCalculateHydrogenBonds"
$m add command -label "Correlation..."       -command "lulCalculateCorrelation"
$m add command -label "Geometry..."          -command "lulMeasureGeometry"
$m add command -label "Mean sqr displ..."    -command "lulMeanSquareDisplacement"
$m add command -label "RDF..."               -command "lulSelectRDFatoms"
$m add command -label "Superimpose..."       -command "lulSuperimposeStructures"
$m add command -label "Surface centroid"     -command "lulUtility::gomCalculateSurfaceCentroid"

# Plot menu ...
menubutton .frame.plot -text "Plot" -menu .frame.plot.menu -underline 0
pack .frame.plot -side left
set m .frame.plot.menu
menu .frame.plot.menu
$m add command -label "Axis..."         -command "lulPlotAxis"
$m add command -label "Clip plane..."   -command "lulContourClipPlaneControl"
$m add command -label "Colour scale..." -command "lulPlotColourScale"
$m add command -label "Contour..."      -command "package require gom::gui::Contour; gom::gui::Contour::Control"
$m add command -label "Cutplane..."     -command "lulContourPlaneControl"
$m add command -label "LDP..."          -command "lulSelectLDPatoms"
$m add command -label "Plumber..."      -command "lulPlumber::DefinePlumber"
$m add command -label "RMSD..."         -command "lulDefineRMSDplot "
$m add command -label "Vector file..."  -command "lulPlotVectorFile "

# Trajectory menu ...
menubutton .frame.trajectory -text "Trajectory" -menu .frame.trajectory.menu -underline 0
pack .frame.trajectory -side left
set m .frame.trajectory.menu
menu .frame.trajectory.menu
#$m add command -label "Control..."  -command {error "this is just a demo: no action has been defined for the \"Open ...\" entry"}
$m add cascade -label "Fill"        -menu    $m.fill
$m add cascade -label "Delete"      -menu    $m.delete
$m add command -label "Main..."     -command {lulTrajControl}
$m add cascade -label "Monitor"     -menu    $m.monitor
$m add command -label "Time series..." -command {lulTimeSeriesManipulation}
$m add command -label "Trace..."       -command "lulDefineTraceAtoms"
$m add separator
$m add command -label "Make movie..."   -command "lulMakeMovie"

menu  $m.monitor
set w $m.monitor
$w add command -label "Distance..." -command "lulMonitorDistance"
$w add command -label "Angle..."    -command "lulMonitorAngle"
$w add command -label "Torsion..."  -command "lulMonitorTorsion"

menu  $m.fill
set w $m.fill
$w add command -label "Distance array" -command "fill distance array"
$w add command -label "Angle array"    -command "fill angle array"
$w add command -label "Torsion array"  -command "fill torsion array"

menu  $m.delete
set w $m.delete
$w add command -label "Distance array" -command "fill -distance"
$w add command -label "Angle array"    -command "fill -angle"
$w add command -label "Torsion array"  -command "fill -torsion"

# Run menu
menubutton .frame.run -text "Run" -menu .frame.run.menu -underline 0
pack .frame.run -side left
set m .frame.run.menu
menu .frame.run.menu
$m add command -label "AutoDock2plt (map)..."             -command {lulRungAutoDock2Plt}
$m add command -label "ContMan (plt)..."                  -command {lulRunContMan}
$m add command -label "Gamess2plt (cube)..."              -command {lulRungGamess2Plt}
$m add command -label "gCube2plt/g94cub2pl (cube)..."     -command {lulRungCube2Plt}
$m add command -label "Jaguar2plt (cube)..."              -command {lulRungJaguar2Plt}
$m add command -label "Join Gamess IRC files..."          -command {lulRunJoinGamessIRCfiles}
$m add command -label "Kont2plt (grid)..."                -command {lulRungKont2Plt}
$m add command -label "Pltfile (conversion)..."           -command {lulRunPltfile}
$m add command -label "Probesurf (Connolly)..."           -command {lulRunProbesurf}
$m add cascade -label "Socket server/client..."           -menu $m.socket
$m add command -label "TMole2plt (grid)..."               -command {lulRunTurboMole2Plt}
$m add command -label "UHBD2plt (map)..."                 -command {lulRunUHBD2Plt}
$m add command -label "Xvibs (conversion)..."             -command {lulRunXvibs}

menu $m.socket
set w $m.socket
$w add command -label "Server control" -command {lulSocketServerControl}
$w add command -label "Client control" -command {lulSocketClientControl}

# Tools menu
menubutton .frame.tools -text "Tools" -menu .frame.tools.menu -underline 0
pack .frame.tools -side left
set m .frame.tools.menu
menu .frame.tools.menu
$m add cascade -label "Plugins" -menu $m.plugins
$m add command -label "AutoDock..."  -command {lulAutoDock::MainControl}

# Plugins menu
set   gomPluginsMenu $m.plugins
menu $gomPluginsMenu

# Help menu ...
menubutton .frame.help -text "Help" -menu .frame.help.menu -underline 0
pack .frame.help -side right
set m .frame.help.menu
menu .frame.help.menu
$m add command -label "About..."         -command {lulAboutDialog $gOpenMolInfo}
$m add command -label "Demo..."          -command {lulShowDemo}

#if {$tcl_platform(platform) == "windows"} {
#$m add command -label "Help..."          -command {exec $gomEnv(GOM_HELP)/gopenmol.html}
#$m add command -label "Menu help..."     -command {exec cmd.exe /c $gomEnv(GOM_HELP)/gOpenMol_begin.html#COMMANDS}
#} else {
$m add command -label "Help..."          -command {htmlShowHelp gopenmol.html}
$m add command -label "Menu help..."     -command {htmlShowHelp gOpenMol_begin.html#COMMANDS}
#}
$m add command -label "Peek version..."  -command {lulUtility::gomPeekCurrentProgramVersion}
$m add cascade -label "Tutorials..."     -menu $m.tutorials 

menu $m.tutorials
set  mm $m.tutorials

$mm add command -label "Scott's Intro..."    -command {htmlShowHelp tutorials/scottanderson1/tutorial.html}


# Prelude ...
tk_menuBar .frame .frame.file .frame.file.import .frame.file.export \
.frame.edit .frame.help 
focus .frame
# ....

########################################################################

if {[file exists [file join $gomEnv(GOM_DATA) images logo.gif]]} { 
    image create photo gomlogo  -file [file join $gomEnv(GOM_DATA) images logo.gif]
    button  .dummy.right.logo   -image gomlogo -command "display 1"
} else {
    button  .dummy.right.logo     -text gOpenMol -command "display 1"
}

pack .dummy.right.logo     -side top -anchor ne

labelframe  .dummy.right.collection -text "Atom control" -padx 2 -pady 2 
#frame  .dummy.right.collection 
pack    .dummy.right.collection -side top -anchor e -pady 4

frame   .dummy.right.collection.buttons -borderwidth 2 -relief raised
pack    .dummy.right.collection.buttons -side top -anchor w 
radiobutton .dummy.right.collection.buttons.rotate    -text "Rotate"  -value 0   \
             -command  {define rotation state on}         \
             -variable gomManipulationState
radiobutton .dummy.right.collection.buttons.translate -text "Translate" -value 1 \
             -command  {define translation state on}      \
             -variable gomManipulationState
pack    .dummy.right.collection.buttons.rotate .dummy.right.collection.buttons.translate -side top -anchor w

  if {[show rotation state]} {
      .dummy.right.collection.buttons.rotate    select
  } else {
      .dummy.right.collection.buttons.translate select
  }

frame   .dummy.right.collection.selstate -borderwidth 2 -relief raised
pack    .dummy.right.collection.selstate -side top -anchor w -pady 3 -fill x

label   .dummy.right.collection.selstate.label -text "Selection:"
pack    .dummy.right.collection.selstate.label -side top  -anchor w
radiobutton .dummy.right.collection.selstate.on  -text "On"  -value 1   \
             -variable gomSelectionState         \
             -command  {define atom selection on}
radiobutton .dummy.right.collection.selstate.off -text "Off" -value 0   \
             -variable gomSelectionState         \
             -command  {define atom selection off}
pack    .dummy.right.collection.selstate.on .dummy.right.collection.selstate.off    \
         -side top -anchor w 
        .dummy.right.collection.selstate.off select

# instant update
frame   .dummy.right.collection.update -borderwidth 2 -relief raised
pack    .dummy.right.collection.update -side top -anchor w -pady 3 -fill x

label   .dummy.right.collection.update.label -text "Instant update:"
pack    .dummy.right.collection.update.label -side top  -anchor w
radiobutton .dummy.right.collection.update.on  -text "On"  -value 1   \
             -variable gomInstantDisplay -command "define window update automatic"     
radiobutton .dummy.right.collection.update.off -text "Off" -value 0   \
             -variable gomInstantDisplay -command "define window update manual"       
pack    .dummy.right.collection.update.on .dummy.right.collection.update.off    \
         -side top -anchor w 

# identify
frame   .dummy.right.collection.identify -borderwidth 2 -relief raised
pack    .dummy.right.collection.identify -side top -anchor w -pady 3 -fill x

label   .dummy.right.collection.identify.label -text "Pick atom(s):"
pack    .dummy.right.collection.identify.label -side top  -anchor w
radiobutton .dummy.right.collection.identify.on  -text "On"  -value 1   \
             -variable gomPickSwitch -command {lulActivateIdentifyAtom 1}     
radiobutton .dummy.right.collection.identify.off -text "Off" -value 0   \
             -variable gomPickSwitch -command {lulActivateIdentifyAtom 0}       
button      .dummy.right.collection.identify.unpickall -text "Unpick all" \
             -command {picking reset;lulInvalidateDisplay}
pack    .dummy.right.collection.identify.on .dummy.right.collection.identify.off    \
         -side top -anchor w 
pack    .dummy.right.collection.identify.unpickall -side top -pady 2

########################################################################

frame .statusline -relief raised -borderwidth 2
pack  .statusline -side bottom

# handle the text to File Size entry
entry .statusline.filesize  -width 22 -relief sunken -bd 2 -state disabled -textvariable gomFileSize
#### PROC
proc lulPutText2Info1 {Text} {
     .statusline.filesize configure -state normal
     .statusline.filesize delete 0 end
	 .statusline.filesize insert 0 $Text
	 .statusline.filesize configure -state disabled
}

# handle the text to Number of Atoms entry
entry .statusline.numatoms  -width 17 -relief sunken -bd 2 -state disabled -textvariable gomNumAtoms
### PROC
proc lulPutText2Info3 {Text} {
     .statusline.numatoms configure -state normal
     .statusline.numatoms delete 0 end
	 .statusline.numatoms insert 0 $Text
     .statusline.numatoms configure -state disabled
}

entry .statusline.timeflies -width 20 -relief sunken -bd 2       \
       -state disabled -textvariable time_flies -foreground blue \
	   -font {Helvetica 12 bold}
pack  .statusline.filesize .statusline.numatoms .statusline.timeflies -side left -expand 1
### PROC
proc lulTimeFliesWrapper { Fvalue } {

set Ivalue [expr int($Fvalue * 10.0)]

switch $Ivalue {
    0  {set Text "|*|"}
    1  {set Text "|*|*|"}
    2  {set Text "|*|*|*|"}
    3  {set Text "|*|*|*|*|"}
    4  {set Text "|*|*|*|*|*|"}
    5  {set Text "|*|*|*|*|*|*|"}
    6  {set Text "|*|*|*|*|*|*|*|"}
    7  {set Text "|*|*|*|*|*|*|*|*|"}
    8  {set Text "|*|*|*|*|*|*|*|*|*|"}
    9  {set Text "|*|*|*|*|*|*|*|*|*|*|"}
    10 {set Text "|*|*|*|*|*|*|*|*|*|*|*|"}
    default {set Text " "}
}

lulPushText2TimeFlies $Text

}
######################### push text to time flies ###############
proc lulPushText2TimeFlies { Text } {

.statusline.timeflies configure -state normal
.statusline.timeflies delete 0 end
.statusline.timeflies insert end $Text
.statusline.timeflies configure -state disabled

# update the widget
update idletasks

}

frame   .commandframe
label   .commandframe.commandlabel  -text "Command: " -pady 2 -padx 2
entry   .commandframe.commandline   -width  40 -relief sunken -bd 2 \
         -textvariable command_line
pack    .commandframe -side bottom -fill x
pack    .commandframe.commandline .commandframe.commandlabel  -side right

bind    .commandframe.commandline <Return> {
    if {[catch $command_line errmsg]} {
	lulErrorStackDialog $errmsg
    }
    lulPushToCommandStack $command_line;
    # protect against "double display"
    if {![string match $command_line* "display"]} {lulInvalidateDisplay}
    set command_line "";
}

bind    .commandframe.commandline <Up>   {lulGetNextFromStack .commandframe.commandline}
bind    .commandframe.commandline <Down> {lulGetPrevFromStack .commandframe.commandline}

#
# Fill the structure list box with an entry:
# if $Action == 0 create new list
#            == 1 append to old list
#             < 0 just delete the list
#
### PROC
proc lulFillStructureListbox {Action Entry} {
     set Size [.dummy.left.structurelist size]

     if {$Action < 0 } {
       if { $Size > 0 } {
         .dummy.left.structurelist delete 1 $Size
          return
       }
       else
          return
     }
     if { ($Action == 0) && ($Size > 0)} {
      .dummy.left.structurelist delete 1 $Size
     }
      .dummy.left.structurelist insert end $Entry
      .dummy.left.structurelist selection set $Size
}

#############################################################
# PROC
proc lulCreateAtomInputEntries { baseFrame entryWidth SegDef ResDef AtmDef baseVar allowMultiple } {
    global gomAtomPickingTargetWidgets

    # Segment
    frame  $baseFrame.segment
    label  $baseFrame.segment.label -text " Segment:" -width 10 -anchor w
    entry  $baseFrame.segment.input -width $entryWidth
    pack   $baseFrame.segment -side top -anchor w
    pack   $baseFrame.segment.label $baseFrame.segment.input -side left

    if { "" != $baseVar } {$baseFrame.segment.input configure -textvariable ${baseVar}Segment}
    if { "" == [$baseFrame.segment.input get] } {$baseFrame.segment.input insert 0 $SegDef}

    # Residue
    frame  $baseFrame.residue
    label  $baseFrame.residue.label -text " Residue:" -width 10 -anchor w
    entry  $baseFrame.residue.input -width $entryWidth
    pack   $baseFrame.residue -side top -anchor w
    pack   $baseFrame.residue.label $baseFrame.residue.input -side left

    if { "" != $baseVar } {$baseFrame.residue.input configure -textvariable ${baseVar}Residue}
    if { "" == [$baseFrame.residue.input get] } {$baseFrame.residue.input insert 0 $ResDef}

    # Atom
    frame  $baseFrame.atom
    label  $baseFrame.atom.label -text " Atom:" -width 10 -anchor w
    entry  $baseFrame.atom.input -width $entryWidth
    pack   $baseFrame.atom -side top -anchor w
    pack   $baseFrame.atom.label $baseFrame.atom.input -side left

    if { "" != $baseVar } {$baseFrame.atom.input configure -textvariable ${baseVar}Atom}
    if { "" == [$baseFrame.atom.input get] } {$baseFrame.atom.input insert 0 $AtmDef}

    set gomAtomPickingTargetWidgets($baseFrame) [list \
        segment.input residue.input atom.input $allowMultiple]
}

#############################################################
# PROC
proc lulModifyTransformationRadiobutton {name1 name2 op} {
     
     global lulTraceTransformationValue

  if {[show transformation type] == "global"} {
      .dummy.left.buttons.bottom.global     select
  } else {
      .dummy.left.buttons.bottom.individual select
  }


}


###############################################################################
#
# dialog widget
#
proc lulAboutDialog InputText {

after idle {.dialog1.msg configure -wraplength 4i}
set i [tk_dialog .dialog1 "gOpenMol Info" $InputText \
info 0 OK ]

#switch $i {
#    0 {puts "You pressed OK"}
#}
}

##################### update background colour ############################
#
# define background colour widget
#
# PROC
proc DefineBGColourRGB {} {

     global env
     global gomControlFont


set w .gomcolourscale
catch {destroy $w}
toplevel    $w
wm title    $w "Colour scale"
wm iconname $w "Colour"

frame $w.buttons -borderwidth 2 -relief raised -bd 2
pack $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
pack $w.buttons.dismiss -side left -expand 1

scan [show bgcolor] "%f %f %f" red green blue

button $w.button -text "Click to choose the colour" -height 2 \
        -command UpdateBGColour
scale  $w.red   -label Red   -from 0 -to 255 -length 8c -orient vertical \
      -command "ColourButton $w"
scale  $w.green -label Green -from 0 -to 255 -length 8c -orient vertical \
      -command "ColourButton $w"
scale  $w.blue  -label Blue  -from 0 -to 255 -length 8c -orient vertical \
      -command "ColourButton $w"

#listbox $w.listbox -width 10 -height 10

pack   $w.button -padx 2 -pady 2
pack   $w.red $w.green $w.blue -side left
#pack   $w.listbox

# put the sliders to the current background colour
$w.red   set   [expr int($red   * 255.)]
$w.green set   [expr int($green * 255.)]
$w.blue  set   [expr int($blue  * 255.)]

#set f [open [file join $gomEnv(GOM_DATA) colour_table.data]]
#while {[gets $f line] >= 0} {
#    $w.listbox insert end $line
#}
}
#######################################################################
# change the colour button background colour according to the sliders
#
proc ColourButton {w a} {
  set colour [format "#%02x%02x%02x" [$w.red get] [$w.green get] [$w.blue get]]
  $w.button config -bg $colour -fg [lulGetVisibleForegroundColour $colour]
}
#
# take the red green blue colour values (as integers) and convert them to
# float and change the background colour
#
proc UpdateBGColour {} {
     set red   [.gomcolourscale.red   get]
     set green [.gomcolourscale.green get]
     set blue  [.gomcolourscale.blue  get]
     define bgcolour "[expr $red/255.] [expr $green/255.] [expr $blue/255.]"
     display
}

#############################################################################
# PROC
proc DefineBGColour {} {

# get current background colour
    scan [show bgcolour] "%f %f %f" red green blue

    set colour [tk_chooseColor -title "Choose background colour" \
         -initialcolor [lulColourFloat2Hex $red $green $blue]]

    if { $colour != "" } {
         define bgcolour "[lulColourHex2Float $colour]"
         display
    }
}

########################## end of update bgcolour #########################

##################### Atom type display ############################
#
# define atom type widget
#
proc lulDefineADType {} {

     global gomHelpFile
     global gomControlFont

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomatomtype
catch {destroy $w}
toplevel $w 
wm title $w "Atom type"
wm iconname $w "Atom type"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulDoDisplayType $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(atomtype)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Atom type display:"
pack   $w.label -side top -anchor w

lulCreateAtomInputEntries $w 30 "*" "*" "*" "gomTypeDisplay" 1

labelframe  $w.frame4 -borderwidth 2 -relief ridge -text "Display state" -padx 2 -pady 2
radiobutton $w.frame4.on  -text "On"  -variable disp_action -value 1
radiobutton $w.frame4.off -text "Off" -variable disp_action -value 0
pack $w.frame4 -side left 
pack $w.frame4.on $w.frame4.off -side top -anchor w
$w.frame4.on select

labelframe  $w.frame5 -text "Display type" -borderwidth 2 -relief ridge -padx 2 -pady 2

radiobutton $w.frame5.stick    -text "Stick"    \
            -variable disp_type -value stick    \
            -command "$w.frame6.cylinder configure -state disabled;\
                      $w.frame6.sphere   configure -state disabled;\
                      $w.frame6.cpkscale configure -state disabled;\
                      lulInvalidateDisplay"
			          
radiobutton $w.frame5.cpk      -text "CPK"      \
            -variable disp_type -value cpk      \
            -command "$w.frame6.cylinder configure -state disabled;\
                      $w.frame6.sphere   configure -state disabled;\
                      $w.frame6.cpkscale configure -state normal;\
                      lulInvalidateDisplay"

radiobutton $w.frame5.licorice -text "Licorice" \
            -variable disp_type -value licorice \
            -command "$w.frame6.cylinder configure -state normal;\
                      $w.frame6.sphere   configure -state normal;\
                      $w.frame6.cpkscale configure -state disabled;\
                      lulInvalidateDisplay"

pack $w.frame5 -side top -anchor e
pack $w.frame5.stick $w.frame5.cpk $w.frame5.licorice -side left
$w.frame5.stick select

labelframe $w.frame6 -text "Display params" -borderwidth 2 -relief ridge -padx 2 -pady 2
label $w.frame6.cylinderlabel -text "Cyl. rad: "
entry $w.frame6.cylinder      -width 10
label $w.frame6.spherelabel   -text "Sph. rad: "
entry $w.frame6.sphere        -width 10
label $w.frame6.cpkscalelabel -text "CPK scale: "
entry $w.frame6.cpkscale      -width 10

pack $w.frame6 -side bottom
pack $w.frame6.cylinderlabel $w.frame6.cylinder -side left
pack $w.frame6.spherelabel   $w.frame6.sphere   -side left
pack $w.frame6.cpkscalelabel $w.frame6.cpkscale -side left

set NumStruct [show molstructures]

if {$NumStruct != 0} {
    $w.frame6.cylinder delete 0 end
    $w.frame6.cylinder insert 0 [format %.6g [show licocylinder]]
    $w.frame6.sphere   delete 0 end
    $w.frame6.sphere   insert 0 [format %.6g [show licosphere]]
    $w.frame6.cpkscale delete 0 end
}

$w.frame6.cpkscale insert 0 "1.0"

# original state ...
$w.frame6.cylinder configure -state disabled
$w.frame6.sphere   configure -state disabled
$w.frame6.cpkscale configure -state disabled


}
#######################################################################
# PROC
# react to display type
#
proc lulDoDisplayType { w } {
  global disp_type
  global disp_action
  global gomTypeDisplaySegment
  global gomTypeDisplayResidue
  global gomTypeDisplayAtom

  if { $disp_type == "stick" } {
     if { $disp_action } {
          set Action "atom display"
     } else {
          set Action "atom -display"
     }
  } else {
     if { $disp_action } {
          set Action "atom  $disp_type"
     } else {
          set Action "atom -$disp_type"
     }
  }

  eval [concat $Action \{$gomTypeDisplaySegment\} \{$gomTypeDisplayResidue\} \{$gomTypeDisplayAtom\}]

  if {$disp_type == "cpk"} {
     set cpkScale [$w.frame6.cpkscale get]
	 eval "atom scale cpk $cpkScale \{$gomTypeDisplaySegment\} \{$gomTypeDisplayResidue\} \{$gomTypeDisplayAtom\}"
  } elseif {$disp_type == "licorice"} {
     define licosphere   [$w.frame6.sphere   get]  
     define licocylinder [$w.frame6.cylinder get]
  }

  lulInvalidateDisplay
}
#
# take the red green blue colour values (as integers) and convert them to
# float and change the background colour
#
proc ApplyAtomTypeDisplay {} {
     display
}
########################## end of update bgcolour #########################

##################### import gOpenMol file ############################

proc lulImportModelFile { w } {
    global gomFileName

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"GOM files"		{.gom}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .gom -title "Read Model File"]]

    if [string compare $file ""] {
      import model "$file"
      set gomFileName "$file"
# change to current directory
      lulChangeDirectory "$file"
    }
}

########################## end of import model #########################

##################### export gOpenMol file ############################

proc lulExportModelFile { w } {
    global gomFileName


# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"GOM files"		{.gom}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .gom -title "Save Model File"]]

    if [string compare $file ""] {
      export model "$file"
      set gomFileName "$file"
# change to current directory
      lulChangeDirectory "$file"
    }
}

#
# if there is a name already available save with that name
# if not ask for a name
#
###############################################################
# PROC
proc lulSaveModelFile { w } {
    global gomFileName

    if {$gomFileName == "" } {
        lulExportModelFile $w
    } else {
        export model "$gomFileName"
    }
}

########################## end of export model #########################


###################################################################
# PROC
proc lulSetSystemTranslationState { } {

   if {[show system translation] == "on"} {
      puts [show system translation]
      define system translation off
      puts [show system translation]
   } elseif {[show system translation] == "off"} {
      define system translation on
   } else {
      gomError "'show system translation' does not return 'on' or 'off'"
   }

   lulInvalidateDisplay
}

##################### import dictionary file ############################

proc lulImportDictionaryFile { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Dictionary files"		{.dic}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .dic -title "Read Dictionary File"]]

    if {$file != ""} {
      puts "Importing dictionary file '$file'"
      import dictionary "$file"
# change to current directory
      lulChangeDirectory "$file"
    }
}

########################## end of import dictionary #########################

##################### import GBASIS file ############################
# PROC
proc lulImportGBasisFile {} {

     global gomHelpFile
     global gomControlFont

set w .gomgbasis
catch {destroy $w}
toplevel $w 
wm title $w "Import GBasis"
wm iconname $w "Import GBasis"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -command "destroy $w" \
        -font "$gomControlFont"
button $w.buttons.apply   -text Apply -command "lulProcessGBasisFile $w" \
        -font "$gomControlFont"
button $w.buttons.help    -text Help -command \
       "htmlShowHelp $gomHelpFile(importgbasis)" \
	    -font "$gomControlFont"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label $w.label -text "Import GBasis file:" -font "$gomControlFont"
pack  $w.label -side top -anchor w

frame       $w.fileinput  -borderwidth 2 -relief raised
#label       $w.fileinput.label    -text "File name: "
labelframe  $w.fileinput.input  -text "File name"  -borderwidth 2 -relief ridge -padx 2 -pady 2
entry       $w.fileinput.input.filename -width 40
button      $w.fileinput.input.browse   -text "Browse..."    \
             -command "lulDoImportGBasisFile $w.fileinput.input.filename" \
             -width 10

pack   $w.fileinput -side top -anchor w
#pack   $w.fileinput.label     -side left -anchor w
pack   $w.fileinput.input     -side left
pack   $w.fileinput.input.filename  -side left
pack   $w.fileinput.input.browse   -padx 2 -side left

frame $w.left  -borderwidth 2 -relief raised
frame $w.right -borderwidth 2 -relief raised

pack $w.left  -side left 
pack $w.right -side right

listbox $w.left.entries -width 40 -height 20 \
         -yscrollcommand "$w.left.scrolly set"
pack    $w.left.entries -side left
scrollbar $w.left.scrolly -command "$w.left.entries yview"
pack      $w.left.scrolly -side right -fill y

text      $w.right.entry  -width 40 -height 20 \
          -yscrollcommand "$w.right.scrolly set"
pack      $w.right.entry  -side left
scrollbar $w.right.scrolly -command "$w.right.entry yview"
pack      $w.right.scrolly -side right -fill y

bind      $w.left.entries <Double-Button-1> {lulShowGBasisEntry [selection get] .gomgbasis}

}
##########################################################################
# PROC
proc lulDoImportGBasisFile { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"GBasis files"		{.data}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .data -title "Read GBasis File"]]

    if [string compare $file ""] {
    $w delete 0 end
    $w insert 0 $file
    }
}

##########################################################################
# PROC
proc lulProcessGBasisFile { w } {

   set FileName [$w.fileinput.input.filename get]

   if {$FileName == ""} return

   import gbasis "$FileName"

# change to current directory
    lulChangeDirectory "$FileName"

   set Entries [show gbasis sets]

   for {set i 1} {$i <= $Entries} {incr i} {
     $w.left.entries inser end "$i [show gbasis tag $i]"
   }

}
########################################################################
# PROC
proc lulShowGBasisEntry {Index w} {

   set Index [lindex $Index 0]

   $w.right.entry delete 0.0 end
   $w.right.entry insert 0.0 [show gbasis entry $Index]
}
########################## end of import GBasis #########################

##################### import vector file ############################

proc lulImportVectorFile { what w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    if {$what == 1} {
      set types {
	  {"Vector files"		{.crd .CRD}		TEXT}
	  {"All files"		*}
      }
    } elseif {$what == 2} {
      set types {
	  {"Vector files"		{.txt .TXT}		TEXT}
	  {"All files"		*}
      }
    } else {
      gomError "unknown 'vector file type' (= $what)"
      return
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .crd -title "Read Vector File"]]

    if [string compare $file ""] {
      if {$what == 1} {
        import vector charmm   "$file"
      } elseif {$what == 2} {
        import vector flatfile "$file"
      } else {
        gomError "unknown vector data file type: '$what'"
        return
      }
# change to current directory
      lulChangeDirectory "$file"
    }
}

########################## end of import dictionary #########################

##################### import cluster file ############################

proc lulImportClusterFile { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Cluster data files"		{.dat .DAT}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .dat -title "Read Cluster File"]]

    if [string compare $file ""] {
      import cluster "$file"
# change to current directory
      lulChangeDirectory "$file"
    }
}

##################### export cluster file ############################

proc lulExportClusterFile { w } {

    set ClusterStatus [show cluster status]

	if {!$ClusterStatus} {
	    lulErrorDialog "ERROR - no cluster data is available to be exported"
		return
    }

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Cluster file"		{.dat}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .dat -title "Save Cluster File"]]

    if [string compare $file ""] {
      export cluster "$file"
# change to current directory
      lulChangeDirectory "$file"
    }
}
########################## end of export cluster #########################

##################### export correlation file ############################
# PROC
proc lulExportCorrelationFile { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Correlation data files"	{.dat .DAT}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .dat -title "Save Correlation File"]]

    if [string compare $file ""] {
      export correlation "$file"
# change to current directory
      lulChangeDirectory "$file"
    }
}
########################## end of export correlation #########################

##################### Atom display ############################
#
# define atom widget
#
proc lulDefineDisplayAtoms {} {

     global gomHelpFile
     global gomControlFont
     global gomAroundColourValue

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

  set gomAroundColourValue  #ffffff

set w .gomatoms
catch {destroy $w}
toplevel $w 
wm title $w "Display Atoms"
wm iconname $w "Display Atoms"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulDoDisplayAtoms $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(atomdisplaymask)"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Display atoms:"
pack   $w.label -side top -anchor w

lulCreateAtomInputEntries $w 30 "*" "*" "*" "gomDisplay" 1

frame       $w.frame4 -relief raised -borderwidth 2
label       $w.frame4.label -text "Display state:"
radiobutton $w.frame4.on  -text "On"  -variable disp_action -value 1
radiobutton $w.frame4.off -text "Off" -variable disp_action -value 0
pack        $w.frame4 -side top -pady 4 -anchor w 
pack        $w.frame4.label $w.frame4.on $w.frame4.off -side top -anchor w
$w.frame4.on select

# selection residue/atom

frame       $w.around6 -relief raised -borderwidth 2
label       $w.around6.label -text "Select by:"
radiobutton $w.around6.residue  -text "Residue"  -value 1 \
             -command {atom selection residue}
radiobutton $w.around6.atom     -text "Atom"     -value 0 \
             -command {atom selection atom}
pack $w.around6 -side top -pady 4
pack $w.around6.label $w.around6.residue $w.around6.atom -side left -anchor w

set SelectionMode [show atom selection]

if       {$SelectionMode == "residue"} {
     $w.around6.residue select
} elseif {$SelectionMode == "atom"} {
     $w.around6.atom    select
} else {
          lulErrorDialog {ERROR: Wrong selection mode. Has to be Residue/Atom!}
	      return
}

# around

frame       $w.around5 -relief raised -borderwidth 2
label       $w.around5.label -text "Use around:"
radiobutton $w.around5.on  -text "On"  -variable around_action -value 1 \
             -command {.gomatoms.glob.segment.input configure -state normal;
                       .gomatoms.glob.residue.input configure -state normal;
                       .gomatoms.glob.atom.input    configure -state normal;
                       .gomatoms.glob.around5.distance configure -state normal;
                       .gomatoms.glob.around6.on   configure -state normal;
                       .gomatoms.glob.around6.off  configure -state normal;
                       .gomatoms.glob.around7.colour configure -state normal}
radiobutton $w.around5.off -text "Off" -variable around_action -value 0 \
             -command {.gomatoms.glob.segment.input configure -state disabled;
                       .gomatoms.glob.residue.input configure -state disabled;
                       .gomatoms.glob.atom.input    configure -state disabled;
                       .gomatoms.glob.around5.distance configure -state disabled;
                       .gomatoms.glob.around6.on   configure -state disabled;
                       .gomatoms.glob.around6.off  configure -state disabled;
                       .gomatoms.glob.around7.colour configure -state disabled}
pack $w.around5 -side top -pady 4
pack $w.around5.label $w.around5.on $w.around5.off -side left -anchor w
$w.around5.off select

frame       $w.glob -relief raised -borderwidth 2
pack        $w.glob -side top

lulCreateAtomInputEntries $w.glob 30 "*" "*" "*" "gomAroundDisplay" 1

$w.glob.segment.input configure -state disabled
$w.glob.residue.input configure -state disabled
$w.glob.atom.input    configure -state disabled

frame  $w.glob.around5
label  $w.glob.around5.label     -text "Distance: " -width 10
entry  $w.glob.around5.distance  -width 10 
pack   $w.glob.around5 -side top -anchor w -pady 4
pack   $w.glob.around5.label $w.glob.around5.distance -side left

$w.glob.around5.distance configure -state disabled

frame       $w.glob.around6 -relief raised -borderwidth 2
label       $w.glob.around6.label -text "Use colour:"
radiobutton $w.glob.around6.on  -text "On"  -variable gomAroundColour -value 1 \
             -command {.gomatoms.glob.around7.colour configure -state normal}
radiobutton $w.glob.around6.off -text "Off" -variable gomAroundColour -value 0 \
             -command {.gomatoms.glob.around7.colour configure -state disabled}
pack $w.glob.around6 -side top -pady 4
pack $w.glob.around6.label $w.glob.around6.on $w.glob.around6.off -side left -anchor w
$w.glob.around6.off select

$w.glob.around6.on   configure -state disabled
$w.glob.around6.off  configure -state disabled

frame   $w.glob.around7 
label   $w.glob.around7.label     -text "Colour: "
button  $w.glob.around7.colour    -text "Click to change colour"    \
	-command "::lulChooseButtonColour $w.glob.around7.colour \
	    ::gomAroundColourValue {Choose colour}" \
	-fg [::lulGetVisibleForegroundColour $gomAroundColourValue] \
	-bg $gomAroundColourValue
pack    $w.glob.around7 -side top -anchor w -pady 4
pack    $w.glob.around7.label $w.glob.around7.colour -side left

$w.glob.around7.colour configure -state disabled

}
###########################################################################
# react to display atom
#
# PROC
proc lulDoDisplayAtoms { w } {

     global  gomDisplaySegment gomDisplayResidue gomDisplayAtom
     global  gomAroundDisplaySegment gomAroundDisplayResidue gomAroundDisplayAtom
     global  disp_action
     global  gomAroundColour
     global  gomAroundColourValue
     global  around_action

     if { $disp_action } {
          set Action "atom display"
     } else {
          set Action "atom -display"
     }

  set gomDisplaySegment [string trim $gomDisplaySegment]
  if {$gomDisplaySegment == ""} {set gomDisplaySegment "*"}

  set gomDisplayResidue [string trim $gomDisplayResidue]
  if {$gomDisplayResidue == ""} {set gomDisplayResidue "*"}
  
  set gomDisplayAtom    [string trim $gomDisplayAtom]
  if {$gomDisplayAtom == ""} {set gomDisplayAtom "*"}

  if {$around_action} {
    set gomAroundDisplaySegment [string trim $gomAroundDisplaySegment]
    set gomAroundDisplayResidue [string trim $gomAroundDisplayResidue]
    set gomAroundDisplayAtom    [string trim $gomAroundDisplayAtom]
    set Distance                [string trim [$w.glob.around5.distance get]]

    if {$Distance == ""} {
          lulErrorDialog {ERROR: Distance value missing!}
	      return
	}

      if {$gomAroundDisplaySegment == ""} {set gomAroundDisplaySegment "*"}
      if {$gomAroundDisplayResidue == ""} {set gomAroundDisplayResidue "*"}
      if {$gomAroundDisplayAtom    == ""} {set gomAroundDisplayAtom    "*"}

      if {$gomAroundColour} {
      eval $Action $gomDisplaySegment $gomDisplayResidue $gomDisplayAtom around $Distance      \
                   $gomAroundDisplaySegment $gomAroundDisplayResidue $gomAroundDisplayAtom  \
                   {[lulColourHex2Float $gomAroundColourValue]}
      } else {
      eval $Action \{$gomDisplaySegment\} \{$gomDisplayResidue\} \{$gomDisplayAtom\} around $Distance \
                   \{$gomAroundDisplaySegment\} \{$gomAroundDisplayResidue\} \{$gomAroundDisplayAtom\}  
      }

  } else {
    eval [concat $Action \{$gomDisplaySegment\} \{$gomDisplayResidue\} \{$gomDisplayAtom\}]
  }

  lulInvalidateDisplay
}

##################### Atom colour display ############################
#
# define atom widget
#
proc lulDefineAtomColour {} {

     global gomAtomColour
     global gomHelpFile
     global gomControlFont
     global gomColByCharge

     package require gom::gui::Widgets
     namespace import -force ::gom::gui::Widgets::ColorButton

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomatomcolor
catch {destroy $w}
toplevel $w 
wm title $w "Define Atom Colour"
wm iconname $w "Define Atom Colour"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -command "destroy $w" \
        -font "$gomControlFont"
button $w.buttons.apply   -text Apply -command "lulDoDefineAtomColour $w" \
        -font "$gomControlFont"
button $w.buttons.help    -text Help -command \
       "htmlShowHelp $gomHelpFile(colourselection)" \
        -font "$gomControlFont"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Define atom colour:"
pack   $w.label -side top -anchor w

lulCreateAtomInputEntries $w 30 "*" "*" "*" "gomDisplay" 1

frame $w.frame4
ColorButton $w.frame4.cbutton -variable gomAtomColour
pack $w.frame4 -anchor w
pack $w.frame4.cbutton -padx 4 -pady 10

frame $w.frame5
button $w.frame5.resetbutton -text "Reset atom colour" \
        -command {reset atomcolour;lulInvalidateDisplay}
pack $w.frame5 -anchor w
pack $w.frame5.resetbutton -padx 4 -pady 5

labelframe  $w.frame6 -text "Colour control" -borderwidth 2 -relief ridge -padx 2 -pady 2
checkbutton $w.frame6.label0  -text "Colour by charge. Min: " -variable gomColByCharge
entry       $w.frame6.min -width 6
label       $w.frame6.label1  -text " Max: "
entry       $w.frame6.max -width 6
pack        $w.frame6 -anchor w -padx 4 -pady 5
pack        $w.frame6.label0 $w.frame6.min $w.frame6.label1 $w.frame6.max -side left  
 
}

###########################################################################
# react to atom colour
#
# PROC
proc lulDoDefineAtomColour { w } {
    global gomDisplaySegment
    global gomDisplayResidue
    global gomDisplayAtom
    global gomAtomColour
    global gomColByCharge
    
    if {!$gomColByCharge} {
	atom colour $gomDisplaySegment $gomDisplayResidue $gomDisplayAtom \
		[lulColourHex2Float $gomAtomColour]
    } else {
	atom colour bycharge $gomDisplaySegment $gomDisplayResidue $gomDisplayAtom \
		[$w.frame6.min get] [$w.frame6.max get]
    }

    lulInvalidateDisplay
}
#########################################################################

######################### update structure name list ###############
#
# PROC
proc lulUpdateStructureList { } {

    global gomNumStructures
    global gomSelectStruct
    global gomStructureNames

    catch {eval destroy [.dummy.left.structures.text window names]}
    array set gomStructureNames {}
    set gomNumStructures [show molstructures]
# select by default
    if { $gomNumStructures > 0 } { select structure $gomNumStructures }
#
    for { set i 1 } { $i <= $gomNumStructures } { incr i } {
        set Name [show atom structure name $i]
        frame .dummy.left.structures.text.entry$i -borderwidth 2 -relief ridge
        pack  .dummy.left.structures.text.entry$i -side top -anchor w
        checkbutton .dummy.left.structures.text.entry$i.structurelist$i \
            -text "($i) $Name" -cursor hand2 \
            -variable gomSelectStruct($i) \
            -command "lulHandleStructureSelection $i"
        pack .dummy.left.structures.text.entry$i.structurelist$i \
            -side top -fill y -anchor w
        set gomStructureNames($i) $Name
        .dummy.left.structures.text window create $i.0 \
            -window .dummy.left.structures.text.entry$i
        .dummy.left.structures.text.entry$i.structurelist$i select
    }

    if { $gomNumStructures == 0 } {
        frame .dummy.left.structures.text.entry -borderwidth 2 -relief ridge
        pack  .dummy.left.structures.text.entry -side top -anchor w
        label .dummy.left.structures.text.entry.nothing \
            -text "*No structure(s) defined*" -bg red
        pack  .dummy.left.structures.text.entry.nothing -side top -anchor w
        .dummy.left.structures.text window create 1.0 \
            -window .dummy.left.structures.text.entry
    }
}
############################################################################
# PROC
proc lulUpdateStructureSelection {what , which} {

  if {$what} {
     .dummy.left.structures.text.entry$which.structurelist$which select
  } else {
     .dummy.left.structures.text.entry$which.structurelist$which deselect
  }
}
############################################################################
# PROC
proc lulHandleStructureSelection {which} {

   global gomSelectStruct

 if {$gomSelectStruct($which)} {
     select  structures $which
 } else {
     select -structures $which
 }

}
############################################################################
# PROC
proc lulTrajControl {} {

     global gomHelpFile
     global gomTrajType
     global gomTrajFileName
     global gomControlFont
     global gomControlFontHuge

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

# reset stuff...
set gomTrajFileName ""

set w .gomtraj
catch {destroy $w}
toplevel $w 
wm title $w "Trajectory Control"
wm iconname $w "Trajectory Control"

#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulDoTrajFileParams  $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(maintrajectory)"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1
#
    set TrajectoryDefined [show trajectory defined]
    
	frame $w.left 

	labelframe $w.right -borderwidth 2 -text "Formats:" -relief ridge -padx 2 -pady 2 

	pack $w.left   -side left   -anchor n -fill both
	pack $w.right  -side right  -anchor e

    frame $w.left.frame1 -borderwidth 2 -relief raised
    set base $w.left.frame1
	label $base.filenamelabel \
		-text {File name: }

	entry $base.filename \
		-textvariable gomTrajFileName -width 40

# make a <return> bind
    bind $base.filename <Return> "lulGetFileData $w"

	button $base.filebrowse \
		-padx 10 \
		-pady 3 \
		-text Browse... -command "lulGetTrajFile $base.filename"

    pack $base -side top -anchor n
	pack $base.filenamelabel -side left
    pack $base.filename -side left -padx 3
	pack $base.filebrowse -side left

    frame $w.left.buttonframe 
    set base $w.left.buttonframe
    button $base.button -text "Import file" -command "lulGetFileData $w" \
            -font "$gomControlFont"
	pack $base -side top -anchor n -fill x
	pack $base.button -pady 3

    set base $w.right

    set TrajType [show trajectory type]
    if {$TrajType != ""} {
       set gomTrajType $TrajType
    }
    
	radiobutton $base.amber \
		-text AMBER \
		-value amber \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.amber     -side top -anchor w

	radiobutton $base.famber \
		-text AMBER(F) \
		-value famber \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector on}

    pack $base.famber     -side top -anchor w

	radiobutton $base.cerius2 \
		-text Cerius2 \
		-value cerius2 \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.cerius2     -side top -anchor w

	radiobutton $base.charmm \
		-text CHARMM \
		-value charmm \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.charmm     -side top -anchor w

	radiobutton $base.discover \
		-text Discover \
		-value discover \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.discover     -side top -anchor w

	radiobutton $base.udlpoly \
		-text UDL_Poly \
		-value udl_poly \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.udlpoly     -side top -anchor w

	radiobutton $base.fdlpoly \
		-text FDL_Poly \
		-value fdl_poly \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector on}

    pack $base.fdlpoly     -side top -anchor w

	radiobutton $base.gromacs \
		-text Gromacs \
		-value gromacs \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.gromacs     -side top -anchor w

	radiobutton $base.gromos \
		-text Gromos \
		-value gromos \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.gromos     -side top -anchor w

	radiobutton $base.gromos96a \
		-text Gromos96A \
		-value gromos96a \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector on}

    pack $base.gromos96a  -side top -anchor w

	radiobutton $base.hyperchem \
		-text HyperChem \
		-value hyperchem \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.hyperchem     -side top -anchor w

	radiobutton $base.mumod \
		-text MUMOD \
		-value mumod \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.mumod     -side top -anchor w

	radiobutton $base.tinker \
		-text TINKER \
		-value tinker \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector on}

    pack $base.tinker    -side top -anchor w

	radiobutton $base.xmol \
		-text XMOL \
		-value xmol \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector on}

    pack $base.xmol     -side top -anchor w

	radiobutton $base.xplor \
		-text XPLOR \
		-value xplor \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.xplor     -side top -anchor w

	radiobutton $base.yasp \
		-text YASP \
		-value yasp \
		-variable gomTrajType \
        -command {lulTouchTrajectoryMethodSelector off}

    pack $base.yasp     -side top -anchor w

#    frame $w.left.frame3 -borderwidth 2 -relief ridge
    labelframe $w.left.frame3 -text "Trajectory info" -borderwidth 2 -relief ridge -padx 2 -pady 2 
    pack       $w.left.frame3 -anchor w -side top 

    frame $w.left.frame3.frame31
    set base $w.left.frame3.frame31
	label $base.frameslabel -text "Num of Frames:" -width 15
	entry $base.frames      -width 20 -state disabled
	pack $base  -anchor w -side top
	pack $base.frameslabel   -side left
	pack $base.frames        

    frame $w.left.frame3.frame32
    set base $w.left.frame3.frame32
    label $base.firstlabel  -text "First frame:" -width 15
	entry $base.first       -width 20
	pack  $base -anchor w -side top
	pack  $base.firstlabel   -side left
    pack  $base.first        

    frame $w.left.frame3.frame33
    set base $w.left.frame3.frame33
    label $base.lastlabel  -text "Last frame:" -width 15
	entry $base.last       -width 20
	pack  $base -anchor w -side top
	pack  $base.lastlabel $base.last  -side left

    frame $w.left.frame3.frame34
    set base $w.left.frame3.frame34
    label $base.steplabel  -text "Step frame(s):" -width 15
	entry $base.step       -width 20
	pack  $base -anchor w -side top
	pack  $base.steplabel $base.step  -side left

    frame $w.left.frame3.frame35 -borderwidth 2 
    pack  $w.left.frame3.frame35 -side top -anchor w -pady 2 -padx 2 -fill x
    set  base $w.left.frame3.frame35
    label   $base.label  -text "Slowdown display (msec): "
    entry   $base.input  -width 10
    button  $base.apply  -text "Apply" \
               -command "gomApplySnail2Loop $base.input"
    pack    $base.label -side left  -anchor w
    pack    $base.input -side left  -anchor w
    pack    $base.apply -side left  -anchor w -padx 3

    if {[info exists gomDisplayLoopSleep]} {
            $base.input delete 0 end
            $base.input inser  0 $gomDisplayLoopSleep
    }

    bind $base.input <Return> "gomApplySnail2Loop $base.input"

#    frame $w.left.frame3.frame35
#    set base $w.left.frame3.frame35
#   label $base.currlabel  -text "Current frame:" -width 15
#	entry $base.current       -width 20
#	pack  $base -anchor w -side top
#	pack  $base.currlabel $base.current  -side left

#    frame $w.left.frame2 -borderwidth 2 -relief ridge
    labelframe $w.left.frame2 -text "Trajectory control" -borderwidth 2 -relief ridge -padx 2 -pady 2
    set base $w.left.frame2
    pack   $base -side left -anchor w

	button $base.first    -text "|<"   -command "trajectory loop first"    \
            -font "$gomControlFontHuge"
	pack   $base.first    -side left  -padx 3 -pady 4 
	button $base.back     -text "<"    -command "trajectory loop backward" \
            -font "$gomControlFontHuge"
	pack   $base.back     -side left  -padx 3 -pady 4
	button $base.stop     -text "\[\]" -command "trajectory loop stop"     \
            -font "$gomControlFontHuge"
   	pack   $base.stop     -side left  -padx 3 -pady 4
	button $base.play     -text "|>"   -command "trajectory loop play"     \
            -font "$gomControlFontHuge"
	pack   $base.play     -side left  -padx 3 -pady 4
	button $base.forward  -text ">"    -command "trajectory loop forward"  \
            -font "$gomControlFontHuge"
	pack   $base.forward  -side left  -padx 3 -pady 4
	button $base.last     -text ">|"   -command "trajectory loop last"     \
            -font "$gomControlFontHuge"
	pack   $base.last     -side left  -padx 3 -pady 4

   if {$TrajectoryDefined} {
      $w.left.frame1.filename delete 0 end
      $w.left.frame1.filename insert 0 [show trajectory filename]

      $w.left.frame3.frame31.frames configure -state normal
      $w.left.frame3.frame31.frames delete 0 end
      $w.left.frame3.frame31.frames insert 0 [show trajectory frames]
      $w.left.frame3.frame31.frames configure -state disabled

      scan [show trajectory display] "%d %d %d" First Last Step

      $w.left.frame3.frame32.first delete 0  end
      $w.left.frame3.frame32.first insert 0 $First

      $w.left.frame3.frame33.last delete 0  end
      $w.left.frame3.frame33.last insert 0 $Last

      $w.left.frame3.frame34.step delete 0  end
      $w.left.frame3.frame34.step insert 0 $Step

#      $w.left.frame3.frame35.current delete 0  end
#     $w.left.frame3.frame35.current insert 0 [show trajectory current]

   }
      frame $w.left.frame4 
      pack  $w.left.frame4 -side top -anchor e -pady 5 -padx 10

#      frame $w.left.frame4.conn -borderwidth 2 -relief ridge
      labelframe $w.left.frame4.conn -text "Atom reconnection" -borderwidth 2 -relief ridge -padx 2 -pady 2
      pack       $w.left.frame4.conn -side top -anchor e -pady 2 -padx 2

#      label   $w.left.frame4.conn.label -text "Atom reconnection:   "
#      pack    $w.left.frame4.conn.label -side top  -anchor w
      radiobutton $w.left.frame4.conn.on  -text "On"  -value 1 -variable gomConnect  \
                  -command  {define atom reconnect on;lulInvalidateDisplay}
      radiobutton $w.left.frame4.conn.off -text "Off" -value 0 -variable gomConnect  \
                  -command  {define atom reconnect off;lulInvalidateDisplay}
      pack    $w.left.frame4.conn.on $w.left.frame4.conn.off          \
                  -side top -anchor w

      if {![show atom reconnect]} {
              $w.left.frame4.conn.off select
      } else {
              $w.left.frame4.conn.on select
      }
#
#      frame $w.left.frame4.hbond -borderwidth 2 -relief ridge
      labelframe $w.left.frame4.hbond -text "Hbond recalculation" -borderwidth 2 -relief ridge -padx 2 -pady 2
      pack       $w.left.frame4.hbond -side top -anchor e -pady 2 -padx 2

#      label   $w.left.frame4.hbond.label -text "Hbond recalculation:   "
#      pack    $w.left.frame4.hbond.label -side top  -anchor w
      radiobutton $w.left.frame4.hbond.on  -text "On"  -value 1 -variable gomHbondConnect  \
                  -command  {define atom hbreconnect on;lulInvalidateDisplay}
      radiobutton $w.left.frame4.hbond.off -text "Off" -value 0 -variable gomHbondConnect  \
                  -command  {define atom hbreconnect off;lulInvalidateDisplay}
      pack    $w.left.frame4.hbond.on $w.left.frame4.hbond.off          \
                  -side top -anchor w

      if {![show atom hbreconnect]} {
              $w.left.frame4.hbond.off select
      } else {
              $w.left.frame4.hbond.on select
      }
#
#      frame $w.left.frame4.method -borderwidth 2 -relief ridge
      labelframe $w.left.frame4.method -text "Trajectory retrieval" -borderwidth 2 -relief ridge -padx 2 -pady 2
      pack       $w.left.frame4.method -side top -anchor w -pady 2 -padx 2 -fill x

#      label   $w.left.frame4.method.label -text "Trajectory retrieval:    "
#      pack    $w.left.frame4.method.label -side top  -anchor w
      radiobutton $w.left.frame4.method.fast  -text "Fast"  -value 1  -variable gomMethod \
                  -command  {trajectory action 0}
      radiobutton $w.left.frame4.method.slow  -text "Slow"  -value 0  -variable gomMethod \
                  -command  {trajectory action 1}
      pack    $w.left.frame4.method.fast $w.left.frame4.method.slow          \
                  -side top -anchor w

      if {![show trajectory action]} {
              $w.left.frame4.method.fast select
      } else {
              $w.left.frame4.method.slow select
      }

# look for state 

    if { $gomTrajType == "xmol"      ||
         $gomTrajType == "tinker"    ||
         $gomTrajType == "gromos96a" ||
         $gomTrajType == "fdl_poly"} {
         eval {lulTouchTrajectoryMethodSelector on}
    } else {
         eval {lulTouchTrajectoryMethodSelector off}
    }

#
# trajectory frame number display

#      frame $w.left.frame5 -borderwidth 2 -relief ridge
      labelframe $w.left.frame5 -text "Display frame nr" -borderwidth 2 -relief ridge -padx 2 -pady 2
      pack       $w.left.frame5 -side top -anchor w -pady 2 -padx 2

#      label   $w.left.frame5.label -text "Display frame number:"
#      pack    $w.left.frame5.label -side top  -anchor w
      radiobutton $w.left.frame5.on  -text "On"  -value 1 -variable gomFrameDisplayID  \
                  -command  {define trajectory fid on; lulInvalidateDisplay}
      radiobutton $w.left.frame5.off -text "Off" -value 0 -variable gomFrameDisplayID  \
                  -command  {define trajectory fid off; lulInvalidateDisplay}
      pack    $w.left.frame5.on $w.left.frame5.off          \
                  -side top -anchor w

      if {![show trajectory fid]} {
              $w.left.frame5.off select
      } else {
              $w.left.frame5.on  select
      }
#


}

############################################################################
# PROC
proc gomApplySnail2Loop { w } {

     global gomDisplayLoopSleep

  set Input [string trim [$w get]]
  if {$Input != ""} {
     set gomDisplayLoopSleep $Input
  } else {
     gomError "you have to supply a value in milliseconds!"
  } 

}
############################################################################
# PROC
proc lulTouchTrajectoryMethodSelector { Value } {

     if {$Value == "on"} {
        .gomtraj.left.frame4.method.fast configure -state normal
        .gomtraj.left.frame4.method.slow configure -state normal
     } elseif { $Value == "off"} {
        .gomtraj.left.frame4.method.fast configure -state disabled
        .gomtraj.left.frame4.method.slow configure -state disabled
     }
}


############################################################################
proc lulGetTrajFile { w } {


    global gomTrajType
	global gomTrajFileName


    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    
	if {$gomTrajType == "amber"} {
       set types {
         {"Amber trajectory"		"*"		       TEXT}
         {"All files"		*}
       }
       set def_type "amber"
	} elseif {$gomTrajType == "famber"} {
       set types {
         {"Amber(f) trajectory"		"*"		       TEXT}
         {"All files"		*}
       }
       set def_type "famber"
    } elseif {$gomTrajType == "cerius2"} {
       set types {
         {"CERIUS2 trajectory"	".trj .TRJ"			TEXT}
         {"All files"		*}
       }
       set def_type "cerius2"
    } elseif {$gomTrajType == "charmm"} {
       set types {
         {"CHARMM trajectory"	".dcd .DCD"			TEXT}
         {"All files"		*}
       }
       set def_type "charmm"
    } elseif {$gomTrajType == "discover"} {
       set types {
         {"DISCOVER trajectory"	"*"         	TEXT}
         {"All files"		*}
       }
       set def_type "discover"
    } elseif {$gomTrajType == "fdl_poly"} {
       set types {
         {"Formatted DL_Poly trajectory"	"*"         	TEXT}
         {"All files"		*}
       }
       set def_type "fdl_poly"
    } elseif {$gomTrajType == "udl_poly"} {
       set types {
         {"Unformatted DL_Poly trajectory"	"*"         	TEXT}
         {"All files"		*}
       }
       set def_type "udl_poly"
    } elseif {$gomTrajType == "gromacs"} {
       set types {
         {"GROMACS trajectory"	".trj .trr .xtc .TRJ .TRR .XTC"         	TEXT}
         {"All files"		*}
       }
       set def_type "gromacs"
    } elseif {$gomTrajType == "gromos"} {
       set types {
         {"GROMOS trajectory"	"*"         	TEXT}
         {"All files"		*}
       }
       set def_type "gromos"
    } elseif {$gomTrajType == "gromos96a"} {
       set types {
         {"GROMOS trajectory"	".dat .DAT"    	TEXT}
         {"All files"		*}
       }
       set def_type "gromos96a"
    } elseif {$gomTrajType == "hyperchem"} {
       set types {
         {"HyperChem trajectory"	".snp .SNP"        	TEXT}
         {"All files"		*}
       }
       set def_type "hyperchem"
    } elseif {$gomTrajType == "mumod"} {
       set types {
         {"MUMOD trajectory"	"*"				TEXT}
         {"All files"		*}
       }
       set def_type "mumod"
    } elseif {$gomTrajType == "tinker"} {
       set types {
         {"TINKER trajectory"	".arc .ARC"				TEXT}
         {"All files"		*}
       }
       set def_type "tinker"
    } elseif {$gomTrajType == "xmol"} {
       set types {
         {"XMOL trajectory"	".xyz .xmol .XYZ .XMOL"	TEXT}
         {"All files"		*}
       }
       set def_type "xmol"
    } elseif {$gomTrajType == "xplor"} {
       set types {
         {"XPLOR trajectory"	".dcd .DCD"      TEXT}
         {"All files"		*}
       }
       set def_type "xplor"
    } elseif {$gomTrajType == "yasp"} {
       set types {
         {"YASP trajectory"	"*"	                 TEXT}
         {"All files"		*}
       }
       set def_type "yasp"
    } 

    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .dic -title "Trajectory file"]]

    if [string compare $file ""] {
    $w delete 0 end
	$w insert 0 $file

    set gomTrajFileName $file
    }
}

############################################################################
# PROC
proc lulGetFileData {w} {

   global gomTrajFileName
   global gomTrajType
   global gomMethod

   trajectory file $gomTrajType "$gomTrajFileName"  

#   set TrajectoryDefined [show trajectory defined]

   $w.left.frame3.frame31.frames configure -state normal
   $w.left.frame3.frame31.frames delete 0 end
   $w.left.frame3.frame31.frames insert 0 [show trajectory frames]
   $w.left.frame3.frame31.frames configure -state disabled

   scan [show trajectory display] "%d %d %d" First Last Step

   $w.left.frame3.frame32.first delete 0  end
   $w.left.frame3.frame32.first insert 0 $First

   $w.left.frame3.frame33.last delete 0  end
   $w.left.frame3.frame33.last insert 0 $Last

   $w.left.frame3.frame34.step delete 0  end
   $w.left.frame3.frame34.step insert 0 $Step

#   $w.left.frame3.frame35.current delete 0  end
#   $w.left.frame3.frame35.current insert 0 [show trajectory current]

# change to current directory
    lulChangeDirectory "$gomTrajFileName"
}

############################################################################
# PROC
proc lulDoTrajFileParams { w } {

 set First [$w.left.frame3.frame32.first get]
 set Last  [$w.left.frame3.frame33.last  get]
 set Step  [$w.left.frame3.frame34.step  get]

 eval trajectory limits $First $Last $Step

 puts "Trajectory display limits changed to:"
 puts "First: $First, Last: $Last, Step: $Step"
}
############################################################################
proc lulUpdateCurrentFrameDisplay { frame } {

   $w.left.frame3.frame35.current delete 0  end
   $w.left.frame3.frame35.current insert 0 [show trajectory current]
}


##################### Superimpose  ############################
#
# define superimpose
#
proc lulSuperimposeStructures {} {

     global gomControlFont
     global gomHelpFile

# return if no molecular systems defined
     if {[show molstructures] < 2} {
          lulErrorDialog {ERROR: structures < 2. Read first in 2 structures!}
	      return
	 }

set w .gomimpose
catch {destroy $w}
toplevel $w 
wm title $w "Superimpose"
wm iconname $w "Superimpose"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulDoSuperimpose $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(simpose)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Superimpose:"
pack   $w.label -side top -anchor w

# structure #1
frame  $w.structure1 -borderwidth 2 -relief ridge
pack   $w.structure1 -side top -anchor w
frame  $w.structure1.topf
pack   $w.structure1.topf -side top -anchor w
label  $w.structure1.topf.label -text "Structure #" -width 10
pack   $w.structure1.topf.label -side left -anchor w
entry  $w.structure1.topf.numb  -width 5
pack   $w.structure1.topf.numb  -side left -anchor w 

$w.structure1.topf.numb delete 0 end
$w.structure1.topf.numb insert  0 "1"  

frame  $w.structure1.frame1 
label  $w.structure1.frame1.segment -text "Segment:" -width 10
entry  $w.structure1.frame1.segment_input -width 30 \
       -textvariable segment_simpose1 
pack   $w.structure1.frame1 -side top -anchor w
pack   $w.structure1.frame1.segment \
       $w.structure1.frame1.segment_input -side left

if { [$w.structure1.frame1.segment_input get] == "" } {
      $w.structure1.frame1.segment_input insert 0 "*"}

frame  $w.structure1.frame2 
label  $w.structure1.frame2.residue -text "Residue:" -width 10
entry  $w.structure1.frame2.residue_input -width 30 \
       -textvariable residue_simpose1 
pack   $w.structure1.frame2 -side top -anchor w
pack   $w.structure1.frame2.residue \
       $w.structure1.frame2.residue_input -side left 

if { [$w.structure1.frame2.residue_input get] == ""} {
      $w.structure1.frame2.residue_input insert 0 "*"}

frame  $w.structure1.frame3 
label  $w.structure1.frame3.atom -text "Atom:" -width 10
entry  $w.structure1.frame3.atom_input -width 30 \
       -textvariable atom_simpose1 
pack   $w.structure1.frame3 -side top -anchor w
pack   $w.structure1.frame3.atom \
       $w.structure1.frame3.atom_input -side left

if { [$w.structure1.frame3.atom_input get] == "" } {
      $w.structure1.frame3.atom_input insert 0 "*"}

# structure #2

frame  $w.structure2 -borderwidth 2 -relief ridge
pack   $w.structure2 -side top -anchor w
frame  $w.structure2.topf
pack   $w.structure2.topf -side top -anchor w
label  $w.structure2.topf.label -text "Structure #" -width 10
pack   $w.structure2.topf.label -side left -anchor w
entry  $w.structure2.topf.numb  -width 5
pack   $w.structure2.topf.numb  -side left -anchor w 

$w.structure2.topf.numb delete 0 end
$w.structure2.topf.numb insert  0 "2"  

frame  $w.structure2.frame1 
label  $w.structure2.frame1.segment -text "Segment:" -width 10
entry  $w.structure2.frame1.segment_input -width 30 \
       -textvariable segment_simpose2 
pack   $w.structure2.frame1 -side top -anchor w
pack   $w.structure2.frame1.segment \
       $w.structure2.frame1.segment_input -side left

if { [$w.structure2.frame1.segment_input get] == "" } {
      $w.structure2.frame1.segment_input insert 0 "*"}

frame  $w.structure2.frame2 
label  $w.structure2.frame2.residue -text "Residue:" -width 10
entry  $w.structure2.frame2.residue_input -width 30 \
       -textvariable residue_simpose2 
pack   $w.structure2.frame2 -side top -anchor w
pack   $w.structure2.frame2.residue \
       $w.structure2.frame2.residue_input -side left 

if { [$w.structure2.frame2.residue_input get] == ""} {
      $w.structure2.frame2.residue_input insert 0 "*"}

frame  $w.structure2.frame3 
label  $w.structure2.frame3.atom -text "Atom:" -width 10
entry  $w.structure2.frame3.atom_input -width 30 \
       -textvariable atom_simpose2 
pack   $w.structure2.frame3 -side top -anchor w
pack   $w.structure2.frame3.atom \
       $w.structure2.frame3.atom_input -side left

if { [$w.structure2.frame3.atom_input get] == "" } {
      $w.structure2.frame3.atom_input insert 0 "*"}

# display results

frame  $w.frame3 -borderwidth 2 -relief ridge
label       $w.frame3.label -text "Display result:"
radiobutton $w.frame3.on  -text "On"  \
            -variable gomSimposeAction -value "on" 
radiobutton $w.frame3.off -text "Off" \
            -variable gomSimposeAction -value "off"
pack $w.frame3 -side top -anchor w 
pack $w.frame3.label -side top -anchor w
pack $w.frame3.on $w.frame3.off -side top -anchor w
$w.frame3.on select


}

###############################################################
proc lulDoSuperimpose { w } {

     global gomSimposeAction

# structure #1
     set Struct1  [string trim [$w.structure1.topf.numb get]]
     set Segment1 [string trim [$w.structure1.frame1.segment_input get]]
	 if {$Segment1 == ""} {set Segment1 "*"}
     set Residue1 [string trim [$w.structure1.frame2.residue_input get]]
	 if {$Residue1 == ""} {set Residue1 "*"}
     set Atom1    [string trim [$w.structure1.frame3.atom_input    get]]
	 if {$Atom1    == ""} {set Atom1 "*"}

# structure #2
     set Struct2  [string trim [$w.structure2.topf.numb get]]
     set Segment2 [string trim [$w.structure2.frame1.segment_input get]]
	 if {$Segment2 == ""} {set Segment2 "*"}
     set Residue2 [string trim [$w.structure2.frame2.residue_input get]]
	 if {$Residue2 == ""} {set Residue2 "*"}
     set Atom2    [string trim [$w.structure2.frame3.atom_input    get]]
	 if {$Atom2    == ""} {set Atom2 "*"}

     if {$Struct1 == $Struct2} {
         gomError "Structure #1 == Structure # 2"
         return
     }

     set Command [concat calculate quatfit          \
	                  $Struct1  \{$Segment1\} \{$Residue1\} \{$Atom1\} \
	                  $Struct2  \{$Segment2\} \{$Residue2\} \{$Atom2\} \
						 $gomSimposeAction]

	 eval $Command
	 lulInvalidateDisplay

}

##################### Trace display ############################
#
# define trace display
#
proc lulDefineTraceAtoms {} {

     global gomControlFont
     global gomAppend
     global gomHelpFile

# return if no molecular systems defined
     if {[show trajectory frames] < 1} {
          lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
	      return
	 }

set w .gomtraceatoms
catch {destroy $w}
toplevel $w 
wm title $w "Trace atoms"
wm iconname $w "Trace atoms"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulDoTraceAction $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(traceatoms)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Define trace atoms:"
pack   $w.label -side top -anchor w

frame       $w.frame8 -borderwidth 2 -relief ridge
button      $w.frame8.delete -text "Delete" \
            -command "ptrace -atoms;lulInvalidateDisplay"

pack        $w.frame8 -side bottom -anchor e
pack        $w.frame8.delete -side left

frame       $w.frame9 -borderwidth 2 -relief ridge
checkbutton $w.frame9.append -text "Append" -variable gomAppend 

pack        $w.frame9 -side bottom -anchor e
pack        $w.frame9.append -side left

$w.frame9.append deselect

frame  $w.frame1 
label  $w.frame1.segment -text "Segment:" -width 10
entry  $w.frame1.segment_input -width 30 -textvariable gomSegmentTrace 
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.segment $w.frame1.segment_input -side left

if { [$w.frame1.segment_input get] == "" } {$w.frame1.segment_input insert 0 "*"}

frame  $w.frame2 
label  $w.frame2.residue -text "Residue:" -width 10
entry  $w.frame2.residue_input -width 30 -textvariable gomResidueTrace 
pack   $w.frame2 -side top -anchor w
pack   $w.frame2.residue $w.frame2.residue_input -side left 

if { [$w.frame2.residue_input get] == ""} {$w.frame2.residue_input insert 0 "*"}

frame  $w.frame3 
label  $w.frame3.atom -text "Atom:" -width 10
entry  $w.frame3.atom_input -width 30 -textvariable gomAtomTrace 
pack   $w.frame3 -side top -anchor w
pack   $w.frame3.atom $w.frame3.atom_input -side left

if { [$w.frame3.atom_input get] == "" } {$w.frame3.atom_input insert 0 "*"}

frame       $w.frame4 -borderwidth 2 -relief ridge
label       $w.frame4.label -text "Display:"
radiobutton $w.frame4.on  -text "On"  \
            -variable gomTraceAction -value 1 \
			-command  "ptrace display on;lulInvalidateDisplay" 
radiobutton $w.frame4.off -text "Off" \
            -variable gomTraceAction -value 0 \
			-command  "ptrace display off;lulInvalidateDisplay" 
pack $w.frame4 -side left 
pack $w.frame4.label $w.frame4.on $w.frame4.off \
     -side top -anchor w
$w.frame4.off select

}

########################################################################
proc lulDoTraceAction { w } {

    global gomAppend

    set Segment [string trim [$w.frame1.segment_input get]]
	if {$Segment == ""} {set Segment "*"}
	set Residue [string trim [$w.frame2.residue_input get]]
	if {$Residue == ""} {set Residue "*"}
	set Atom    [string trim [$w.frame3.atom_input    get]]
	if {$Atom    == ""} {set Atom    "*"}

    if {$gomAppend} {
	     set Append "append"
	} else {
	     set Append ""
    }

    set Command [concat ptrace atoms \{$Segment\} \{$Residue\} \{$Atom\} $Append]
	eval $Command
	lulInvalidateDisplay

}

####################################################################
proc lulTraceDisplayState { state } {

    if {$state} {
       if {[show traces]} {ptrace display on}
    } else {
	   ptrace display off
	}
}

####################################################################
# 
# identify an atom
#
# PROC
proc lulIdentifyAtom { } {

     global gomControlFont
     global gomHelpFile
     global gomPickSwitch


# return if no molecular systems defined
     if {[show molstructures] < 1} {
	      return
     }

set w .gomidentifyatom
catch {destroy $w}
toplevel $w 
wm title $w "Identify atom"
wm iconname $w "Identify atom"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command \
         "eval define atom identify off;destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command {gomPrint "Identify atom is disabled"} -state disabled
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(identifyatom)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

frame $w.frame1              -borderwidth 2 -relief ridge
label $w.frame1.segmentlabel -text "Segment: "
entry $w.frame1.segment      -width 10
label $w.frame1.residuelabel -text "Residue: "
entry $w.frame1.residue      -width 10
label $w.frame1.atomlabel    -text "Atom: "
entry $w.frame1.atom         -width 10

pack $w.frame1 -side top
pack $w.frame1.segmentlabel $w.frame1.segment -side left
pack $w.frame1.residuelabel $w.frame1.residue -side left
pack $w.frame1.atomlabel    $w.frame1.atom    -side left

frame $w.frame4 -borderwidth 2 -relief ridge
label $w.frame4.label -text "Gaussian basis set tag:"
entry $w.frame4.gbasis       -width 30

pack  $w.frame4 -side top          -anchor w
pack  $w.frame4.label -side top    -anchor w
pack  $w.frame4.gbasis -side top   -anchor w

frame $w.frame5 -borderwidth 2 -relief ridge
label $w.frame5.label -text "Nuclear charge:"
entry $w.frame5.nuccharge -width 15
pack  $w.frame5 -side right -anchor n
pack  $w.frame5.label -side top -anchor w
pack  $w.frame5.nuccharge -side top -anchor w

frame $w.frame3 -borderwidth 2 -relief ridge
label $w.frame3.xlabel -text "X-Coord:"
entry $w.frame3.x      -width 15
label $w.frame3.ylabel -text "Y-Coord:"
entry $w.frame3.y      -width 15
label $w.frame3.zlabel -text "Z-Coord:"
entry $w.frame3.z      -width 15
pack  $w.frame3 -side left
pack  $w.frame3.xlabel -side top -anchor w
pack  $w.frame3.x -side top
pack  $w.frame3.ylabel -side top -anchor w
pack  $w.frame3.y -side top
pack  $w.frame3.zlabel -side top -anchor w
pack  $w.frame3.z -side top

frame $w.frame2 -borderwidth 2 -relief ridge
label $w.frame2.label -text "Identify atom(s):"
radiobutton $w.frame2.on  -text "On"  -variable gomPickSwitch -value 1 \
            -command {lulActivateIdentifyAtom 1}
radiobutton $w.frame2.off -text "Off" -variable gomPickSwitch -value 0 \
            -command {lulActivateIdentifyAtom 0}

pack  $w.frame2 -side right -anchor se 
pack  $w.frame2.label -side top
pack  $w.frame2.on    -side top
pack  $w.frame2.off   -side top

set State [show atom identify]

  if {$State} {
     $w.frame2.on select
  } else {
     $w.frame2.off select
  }
}

####################################################################
# 
# PROC
proc lulActivateIdentifyAtom { value } {

  if {$value} {
      eval "define atom identify on"
  } else {
      eval "define atom identify off"
  }
}

####################################################################
# 
# PROC
proc lulPassAtomID2GUI { sid aid } {

# this is taken from the atom identify widget
   set w .gomidentifyatom
   set status [winfo exists $w]
   if {!$status} {
#          lulErrorDialog {ERROR: Identify atom widget does not exist. Open the Atom identify widget!}
	      return      
   }

   set SegName [show atom segment $aid $sid]
   $w.frame1.segment delete 0 end
   $w.frame1.segment insert 0 $SegName

   set ResName [show atom residue $aid $sid]
   $w.frame1.residue delete 0 end
   $w.frame1.residue insert 0 $ResName

   set AtmName [show atom atom $aid $sid]
   $w.frame1.atom    delete 0 end
   $w.frame1.atom    insert 0 $AtmName

   scan [show atom coord $aid $sid] "%f %f %f"  Xc Yc Zc

   $w.frame3.x    delete 0 end
   $w.frame3.x    insert 0 $Xc

   $w.frame3.y    delete 0 end
   $w.frame3.y    insert 0 $Yc

   $w.frame3.z    delete 0 end
   $w.frame3.z    insert 0 $Zc

   set GBasis [show atom gbasis $aid $sid]

   $w.frame4.gbasis    delete 0 end
   $w.frame4.gbasis    insert 0 $GBasis

   set NuclearCharge [show atom nuclear $aid $sid]

   $w.frame5.nuccharge delete 0 end
   $w.frame5.nuccharge insert 0 $NuclearCharge

# change window title line to reflect atom number
   if {$AtmName > 0 } {
      wm title $w "Identify atom (#$sid:$aid)"
      }

}

##################### hardcopy ############################
#
# define hardcopy
#
proc lulMakeHardcopy {} {

     global gomControlFont
     global gomHardcopyType
     global gomHardcopyExt
     global gomHelpFile
     global gomPostOrient

set w .gomhardcopy
catch {destroy $w}
toplevel $w 
wm title $w "Hardcopy"
wm iconname $w "Hardcopy"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulApplyHardcopy $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(hardcopy)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Hardcopy:"
pack   $w.label -side top -anchor w

frame        $w.frame1 -borderwidth 2 -relief ridge 
label        $w.frame1.label -text "File formats: "
radiobutton  $w.frame1.bmp -text "BMP"        -variable gomHardcopyType -value bmp \
              -command "lulPostScriptOrientationControl $w off"
radiobutton  $w.frame1.jpg -text "JPG"        -variable gomHardcopyType -value jpeg \
              -command "lulPostScriptOrientationControl $w off"
radiobutton  $w.frame1.ps  -text "PostScript" -variable gomHardcopyType -value postscript \
              -command "lulPostScriptOrientationControl $w on"
radiobutton  $w.frame1.rgb -text "RGB"        -variable gomHardcopyType -value rgb \
              -command "lulPostScriptOrientationControl $w off"
radiobutton  $w.frame1.tga -text "TGA"        -variable gomHardcopyType -value tga \
              -command "lulPostScriptOrientationControl $w off"
radiobutton  $w.frame1.xwd -text "XWD"        -variable gomHardcopyType -value xwd \
              -command "lulPostScriptOrientationControl $w off"
               
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.bmp $w.frame1.jpg $w.frame1.ps $w.frame1.rgb $w.frame1.tga $w.frame1.xwd -side left

frame  $w.frame2
label  $w.frame2.label     -text "File name: "
entry  $w.frame2.filename  -width 40
button $w.frame2.browse    -text "Browse..." -command "lulDoHardcopyAction $w"

pack   $w.frame2 -side top -pady 6
pack   $w.frame2.label $w.frame2.filename -side left
pack   $w.frame2.browse -padx 3 -side left
bind   $w.frame2.filename <Return> "eval lulApplyHardcopy $w"

$w.frame1.ps select

frame       $w.frame3 -borderwidth 2 -relief ridge
pack        $w.frame3 -side top -anchor w -pady 4
label       $w.frame3.orientation -text "Orientation:"
pack        $w.frame3.orientation -side top -anchor w
radiobutton $w.frame3.portrait  -text "Portrait"  -variable gomPostOrient -value 0
radiobutton $w.frame3.landscape -text "Landscape" -variable gomPostOrient -value 1
pack        $w.frame3.portrait $w.frame3.landscape -side top -anchor w

$w.frame3.portrait select

frame       $w.frame4 -borderwidth 2 -relief ridge
pack        $w.frame4 -side top -anchor w -pady 4
label       $w.frame4.label -text "Choose window:"
pack        $w.frame4.label -side top -anchor w
pack        $w.frame4.label -side top -anchor w

for {set i 1} {$i <= [show window defined]} {incr i} { 
    radiobutton $w.frame4.win$i  -text "Window: ($i)"  -variable gomFocusWindow -value $i
    pack        $w.frame4.win$i  -side top -anchor w
}

$w.frame4.win1 select
 
}
########################################################################
# PROC
proc lulPostScriptOrientationControl { w action} {

  if {$action == "on"} {
      $w.frame3.portrait  configure -state normal
      $w.frame3.landscape configure -state normal
  } elseif {$action == "off"} {
      $w.frame3.portrait  configure -state disabled
      $w.frame3.landscape configure -state disabled
  }
}
########################################################################
# PROC
proc lulDoHardcopyAction { w } {

     global gomHardcopyType
	 global gomHardcopyExt
     global gomPostOrient

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    if {$gomHardcopyType == "bmp"} {
       set types {
         {"BMP files"		".bmp .BMP"		       TEXT}
         {"All files"		*}
       }
    } elseif {$gomHardcopyType == "jpeg"} {
       set types {
         {"PostScript files"	".jpg .JPG"		TEXT}
         {"All files"		*}
       }
    } elseif {$gomHardcopyType == "postscript"} {
       set types {
         {"PostScript files"	".ps .PS"		TEXT}
         {"All files"		*}
       }
    } elseif {$gomHardcopyType == "rgb"} {
       set types {
         {"SGI RGB files"	".rgb .RGB"         	TEXT}
         {"All files"		*}
       }
    } elseif {$gomHardcopyType == "tga"} {
       set types {
         {"XWD files"	".tga .TGA"         	TEXT}
         {"All files"		*}
       }
    } elseif {$gomHardcopyType == "xwd"} {
       set types {
         {"XWD files"	".xwd .XWD"         	TEXT}
         {"All files"		*}
       }
    }

    set file [string trim [tk_getSaveFile -filetypes $types -parent $w   \
         -defaultextension $gomHardcopyExt($gomHardcopyType)\
		 -title "Save Hardcopy File"]]

    if [string compare $file ""] {
    $w.frame2.filename delete 0 end
	$w.frame2.filename insert 0 $file
    }

}
###########################################################################
# PROC
proc lulApplyHardcopy { w } {

     global gomHardcopyType
     global gomPostOrient
     global gomFocusWindow

     set File [$w.frame2.filename get]

     if {$gomPostOrient} {
         set Extra "\{-l\}"
     } else {
         set Extra ""
     }

	 if {$File != ""} {
	   puts "Saving file: $File in format $gomHardcopyType ..."
       eval hardcopy $gomFocusWindow $gomHardcopyType \"$File\" $Extra
	   puts "Done!"
# change to current directory
       lulChangeDirectory "$File"
     }
}

####################################################################
# 
# select atoms ...
#
proc lulSelectAtoms { } {

     global gomHelpFile
     global gomStructure
     global gomControlFont
     global gomStructure

# return if no molecular systems defined
     if {[show molstructures] < 1} {
	      gomError "no structures available, import first one"
          return
	 }

set w .gomselectatoms
catch {destroy $w}
toplevel $w 
wm title $w "Select atoms"
wm iconname $w "Select"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplySelection $w"
button $w.buttons.unselect -text "Unselect all" -font "$gomControlFont" \
        -command "select -atoms;lulInvalidateDisplay"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(selectatom)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.unselect $w.buttons.help -side left -expand 1 -padx 4

label  $w.label -text "Select atoms in:" -font {Arial 12 bold}
pack   $w.label -side top -anchor w -pady 4

set  NumStruct [show molstructures]

frame $w.frame0 -borderwidth 2 -relief ridge
pack $w.frame0 -side top

for { set i 0} { $i < $NumStruct} {incr i} {

radiobutton $w.frame0.structure$i -text "Structure nr: [expr $i + 1]" \
            -value [expr $i + 1] -variable gomStructure

pack $w.frame0.structure$i -side top
}

set gomStructure 1

lulCreateAtomInputEntries $w 30 "*" "*" "*" "" 1

}

####################################################################
# PROC
proc lulApplySelection { w } {

     global gomStructure

     set Segment [string trim [$w.segment.input get]]
     set Residue [string trim [$w.residue.input get]]
     set Atom    [string trim [$w.atom.input    get]]

     select atoms $Segment $Residue $Atom $gomStructure

     lulInvalidateDisplay
}
####################################################################
# 
# define cell ...
#
# PROC
proc lulDefineCell { } {

     global gomHelpFile
     global gomSwitch
     global gomCellColour
     global gomControlFont

# return if no molecular systems defined
     set  NumStruct [show molstructures]
     if {$NumStruct < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomdefinecell
catch {destroy $w}
toplevel $w 
wm title $w "Define Cell"
wm iconname $w "Cell"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyDefineCell $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         {htmlShowHelp $gomHelpFile(definecell)}
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

scan [show cell dimension] "%e %e %e" a b c

# cell sides
frame $w.frame1              -borderwidth 2 -relief ridge
pack  $w.frame1 -side top    -anchor w

label $w.frame1.labela       -text "a:" -width 10
entry $w.frame1.a            -width 12
pack  $w.frame1.labela $w.frame1.a -side left
# default
$w.frame1.a insert 0 [format "%.4e" $a]

label $w.frame1.labelb       -text "b:" -width 10
entry $w.frame1.b            -width 12
pack  $w.frame1.labelb $w.frame1.b -side left
# default
$w.frame1.b insert 0 [format "%.4e" $b]

label $w.frame1.labelc       -text "c:" -width 10
entry $w.frame1.c            -width 12
pack  $w.frame1.labelc $w.frame1.c    -side left
# default
$w.frame1.c insert 0 [format "%.4e" $c]

# cell translation
frame $w.frame5              -borderwidth 2 -relief ridge
pack  $w.frame5 -side top    -anchor w

label $w.frame5.labela       -text "X trans:" -width 10
entry $w.frame5.xtrans       -width 12
pack  $w.frame5.labela $w.frame5.xtrans -side left
# default
$w.frame5.xtrans insert 0 [format "%.4e" [lindex [show cell trans] 0]]

label $w.frame5.labelb       -text "Y trans:" -width 10
entry $w.frame5.ytrans            -width 12
pack  $w.frame5.labelb $w.frame5.ytrans -side left
# default
$w.frame5.ytrans insert 0 [format "%.4e" [lindex [show cell trans] 1]]

label $w.frame5.labelc       -text "Z trans:" -width 10
entry $w.frame5.ztrans            -width 12
pack  $w.frame5.labelc $w.frame5.ztrans    -side left
# default
$w.frame5.ztrans insert 0 [format "%.4e" [lindex [show cell trans] 2]]

# cell angles
frame $w.frame2              -borderwidth 2 -relief ridge
pack  $w.frame2 -side top    -anchor w

label $w.frame2.labelalpha   -text "alpha:" -width 10
entry $w.frame2.alpha        -width 12
pack  $w.frame2.labelalpha $w.frame2.alpha -side left

$w.frame2.alpha insert 0 "90.0"
$w.frame2.alpha configure  -state disabled


label $w.frame2.labelbeta    -text "beta:" -width 10
entry $w.frame2.beta         -width 12
pack  $w.frame2.labelbeta $w.frame2.beta -side left

$w.frame2.beta insert 0 "90.0"
$w.frame2.beta configure  -state disabled

label $w.frame2.labelgamma   -text "gamma:" -width 10
entry $w.frame2.gamma        -width 12
pack  $w.frame2.labelgamma $w.frame2.gamma    -side left

$w.frame2.gamma insert 0 "90.0"
$w.frame2.gamma configure  -state disabled

frame $w.frame3 -borderwidth 2 -relief ridge
label $w.frame3.label -text "Display state:"
radiobutton $w.frame3.on  -text "On"  -variable gomSwitch  -value 1 \
            -command "plot cell on;lulInvalidateDisplay"
radiobutton $w.frame3.off -text "Off" -variable gomSwitch  -value 0 \
            -command "plot cell off;lulInvalidateDisplay"
pack $w.frame3 -side top    -anchor w
pack $w.frame3.label  -side top
pack $w.frame3.on  -side top
pack $w.frame3.off -side top
$w.frame3.off select

scan [show cell colour] "%f %f %f" red green blue
set gomCellColour [lulColourFloat2Hex $red $green $blue]

frame  $w.frame4 -borderwidth 2 -relief ridge
button $w.frame4.colour -text "Colour..." -bg $gomCellColour \
       -command "lulDefineCellColour $w"
label  $w.frame4.linewlabel -text "Line width: "
entry  $w.frame4.linewidth  -width 5
pack   $w.frame4 -side left
pack   $w.frame4.colour -side left
pack   $w.frame4.linewlabel $w.frame4.linewidth -side left

$w.frame4.linewidth insert 0 [show cell linewidth]

}

####################################################################
# 
# PROC
proc lulDefineCellColour { w } {

   global gomCellColour

#   set StartColour [lulColourFloat2Hex $gomCellColour]

   set colour [tk_chooseColor -title "Choose Cell Colour" \
              -parent $w  \
			  -initialcolor $gomCellColour]

   if { $colour != "" } {
       $w.frame4.colour configure -bg $colour
       set RGBcolour [lulColourHex2Float $colour]
       eval "define cell colour {$RGBcolour}"
# save colour for later use
	   set gomCellColour $colour
   }

}

####################################################################
# 
# PROC
proc lulApplyDefineCell { w } {

     global gomSwitch
     global gomCellColour

     set a [$w.frame1.a get]
     set b [$w.frame1.b get]
     set c [$w.frame1.c get]

     set Command "define cell dimensions $a $b $c"

	 eval $Command

     set alpha     [$w.frame2.alpha get]
     set beta      [$w.frame2.beta  get]
     set gamma     [$w.frame2.gamma get]
     set linewidth [$w.frame4.linewidth get]

     set Command "define cell angles $alpha $beta $gamma"

	 eval $Command

     set xtrans     [$w.frame5.xtrans get]
     set ytrans     [$w.frame5.ytrans get]
     set ztrans     [$w.frame5.ztrans get]

     set Command "define cell trans $xtrans $ytrans $ztrans"

	 eval $Command

     if {$gomSwitch} {
	    eval "plot cell on"
     } else {
	    eval "plot cell off"
     }

     set colour [lulColourHex2Float $gomCellColour]

     eval "define cell colour {$colour}"

     eval "define cell linewidth $linewidth"

     lulInvalidateDisplay
}

####################################################################
# 
# select LDP atoms ...
#
proc lulSelectLDPatoms { } {

     global gomHelpFile
     global gomStructure
     global gomSwitch
     global gomControlFont

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomselectldpatoms
catch {destroy $w}
toplevel $w 
wm title $w "Select LDP atoms"
wm iconname $w "Select"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyLDPselection $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(plotldp)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1 -padx 3

set  NumStruct [show molstructures]

#frame $w.frame0 -borderwidth 2 -relief ridge
#pack $w.frame0 -side top
#
#for { set i 0} { $i < $NumStruct} {incr i} {
#
#radiobutton $w.frame0.structure$i -text "Structure nr: [expr $i + 1]" \
#            -value $i -variable gomStructure
#
#pack $w.frame0.structure$i -side left
#}

frame $w.frame1              -borderwidth 2 -relief ridge
label $w.frame1.segmentlabel -text "Segment:" -width 10
entry $w.frame1.segment      -width 15
pack  $w.frame1 -side top    -anchor w
pack  $w.frame1.segmentlabel $w.frame1.segment -side left
# default
$w.frame1.segment insert 0 "*"

frame $w.frame2              -borderwidth 2 -relief ridge
label $w.frame2.residuelabel -text "Residue:" -width 10
entry $w.frame2.residue      -width 15
pack  $w.frame2 -side top    -anchor w
pack  $w.frame2.residuelabel $w.frame2.residue -side left
# default
$w.frame2.residue insert 0 "*"

frame $w.frame3              -borderwidth 2 -relief ridge
label $w.frame3.atomlabel    -text "Atom:" -width 10
entry $w.frame3.atom         -width 15
pack  $w.frame3 -side top    -anchor w
pack  $w.frame3.atomlabel    $w.frame3.atom    -side left
# default
$w.frame3.atom insert 0 "CA"

set Segment [$w.frame1.segment get]
set Residue [$w.frame2.residue get]
set Atom    [$w.frame3.atom    get]

frame       $w.frame4 -borderwidth 2 -relief ridge
label       $w.frame4.label -text "LDP plot:"
radiobutton $w.frame4.on    -text "On"  -variable gomSwitch -value 1
radiobutton $w.frame4.off   -text "Off"   \
             -variable gomSwitch -value 0 \
             -command {plot -arrow}

pack        $w.frame4 -side top -anchor w -pady 4
pack        $w.frame4.label $w.frame4.on $w.frame4.off -side top

  if {[show ldpstatus] == "on"} {
       $w.frame4.on  select
  } else {
       $w.frame4.off select
  }
}

####################################################################
# PROC
proc lulApplyLDPselection { w } {

     global gomStructure
     global gomSwitch

     set Segment [string trim [$w.frame1.segment get]]
     set Residue [string trim [$w.frame2.residue get]]
     set Atom    [string trim [$w.frame3.atom    get]]

     eval "plot ldp atoms \{$Segment\} \{$Residue\} \{$Atom\}"
#     puts "plot ldp atoms \{$Segment\} \{$Residue\} \{$Atom\}"

     if { $gomSwitch } {
	      eval "plot ldp on"
     } else {
	      eval "plot ldp off"
     }

     lulInvalidateDisplay
}

##################### export probesurf file ############################
# PROC
proc lulExportProbesurfFile {} {

     global gomHelpFile
     global gomControlFont
     global gomProbeStructure


# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomexportprobesurf
catch {destroy $w}
toplevel $w 
wm title $w "Export probesurf inp"
wm iconname $w "Export probesurf inp"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulSaveProbesurfFile $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(writeinput)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

set Size     [show sizeofsystem]
set List ""
for {set i 1} { $i <= [show molstruct]} {incr i} {
 lappend List "$i"
}

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "File name: "
entry  $w.frame0.filename  -width 30
button $w.frame0.browse    -text "Browse..." -command "lulDoExportProbesurfFile $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse  -pady 3 -side left

frame  $w.frame8
pack   $w.frame8 -side top -anchor w -pady 3
label  $w.frame8.label -text "Structure # "
eval tk_optionMenu $w.frame8.struct gomProbeStructure $List
pack   $w.frame8.label $w.frame8.struct -side left -anchor w -pady 3

frame  $w.frame1
pack   $w.frame1 -side top -anchor w -pady 3

label  $w.frame1.label1 -text "X min: " -width 8 
entry  $w.frame1.xmin   -width 15
$w.frame1.xmin insert 0 "-$Size"

label  $w.frame1.label2 -text "X max: " -width 8
entry  $w.frame1.xmax   -width 15
$w.frame1.xmax insert 0 "$Size"

pack   $w.frame1.label1 $w.frame1.xmin  -side left
pack   $w.frame1.label2 $w.frame1.xmax  -side left

frame  $w.frame2
pack   $w.frame2 -side top -anchor w

label  $w.frame2.label1 -text "Y min: " -width 8
entry  $w.frame2.ymin   -width 15
$w.frame2.ymin insert 0 "-$Size"

label  $w.frame2.label2 -text "Y max: " -width 8
entry  $w.frame2.ymax   -width 15
$w.frame2.ymax insert 0 "$Size"

pack   $w.frame2.label1 $w.frame2.ymin  -side left
pack   $w.frame2.label2 $w.frame2.ymax  -side left

frame  $w.frame3
pack   $w.frame3 -side top -anchor w -pady 3

label  $w.frame3.label1 -text "Z min: " -width 8
entry  $w.frame3.zmin   -width 15
$w.frame3.zmin insert 0 "-$Size"

label  $w.frame3.label2 -text "Z max: " -width 8
entry  $w.frame3.zmax   -width 15
$w.frame3.zmax insert 0 "$Size"

pack   $w.frame3.label1 $w.frame3.zmin  -side left
pack   $w.frame3.label2 $w.frame3.zmax  -side left

frame  $w.frame4
pack   $w.frame4 -side top -anchor w
label  $w.frame4.label1 -text "X points" -width 9
entry  $w.frame4.xpts   -width 10
$w.frame4.xpts insert 0 "60"
pack   $w.frame4.label1 $w.frame4.xpts -side left

frame  $w.frame5
pack   $w.frame5 -side top -anchor w
label  $w.frame5.label2 -text "Y points" -width 9
entry  $w.frame5.ypts   -width 10
$w.frame5.ypts insert 0 "60"
pack   $w.frame5.label2 $w.frame5.ypts -side left


frame  $w.frame6
pack   $w.frame6 -side top -anchor w
label  $w.frame6.label3 -text "Z points" -width 9
entry  $w.frame6.zpts   -width 10
$w.frame6.zpts insert 0 "60"
pack   $w.frame6.label3 $w.frame6.zpts -side left

frame  $w.frame7
label  $w.frame7.label     -text "Probe radius: " -width 14
entry  $w.frame7.proberad  -width 10
$w.frame7.proberad insert 0 "2.0"
pack   $w.frame7 -side top -anchor w -pady 4
pack   $w.frame7.label $w.frame7.proberad -side left
 
}

##################### export probesurf file ############################
#
# PROC
proc lulDoExportProbesurfFile { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Probesurf files"		{.inp}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .inp -title "Save Probesurf File"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 $file
    }
}

#########################################################################
# if there is a name already available save with that name
# if not ask for a name
#
# PROC
proc lulSaveProbesurfFile { w } {

     global gomProbeStructure

set FileName [$w.frame0.filename get]
set Size     [show sizeofsystem]

    if {[string trim $FileName] != "" } {

              set Xmin     [$w.frame1.xmin get]
              if {[string trim $Xmin] == ""} {
			      set Xmin "-$Size"
				  }
              set Xmax     [$w.frame1.xmax get]
              if {[string trim $Xmax] == ""} {
			      set Xmax $Size
				  }

              set Ymin     [$w.frame2.ymin get]
              if {[string trim $Ymin] == ""} {
			      set Ymin "-$Size"
				  }

              set Ymax     [$w.frame2.ymax get]
              if {[string trim $Ymax] == ""} {
			      set Ymax $Size
				  }

              set Zmin     [$w.frame3.zmin get]
              if {[string trim $Zmin] == ""} {
			      set Zmin "-$Size"
				  }
              set Zmax     [$w.frame3.zmax get]
              if {[string trim $Zmax] == ""} {
			      set Zmax $Size
				  }

              set Xpts     [$w.frame4.xpts get]
              if {[string trim $Xpts] == ""} {
			      set Xpts 60
				  }
              set Ypts     [$w.frame5.ypts get]
              if {[string trim $Ypts] == ""} {
			      set Ypts 60
				  }
              set Zpts     [$w.frame6.zpts get]
              if {[string trim $Zpts] == ""} {
			      set Zpts 60
				  }

              set ProbeRad [$w.frame7.proberad get]
              if {[string trim $ProbeRad] == ""} {
			      set ProbeRad 2.0
				  }

    eval "export input probesurf $gomProbeStructure \"$FileName\" $Xmin $Xmax $Ymin $Ymax $Zmin $Zmax $Xpts $Ypts $Zpts $ProbeRad"

# change to current directory
    lulChangeDirectory "$FileName"
    }
}

########################## end of export model #########################

##################### run gcube2plt file ############################
# PROC
proc lulRungCube2Plt {} {

     global gomHelpFile
     global gomControlFont
     global gomCubeProg
     global gomGaussianCubeDataCB
     global gomGaussianCubeProg


# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrungcube2plt
catch {destroy $w}
toplevel $w 
wm title $w "Run gcube2plt/g94cub2pl"
wm iconname $w "Run gcube2plt/g94cub2pl"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRungCube2pltCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(rungcube2plt)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Convert a Gaussian94/98 cube file to a plt file:" -font $gomControlFont
pack   $w.label  -side top -anchor w

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRungCube2pltInput $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRungCube2pltOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left

frame  $w.frame4
pack   $w.frame4 -side top -anchor w
label  $w.frame4.label -text "Used: "
radiobutton $w.frame4.gaussian95  -text "Gaussian94"  -variable gomGaussianCubeProg  -value 0 \
            -command ""
radiobutton $w.frame4.gaussian98  -text "Gaussian98"  -variable gomGaussianCubeProg  -value 1 \
            -command "" 
pack $w.frame4.label       -side left 
pack $w.frame4.gaussian95  -side left -padx 3
pack $w.frame4.gaussian98  -side left

$w.frame4.gaussian98 select

frame  $w.frame5
pack   $w.frame5 -side top -anchor w
label       $w.frame5.label -text "Calculate from Cube file: "
checkbutton $w.frame5.orbdens     -text "Orbital/density" \
             -variable gomGaussianCubeDataCB(1) -command "set gomCubeProg 1;gomGaussianCubeAction $w"
checkbutton $w.frame5.gradient     -text "Gradient" \
             -variable gomGaussianCubeDataCB(2) -command "set gomCubeProg 1;gomGaussianCubeAction $w"
checkbutton $w.frame5.gradnorm     -text "Gradient norm" \
             -variable gomGaussianCubeDataCB(3) -command "set gomCubeProg 1;gomGaussianCubeAction $w"
checkbutton $w.frame5.laplacian     -text "Laplacian" \
             -variable gomGaussianCubeDataCB(4) -command "set gomCubeProg 1;gomGaussianCubeAction $w"
pack   $w.frame5.label      -side left -anchor w
pack   $w.frame5.orbdens    -side left -anchor w
pack   $w.frame5.gradient   -side left -anchor w
pack   $w.frame5.gradnorm   -side left -anchor w
pack   $w.frame5.laplacian  -side left -anchor w

$w.frame5.orbdens select
$w.frame5.orbdens configure -state disable

frame  $w.frame2
pack   $w.frame2 -side top -anchor w 

label  $w.frame2.orbtext   -text "Orbital nr: "
entry  $w.frame2.orbnumber -width 10
pack   $w.frame2.orbtext   -side left
pack   $w.frame2.orbnumber -side left

frame  $w.frame6
pack   $w.frame6 -side top -anchor w 
label  $w.frame6.program   -text "Run program: "
pack   $w.frame6.program   -side left

radiobutton $w.frame6.gcube2plt  -text "gcube2plt"  -variable gomCubeProg  -value 1 \
            -command "$w.frame2.orbnumber configure -state normal;\
                      $w.frame5.orbdens select; $w.frame5.orbdens configure -state disable;\
                      $w.frame5.gradient configure  -state normal; \
                      $w.frame5.gradnorm configure  -state normal; \
                      $w.frame5.laplacian configure -state normal"
radiobutton $w.frame6.g94cub2pl  -text "g94cub2pl"  -variable gomCubeProg  -value 0 \
            -command "$w.frame2.orbnumber configure -state disabled;\
                      $w.frame5.orbdens select; $w.frame5.orbdens configure -state disable;\
                      $w.frame5.gradient  deselect; $w.frame5.gradient configure  -state disable; \
                      $w.frame5.gradnorm  deselect; $w.frame5.gradnorm configure  -state disable; \
                      $w.frame5.laplacian deselect; $w.frame5.laplacian configure -state disable"

 
pack $w.frame6.gcube2plt  -side left -padx 3
pack $w.frame6.g94cub2pl  -side left
set gomCubeProg     1
$w.frame6.gcube2plt select

if {$gomCubeProg} {
$w.frame2.orbnumber configure -state normal
} else {
$w.frame2.orbnumber configure -state disabled
}

frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}


#########################################################################
# PROC
proc gomGaussianCubeAction { w } {

     global gomGaussianCubeDataCB

if {$gomGaussianCubeDataCB(1) || $gomGaussianCubeDataCB(2) || $gomGaussianCubeDataCB(3)
    || $gomGaussianCubeDataCB(4)} {
    $w.frame2.orbnumber configure -state disabled
    return;
} else {
    $w.frame2.orbnumber configure -state normal
    return;
}


}
##################### lulDoRungcube2pltInput ############################
#
# PROC
proc lulDoRungCube2pltInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Cube data"		{.data .cub .cube .DATA .CUB .CUBE}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .data -title "Input Cube File"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRungcube2pltOutput ############################
#
# PROC
proc lulDoRungCube2pltOutput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.plt}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Save Plot File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRungcube2pltInput ############################
#
# PROC
proc lulRungCube2pltCommand { w } {

     global env
     global gomCubeProg
	 global gomEnv
     global gomGaussianCubeDataCB
     global gomGaussianCubeProg

   set InputFile  [string trim [$w.frame0.filename get]]
   if {$InputFile == ""} {
        gomError "Input file name missing!"
        return
   }
   set OutputFile [string trim [$w.frame1.filename get]]
   if {$OutputFile == ""} {
       set OutputFile [file rootname $InputFile]
       set OutputFile "[file join [file dirname $InputFile] $OutputFile.plt]"
       $w.frame1.filename delete 0 end
       $w.frame1.filename insert 0 "$OutputFile"
       update idletasks
   }
         
# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

set OrbNumber  [string trim [$w.frame2.orbnumber get]]

if {$gomCubeProg} {
    if {$OrbNumber != ""} {
      if {$gomGaussianCubeProg} {
        set Program  "gcube2plt -tg98 -i$InputFile -o$OutputFile -m$OrbNumber"
      } else {
        set Program  "gcube2plt -tg94 -i$InputFile -o$OutputFile -m$OrbNumber"
      }
    } else {
      set Types ""
      if {$gomGaussianCubeDataCB(2)} {
          set Types "g"
      }
      if {$gomGaussianCubeDataCB(3)} {
          set Types [append Types n]
      }
      if {$gomGaussianCubeDataCB(4)} {
          set Types [append Types l]
      }

      if {$gomGaussianCubeProg} {
        set Program "gcube2plt -tg98  -d$Types  -i$InputFile -o$OutputFile"
      } else {
        set Program "gcube2plt -tg94  -d$Types  -i$InputFile -o$OutputFile"
      }
    }
} else {
    set Program    "g94cub2pl $InputFile $OutputFile"
    set OrbNumber  ""
}

set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

  catch {close $f}
}

##################### run probesurf file ############################
# PROC
proc lulRunProbesurf {} {

     global gomHelpFile
     global gomControlFont
     global gomRunProbesurfRunState
     global gomRunProbesurfRMethod

# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrunprobesurf
catch {destroy $w}
toplevel $w 
wm title $w "Run ProbeSurf"
wm iconname $w "Run probesurf"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRunProbesurfCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(runprobesurf)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRunProbesurfInput $w"
radiobutton $w.frame0.foreground -text "Calc in foreground: " -variable gomRunProbesurfRunState -value fg
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left -padx 5
pack   $w.frame0.foreground -padx 4

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRunProbesurfOutput $w"
radiobutton $w.frame1.background -text "Calc in background: " -variable gomRunProbesurfRunState -value bg
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left -padx 5
pack   $w.frame1.background -padx 4

$w.frame0.foreground select

frame  $w.frame4
pack   $w.frame4 -side top -anchor w
label       $w.frame4.text -text "Distance criteria: " -width 18
pack        $w.frame4.text -side left -anchor w
radiobutton $w.frame4.r  -text "Direct (r)"     -variable gomRunProbesurfRMethod -value 1
radiobutton $w.frame4.r2 -text "Squared (r**2)" -variable gomRunProbesurfRMethod -value 2
pack        $w.frame4.r  -padx 4 -side left -anchor w
pack        $w.frame4.r2 -padx 4 -side left -anchor w

$w.frame4.r select

frame  $w.frame3
pack   $w.frame3 -side top -anchor center

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRunProbesurfInput ############################
#
# PROC
proc lulDoRunProbesurfInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Probesurf input"		{.inp}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .inp -title "Input Probesurf File"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRunProbesurfOutput ############################
#
# PROC
proc lulDoRunProbesurfOutput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.plt}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Save Plot File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRunProbesurfInput ############################
#
# PROC
proc lulRunProbesurfCommand { w } {

     global env
     global lulRunProbesurfLineCount
     global gomRunProbesurfRunState
     global gomRunProbesurfRMethod
     global gomRunProbesurfFileName
     global gomEnv

set InputFile  [string trim [$w.frame0.filename get]]
set OutputFile [string trim [$w.frame1.filename get]]
# check to see that output is not a directory
  if {[file isdirectory "$OutputFile"]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

if {$gomRunProbesurfRMethod == 1} {
    set Method 1
} elseif {$gomRunProbesurfRMethod == 1} {
    set Method 2
} else {
    gomError "wrong method value (has to be 1 for 'r' or 2 for 'r**2'"
    return 
}

if {[string trim $OutputFile] == ""} {
    set OutputFile "[file join [file dirname $InputFile] probesurf.plt]"
    $w.frame1.filename delete 0 end
    $w.frame1.filename insert 0 "$OutputFile"
    update idletasks
}
    set gomRunProbesurfFileName  $OutputFile

    set Program    "probsurf"

set Command [file join $gomEnv(GOM_BIN) $Program]

if {$gomRunProbesurfRunState == "fg"} {
  $w.frame3.text insert 1.0 "Please stand by running job in foreground...\n"
  update idletasks
  set i 1
  set f [open "|$Command -m$Method -o$OutputFile \< $InputFile " r]
  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }
  close $f

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulMessageDialog "Your foreground probesurf job has finished!"
	      return
	 }

  lulProbesurfEOFMessageDialog "Your foreground probesurf job has finished!\nDo you want to read in the contour file?"
} elseif {$gomRunProbesurfRunState == "bg"} {
  set lulRunProbesurfLineCount 1
  $w.frame3.text insert 1.0 "Please stand by running job in background...\n"
  update idletasks
  set f [open "|$Command -m$Method -o$OutputFile \< $InputFile " r]
  fileevent $f readable [list lulRunProbesurfCommandReader $f]
  fconfigure $f -buffering none -blocking 0
} else {
  gomError "unknow run state (only 'fg' (foreground) and 'bg' (background) allowed)"
  return
}

}
###############################################################################
#
# special probesurf dialog widget
#
proc lulProbesurfEOFMessageDialog {InputText} {

     global gomRunProbesurfFileName

after idle {.gomdialogps.msg configure -wraplength 4i}
set i [tk_dialog .gomdialogps "gOpenMol Question" $InputText \
question 0 "Yes" "No" ]

switch $i {
    0 {contour file $gomRunProbesurfFileName;lulContourControl}
    1 {gomPrint "Countour file $gomRunProbesurfFileName is now available"}
}

}

##################### lulRunProbesurfCommandReader ############################
#
# PROC
proc lulRunProbesurfCommandReader { pipe } {

   global lulRunProbesurfLineCount

   if [eof $pipe] {
      catch "close $pipe"
      lulMessageDialog "Your background probesurf job has finished!"

# return if no molecular systems defined
      if {[show molstructures] < 1} {
           return
      }

      lulProbesurfEOFMessageDialog "Do you want to read in the contour file?" 
      return
   }

    gets $pipe Text
    if {[winfo exists  .gomrunprobesurf]} {
      .gomrunprobesurf.frame3.text insert $lulRunProbesurfLineCount.0 "$Text\n"
      .gomrunprobesurf.frame3.text see end
       update idletasks
    } else {
       close $pipe
       return
    }
    puts $Text
    incr lulRunProbesurfLineCount
}
##################### run pltfile ############################
# PROC
proc lulRunPltfile {} {

     global gomHelpFile
     global gomSwitch
     global gomControlFont


# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrunpltfile
catch {destroy $w}
toplevel $w 
wm title $w "Run pltfile"
wm iconname $w "Run pltfile"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRunPltfileCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(runpltfile)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRunPltfileInput $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRunPltfileOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left

frame  $w.frame2 -borderwidth 2 -relief ridge
pack   $w.frame2 -side top -anchor w

radiobutton  $w.frame2.uf    -text "Unformatted ==> Formatted"    \
             -variable gomSwitch                                  \
			 -value 1

radiobutton  $w.frame2.fu    -text "Formatted   ==> Unformatted"  \
             -variable gomSwitch                                  \
			 -value 0

pack         $w.frame2.uf $w.frame2.fu -side top -anchor w

$w.frame2.uf select

frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRunPltfileInput ############################
#
# PROC
proc lulDoRunPltfileInput { w } {

     global gomSwitch

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    if {$gomSwitch} {
      set types {
	  {"Pltfile input"		{.plt}		TEXT}
	  {"All files"		*}
      }
	} else {
      set types {
	  {"Pltfile input"		{.txt}		TEXT}
	  {"All files"		*}
      }
	}
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension * -title "Input Pltfile File"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRunPltfileOutput ############################
#
# PROC
proc lulDoRunPltfileOutput { w } {

     global gomSwitch

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    if {$gomSwitch} {
      set types {
        {"Plotfile output"		{.txt}		TEXT}
	    {"All files"		*}
        }
    } else {
      set types {
        {"Plotfile output"		{.plt}		TEXT}
	    {"All files"		*}
        }
	}
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension * -title "Save Pltfile File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRunPltfile command ############################
#
# PROC
proc lulRunPltfileCommand { w } {

     global env
	 global gomSwitch
	 global gomEnv

set InputFile  [string trim [$w.frame0.filename get]]
if {$InputFile == ""} {
   $w.frame3.text insert 0.0 "Input file name missing"
   return
   }

set OutputFile [string trim [$w.frame1.filename get]]
if {$OutputFile == ""} {
   $w.frame3.text insert 0.0 "Output file name missing"
   return
   }
# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

if {$gomSwitch != 0} {
  set Program    "pltfile -uf -i$InputFile -o$OutputFile"
} else {
  set Program    "pltfile -fu -i$InputFile -o$OutputFile"
}

set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }
close $f
}


##################### run kont2plt file ############################
# PROC
proc lulRungKont2Plt {} {

     global gomHelpFile
     global gomControlFont
     global gomOutPltFileType

# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrungkont2plt
catch {destroy $w}
toplevel $w 
wm title $w "Run kont2plt"
wm iconname $w "kont2plt"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRungKont2pltCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(runkont2plt)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Convert a Grid data file to a (formatted/unformatted) plt file:" -font $gomControlFont
pack   $w.label  -side top -anchor w

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRungKont2pltInput $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRungKont2pltOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left

frame  $w.frame2
pack   $w.frame2 -side top -anchor w

label       $w.frame2.label -text "Output plt file type: "
pack        $w.frame2.label -side left -anchor w
radiobutton $w.frame2.unformatted  -text "Unformatted"  -variable gomOutPltFileType  -value 2
radiobutton $w.frame2.formatted    -text "Formatted"    -variable gomOutPltFileType  -value 1 
pack $w.frame2.unformatted  -side left -padx 3
pack $w.frame2.formatted    -side left
$w.frame2.unformatted select


frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRungKont2pltInput ############################
#
# PROC
proc lulDoRungKont2pltInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Grid data"		{.grd .GRD}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .data -title "Input Grid File"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRungKont2pltOutput ############################
#
# PROC
proc lulDoRungKont2pltOutput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.plt}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Save Plot File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRungKont2pltInput ############################
#
# PROC
proc lulRungKont2pltCommand { w } {

     global env
     global gomOutPltFileType
	 global gomEnv

set InputFile  [string trim [$w.frame0.filename get]]
set OutputFile [string trim [$w.frame1.filename get]]

if {$InputFile == "" || $OutputFile == ""} {

   gomError "either input: '$InputFile' or output: '$OutputFile' file name is missing"
   return
}

# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

  set TempFileName "kont2plt[clock seconds].tmp"

# write it out in temp file ...
  set f [open $TempFileName w]
  if {$f == ""} {
    gomError "can't open '$TempFileName' file"
    return
  }
  puts $f "$InputFile"
# output file type unformatted/formatted
  if {$gomOutPltFileType == 2} {
    puts $f $gomOutPltFileType
  } elseif {$gomOutPltFileType == 1} {
    puts $f $gomOutPltFileType
  } else {
    gomError "output file type is unknown: '$gomOutPltFileType' (should be 2 = unform , 1 = form)"
    close $f
    return
  }

  puts $f "$OutputFile"
  close $f
# done

  set Program    "kont2plt \< $TempFileName"

set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

  catch {close $f}
# delete temp file
  file delete $TempFileName

}

##################### run autodock2plt file ############################
# PROC
proc lulRungAutoDock2Plt {} {

     global gomHelpFile
     global gomControlFont


# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrungautodock2plt
catch {destroy $w}
toplevel $w 
wm title $w "Run AutoDock2plt"
wm iconname $w "AutoDock2plt"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRunAutoDock2pltCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(runautodock2plt)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Convert an AutoDock map file to a plt file:" -font $gomControlFont
pack   $w.label  -side top -anchor w

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRungAutoDock2pltInput $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRungAutoDock2pltOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left

frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRungAutoDock2pltInput ############################
#
# PROC
proc lulDoRungAutoDock2pltInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"AutoDock map file"		{.map .MAP}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .map -title "Input AutoDock map File"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRungAutoDock2pltOutput ############################
#
# PROC
proc lulDoRungAutoDock2pltOutput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.plt}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Save Plot File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRungAutoDock2pltInput ############################
#
# PROC
proc lulRunAutoDock2pltCommand { w } {

     global env
	 global gomEnv

set InputFile  [string trim [$w.frame0.filename get]]
set OutputFile [string trim [$w.frame1.filename get]]

if {$InputFile == "" || $OutputFile == ""} {

   gomError "either input: '$InputFile' or output: '$OutputFile' file name is missing"
   return
}

# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

  set TempFileName "autodock2plt[clock seconds].tmp"

# write it out in temp file ...
  set f [open $TempFileName w]
  if {$f == ""} {
    gomError "can't open '$TempFileName' file"
    return
  }
  puts $f "$InputFile"
  puts $f "$OutputFile"
  close $f
# done

  set Program    "autodock2plt \< $TempFileName"

set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

  catch {close $f}
# delete temp file
  file delete $TempFileName

}

##################### run contman file ############################
# PROC
proc lulRunContMan {} {

     global gomHelpFile
     global gomControlFont
     global gomOutFileAction

# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomruncontman
catch {destroy $w}
toplevel $w 
wm title $w "Run ContMan"
wm iconname $w "contman"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRunContManCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(runcontman)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Contour file manipulation:" -font $gomControlFont
pack   $w.label  -side top -anchor w

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file #1 name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRunContManInput 1 $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left

frame  $w.frame4
pack   $w.frame4 -side top -anchor w

label  $w.frame4.label     -text "Input file #2 name:" -width 18
entry  $w.frame4.filename  -width 40
button $w.frame4.browse    -text "Browse..." -command "lulDoRunContManInput 2 $w"
pack   $w.frame4 -side top -anchor w
pack   $w.frame4.label $w.frame4.filename $w.frame4.browse -side left

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRunContManOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left

frame  $w.frame2
pack   $w.frame2 -side top -anchor w

label       $w.frame2.label -text "Action type (-/+): "
pack        $w.frame2.label -side left -anchor w
radiobutton $w.frame2.minus  -text "Minus (-)" -variable gomOutFileAction  -value 0
radiobutton $w.frame2.plus   -text "Plus (+)"  -variable gomOutFileAction  -value 1 
pack $w.frame2.minus  -side left -padx 3
pack $w.frame2.plus   -side left
$w.frame2.minus select


frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRunContManInput ############################
#
# PROC
proc lulDoRunContManInput { Input w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.plt .PLT}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Input #$Input Plot File"]]

    if [string compare $file ""] {

      if {$Input == 1} {
        $w.frame0.filename delete 0 end
	    $w.frame0.filename insert 0 $file
	    $w.frame0.filename xview end
      } elseif {$Input == 2} {
        $w.frame4.filename delete 0 end
	    $w.frame4.filename insert 0 "$file"
	    $w.frame4.filename xview end
      } else {
        gomError "wrong input file index '$Input'! Has to be 1 or 2"
        return
      }

# change to current directory
      lulChangeDirectory "$file"

    }
}

##################### lulDoRunContManOutput ############################
#
# PROC
proc lulDoRunContManOutput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.plt}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Save Plot File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"
	$w.frame1.filename xview end

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRunContManInput ############################
#
# PROC
proc lulRunContManCommand { w } {

     global env
     global gomOutFileAction
	 global gomEnv

set InputFile1  [string trim [$w.frame0.filename get]]
set InputFile2  [string trim [$w.frame4.filename get]]
set OutputFile [string trim [$w.frame1.filename get]]

if {$InputFile1 == "" || $InputFile2 == "" || $OutputFile == ""} {

   gomError "either input1: '$InputFile1', input2: '$InputFile2'  or output: '$OutputFile' file name is missing"
   return
}

# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

  if {$gomOutFileAction == 0} {
     set Program    "contman -i$InputFile1-$InputFile2 -o$OutputFile"
  } elseif {$gomOutFileAction == 1} {
     set Program    "contman -i$InputFile1+$InputFile2 -o$OutputFile"
  } else {
     gomError "wrong action type: '$gomOutFileAction'. Allowed values are 0 or 1"
     return
  }

set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

  catch {close $f}

}

##################### run uhbd2plt file ############################
# PROC
proc lulRunUHBD2Plt {} {

     global gomHelpFile
     global gomControlFont
     global gomOutUHBDFileAction

# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrunuhbd2plt
catch {destroy $w}
toplevel $w 
wm title $w "Run UHBD2Plt"
wm iconname $w "UHBD2Plt"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRunUHBD2PltCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(runuhbd2plt)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Run UHBD2Plt program:" -font $gomControlFont
pack   $w.label  -side top -anchor w

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file  name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRunUHBD2PltInput $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRunUHBD2PltOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left

frame  $w.frame2
pack   $w.frame2 -side top -anchor w

label       $w.frame2.label -text "Input file is: "
pack        $w.frame2.label -side left -anchor w
radiobutton $w.frame2.formatted   -text "Formatted"    -variable gomOutUHBDFileAction  -value 0
radiobutton $w.frame2.unformatted -text "Unformatted"  -variable gomOutUHBDFileAction  -value 1 
pack $w.frame2.formatted   -side left -padx 3
pack $w.frame2.unformatted -side left
$w.frame2.formatted select

frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRunUHBD2PltInput ############################
#
# PROC
proc lulDoRunUHBD2PltInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"UHBD grid file"		{.grd .GRD}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .grd -title "Input GRID File"]]

    if [string compare $file ""] {

        $w.frame0.filename delete 0 end
	    $w.frame0.filename insert 0 "$file"
	    $w.frame0.filename xview end

# change to current directory
      lulChangeDirectory "$file"

    }
}

##################### lulDoRunUHBD2PltOutput ############################
#
# PROC
proc lulDoRunUHBD2PltOutput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.plt}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Save Plot File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 $file
	$w.frame1.filename xview end

# change to current directory
    lulChangeDirectory $file

    }
}

##################### lulDoRunUHBD2PltInput ############################
#
# PROC
proc lulRunUHBD2PltCommand { w } {

     global env
     global gomOutUHBDFileAction
	 global gomEnv

set InputFile  [string trim [$w.frame0.filename get]]
set OutputFile [string trim [$w.frame1.filename get]]

if {$InputFile == "" || $OutputFile == ""} {

   gomError "either input: '$InputFile' or output: '$OutputFile' file name is missing"
   return
}

# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

  if {$gomOutUHBDFileAction == 0} {
     set Program    "gridasc2plt"
  } elseif {$gomOutUHBDFileAction == 1} {
     set Program    "gridbin2plt"
  } else {
     gomError "wrong action type: '$gomOutFileAction'. Allowed are 0 or 1"
     return
  }

  set TempFileName "$Program[clock seconds].tmp"

# write it out in temp file ...
  set f [open $TempFileName w]
  if {$f == ""} {
    gomError "can't open '$TempFileName' file"
    return
  }
  puts $f "$InputFile"
  puts $f "$OutputFile"
  close $f
# done

set Program "$Program \< $TempFileName"
set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

  catch {close $f}
# delete temp file
  file delete $TempFileName

}

##################### run turbomole2plt file ############################
# PROC
proc lulRunTurboMole2Plt {} {

     global gomHelpFile
     global gomControlFont


# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrungturbomole2plt
catch {destroy $w}
toplevel $w 
wm title $w "Run TurboMole2plt"
wm iconname $w "TurboMole2plt"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRungTurboMole2pltCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(runturbomole2plt)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Convert an TurboMole grid file to a plt file:" -font $gomControlFont
pack   $w.label  -side top -anchor w

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRungTurboMole2pltInput $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRungTurboMole2pltOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left

frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRungTurboMole2pltInput ############################
#
# PROC
proc lulDoRungTurboMole2pltInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"TurboMole grid file"		{.dat *}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .map -title "Input TurboMole grid File"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 $file

# change to current directory
    lulChangeDirectory $file

    }
}

##################### lulDoRungTurboMole2pltOutput ############################
#
# PROC
proc lulDoRungTurboMole2pltOutput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.plt}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Save Plot File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRungTurboMole2pltInput ############################
#
# PROC
proc lulRungTurboMole2pltCommand { w } {

     global env
	 global gomEnv

set InputFile  [string trim [$w.frame0.filename get]]
set OutputFile [string trim [$w.frame1.filename get]]

if {$InputFile == "" || $OutputFile == ""} {

   gomError "either input: '$InputFile' or output: '$OutputFile' file name is missing"
   return
}

# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

  set Program    "tmole2plt $InputFile $OutputFile"

set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

  catch {close $f}

}

##################### run gamess2plt file ############################
# PROC
proc lulRungGamess2Plt {} {

     global gomHelpFile
     global gomControlFont


# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrungamess2plt
catch {destroy $w}
toplevel $w 
wm title $w "Run Gamess2plt"
wm iconname $w "gamess2plt"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRunGamess2pltCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(rungamess2plt)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Coonvert a Gamess 'cube' PUNCH file to a plt file:" -font $gomControlFont
pack   $w.label  -side top -anchor w

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRungGamess2pltInput $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRungGamess2pltOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left

frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRungGamess2pltInput ############################
#
# PROC
proc lulDoRungGamess2pltInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Gamess PUNCH file"		{*}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension * -title "Input Gamess PUNCH File"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRungGamess2pltOutput ############################
#
# PROC
proc lulDoRungGamess2pltOutput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.plt}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Save Plot File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRungGamess2pltInput ############################
#
# PROC
proc lulRunGamess2pltCommand { w } {

     global env
	 global gomEnv

set InputFile  [string trim [$w.frame0.filename get]]
set OutputFile [string trim [$w.frame1.filename get]]

if {$InputFile == ""} {

   gomError "input: '$InputFile' file name is missing"
   return
}

# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

  if {$OutputFile == ""} {
     set Program    "gamess2plt  -i$InputFile"
  } else {
     set Program    "gamess2plt  -i$InputFile -o$OutputFile"
  }

set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

  catch {close $f}

}

##################### run join gamess irc files ############################
# PROC
proc lulRunJoinGamessIRCfiles {} {

     global gomHelpFile
     global gomSwitch
     global gomControlFont


# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrunpltfile
catch {destroy $w}
toplevel $w 
wm title $w "Run Join Gamess IRC Files"
wm iconname $w "Run JIRC"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRunJoinGamessIRCfilesCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(joingamessirc)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name #1:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRunJoinGIRCInput1 $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left -padx 4

frame  $w.frame4
pack   $w.frame4 -side top -anchor w

label  $w.frame4.label     -text "Input file name #2:" -width 18
entry  $w.frame4.filename  -width 40
button $w.frame4.browse    -text "Browse..." -command "lulDoRunJoinGIRCInput2 $w"
pack   $w.frame4 -side top -anchor w
pack   $w.frame4.label $w.frame4.filename $w.frame4.browse -side left -padx 4

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRunRunJoinGIRCOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left -padx 4

}

##################### lulDoRunJoinGIRCInput1 ############################
#
# PROC
proc lulDoRunJoinGIRCInput1 { w } {


    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
      set types {
	  {"Pltfile input"		{.irc .IRC}		TEXT}
	  {"All files"		*}
      }

    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension * -title "Input IRC File #1"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}
##################### lulDoRunJoinGIRCInput2 ############################
#
# PROC
proc lulDoRunJoinGIRCInput2 { w } {

     global gomSwitch

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
      set types {
	  {"Pltfile input"		{.irc .IRC}		TEXT}
	  {"All files"		*}
      }

    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension * -title "Input IRC File #1"]]

    if [string compare $file ""] {
    $w.frame4.filename delete 0 end
	$w.frame4.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRunJoinGIRCOutput ############################
#
# PROC
proc lulDoRunRunJoinGIRCOutput { w } {


    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
      set types {
        {"Plotfile output"		{.irc .IRC}		TEXT}
	    {"All files"		*}
        }

    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension * -title "Save IRC File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRunPltfile command ############################
#
# PROC
proc lulRunJoinGamessIRCfilesCommand { w } {

     global env
	 global gomEnv

set InputFile1  [string trim [$w.frame0.filename get]]
if {$InputFile1 == ""} {
   gomError "Input file name #1 missing"
   return
   }

set InputFile2  [string trim [$w.frame4.filename get]]
if {$InputFile2 == ""} {
   gomError "Input file name #2 missing"
   return
   }

set OutputFile [string trim [$w.frame1.filename get]]
if {$OutputFile == ""} {
   gomError "Output file name missing"
   return
   }
# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

  lulUtility::join-gamess-irc $InputFile1 $InputFile2 $OutputFile

}

##################### run Jaguar2plt file ############################
# PROC
proc lulRungJaguar2Plt {} {

     global gomHelpFile
     global gomControlFont
     global gomCubeProg


# return if no molecular systems defined
#     if {[show molstructures] < 1} {
#          lulErrorDialog {ERROR: no structure available. Read a structure first!}
#	      return
#	 }

set w .gomrunjaguar2plt
catch {destroy $w}
toplevel $w 
wm title $w "Run Jaguar2plt"
wm iconname $w "Run Jaguar2Plt"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRunJaguar2pltCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(runjaguar2plt)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Convert a Jaguar plot file to a gOpenMol plt file:" -font $gomControlFont
pack   $w.label  -side top -anchor w

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRunJaguar2pltInput $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Output file name:" -width 18
entry  $w.frame1.filename  -width 40
button $w.frame1.browse    -text "Browse..." -command "lulDoRunJaguar2pltOutput $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.filename $w.frame1.browse -side left

frame  $w.frame2
pack   $w.frame2 -side top -anchor w

frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRungJaguar2pltInput ############################
#
# PROC
proc lulDoRunJaguar2pltInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Jaguar plot file"		{.plt *.PLT}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Input Jaguar File"]]

    if [string compare $file ""] {
    $w.frame0.filename delete 0 end
	$w.frame0.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRunJaguar2pltOutput ############################
#
# PROC
proc lulDoRunJaguar2pltOutput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"gOpenMol Plot file"		{.plt}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Save Plot File"]]

    if [string compare $file ""] {
    $w.frame1.filename delete 0 end
	$w.frame1.filename insert 0 "$file"

# change to current directory
    lulChangeDirectory "$file"

    }
}

##################### lulDoRunJaguar2pltInput ############################
#
# PROC
proc lulRunJaguar2pltCommand { w } {

     global env
     global gomCubeProg
	 global gomEnv

set InputFile  [string trim [$w.frame0.filename get]]
set OutputFile [string trim [$w.frame1.filename get]]
# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

if {[string trim $OutputFile] == ""} {
 set Program    "jaguar2plt -i$InputFile"
} else {
 set Program    "jaguar2plt -i$InputFile -o$OutputFile"
}

set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

  catch {close $f}
}

##################### monitor distance ############################
#
# monitor distance
#
# PROC
proc lulMonitorDistance {} {

     global gomHelpFile
     global gomSwitchDistMon
     global gomAppendDistList
     global gomControlFont

# return if no molecular systems defined
#     if {[show trajectory frames] < 1} {
#          lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
#	      return
#	 }

set w .gommonitordistance
catch {destroy $w}
toplevel $w 
wm title $w "Monitor distance"
wm iconname $w "Monitor distance"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss    -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply      -font "$gomControlFont" \
        -command "lulDoMonitorDistance $w"
button $w.buttons.help    -text Help       -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(monitordistance)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Monitor distance:"
pack   $w.label -side top -anchor w

frame  $w.atoms
pack   $w.atoms -side top -anchor w

# set #1,#2
foreach set {1 2} {
    frame  $w.atoms.set$set -borderwidth 2 -relief ridge
    pack   $w.atoms.set$set -side left
    label  $w.atoms.set$set.label -text "Atom(s) #$set:"
    pack   $w.atoms.set$set.label -side top -anchor w
    lulCreateAtomInputEntries $w.atoms.set$set 15 "*" "*" "*" "" 1
}

frame       $w.properties -borderwidth 2 -relief ridge
pack        $w.properties -side top -anchor w -fill x

button      $w.properties.color -text "Line colour..." -bg #ffffff \
             -command "lulMonitorDistanceLineColor $w"
pack        $w.properties.color -side left -anchor w -pady 3 -padx 3

menubutton  $w.properties.type  -text "Line type"   \
             -menu $w.properties.type.menu          \
             -borderwidth 2 -relief raised
menu        $w.properties.type.menu
set m       $w.properties.type.menu
$m add radiobutton -label "* * * *"    -variable gomDistMonLineType -value 1
$m add radiobutton -label "- - - -"    -variable gomDistMonLineType -value 2
$m add radiobutton -label "* - - *"    -variable gomDistMonLineType -value 3
$m add radiobutton -label "-------"    -variable gomDistMonLineType -value 4
pack        $w.properties.type  -side left -anchor w
$m invoke 1

frame       $w.options 
pack        $w.options -side top -anchor w -fill x

frame       $w.options.left
pack        $w.options.left   -side left -anchor w

frame       $w.options.right 
pack        $w.options.right  -side right -anchor e

frame       $w.options.left.frame8 -borderwidth 2 -relief ridge
button      $w.options.left.frame8.delete -text "Delete all" \
            -command "lulDeleteMonitorDistance $w;lulInvalidateDisplay"

pack        $w.options.left.frame8 -side bottom -anchor w
pack        $w.options.left.frame8.delete -side left

frame       $w.options.left.frame9 -borderwidth 2 -relief ridge
checkbutton $w.options.left.frame9.append -text "Append"    \
             -variable gomAppendDistList \
              -onvalue 1 -offvalue 0

pack        $w.options.left.frame9 -side bottom -anchor w
pack        $w.options.left.frame9.append -side left

$w.options.left.frame9.append select

frame       $w.options.right.frame4 -borderwidth 2 -relief ridge
label       $w.options.right.frame4.label -text "Display:"
radiobutton $w.options.right.frame4.on  -text "On"  \
            -variable gomSwitchDistMon -value 1 \
			-command  {eval monitor display distance on;lulInvalidateDisplay}  
radiobutton $w.options.right.frame4.off -text "Off" \
            -variable gomSwitchDistMon -value 0 \
			-command  {eval monitor display distance off;lulInvalidateDisplay}
pack $w.options.right.frame4 -side top -anchor e
pack $w.options.right.frame4.label \
     $w.options.right.frame4.on    \
	 $w.options.right.frame4.off \
     -side top -anchor w

# set the display state
set DispState [show monitor distance state]

if {$DispState == 1} {
    $w.options.right.frame4.on select
  } else {
    $w.options.right.frame4.off select
  }

# 
frame $w.distmon -borderwidth 2 -relief ridge 
pack  $w.distmon -side top 
text  $w.distmon.text \
      -yscrollcommand "$w.distmon.scrollbar set" \
	  -width 60 -height 10
scrollbar $w.distmon.scrollbar -command "$w.distmon.text yview"
pack  $w.distmon.scrollbar -side right -fill y
pack  $w.distmon.text -side left -fill both

set DefDist [string trim [show monitor distance list]]

if {[llength $DefDist] < 2} return

set DistValues [lindex $DefDist 0]

        for {set i 1} {$i <= $DistValues} {incr i} {

        set Type  [string trim [show monitor distance type $i]]
        set Color [show monitor distance colour $i]
        scan $Color "%f %f %f" Red Green Blue

        frame  $w.distmon.text.frame$i
		entry  $w.distmon.text.frame$i.entry$i -width 34

        set AtomIndex   [lindex $DefDist $i]
		set AtomIndex1  [lindex $AtomIndex 0]
		set AtomIndex2  [lindex $AtomIndex 1]

        set AtomInfo "\[$i\]\
		              ([show atom segmentname $AtomIndex1 1]:\
		               [show atom residuename $AtomIndex1 1]:\
					   [show atom atomname    $AtomIndex1 1]($AtomIndex1)) - \
					  ([show atom segmentname $AtomIndex2 1]:\
					   [show atom residuename $AtomIndex2 1]:\
					   [show atom atomname    $AtomIndex2 1]($AtomIndex2)) "
        $w.distmon.text.frame$i.entry$i delete 0 end
		$w.distmon.text.frame$i.entry$i insert 0 $AtomInfo

        menubutton  $w.distmon.text.frame$i.type$i  -text "Line type"   \
                    -menu $w.distmon.text.frame$i.type$i.menu                       \
					-cursor hand2                                       \
                    -borderwidth 2 -relief raised -padx 2
        menu        $w.distmon.text.frame$i.type$i.menu
        set m       $w.distmon.text.frame$i.type$i.menu
        $m add radiobutton -label "* * * *"   -variable gomDistMonLineType -value 1 \
           -command "monitor distance type 1 $i;lulInvalidateDisplay"
        $m add radiobutton -label "- - - -"   -variable gomDistMonLineType -value 2 \
           -command "monitor distance type 2 $i;lulInvalidateDisplay"
        $m add radiobutton -label "* - - *"   -variable gomDistMonLineType -value 3 \
           -command "monitor distance type 3 $i;lulInvalidateDisplay"
        $m add radiobutton -label "-------"   -variable gomDistMonLineType -value 4 \
           -command "monitor distance type 4 $i;lulInvalidateDisplay"
        $m invoke [show monitor distance type $i]

		button $w.distmon.text.frame$i.color$i  -text "Colour"    \
		       -cursor hand2                                      \
			   -bg [lulColourFloat2Hex $Red $Green $Blue] -padx 2 \
			   -command "lulChangeMonitorDistanceLineColor $i $w.distmon.text.frame$i.color$i;lulInvalidateDisplay"  

        button $w.distmon.text.frame$i.edit$i -text "Edit"        \
               -cursor hand2                                      \
               -command "lulMonitorEditDistance $i"
		pack   $w.distmon.text.frame$i.entry$i -side left -padx 3
		pack   $w.distmon.text.frame$i.type$i  -side left
		pack   $w.distmon.text.frame$i.color$i -side left
        pack   $w.distmon.text.frame$i.edit$i  -side left
        $w.distmon.text window create $i.0 -window $w.distmon.text.frame$i
        }
}
###############################################################
# PROC
proc lulDoMonitorDistance { w } {

     global gomAppendDistList
     global gomDistMonLineType

   foreach set {1 2} {
      set Seg$set [string trim [$w.atoms.set$set.segment.input get]]
      set Res$set [string trim [$w.atoms.set$set.residue.input get]]
      set Atm$set [string trim [$w.atoms.set$set.atom.input get]]
   }

#   set LineType    "0x0101"
   set LineType    $gomDistMonLineType
   set LineColour  [$w.properties.color cget -bg]
#
# check if "append" is disabled ...
   if {!$gomAppendDistList} {
       puts "Reseting distance monitoring..."
       lulDeleteMonitorDistance $w
   }

#   puts "Distance monitoring on: $Seg1 $Res1 $Atm1 $Seg2 $Res2 $Atm2 $LineType {[lulColourHex2Float $LineColour]}"
   monitor distance $Seg1 $Res1 $Atm1 $Seg2 $Res2 $Atm2 $LineType [lulColourHex2Float $LineColour]

# update the list into widget ...
   set DefDist [string trim [show monitor distance list]]
   set i [lindex $DefDist 0]

        set Color [show monitor distance colour $i]
        scan $Color "%f %f %f" Red Green Blue

        frame  $w.distmon.text.frame$i
		entry  $w.distmon.text.frame$i.entry$i -width 34

        set AtomIndex   [lindex $DefDist $i]
		set AtomIndex1  [lindex $AtomIndex 0]
		set AtomIndex2  [lindex $AtomIndex 1]

        set AtomInfo "\[$i\]\
		              ([show atom segmentname $AtomIndex1 1]:\
		               [show atom residuename $AtomIndex1 1]:\
					   [show atom atomname    $AtomIndex1 1]($AtomIndex1)) - \
					  ([show atom segmentname $AtomIndex2 1]:\
					   [show atom residuename $AtomIndex2 1]:\
					   [show atom atomname    $AtomIndex2 1]($AtomIndex2)) "
        $w.distmon.text.frame$i.entry$i delete 0 end
		$w.distmon.text.frame$i.entry$i insert 0 $AtomInfo

        menubutton  $w.distmon.text.frame$i.type$i  -text "Line type"   \
                    -menu $w.distmon.text.frame$i.type$i.menu                       \
					-cursor hand2                                       \
                    -borderwidth 2 -relief raised  -padx 2
        menu        $w.distmon.text.frame$i.type$i.menu
        set m       $w.distmon.text.frame$i.type$i.menu
        $m add radiobutton -label "* * * *"   -variable gomDistMonLineType -value 1 \
           -command "monitor distance type 1 $i;lulInvalidateDisplay"
        $m add radiobutton -label "- - - -"   -variable gomDistMonLineType -value 2 \
           -command "monitor distance type 2 $i;lulInvalidateDisplay"
        $m add radiobutton -label "* - - *"   -variable gomDistMonLineType -value 3 \
           -command "monitor distance type 3 $i;lulInvalidateDisplay"
        $m add radiobutton -label "-------"   -variable gomDistMonLineType -value 4 \
           -command "monitor distance type 4 $i;lulInvalidateDisplay"
        $m invoke [show monitor distance type $i]

		button $w.distmon.text.frame$i.color$i  -text "Colour"    \
		       -cursor hand2                                      \
			   -bg [lulColourFloat2Hex $Red $Green $Blue] -padx 2 \
			   -command "lulChangeMonitorDistanceLineColor $i $w.distmon.text.frame$i.color$i;lulInvalidateDisplay"  

        button $w.distmon.text.frame$i.edit$i -text "Edit"        \
               -cursor hand2                                      \
               -command "lulMonitorEditDistance $i"

		pack   $w.distmon.text.frame$i.entry$i -side left
		pack   $w.distmon.text.frame$i.type$i  -side left -padx 3
		pack   $w.distmon.text.frame$i.color$i -side left
        pack   $w.distmon.text.frame$i.edit$i  -side left
        $w.distmon.text window create $i.0 -window $w.distmon.text.frame$i

	lulInvalidateDisplay
}

###############################################################
# PROC
proc lulDeleteMonitorDistance { w } {

   set Loop [lindex [string trim [show monitor distance list]] 0]

   monitor -distance

   set DispState [show monitor distance state]

   for {set i 1} {$i <= $Loop} {incr i} {

       destroy $w.distmon.text.frame$i
   }

   if {$DispState == 1} {
       $w.options.right.frame4.on select
      } else {
       $w.options.right.frame4.off select
   }
}
###############################################################
# PROC
proc lulCleanMonitorDistanceWidget { } {

   set Loop [lindex [string trim [show monitor distance list]] 0]

   set DispState [show monitor distance state]

   for {set i 1} {$i <= $Loop} {incr i} {

       destroy .gommonitordistance.distmon.text.frame$i
   }

}


###############################################################
# PROC
proc lulMonitorDistanceLineColor { w } {

    set colour [tk_chooseColor -title "Choose line colour" -parent $w \
         -initialcolor [$w.properties.color cget -bg]]

    if { $colour != "" } {
        $w.properties.color configure -bg $colour
    }
}

###############################################################
# PROC
proc lulChangeMonitorDistanceLineColor { i w } {

    set colour [tk_chooseColor -title "Change line colour" -parent $w \
         -initialcolor [$w cget -bg]]

    if { $colour != "" } {
        $w configure -bg $colour
        eval monitor distance colour {[lulColourHex2Float $colour]} $i
    }
}

##################### monitor angle ############################
#
# monitor angle
#
# PROC
proc lulMonitorAngle {} {

     global gomHelpFile
     global gomSwitchAngMon
     global gomAppendDistList
     global gomControlFont

# return if no molecular systems defined
#     if {[show trajectory frames] < 1} {
#          lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
#	      return
#	 }

set w .gommonitorangle
catch {destroy $w}
toplevel $w 
wm title $w "Monitor angle"
wm iconname $w "Monitor angle"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulDoMonitorAngle $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(monitorangle)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Monitor angle:"
pack   $w.label -side top -anchor w

frame  $w.atoms
pack   $w.atoms -side top -anchor w

# set #1,#2,#3
foreach set {1 2 3} {
    frame  $w.atoms.set$set -borderwidth 2 -relief ridge
    pack   $w.atoms.set$set -side left
    label  $w.atoms.set$set.label -text "Atom(s) #$set:"
    pack   $w.atoms.set$set.label -side top -anchor w
    lulCreateAtomInputEntries $w.atoms.set$set 10 "*" "*" "*" "" 1
}

frame       $w.properties -borderwidth 2 -relief ridge
pack        $w.properties -side top -anchor w -fill x

button      $w.properties.color -text "Line colour..." -bg #ffffff \
             -command "lulMonitorAngleLineColor $w"
pack        $w.properties.color -side left -anchor w -padx 3

menubutton  $w.properties.type  -text "Line type"   \
             -menu $w.properties.type.menu          \
             -borderwidth 2 -relief raised
menu        $w.properties.type.menu
set m       $w.properties.type.menu
$m add radiobutton -label "* * * *"    -variable gomAngMonLineType -value 1
$m add radiobutton -label "- - - -"    -variable gomAngMonLineType -value 2
$m add radiobutton -label "* - - *"    -variable gomAngMonLineType -value 3
$m add radiobutton -label "-------"    -variable gomAngMonLineType -value 4
pack        $w.properties.type  -side left -anchor w
$m invoke 2

frame       $w.options 
pack        $w.options -side top -anchor w -fill x

frame       $w.options.left
pack        $w.options.left   -side left -anchor w

frame       $w.options.right 
pack        $w.options.right  -side right -anchor e

frame       $w.options.left.frame8 -borderwidth 2 -relief ridge
button      $w.options.left.frame8.delete -text "Delete" \
            -command "lulDeleteMonitorAngle $w;lulInvalidateDisplay"

pack        $w.options.left.frame8 -side bottom -anchor w
pack        $w.options.left.frame8.delete -side left

frame       $w.options.left.frame9 -borderwidth 2 -relief ridge
checkbutton $w.options.left.frame9.append -text "Append"    \
             -variable gomAppendAngList \
              -onvalue 1 -offvalue 0

pack        $w.options.left.frame9 -side bottom -anchor w
pack        $w.options.left.frame9.append -side left

$w.options.left.frame9.append select

frame       $w.options.right.frame4 -borderwidth 2 -relief ridge
label       $w.options.right.frame4.label -text "Display:"
radiobutton $w.options.right.frame4.on  -text "On"  \
            -variable gomSwitchAngMon -value 1 \
			-command  {eval monitor display angle on;lulInvalidateDisplay}  
radiobutton $w.options.right.frame4.off -text "Off" \
            -variable gomSwitchAngMon -value 0 \
			-command  {eval monitor display angle off;lulInvalidateDisplay}
pack $w.options.right.frame4 -side top -anchor e
pack $w.options.right.frame4.label \
     $w.options.right.frame4.on    \
	 $w.options.right.frame4.off \
     -side top -anchor w

# set the display state
set DispState [show monitor angle state]

if {$DispState == 1} {
    $w.options.right.frame4.on select
  } else {
    $w.options.right.frame4.off select
  }

# 
frame $w.angmon -borderwidth 2 -relief ridge 
pack  $w.angmon -side top 
text  $w.angmon.text \
      -yscrollcommand "$w.angmon.scrollbar set" \
	  -width 60 -height 10
scrollbar $w.angmon.scrollbar -command "$w.angmon.text yview"
pack  $w.angmon.scrollbar -side right -fill y
pack  $w.angmon.text -side left -fill both

set DefAng [string trim [show monitor angle list]]

if {[llength $DefAng] < 2} return

set AngValues [lindex $DefAng 0]

        for {set i 1} {$i <= $AngValues} {incr i} {

        set Color [show monitor angle colour $i]
        scan $Color "%f %f %f" Red Green Blue

        frame  $w.angmon.text.frame$i
		entry  $w.angmon.text.frame$i.entry$i -width 34

        set AtomIndex   [lindex $DefAng $i]
		set AtomIndex1  [lindex $AtomIndex 0]
		set AtomIndex2  [lindex $AtomIndex 1]
		set AtomIndex3  [lindex $AtomIndex 2]

        set AtomInfo "\[$i\]\
		              ([show atom segmentname $AtomIndex1 1]:\
		               [show atom residuename $AtomIndex1 1]:\
					   [show atom atomname    $AtomIndex1 1]($AtomIndex1)) - \
					  ([show atom segmentname $AtomIndex2 1]:\
					   [show atom residuename $AtomIndex2 1]:\
					   [show atom atomname    $AtomIndex2 1]($AtomIndex2)) - \
					  ([show atom segmentname $AtomIndex3 1]:\
					   [show atom residuename $AtomIndex3 1]:\
					   [show atom atomname    $AtomIndex3 1]($AtomIndex3))"
        $w.angmon.text.frame$i.entry$i delete 0 end
		$w.angmon.text.frame$i.entry$i insert 0 $AtomInfo

        menubutton  $w.angmon.text.frame$i.type$i  -text "Line type"   \
                    -menu $w.angmon.text.frame$i.type$i.menu                       \
					-cursor hand2                                       \
                    -borderwidth 2 -relief raised -padx 2
        menu        $w.angmon.text.frame$i.type$i.menu
        set m       $w.angmon.text.frame$i.type$i.menu
        $m add radiobutton -label "* * * *"   -variable gomAngMonLineType -value 1 \
           -command "monitor angle type 1 $i;lulInvalidateDisplay"
        $m add radiobutton -label "- - - -"   -variable gomAngMonLineType -value 2 \
           -command "monitor angle type 2 $i;lulInvalidateDisplay"
        $m add radiobutton -label "* - - *"   -variable gomAngMonLineType -value 3 \
           -command "monitor angle type 3 $i;lulInvalidateDisplay"
        $m add radiobutton -label "-------"   -variable gomAngMonLineType -value 4 \
           -command "monitor angle type 4 $i;lulInvalidateDisplay"
        $m invoke [show monitor angle type $i]

		button $w.angmon.text.frame$i.color$i  -text "Colour"    \
		       -cursor hand2                                      \
			   -bg [lulColourFloat2Hex $Red $Green $Blue] -padx 2 \
			   -command "lulChangeMonitorAngleLineColor $i $w.angmon.text.frame$i.color$i;lulInvalidateDisplay"  

        button $w.angmon.text.frame$i.edit$i -text "Edit"        \
               -cursor hand2                                      \
               -command "lulMonitorEditAngle $i"
		pack   $w.angmon.text.frame$i.entry$i -side left
		pack   $w.angmon.text.frame$i.type$i  -side left -padx 3
		pack   $w.angmon.text.frame$i.color$i -side left
        pack   $w.angmon.text.frame$i.edit$i  -side left
        $w.angmon.text window create $i.0 -window $w.angmon.text.frame$i
        }
}
###############################################################
# PROC
proc lulDoMonitorAngle { w } {

     global gomAppendAngList
     global gomAngMonLineType

   foreach set {1 2 3} {
      set Seg$set [string trim [$w.atoms.set$set.segment.input get]]
      set Res$set [string trim [$w.atoms.set$set.residue.input get]]
      set Atm$set [string trim [$w.atoms.set$set.atom.input get]]
   }

#   set LineType    "0x0101"
   set LineType    $gomAngMonLineType
   set LineColour  [$w.properties.color cget -bg]
#
# check if "append" is disabled ...
   if {!$gomAppendAngList} {
       puts "Reseting angle monitoring..."
       lulDeleteMonitorAngle $w
   }

   monitor angle $Seg1 $Res1 $Atm1 $Seg2 $Res2 $Atm2 $Seg3 $Res3 $Atm3 $LineType [lulColourHex2Float $LineColour]

# update the list into widget ...
   set DefAng [string trim [show monitor angle list]]
   set i [lindex $DefAng 0]

        set Color [show monitor angle colour $i]
        scan $Color "%f %f %f" Red Green Blue

        frame  $w.angmon.text.frame$i
		entry  $w.angmon.text.frame$i.entry$i -width 34

        set AtomIndex   [lindex $DefAng $i]
		set AtomIndex1  [lindex $AtomIndex 0]
		set AtomIndex2  [lindex $AtomIndex 1]
		set AtomIndex3  [lindex $AtomIndex 2]

        set AtomInfo "\[$i\]\
		              ([show atom segmentname $AtomIndex1 1]:\
		               [show atom residuename $AtomIndex1 1]:\
					   [show atom atomname    $AtomIndex1 1]($AtomIndex1)) - \
					  ([show atom segmentname $AtomIndex2 1]:\
					   [show atom residuename $AtomIndex2 1]:\
					   [show atom atomname    $AtomIndex2 1]($AtomIndex2)) - \
					  ([show atom segmentname $AtomIndex3 1]:\
					   [show atom residuename $AtomIndex3 1]:\
					   [show atom atomname    $AtomIndex3 1]($AtomIndex3))"
        $w.angmon.text.frame$i.entry$i delete 0 end
		$w.angmon.text.frame$i.entry$i insert 0 $AtomInfo

        menubutton  $w.angmon.text.frame$i.type$i  -text "Line type"   \
                    -menu $w.angmon.text.frame$i.type$i.menu                       \
					-cursor hand2                                       \
                    -borderwidth 2 -relief raised  -padx 2
        menu        $w.angmon.text.frame$i.type$i.menu
        set m       $w.angmon.text.frame$i.type$i.menu
        $m add radiobutton -label "* * * *"   -variable gomAngMonLineType -value 1 \
           -command "monitor angle type 1 $i;lulInvalidateDisplay"
        $m add radiobutton -label "- - - -"   -variable gomAngMonLineType -value 2 \
           -command "monitor angle type 2 $i;lulInvalidateDisplay"
        $m add radiobutton -label "* - - *"   -variable gomAngMonLineType -value 3 \
           -command "monitor angle type 3 $i;lulInvalidateDisplay"
        $m add radiobutton -label "-------"   -variable gomAngMonLineType -value 4 \
           -command "monitor angle type 4 $i;lulInvalidateDisplay"
        $m invoke [show monitor angle type $i]

		button $w.angmon.text.frame$i.color$i  -text "Colour"    \
		       -cursor hand2                                      \
			   -bg [lulColourFloat2Hex $Red $Green $Blue] -padx 2 \
			   -command "lulChangeMonitorAngleLineColor $i $w.angmon.text.frame$i.color$i;lulInvalidateDisplay"  

        button $w.angmon.text.frame$i.edit$i -text "Edit"        \
               -cursor hand2                                      \
               -command "lulMonitorEditAngle $i"
		pack   $w.angmon.text.frame$i.entry$i -side left
		pack   $w.angmon.text.frame$i.type$i  -side left -padx 3
		pack   $w.angmon.text.frame$i.color$i -side left
        pack   $w.angmon.text.frame$i.edit$i  -side left
        $w.angmon.text window create $i.0 -window $w.angmon.text.frame$i

	lulInvalidateDisplay
}

###############################################################
# PROC
proc lulDeleteMonitorAngle { w } {

   set Loop [lindex [string trim [show monitor angle list]] 0]

   monitor -angle

   set DispState [show monitor angle state]

   for {set i 1} {$i <= $Loop} {incr i} {

       destroy $w.angmon.text.frame$i
   }

   if {$DispState == 1} {
       $w.options.right.frame4.on select
      } else {
       $w.options.right.frame4.off select
   }
}
###############################################################
# PROC
proc lulCleanMonitorAngleWidget { } {

   set Loop [lindex [string trim [show monitor angle list]] 0]

   set DispState [show monitor angle state]

   for {set i 1} {$i <= $Loop} {incr i} {

       destroy .gommonitorangle.angmon.text.frame$i
   }

}

###############################################################
# PROC
proc lulMonitorAngleLineColor { w } {

    set colour [tk_chooseColor -title "Choose line colour" -parent $w \
         -initialcolor [$w.properties.color cget -bg]]

    if { $colour != "" } {
        $w.properties.color configure -bg $colour
    }
}

###############################################################
# PROC
proc lulChangeMonitorAngleLineColor { i w } {

    set colour [tk_chooseColor -title "Change line colour" -parent $w \
         -initialcolor [$w cget -bg]]

    if { $colour != "" } {
        $w configure -bg $colour
        eval monitor angle colour {[lulColourHex2Float $colour]} $i
    }
}

##################### monitor torsion ############################
#
# monitor torsion
#
# PROC
proc lulMonitorTorsion {} {

     global gomHelpFile
     global gomSwitchTorsMon
     global gomAppendDistList
     global gomControlFont

# return if no molecular systems defined
#     if {[show trajectory frames] < 1} {
#          lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
#	      return
#	 }

set w .gommonitortorsion
catch {destroy $w}
toplevel $w 
wm title $w "Monitor torsion"
wm iconname $w "Monitor torsion"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulDoMonitorTorsion $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(monitortorsion)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Monitor torsion:"
pack   $w.label -side top -anchor w

frame  $w.atoms
pack   $w.atoms -side top -anchor w

# set #1,#2,#3,#4
foreach set {1 2 3 4} {
    frame  $w.atoms.set$set -borderwidth 2 -relief ridge
    pack   $w.atoms.set$set -side left
    label  $w.atoms.set$set.label -text "Atom(s) #$set:"
    pack   $w.atoms.set$set.label -side top -anchor w
    lulCreateAtomInputEntries $w.atoms.set$set 7 "*" "*" "*" "" 1
}

frame       $w.properties -borderwidth 2 -relief ridge
pack        $w.properties -side top -anchor w -fill x

button      $w.properties.color -text "Line colour..." -bg #ffffff \
             -command "lulMonitorTorsionLineColor $w"
pack        $w.properties.color -side left -anchor w -padx 3

menubutton  $w.properties.type  -text "Line type"   \
             -menu $w.properties.type.menu          \
             -borderwidth 2 -relief raised
menu        $w.properties.type.menu
set m       $w.properties.type.menu
$m add radiobutton -label "* * * *"    -variable gomTorsMonLineType -value 1
$m add radiobutton -label "- - - -"    -variable gomTorsMonLineType -value 2
$m add radiobutton -label "* - - *"    -variable gomTorsMonLineType -value 3
$m add radiobutton -label "-------"    -variable gomTorsMonLineType -value 4
pack        $w.properties.type  -side left -anchor w
$m invoke 3

frame       $w.options 
pack        $w.options -side top -anchor w -fill x

frame       $w.options.left
pack        $w.options.left   -side left -anchor w

frame       $w.options.right 
pack        $w.options.right  -side right -anchor e

frame       $w.options.left.frame8 -borderwidth 2 -relief ridge
button      $w.options.left.frame8.delete -text "Delete" \
            -command "lulDeleteMonitorTorsion $w;lulInvalidateDisplay"

pack        $w.options.left.frame8 -side bottom -anchor w
pack        $w.options.left.frame8.delete -side left

frame       $w.options.left.frame9 -borderwidth 2 -relief ridge
checkbutton $w.options.left.frame9.append -text "Append"    \
             -variable gomAppendTorsList \
              -onvalue 1 -offvalue 0

pack        $w.options.left.frame9 -side bottom -anchor w
pack        $w.options.left.frame9.append -side left

$w.options.left.frame9.append select

frame       $w.options.right.frame4 -borderwidth 2 -relief ridge
label       $w.options.right.frame4.label -text "Display:"
radiobutton $w.options.right.frame4.on  -text "On"  \
            -variable gomSwitchTorsMon -value 1 \
			-command  {eval monitor display torsion on;lulInvalidateDisplay}  
radiobutton $w.options.right.frame4.off -text "Off" \
            -variable gomSwitchTorsMon -value 0 \
			-command  {eval monitor display torsion off;lulInvalidateDisplay}
pack $w.options.right.frame4 -side top -anchor e
pack $w.options.right.frame4.label \
     $w.options.right.frame4.on    \
	 $w.options.right.frame4.off \
     -side top -anchor w

# set the display state
set DispState [show monitor torsion state]

if {$DispState == 1} {
    $w.options.right.frame4.on select
  } else {
    $w.options.right.frame4.off select
  }

# 
frame $w.torsmon -borderwidth 2 -relief ridge 
pack  $w.torsmon -side top 
text  $w.torsmon.text \
      -yscrollcommand "$w.torsmon.scrollbar set" \
	  -width 60 -height 10
scrollbar $w.torsmon.scrollbar -command "$w.torsmon.text yview"
pack  $w.torsmon.scrollbar -side right -fill y
pack  $w.torsmon.text -side left -fill both

set DefTors [string trim [show monitor torsion list]]

if {[llength $DefTors] < 2} return

set TorsValues [lindex $DefTors 0]

        for {set i 1} {$i <= $TorsValues} {incr i} {

        set Color [show monitor torsion colour $i]
        scan $Color "%f %f %f" Red Green Blue

        frame  $w.torsmon.text.frame$i
		entry  $w.torsmon.text.frame$i.entry$i -width 34

        set AtomIndex   [lindex $DefTors $i]
		set AtomIndex1  [lindex $AtomIndex 0]
		set AtomIndex2  [lindex $AtomIndex 1]
		set AtomIndex3  [lindex $AtomIndex 2]
		set AtomIndex4  [lindex $AtomIndex 3]

        set AtomInfo "\[$i\]\
		              ([show atom segmentname $AtomIndex1 1]:\
		               [show atom residuename $AtomIndex1 1]:\
					   [show atom atomname    $AtomIndex1 1]($AtomIndex1)) - \
					  ([show atom segmentname $AtomIndex2 1]:\
					   [show atom residuename $AtomIndex2 1]:\
					   [show atom atomname    $AtomIndex2 1]($AtomIndex2)) - \
					  ([show atom segmentname $AtomIndex3 1]:\
					   [show atom residuename $AtomIndex3 1]:\
					   [show atom atomname    $AtomIndex3 1]($AtomIndex3)) - \
					  ([show atom segmentname $AtomIndex4 1]:\
					   [show atom residuename $AtomIndex4 1]:\
					   [show atom atomname    $AtomIndex4 1]($AtomIndex4))"
        $w.torsmon.text.frame$i.entry$i delete 0 end
		$w.torsmon.text.frame$i.entry$i insert 0 $AtomInfo

        menubutton  $w.torsmon.text.frame$i.type$i  -text "Line type"   \
                    -menu $w.torsmon.text.frame$i.type$i.menu                       \
					-cursor hand2                                       \
                    -borderwidth 2 -relief raised -padx 2
        menu        $w.torsmon.text.frame$i.type$i.menu
        set m       $w.torsmon.text.frame$i.type$i.menu
        $m add radiobutton -label "* * * *"   -variable gomTorsMonLineType -value 1 \
           -command "monitor torsion type 1 $i;lulInvalidateDisplay"
        $m add radiobutton -label "- - - -"   -variable gomTorsMonLineType -value 2 \
           -command "monitor torsion type 2 $i;lulInvalidateDisplay"
        $m add radiobutton -label "* - - *"   -variable gomTorsMonLineType -value 3 \
           -command "monitor torsion type 3 $i;lulInvalidateDisplay"
        $m add radiobutton -label "-------"   -variable gomTorsMonLineType -value 4 \
           -command "monitor torsion type 4 $i;lulInvalidateDisplay"
        $m invoke [show monitor torsion type $i]

		button $w.torsmon.text.frame$i.color$i  -text "Colour"    \
		       -cursor hand2                                      \
			   -bg [lulColourFloat2Hex $Red $Green $Blue] -padx 2 \
			   -command "lulChangeMonitorTorsionLineColor $i $w.torsmon.text.frame$i.color$i;lulInvalidateDisplay"  

        button $w.torsmon.text.frame$i.edit$i -text "Edit"        \
               -cursor hand2                                      \
               -command "lulMonitorEditTorsion $i"
		pack   $w.torsmon.text.frame$i.entry$i -side left
		pack   $w.torsmon.text.frame$i.type$i  -side left -padx 3
		pack   $w.torsmon.text.frame$i.color$i -side left
        pack   $w.torsmon.text.frame$i.edit$i  -side left
        $w.torsmon.text window create $i.0 -window $w.torsmon.text.frame$i
        }
}
###############################################################
# PROC
proc lulDoMonitorTorsion { w } {

     global gomAppendTorsList
     global gomTorsMonLineType

   foreach set {1 2 3 4} {
      set Seg$set [string trim [$w.atoms.set$set.segment.input get]]
      set Res$set [string trim [$w.atoms.set$set.residue.input get]]
      set Atm$set [string trim [$w.atoms.set$set.atom.input get]]
   }

#   set LineType    "0x0101"
   set LineType    $gomTorsMonLineType
   set LineColour  [$w.properties.color cget -bg]
#
# check if "append" is disabled ...
   if {!$gomAppendTorsList} {
       puts "Reseting torsion monitoring..."
       lulDeleteMonitorTorsion $w
   }

#   puts "Torsion monitoring on: $Seg1 $Res1 $Atm1 $Seg2 $Res2 $Atm2 $Seg3 $Res3 $Atm3 $Seg4 $Res4 $Atm4 $LineType {[lulColourHex2Float $LineColour]}"
   monitor torsion $Seg1 $Res1 $Atm1 $Seg2 $Res2 $Atm2 $Seg3 $Res3 $Atm3 $Seg4 $Res4 $Atm4 $LineType [lulColourHex2Float $LineColour]

# update the list into widget ...
   set DefTors [string trim [show monitor torsion list]]
   set i [lindex $DefTors 0]

        set Color [show monitor torsion colour $i]
        scan $Color "%f %f %f" Red Green Blue

        frame  $w.torsmon.text.frame$i
		entry  $w.torsmon.text.frame$i.entry$i -width 34

        set AtomIndex   [lindex $DefTors $i]
		set AtomIndex1  [lindex $AtomIndex 0]
		set AtomIndex2  [lindex $AtomIndex 1]
		set AtomIndex3  [lindex $AtomIndex 2]
		set AtomIndex4  [lindex $AtomIndex 3]

        set AtomInfo "\[$i\]\
		              ([show atom segmentname $AtomIndex1 1]:\
		               [show atom residuename $AtomIndex1 1]:\
					   [show atom atomname    $AtomIndex1 1]($AtomIndex1)) - \
					  ([show atom segmentname $AtomIndex2 1]:\
					   [show atom residuename $AtomIndex2 1]:\
					   [show atom atomname    $AtomIndex2 1]($AtomIndex2)) - \
					  ([show atom segmentname $AtomIndex3 1]:\
					   [show atom residuename $AtomIndex3 1]:\
					   [show atom atomname    $AtomIndex3 1]($AtomIndex3)) - \
					  ([show atom segmentname $AtomIndex4 1]:\
					   [show atom residuename $AtomIndex4 1]:\
					   [show atom atomname    $AtomIndex4 1]($AtomIndex4))"
        $w.torsmon.text.frame$i.entry$i delete 0 end
		$w.torsmon.text.frame$i.entry$i insert 0 $AtomInfo

        menubutton  $w.torsmon.text.frame$i.type$i  -text "Line type"   \
                    -menu $w.torsmon.text.frame$i.type$i.menu                       \
					-cursor hand2                                       \
                    -borderwidth 2 -relief raised  -padx 2
        menu        $w.torsmon.text.frame$i.type$i.menu
        set m       $w.torsmon.text.frame$i.type$i.menu
        $m add radiobutton -label "* * * *"   -variable gomTorsMonLineType -value 1 \
           -command "monitor torsion type 1 $i;lulInvalidateDisplay"
        $m add radiobutton -label "- - - -"   -variable gomTorsMonLineType -value 2 \
           -command "monitor torsion type 2 $i;lulInvalidateDisplay"
        $m add radiobutton -label "* - - *"   -variable gomTorsMonLineType -value 3 \
           -command "monitor torsion type 3 $i;lulInvalidateDisplay"
        $m add radiobutton -label "-------"   -variable gomTorsMonLineType -value 4 \
           -command "monitor torsion type 4 $i;lulInvalidateDisplay"
        $m invoke [show monitor torsion type $i]

		button $w.torsmon.text.frame$i.color$i  -text "Colour"    \
		       -cursor hand2                                      \
			   -bg [lulColourFloat2Hex $Red $Green $Blue] -padx 2 \
			   -command "lulChangeMonitorTorsionLineColor $i $w.torsmon.text.frame$i.color$i;lulInvalidateDisplay"  

        button $w.torsmon.text.frame$i.edit$i -text "Edit"        \
               -cursor hand2                                      \
               -command "lulMonitorEditTorsion $i"
		pack   $w.torsmon.text.frame$i.entry$i -side left
		pack   $w.torsmon.text.frame$i.type$i  -side left -padx 3
		pack   $w.torsmon.text.frame$i.color$i -side left
        pack   $w.torsmon.text.frame$i.edit$i  -side left
        $w.torsmon.text window create $i.0 -window $w.torsmon.text.frame$i

	lulInvalidateDisplay
}

###############################################################
# PROC
proc lulDeleteMonitorTorsion { w } {

   set Loop [lindex [string trim [show monitor torsion list]] 0]

   monitor -torsion

   set DispState [show monitor torsion state]

   for {set i 1} {$i <= $Loop} {incr i} {

       destroy $w.torsmon.text.frame$i
   }

   if {$DispState == 1} {
       $w.options.right.frame4.on select
      } else {
       $w.options.right.frame4.off select
   }
}
###############################################################
# PROC
proc lulCleanMonitorTorsionWidget { w } {

   set Loop [lindex [string trim [show monitor torsion list]] 0]

   set DispState [show monitor torsion state]

   for {set i 1} {$i <= $Loop} {incr i} {

       destroy .gommonitortorsion.torsmon.text.frame$i
   }

}

###############################################################
# PROC
proc lulMonitorTorsionLineColor { w } {

    set colour [tk_chooseColor -title "Choose line colour" -parent $w \
         -initialcolor [$w.properties.color cget -bg]]

    if { $colour != "" } {
        $w.properties.color configure -bg $colour
    }
}

###############################################################
# PROC
proc lulChangeMonitorTorsionLineColor { i w } {

    set colour [tk_chooseColor -title "Change line colour" -parent $w \
         -initialcolor [$w cget -bg]]

    if { $colour != "" } {
        $w configure -bg $colour
        eval monitor torsion colour {[lulColourHex2Float $colour]} $i
    }
}



##################### edit molecule ############################
#
# edit molecule
#
# PROC
proc lulEditMolecule {} {

     global gomHelpFile
     global gomControlFont

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomeditmol
catch {destroy $w}
toplevel $w 
wm title $w "Edit molecule"
wm iconname $w "Edit molecule"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -state disabled
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(editmol)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Edit molecule:"
pack   $w.label -side top -anchor w

frame  $w.atoms
pack   $w.atoms -side top -anchor w

# set #1
frame  $w.atoms.set1 -borderwidth 2 -relief ridge
pack   $w.atoms.set1 -side left
label  $w.atoms.set1.label1 -text "Atom #1:"
pack   $w.atoms.set1.label1 -side top -anchor w
lulCreateAtomInputEntries $w.atoms.set1 15 "*" "*" "*" "" 1

# set #2
frame  $w.atoms.set2 -borderwidth 2 -relief ridge
pack   $w.atoms.set2 -side left
label  $w.atoms.set2.label1 -text "Atom #2:"
pack   $w.atoms.set2.label1 -side top -anchor w
lulCreateAtomInputEntries $w.atoms.set2 15 "*" "*" "*" "" 1

frame       $w.options 
pack        $w.options -side top -anchor w -fill x

frame       $w.options.left
pack        $w.options.left   -side left -anchor w

frame       $w.options.left.bond -borderwidth 2 -relief ridge
pack        $w.options.left.bond -side left

button      $w.options.left.bond.create -text "Create bond(s)" \
             -command "lulCreateBond $w bond;lulInvalidateDisplay" -width 12
button      $w.options.left.bond.break  -text "Break bond(s)"  \
             -command "lulBreakBond  $w bond;lulInvalidateDisplay" -width 12
pack        $w.options.left.bond.create -side top -padx 2 -pady 3
pack        $w.options.left.bond.break  -side top -padx 2 -pady 3

frame       $w.options.left.hbond -borderwidth 2 -relief ridge
pack        $w.options.left.hbond -side left

button      $w.options.left.hbond.create -text "Create hydrogen bond(s)" \
             -command "lulCreateBond $w hbond;lulInvalidateDisplay" -width 18
button      $w.options.left.hbond.break  -text "Break hydrogen bond(s)"  \
             -command "lulBreakBond  $w hbond;lulInvalidateDisplay" -width 18
pack        $w.options.left.hbond.create -side top -padx 2 -pady 3
pack        $w.options.left.hbond.break  -side top -padx 2 -pady 3

frame       $w.options.right 
pack        $w.options.right  -side right -anchor e

frame       $w.ssbond -borderwidth 2 -relief ridge
pack        $w.ssbond -side top -anchor w -pady 5

label       $w.ssbond.text   -text "Calculate S-S bonds (Angstrom): " 
entry       $w.ssbond.range  -width 10
button      $w.ssbond.doit   -text "Calculate" \
             -command "find ssbonds [$w.ssbond.range get] all;\
             lulInvalidateDisplay"

pack        $w.ssbond.text $w.ssbond.range \
            $w.ssbond.doit -side left -padx 4

$w.ssbond.range insert 0 "2.12132"

}

###############################################################
# PROC
proc lulCreateBond { w type } {

   set Seg1 [string trim [$w.atoms.set1.segment.input get]]
   set Res1 [string trim [$w.atoms.set1.residue.input get]]
   set Atm1 [string trim [$w.atoms.set1.atom.input get]]

   set Seg2 [string trim [$w.atoms.set2.segment.input get]]
   set Res2 [string trim [$w.atoms.set2.residue.input get]]
   set Atm2 [string trim [$w.atoms.set2.atom.input get]]

   if {($Seg2 == ""  && $Res2 == ""  && $Atm2 == "") ||
       ($Seg2 == "*" && $Res2 == "*" && $Atm2 == "*")} {
       if {$type == "bond"} {
          calculate connect  * * *
       } else {
          calculate hbond    * * *
       }
   } else {
       edit $type create $Seg1 $Res1 $Atm1 $Seg2 $Res2 $Atm2
   }
}
###############################################################
# PROC
proc lulBreakBond { w type } {

   set Seg1 [string trim [$w.atoms.set1.segment.input get]]
   set Res1 [string trim [$w.atoms.set1.residue.input get]]
   set Atm1 [string trim [$w.atoms.set1.atom.input get]]

   set Seg2 [string trim [$w.atoms.set2.segment.input get]]
   set Res2 [string trim [$w.atoms.set2.residue.input get]]
   set Atm2 [string trim [$w.atoms.set2.atom.input get]]

   edit $type break $Seg1 $Res1 $Atm1 $Seg2 $Res2 $Atm2
}

##################### edit display ############################
#
# edit display
#
# PROC
proc lulEditDisplay {} {

     global gomHelpFile
     global gomControlFont

set w .gomeditdisplay
catch {destroy $w}
toplevel $w 
wm title $w "Edit display"
wm iconname $w "Edit display"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulDoEditDisplay $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(editdisplay)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Edit Display:"
pack   $w.label -side top -anchor w

frame  $w.frame1 -borderwidth 2 -relief ridge
pack   $w.frame1 -side top -fill x

label  $w.frame1.label -text "Near plane distance"
pack   $w.frame1.label -side top -anchor w

frame  $w.frame1.step 
pack   $w.frame1.step -side top -anchor w
label  $w.frame1.step.label  -text "Step size: "
pack   $w.frame1.step.label  -side left 
entry  $w.frame1.step.value  -width 10
pack   $w.frame1.step.value  -side left

set StepValue [show nearplane step]
$w.frame1.step.value insert end $StepValue

frame  $w.frame1.nearplane 
pack   $w.frame1.nearplane -side top -pady 3
button $w.frame1.nearplane.dec -text "<-\]" -command "lulDecreaseNearplane $w"
pack   $w.frame1.nearplane.dec -side left -padx 4 -pady 3
entry  $w.frame1.nearplane.value -width 15
pack   $w.frame1.nearplane.value -side left -padx 4 -pady 3
button $w.frame1.nearplane.inc -text "\[+>" -command "lulIncreaseNearplane $w"
pack   $w.frame1.nearplane.inc -side left -padx 4 -pady 3

set Value [show nearplane value]
$w.frame1.nearplane.value insert end $Value

frame  $w.frame1.update 
pack   $w.frame1.update -side top -anchor w
checkbutton $w.frame1.update.status -text "Instant update"         \
            -onvalue 1 -offvalue 0 -variable gomNearPlaneUpdate
pack   $w.frame1.update.status -side top

$w.frame1.update.status deselect

}

##################################################################
# PROC
proc lulDecreaseNearplane { w } {

	 global gomNearPlaneUpdate

set StepValue [$w.frame1.step.value      get]
set NearValue [$w.frame1.nearplane.value get]

set Value [expr $NearValue - $StepValue]
if {$Value < 0.1} {set Value 0.1}

$w.frame1.step.value delete 0 end
$w.frame1.step.value insert end $StepValue

$w.frame1.nearplane.value delete 0 end
$w.frame1.nearplane.value insert end $Value

eval "define nearplane step  $StepValue"
eval "define nearplane value $Value"

if {$gomNearPlaneUpdate} {
   display
} else {
   lulInvalidateDisplay
}
}

##################################################################
# PROC
proc lulIncreaseNearplane { w } {

	 global gomNearPlaneUpdate

set StepValue [$w.frame1.step.value      get]
set NearValue [$w.frame1.nearplane.value get]

set Value [expr $NearValue + $StepValue]
if {$Value < 0.1} {set Value 0.1}

$w.frame1.step.value delete 0 end
$w.frame1.step.value insert end $StepValue

$w.frame1.nearplane.value delete 0 end
$w.frame1.nearplane.value insert end $Value

eval "define nearplane step  $StepValue"
eval "define nearplane value $Value"

if {$gomNearPlaneUpdate} {
   display
} else {
   lulInvalidateDisplay
}
}

##################################################################
# PROC
proc lulDoEditDisplay { w } {

display

}

##################### copy timeseries ############################
#
# copy timeseries
#
# PROC
proc lulCopyTimeSeries {} {

     global gomHelpFile
     global gomControlFont


set w .gomcopytimeseries
catch {destroy $w}
toplevel $w 
wm title $w "Copy timeseries"
wm iconname $w "Copy timeseries"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss   -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply     -font "$gomControlFont" \
        -command "lulDoCopyTimeSeries $w"
button $w.buttons.help    -text Help      -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(copytimeseries)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Copy TimeSeries:"
pack   $w.label -side top -anchor w

frame  $w.frame1 -borderwidth 2 -relief ridge -width 30
pack   $w.frame1 -side top -fill x

set DistList [show monitor distance timeseries]

if {$DistList} {
  radiobutton $w.frame1.distlist -text "Distance list #: " \
               -value 0 -variable gomTimeList 
  pack        $w.frame1.distlist -side left -anchor w 
  entry       $w.frame1.value -width 5 
  pack        $w.frame1.value -side right -anchor e
} else {
  label       $w.frame1.label -text "No distance time series are defined"
  pack        $w.frame1.label -side top -anchor w
  }

frame  $w.frame2 -borderwidth 2 -relief ridge
pack   $w.frame2 -side top -fill x

set AngList [show monitor angle timeseries]

if {$AngList} {
  radiobutton $w.frame2.anglist -text "Angle list #: " \
               -value 1 -variable gomTimeList
  pack        $w.frame2.anglist -side left -anchor w 
  entry       $w.frame2.value -width 5
  pack        $w.frame2.value -side right -anchor e
} else {
  label       $w.frame2.label -text "No angle time series are defined"
  pack        $w.frame2.label -side top -anchor w
  }

frame  $w.frame3 -borderwidth 2 -relief ridge
pack   $w.frame3 -side top -fill x

set TorsList [show monitor torsion timeseries]

if {$TorsList} {
  radiobutton $w.frame3.torslist -text "Torsion list #: " \
               -value 2 -variable gomTimeList
  pack        $w.frame3.torslist -side left -anchor w 
  entry       $w.frame3.value -width 5
  pack        $w.frame3.value -side right -anchor e
} else {
  label       $w.frame3.label -text "No torsion time series are defined"
  pack        $w.frame3.label -side top -anchor w
  }

  if {!$DistList && !$AngList && !$TorsList} {
     $w.buttons.apply configure -state disabled
  }
}

#############################################################
# PROC
proc lulDoCopyTimeSeries { w } {

     global gomTimeList

  switch $gomTimeList {
      0 {copy timeseries distance [$w.frame1.value get]}
      1 {copy timeseries angle    [$w.frame2.value get]}
      2 {copy timeseries torsion  [$w.frame3.value get]}
	  default {lulErrorDialog "Please choose by pressing a radiobutton and giving a time series number"}
  }
}

##################### calculate correlation ############################
#
# calculate correlation
#
# PROC
proc lulCalculateCorrelation {} {

     global gomHelpFile
     global gomControlFont

# return if no molecular systems defined
#     if {[show trajectory frames] < 1} {
#          lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
#	      return
#	 }

set w .gomcalculatecorrelation
catch {destroy $w}
toplevel $w 
wm title $w "Calculate correlation"
wm iconname $w "Calculate correlation"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss   -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply     -font "$gomControlFont" \
        -command "lulDoCalculateCorrelation $w"
button $w.buttons.help    -text Help      -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(calccorrelation)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Calculate correlation:"
pack   $w.label -side top -anchor w

frame  $w.frame0 -borderwidth 2 -relief ridge
pack   $w.frame0 -side top -anchor w

radiobutton $w.frame0.distance -text "Distance" -variable gomCorrType \
            -value 0
radiobutton $w.frame0.angle    -text "Angle"    -variable gomCorrType \
            -value 1
radiobutton $w.frame0.torsion  -text "Torsion"  -variable gomCorrType \
            -value 2

pack  $w.frame0.distance -side left
pack  $w.frame0.angle    -side left
pack  $w.frame0.torsion  -side left


frame  $w.frame1 -borderwidth 2 
pack   $w.frame1 -side top -anchor w

frame  $w.frame1.frame2 -relief ridge
pack   $w.frame1.frame2 -side top -anchor w

label  $w.frame1.frame2.label1 -text "(1) Time series #: "
pack   $w.frame1.frame2.label1 -side left
entry  $w.frame1.frame2.entry1 -width 10
pack   $w.frame1.frame2.entry1 -side left

frame  $w.frame1.frame3 -relief ridge
pack   $w.frame1.frame3 -side top -anchor w

label  $w.frame1.frame3.label2 -text "(2) Time series #: "
pack   $w.frame1.frame3.label2 -side left
entry  $w.frame1.frame3.entry2 -width 10
pack   $w.frame1.frame3.entry2 -side left

}

#############################################################
# PROC
proc lulDoCalculateCorrelation { w } {

     global gomCorrType

  switch $gomCorrType {
      0 {calculate correlation distance [$w.frame1.frame2.entry1 get] [$w.frame1.frame3.entry2 get]}
      1 {calculate correlation angle    [$w.frame1.frame2.entry1 get] [$w.frame1.frame3.entry2 get]}
      2 {calculate correlation torsion  [$w.frame1.frame2.entry1 get] [$w.frame1.frame3.entry2 get]}
	  default {lulErrorDialog "Please choose the correlation type (distance/angle/torsion)"}
  }
}

##################### manipulate time series  ############################
#
# manipulate time series
#
# PROC
proc lulManipulateTimeSeries {} {

     global gomHelpFile
     global gomControlFont

# return if no molecular systems defined
#     if {[show trajectory frames] < 1} {
#          lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
#	      return
#	 }

set w .gommanipulatetimeseries
catch {destroy $w}
toplevel $w 
wm title $w "Manipulate time series"
wm iconname $w "Manipulate time series"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss   -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply     -font "$gomControlFont" \
        -command "lulDoManipulateTimeSeries $w"
button $w.buttons.help    -text Help      -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(mantimeseries)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Manipulate time series:"
pack   $w.label -side top -anchor w

frame  $w.frame0 -borderwidth 2 -relief ridge
pack   $w.frame0 -side top -anchor w

radiobutton $w.frame0.distance -text "Distance" -variable gomTimeSeriesType \
            -value 0
radiobutton $w.frame0.angle    -text "Angle"    -variable gomTimeSeriesType \
            -value 1
radiobutton $w.frame0.torsion  -text "Torsion"  -variable gomTimeSeriesType \
            -value 2

pack  $w.frame0.distance -side left
pack  $w.frame0.angle    -side left
pack  $w.frame0.torsion  -side left


frame  $w.frame1 -borderwidth 2 
pack   $w.frame1 -side top -anchor w

frame  $w.frame1.frame2 -relief ridge
pack   $w.frame1.frame2 -side top -anchor w

label  $w.frame1.frame2.label1 -text "Time series #: "
pack   $w.frame1.frame2.label1 -side left
entry  $w.frame1.frame2.entry1 -width 10
pack   $w.frame1.frame2.entry1 -side left

}

#############################################################
# PROC
proc lulDoManipulateTimeSeries { w } {

     global gomTimeSeriesType

  switch $gomTimeSeriesType {
      0 {calculate correlation distance [$w.frame1.frame2.entry1 get] [$w.frame1.frame3.entry2 get]}
      1 {calculate correlation angle    [$w.frame1.frame2.entry1 get] [$w.frame1.frame3.entry2 get]}
      2 {calculate correlation torsion  [$w.frame1.frame2.entry1 get] [$w.frame1.frame3.entry2 get]}
	  default {lulErrorDialog "Please choose the correlation type (distance/angle/torsion)"}
  }
}

####################################################################
# 
# select cluster atoms ...
#
proc lulSelectClusterAtoms { } {

     global gomHelpFile
     global gomStructure
     global gomSwitch
     global gomControlFont


set w .gomselectclusteratoms
catch {destroy $w}
toplevel $w 
wm title $w "Select cluster atoms"
wm iconname $w "Select cluster atoms"
wm geometry $w 250x250

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyClusterSelection $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(selclusteratom)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

set  NumStruct [show molstructures]

if {$NumStruct} {
    $w.buttons.apply configure -state normal
} else {
    $w.buttons.apply configure -state disabled
}

label $w.label -text "Select clustering atoms:"
pack  $w.label -side top -anchor w

frame $w.frame1              -borderwidth 2 -relief ridge
label $w.frame1.segmentlabel -text "Segment:" -width 10
entry $w.frame1.segment      -width 15
pack  $w.frame1 -side top    -anchor w
pack  $w.frame1.segmentlabel $w.frame1.segment -side left

frame $w.frame2              -borderwidth 2 -relief ridge
label $w.frame2.residuelabel -text "Residue:" -width 10
entry $w.frame2.residue      -width 15
pack  $w.frame2 -side top    -anchor w
pack  $w.frame2.residuelabel $w.frame2.residue -side left

frame $w.frame3              -borderwidth 2 -relief ridge
label $w.frame3.atomlabel    -text "Atom:" -width 10
entry $w.frame3.atom         -width 15
pack  $w.frame3 -side top    -anchor w
pack  $w.frame3.atomlabel    $w.frame3.atom    -side left

frame       $w.frame4
pack        $w.frame4 -side top -anchor w -fill both

frame       $w.frame4.left -borderwidth 2 -relief ridge
label       $w.frame4.left.label -text "Cluster plot:"
radiobutton $w.frame4.left.on    -text "On"    -value 0 \
            -command "plot cluster on;lulInvalidateDisplay"
radiobutton $w.frame4.left.off   -text "Off"   -value 1 \
            -command "plot cluster off;lulInvalidateDisplay"

pack        $w.frame4.left -side left -anchor w
pack        $w.frame4.left.label $w.frame4.left.on $w.frame4.left.off -side top

            $w.frame4.left.off select

frame       $w.frame4.right -borderwidth 2 -relief ridge
pack        $w.frame4.right -side right

button      $w.frame4.right.delete -text "Delete cluster data" \
             -command {calculate -cluster;lulInvalidateDisplay}
pack        $w.frame4.right.delete -side top
}

####################################################################
# PROC
proc lulApplyClusterSelection { w } {

     global gomStructure

     set Segment [string trim [$w.frame1.segment get]]
     set Residue [string trim [$w.frame2.residue get]]
     set Atom    [string trim [$w.frame3.atom    get]]

     eval "calculate cluster \{$Segment\} \{$Residue\} \{$Atom\}"

     lulInvalidateDisplay
}

####################################################################
# 
# select RDF atoms ...
#
proc lulSelectRDFatoms { } {

     global gomHelpFile
     global gomStructure
     global gomSwitch
     global gomControlFont


set w .gomselectrdfatoms
catch {destroy $w}
toplevel $w 
wm title $w "Select RDF atoms"
wm iconname $w "Select RDF atoms"
wm geometry $w 
#280x450

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss   -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply     -font "$gomControlFont" \
        -command "lulApplyRDFselection $w"
button $w.buttons.help    -text Help      -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(selrdfatom)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

set  NumStruct [show molstructures]

if {$NumStruct} {
    $w.buttons.apply configure -state normal
} else {
    $w.buttons.apply configure -state disabled
}

label $w.label -text "Select RDF atoms:" -font {Arial 12 bold}
pack  $w.label -side top -anchor w -pady 3

label $w.atom1  -text "From atom(s):"
pack  $w.atom1  -side top -anchor w

frame $w.frame1              -borderwidth 2 -relief ridge
label $w.frame1.segmentlabel -text "Segment:" -width 10
entry $w.frame1.segment      -width 15
pack  $w.frame1 -side top    -anchor e
pack  $w.frame1.segmentlabel $w.frame1.segment -side left

frame $w.frame2              -borderwidth 2 -relief ridge
label $w.frame2.residuelabel -text "Residue:" -width 10
entry $w.frame2.residue      -width 15
pack  $w.frame2 -side top    -anchor e
pack  $w.frame2.residuelabel $w.frame2.residue -side left

frame $w.frame3              -borderwidth 2 -relief ridge
label $w.frame3.atomlabel    -text "Atom:" -width 10
entry $w.frame3.atom         -width 15
pack  $w.frame3 -side top    -anchor e
pack  $w.frame3.atomlabel    $w.frame3.atom    -side left

label $w.atom2  -text "To atom(s):"
pack  $w.atom2  -side top -anchor w

frame $w.frame11              -borderwidth 2 -relief ridge
label $w.frame11.segmentlabel -text "Segment:" -width 10
entry $w.frame11.segment      -width 15
pack  $w.frame11 -side top    -anchor e
pack  $w.frame11.segmentlabel $w.frame11.segment -side left

frame $w.frame12              -borderwidth 2 -relief ridge
label $w.frame12.residuelabel -text "Residue:" -width 10
entry $w.frame12.residue      -width 15
pack  $w.frame12 -side top    -anchor e
pack  $w.frame12.residuelabel $w.frame12.residue -side left

frame $w.frame13              -borderwidth 2 -relief ridge
label $w.frame13.atomlabel    -text "Atom:" -width 10
entry $w.frame13.atom         -width 15
pack  $w.frame13 -side top    -anchor e
pack  $w.frame13.atomlabel    $w.frame13.atom    -side left

frame       $w.frame4
pack        $w.frame4 -side top -anchor w -fill both

frame       $w.frame4.left -borderwidth 2 -relief ridge
pack        $w.frame4.left -side top -anchor w

frame       $w.frame4.left.first
pack        $w.frame4.left.first -side top -anchor w

label       $w.frame4.left.first.label -text "Rcut:" -width 10
entry       $w.frame4.left.first.rcut -width 10
pack        $w.frame4.left.first.label $w.frame4.left.first.rcut -side left
# default
$w.frame4.left.first.rcut insert 0 "10."

frame       $w.frame4.left.second
pack        $w.frame4.left.second -side top -anchor w

label       $w.frame4.left.second.label -text "Nbins:" -width 10
entry       $w.frame4.left.second.nbin -width 10
pack        $w.frame4.left.second.label $w.frame4.left.second.nbin -side left
# default
$w.frame4.left.second.nbin insert 0 "101"

scan [show cell dimension] "%e %e %e" a b c

frame       $w.frame5 -borderwidth 2 -relief ridge
pack        $w.frame5 -side top -anchor w -pady 4

label       $w.frame5.label -text "Cell dimension:"
pack        $w.frame5.label -side top -anchor w

frame       $w.frame5.cell 
pack        $w.frame5.cell -side top -anchor w

label       $w.frame5.cell.alabel -text "a: "
entry       $w.frame5.cell.a      -width 11
pack        $w.frame5.cell.alabel $w.frame5.cell.a -side left -anchor w
# default
$w.frame5.cell.a inser 0 [format "%.5e" $a]

label       $w.frame5.cell.blabel -text "b: "
entry       $w.frame5.cell.b      -width 11
pack        $w.frame5.cell.blabel $w.frame5.cell.b -side left -anchor w
# default
$w.frame5.cell.b inser 0 [format "%.5e" $b]

label       $w.frame5.cell.clabel -text "c: "
entry       $w.frame5.cell.c      -width 11
pack        $w.frame5.cell.clabel $w.frame5.cell.c -side left -anchor w
# default
$w.frame5.cell.c inser 0 [format "%.5e" $c]

set NumObs [show rdf status]

  frame       $w.frame6 -borderwidth 2 -relief ridge
  pack        $w.frame6 -side top

if {$NumObs} {

  frame       $w.frame6.status
  pack        $w.frame6.status -side top -anchor w
  label       $w.frame6.status.label -text "Number of sets: "
  pack        $w.frame6.status.label -side left -anchor w
  entry       $w.frame6.status.value -width 10
  pack        $w.frame6.status.value -side left -anchor w
  $w.frame6.status.value insert 0 $NumObs

  set Average [show rdf average]
  if {!$Average} {
     button      $w.frame6.mean -text "Calc average"             \
                  -command {calculate rdfmean;                   \
                   set w .gomselectrdfatoms;                        \
                   $w.frame6.status.value delete 0 end;          \
                   $w.frame6.status.value insert 0 "** aver **"; \
                   $w.frame6.mean configure -state disabled}
     pack        $w.frame6.mean -side top
  } 
}

  frame       $w.frame7 -borderwidth 2 -relief ridge
  pack        $w.frame7 -side top

     if {[show trajectory frames] > 0} {
       button      $w.frame7.button -text "Calculate average RDF"            \
                    -background green -command "lulApplyRDFFrameLooping $w"
       pack        $w.frame7.button -side top -pady 3
	 }

frame       $w.delete -borderwidth 2 -relief ridge
pack        $w.delete -side top -pady 3

button      $w.delete.button -text "Delete RDF data"            \
             -command "calculate -rdf;destroy $w.frame6.status; \
                       destroy $w.frame6.mean"
pack        $w.delete.button -side top -pady 3

}

####################################################################
# PROC
proc lulApplyRDFFrameLooping { w } {

set Frames [show trajectory frames]
     if {[winfo exists $w.frame7.status]} {
       $w.frame7.status delete 0 end
       $w.frame7.status insert 0 "Starting ..."
	   } else {
       entry       $w.frame7.status -width 40
       pack        $w.frame7.status -side top
	 }

calculate -rdf

     set Segment1 [string trim [$w.frame1.segment get]]
     set Residue1 [string trim [$w.frame2.residue get]]
     set Atom1    [string trim [$w.frame3.atom    get]]

     set Segment2 [string trim [$w.frame11.segment get]]
     set Residue2 [string trim [$w.frame12.residue get]]
     set Atom2    [string trim [$w.frame13.atom    get]]

     set Rcut     [$w.frame4.left.first.rcut  get]
     set Nbin     [$w.frame4.left.second.nbin get]

     set a        [$w.frame5.cell.a get]
     set b        [$w.frame5.cell.b get]
     set c        [$w.frame5.cell.c get]

     eval define cell dimension $a $b $c

    set Limits [show trajectory display]
	set first [lindex $Limits 0]
	set last  [lindex $Limits 1]
	set step  [lindex $Limits 2]

	for {set i $first} {$i <= $last} {set i [expr $i + 1]} {
       import coord frame $i
       $w.frame7.status delete 0 end
       $w.frame7.status insert 0 "Processing frame # $i ..."
       update idletasks
       eval "calculate rdf \{$Segment1\} \{$Residue1\} \{$Atom1\} \{$Segment2\} \{$Residue2\} \{$Atom2\} $Rcut $Nbin"
    }

    calculate rdfmean
       $w.frame7.status delete 0 end
       $w.frame7.status insert 0 "Average for the '$Frames' frames calculated!"

    lulInvalidateDisplay
}

####################################################################
# PROC
proc lulApplyRDFselection { w } {

     global gomStructure

     set Segment1 [string trim [$w.frame1.segment get]]
     set Residue1 [string trim [$w.frame2.residue get]]
     set Atom1    [string trim [$w.frame3.atom    get]]

     set Segment2 [string trim [$w.frame11.segment get]]
     set Residue2 [string trim [$w.frame12.residue get]]
     set Atom2    [string trim [$w.frame13.atom    get]]

     set Rcut     [$w.frame4.left.first.rcut  get]
     set Nbin     [$w.frame4.left.second.nbin get]

     set a        [$w.frame5.cell.a get]
     set b        [$w.frame5.cell.b get]
     set c        [$w.frame5.cell.c get]

     eval define cell dimension $a $b $c

     eval "calculate rdf \{$Segment1\} \{$Residue1\} \{$Atom1\} \{$Segment2\} \{$Residue2\} \{$Atom2\} $Rcut $Nbin"

     set NumObs  [show rdf status]

     if  {$NumObs == 1} {

      if {[winfo exists $w.frame6.mean] && [winfo exists $w.frame6.status]} { 
       $w.frame6.mean configure -state normal
       $w.frame6.status.value delete 0 end
       $w.frame6.status.value insert 0 $NumObs
      } else {
       frame       $w.frame6.status
       pack        $w.frame6.status -side top -anchor w
       label       $w.frame6.status.label -text "Number of sets: "
       pack        $w.frame6.status.label -side left -anchor w
       entry       $w.frame6.status.value -width 10
       pack        $w.frame6.status.value -side left -anchor w
       $w.frame6.status.value insert 0 $NumObs

       set Average [show rdf average]

       if {!$Average} {
          button      $w.frame6.mean -text "Calc average"             \
                       -command {calculate rdfmean;                   \
                        set w .gomselectrdfatoms;                        \
                        $w.frame6.status.value delete 0 end;          \
                        $w.frame6.status.value insert 0 "** aver **"; \
                        $w.frame6.mean configure -state disabled}
          pack        $w.frame6.mean -side top
       }
      }
     }

     set value [show rdf status]

     $w.frame6.status.value delete 0 end
     $w.frame6.status.value insert 0 $value

     lulInvalidateDisplay
}

##################### export rdf file ############################
# PROC
proc lulExportRDFfile { w } {

   if {![show rdf observations]} {
       lulErrorDialog "ERROR - no RDF information is available"
       return
	   }

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"RDF data files"	{.dat .DAT}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .dat -title "Save RDF File"]]

    if [string compare $file ""] {
      export rdf "$file"

# change to current directory
      lulChangeDirectory "$file"

    }
}

####################################################################
# 
# edit material properties ...
#
proc lulEditMaterialProperties { } {

     global gomHelpFile
     global gomStructure
     global gomSwitch
     global gomControlFont
	 global gomScale
     global gomInstantDisplay

set w .gomeditmaterialproperties
catch {destroy $w}
toplevel $w 
wm title $w "Edit material properties"
wm iconname $w "Editmaterial properties"
wm geometry $w 
#280x450

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss   -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply     -font "$gomControlFont" \
        -command "lulApplyEditMaterialProperties $w"
button $w.buttons.help    -text Help      -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(editmatprop)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

frame        $w.matspec  -borderwidth 2 -relief ridge
label        $w.matspec.label -text "Material specular"  -font {Arial 12 bold} 
pack         $w.matspec       -side top -anchor w
pack         $w.matspec.label -side top -anchor w

#red
entry        $w.matspec.entryred -width 10
scale        $w.matspec.scalered -from 0.0 -to 1.0 -length 240 -variable gomMatSpecScaleRed \
             -orient horizontal -label "Material specular RED" -tickinterval 0.2               \
			 -showvalue true -fg red -digits 3 -resolution 0.01 \
			 -command "lulShowMaterialSpecularValueRed $w.matspec.entryred"
pack         $w.matspec.entryred -side top
pack         $w.matspec.scalered -side top -fill x

bind         $w.matspec.entryred <Return> "lulApplyEditMaterialProperties $w"

#green
entry        $w.matspec.entrygreen -width 10
scale        $w.matspec.scalegreen -from 0.0 -to 1.0 -length 240 -variable gomMatSpecScaleGreen \
             -orient horizontal -label "Material specular GREEN" -tickinterval 0.2              \
			 -showvalue true -fg green -digits 3 -resolution 0.01 \
			 -command "lulShowMaterialSpecularValueGreen $w.matspec.entrygreen"
pack         $w.matspec.entrygreen -side top
pack         $w.matspec.scalegreen -side top -fill x

bind         $w.matspec.entrygreen <Return> "lulApplyEditMaterialProperties $w"

#blue
entry        $w.matspec.entryblue -width 10
scale        $w.matspec.scaleblue -from 0.0 -to 1.0 -length 240 -variable gomMatSpecScaleBlue \
             -orient horizontal -label "Material specular BLUE" -tickinterval 0.2               \
			 -showvalue true -fg blue -digits 3 -resolution 0.01 \
			 -command "lulShowMaterialSpecularValueBlue $w.matspec.entryblue"
pack         $w.matspec.entryblue -side top
pack         $w.matspec.scaleblue -side top -fill x

bind         $w.matspec.entryblue <Return> "lulApplyEditMaterialProperties $w"

set          MatSpec [show material specular red]

             $w.matspec.scalered set $MatSpec
             $w.matspec.entryred delete 0 end
             $w.matspec.entryred insert 0 $MatSpec

set          MatSpec [show material specular green]

             $w.matspec.scalegreen set $MatSpec
             $w.matspec.entrygreen delete 0 end
             $w.matspec.entrygreen insert 0 $MatSpec

set          MatSpec [show material specular blue]

             $w.matspec.scaleblue set $MatSpec
             $w.matspec.entryblue delete 0 end
             $w.matspec.entryblue insert 0 $MatSpec

#shininess
entry        $w.matspec.entryshininess -width 10
scale        $w.matspec.scaleshininess -from 0.0 -to 128.0 -length 240 -variable gomMatScaleShininess   \
             -orient horizontal -label "Material Shininess" -tickinterval 20                   \
			 -showvalue true  -digits 4 -resolution 0.1 \
			 -command "lulShowMaterialShininess $w.matspec.entryshininess"
pack         $w.matspec.entryshininess -side top
pack         $w.matspec.scaleshininess -side top -fill x

bind         $w.matspec.entryshininess <Return> "lulApplyEditMaterialProperties $w"

set          MatShin [show material shininess]

             $w.matspec.scaleshininess set $MatShin
             $w.matspec.entryshininess delete 0 end
             $w.matspec.entryshininess insert 0 $MatSpec

frame        $w.control -borderwidth 2 -relief raised
label        $w.control.label -text "Continuous display: "
pack         $w.control -side top -anchor w -padx 4 -pady 4
pack         $w.control.label -side left -anchor w

radiobutton  $w.control.on  -text "On"  -value 1 \
              -command {set gomInstantDisplay 1} 
radiobutton  $w.control.off -text "Off" -value 0 \
              -command {set gomInstantDisplay 0} 
pack         $w.control.on  -side left -anchor w
pack         $w.control.off -side left -anchor w

# pick radiobutton as default
  if {$gomInstantDisplay} {
      $w.control.on  select
  } else {
      $w.control.off select
  }

}

#############################################################
# PROC
proc lulShowMaterialSpecularValueRed { w value } {

     global gomMatSpecScaleRed
     global gomInstantDisplay

     set MatSpec $value

     $w delete 0 end
     $w insert 0 [format "%.4f" $MatSpec]

# check for instant update
     if {$gomInstantDisplay > 0} {
	 
     set MatSpecRed   [$w   get]

	 if {$MatSpecRed != ""} {
	     eval define material specular red $MatSpecRed
     }
	 
	 display}

}

#############################################################
# PROC
proc lulShowMaterialSpecularValueGreen { w value } {

     global gomMatSpecScaleGreen
     global gomInstantDisplay

     set MatSpec $value

     $w delete 0 end
     $w insert 0 [format "%.4f" $MatSpec]

# check for instant update
     if {$gomInstantDisplay > 0} {

     set MatSpecGreen [$w get]

	 if {$MatSpecGreen != ""} {
	     eval define material specular green $MatSpecGreen
     }

     set MatSpecGreen  [$w get]

	 if {$MatSpecGreen != ""} {
	     eval define material specular green $MatSpecGreen
     }
	 
	 display}

}

#############################################################
# PROC
proc lulShowMaterialSpecularValueBlue { w value } {

     global gomMatSpecScaleBlue
     global gomInstantDisplay

     set MatSpec $value

     $w delete 0 end
     $w insert 0 [format "%.4f" $MatSpec]

# check for instant update
     if {$gomInstantDisplay > 0} {

     set MatSpecBlue  [$w get]

	 if {$MatSpecBlue != ""} {
	     eval define material specular blue $MatSpecBlue
     }
	 
	 display}

}

#############################################################
# PROC
proc lulShowMaterialShininess { w value } {

     global gomMatSpecShininess
     global gomInstantDisplay

     set MatShin $value

     $w delete 0 end
     $w insert 0 [format "%.4f" $MatShin]

# check for instant update
     if {$gomInstantDisplay > 0} {
	 
     set MatShininess   [$w   get]

	 if {$MatShininess != ""} {
	     eval define material shininess $MatShininess
     }

	 display}

}

#############################################################
# PROC
proc lulApplyEditMaterialProperties { w } {

     set MatSpecRed   [$w.matspec.entryred   get]

	 if {$MatSpecRed != ""} {
	     eval define material specular red $MatSpecRed
         set  MatSpec $MatSpecRed
         $w.matspec.scalered set $MatSpec
     }

     set MatSpecGreen [$w.matspec.entrygreen get]

	 if {$MatSpecGreen != ""} {
	     eval define material specular green $MatSpecGreen
         set  MatSpec $MatSpecGreen
         $w.matspec.scalegreen set $MatSpec
     }

     set MatSpecBlue  [$w.matspec.entryblue  get]

	 if {$MatSpecBlue != ""} {
	     eval define material specular blue $MatSpecBlue
         set  MatSpec $MatSpecBlue
         $w.matspec.scaleblue set $MatSpec
     }

     set MatShininess  [$w.matspec.entryshininess  get]

	 if {$MatShininess != ""} {
	     eval define material shininess $MatShininess
         set  MatSpec $MatShininess
         $w.matspec.scaleshininess set $MatSpec
     }

     lulInvalidateDisplay
}

####################################################################
# 
# rotation/translation control ...
#
proc lulRotaTransControl { } {

     global gomHelpFile
     global gomStructure
     global gomSwitch
     global gomControlFont
	 global gomScale
     global gomDefaultRotaStep
	 global gomDefaultTranStep

set w .gomrotatranscontrol
catch {destroy $w}
toplevel $w 
wm title $w "Rotation/translation control"
wm iconname $w "Rotation/translation control"
wm geometry $w 
#280x450

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss   -font "$gomControlFont" \
        -command "destroy $w"
#button $w.buttons.apply   -text Apply     -font "$gomControlFont" \
#        -command "lulApplyEditMaterialProperties $w"
button $w.buttons.help    -text Help      -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(rottranscontrol)"
pack   $w.buttons.dismiss $w.buttons.help -side left -expand 1

label  $w.label  -text "Rotation/Translation control" -font {Arial 12 bold}
pack   $w.label  -side top -anchor w -padx 4 -pady 4

frame  $w.frame0 -borderwidth 2 -relief ridge
pack   $w.frame0 -side top
label  $w.frame0.label -text "Rotation:"
pack   $w.frame0.label -side top

### ROTATIONS 
frame  $w.frame0.rotax
pack   $w.frame0.rotax -side top -anchor w
label  $w.frame0.rotax.label -text "X-axis \[degrees\]"
pack   $w.frame0.rotax.label -side top

button $w.frame0.rotax.buttonm -text "<-\]" -command \
        {if {$gomApplyType} {rotate display          \
		[expr -abs([string trim [.gomrotatranscontrol.frame0.rotax.step get]])]} \
		else {rotate selection                       \
		[expr -abs([string trim [.gomrotatranscontrol.frame0.rotax.step get]])]};display 1}
entry  $w.frame0.rotax.step    -width 10
button $w.frame0.rotax.buttonp -text "\[+>" -command \
        {if {$gomApplyType} {rotate display          \
		[expr  abs([string trim [.gomrotatranscontrol.frame0.rotax.step get]])]} \
		else {rotate selection                       \
		[expr  abs([string trim [.gomrotatranscontrol.frame0.rotax.step get]])]};display 1}
pack   $w.frame0.rotax.buttonm $w.frame0.rotax.step $w.frame0.rotax.buttonp \
        -side left
# default x rotation step
$w.frame0.rotax.step insert 0 $gomDefaultRotaStep

frame  $w.frame0.rotay
pack   $w.frame0.rotay -side top -anchor w
label  $w.frame0.rotay.label -text "Y-axis \[degrees\]"
pack   $w.frame0.rotay.label -side top

button $w.frame0.rotay.buttonm -text "<-\]" -command \
        {if {$gomApplyType} {rotate display 0.0      \
		[expr -abs([string trim [.gomrotatranscontrol.frame0.rotay.step get]])]} \
		else {rotate selection              0.0      \
		[expr -abs([string trim [.gomrotatranscontrol.frame0.rotay.step get]])]};display 1}
entry  $w.frame0.rotay.step    -width 10
button $w.frame0.rotay.buttonp -text "\[+>" -command \
        {if {$gomApplyType} {rotate display 0.0      \
		[expr  abs([string trim [.gomrotatranscontrol.frame0.rotay.step get]])]} \
		else {rotate selection              0.0      \
		[expr  abs([string trim [.gomrotatranscontrol.frame0.rotay.step get]])]};display 1}
pack   $w.frame0.rotay.buttonm $w.frame0.rotay.step $w.frame0.rotay.buttonp \
        -side left
# default y rotation step
$w.frame0.rotay.step insert 0 $gomDefaultRotaStep

frame  $w.frame0.rotaz
pack   $w.frame0.rotaz -side top -anchor w
label  $w.frame0.rotaz.label -text "Z-axis \[degrees\]"
pack   $w.frame0.rotaz.label -side top

button $w.frame0.rotaz.buttonm -text "<-\]" -command \
        {if {$gomApplyType} {rotate display 0.0 0.0  \
		[expr -abs([string trim [.gomrotatranscontrol.frame0.rotaz.step get]])]} \
		else {rotate selection              0.0 0.0  \
		[expr -abs([string trim [.gomrotatranscontrol.frame0.rotaz.step get]])]};display 1}
entry  $w.frame0.rotaz.step    -width 10
button $w.frame0.rotaz.buttonp -text "\[+>" -command \
        {if {$gomApplyType} {rotate display 0.0 0.0  \
		[expr  abs([string trim [.gomrotatranscontrol.frame0.rotaz.step get]])]} \
		else {rotate selection              0.0 0.0  \
		[expr  abs([string trim [.gomrotatranscontrol.frame0.rotaz.step get]])]};display 1}
pack   $w.frame0.rotaz.buttonm $w.frame0.rotaz.step $w.frame0.rotaz.buttonp \
        -side left
# default z rotation step
$w.frame0.rotaz.step insert 0 $gomDefaultRotaStep

### TRANSLATION

frame  $w.frame1 -borderwidth 2 -relief ridge
pack   $w.frame1 -side top -pady 5
label  $w.frame1.label -text "Translation:"
pack   $w.frame1.label -side top

frame  $w.frame1.translx
pack   $w.frame1.translx -side top -anchor w
label  $w.frame1.translx.label -text "X-axis \[A\]"
pack   $w.frame1.translx.label -side top

button $w.frame1.translx.buttonm -text "<-\]" -command \
        {if {$gomApplyType} {translate display         \
		[expr -abs([string trim [.gomrotatranscontrol.frame1.translx.step get]])]}  \
		else {translate selection                      \
		[expr -abs([string trim [.gomrotatranscontrol.frame1.translx.step get]])]};display}
entry  $w.frame1.translx.step    -width 10
button $w.frame1.translx.buttonp -text "\[+>" -command \
        {if {$gomApplyType} {translate display         \
		[expr  abs([string trim [.gomrotatranscontrol.frame1.translx.step get]])]}  \
		else {translate selection                      \
		[expr  abs([string trim [.gomrotatranscontrol.frame1.translx.step get]])]};display}
pack   $w.frame1.translx.buttonm $w.frame1.translx.step $w.frame1.translx.buttonp \
        -side left
# default x translation step
$w.frame1.translx.step insert 0 $gomDefaultTranStep

frame  $w.frame1.transly
pack   $w.frame1.transly -side top -anchor w
label  $w.frame1.transly.label -text "Y-axis \[A\]"
pack   $w.frame1.transly.label -side top

button $w.frame1.transly.buttonm -text "<-\]" -command \
        {if {$gomApplyType} {translate display 0.0     \
		[expr -abs([string trim [.gomrotatranscontrol.frame1.transly.step get]])]}  \
		else {translate selection 0.0                  \
		[expr -abs([string trim [.gomrotatranscontrol.frame1.transly.step get]])]};display}
entry  $w.frame1.transly.step    -width 10
button $w.frame1.transly.buttonp -text "\[+>" -command \
        {if {$gomApplyType} {translate display 0.0     \
		[expr  abs([string trim [.gomrotatranscontrol.frame1.transly.step get]])]}  \
		else {translate selection 0.0                  \
		[expr  abs([string trim [.gomrotatranscontrol.frame1.transly.step get]])]};display}
pack   $w.frame1.transly.buttonm $w.frame1.transly.step $w.frame1.transly.buttonp \
        -side left
# default y translation step
$w.frame1.transly.step insert 0 $gomDefaultTranStep

frame  $w.frame1.translz
pack   $w.frame1.translz -side top -anchor w
label  $w.frame1.translz.label -text "Z-axis \[A\]"
pack   $w.frame1.translz.label -side top

button $w.frame1.translz.buttonm -text "<-\]" -command \
        {if {$gomApplyType} {translate display 0.0 0.0 \
		[expr -abs([string trim [.gomrotatranscontrol.frame1.translz.step get]])]}  \
		else {translate selection 0.0 0.0              \
		[expr -abs([string trim [.gomrotatranscontrol.frame1.translz.step get]])]};display}
entry  $w.frame1.translz.step    -width 10
button $w.frame1.translz.buttonp -text "\[+>" -command \
        {if {$gomApplyType} {translate display 0.0 0.0 \
		[expr  abs([string trim [.gomrotatranscontrol.frame1.translz.step get]])]}  \
		else {translate selection 0.0 0.0              \
		[expr  abs([string trim [.gomrotatranscontrol.frame1.translz.step get]])]};display}
pack   $w.frame1.translz.buttonm $w.frame1.translz.step $w.frame1.translz.buttonp \
        -side left
# default z translation step
$w.frame1.translz.step insert 0 $gomDefaultTranStep

frame  $w.frame2 -borderwidth 2 -relief ridge
pack   $w.frame2 -side top -anchor w -padx 4 -pady 4

label  $w.frame2.label -text "Apply on:"
pack   $w.frame2.label -side left -anchor w

radiobutton  $w.frame2.display  -text "Display"    -value 1 \
              -variable gomApplyType \
              -command {define atom selection off;\
                        .dummy.right.collection.selstate.off select}
radiobutton  $w.frame2.selection -text "Selection" -value 0 \
              -variable gomApplyType \
              -command {define atom selection on;\
                        .dummy.right.collection.selstate.on select}
pack         $w.frame2.display  -side left -anchor w
pack         $w.frame2.selection -side left -anchor w

#default
if {[show selection] == "off"} {
$w.frame2.display select
} else {
$w.frame2.selection select
}
}

####################################################################
# 
# mean square displacement ...
#
proc lulMeanSquareDisplacement { } {

     global gomHelpFile
     global gomStructure
     global gomControlFont

# return if no molecular systems defined
     if {[show trajectory frames] < 1} {
          lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
	      return
	 }

set w .gommeansquaredisplacement
catch {destroy $w}
toplevel $w 
wm title $w "Mean square displacement"
wm iconname $w "Mean square displacement"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyMeanSquareDisplacement $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(meansqrdisp)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1 -padx 4

label  $w.label -text "Mean square displacement:" -font {Arial 12 bold}
pack   $w.label -side top -anchor w -pady 4

frame $w.frame1              -borderwidth 2 -relief ridge
label $w.frame1.segmentlabel -text "Segment:" -width 10
entry $w.frame1.segment      -width 15
pack  $w.frame1 -side top    -anchor w
pack  $w.frame1.segmentlabel $w.frame1.segment -side left

frame $w.frame2              -borderwidth 2 -relief ridge
label $w.frame2.residuelabel -text "Residue:" -width 10
entry $w.frame2.residue      -width 15
pack  $w.frame2 -side top    -anchor w
pack  $w.frame2.residuelabel $w.frame2.residue -side left

frame $w.frame3              -borderwidth 2 -relief ridge
label $w.frame3.atomlabel    -text "Atom:" -width 10
entry $w.frame3.atom         -width 15
pack  $w.frame3 -side top    -anchor w
pack  $w.frame3.atomlabel    $w.frame3.atom    -side left

}

####################################################################
# PROC
proc lulApplyMeanSquareDisplacement { w } {

     global gomStructure

     set Segment [string trim [$w.frame1.segment get]]
     set Residue [string trim [$w.frame2.residue get]]
     set Atom    [string trim [$w.frame3.atom    get]]

     eval "calculate msdisplacement atom \{$Segment\} \{$Residue\} \{$Atom\}"
}

##################### measure ############################
#
# measure
#
# PROC
proc lulMeasureGeometry {} {

    global gomHelpFile
    global gomSwitchTorsMon
    global gomAppendDistList
    global gomControlFont

# return if no molecular systems defined
    set  NumStruct [show molstructures]
    if {$NumStruct < 1} {
	lulErrorDialog {ERROR: no structure available. Read a structure first!}
	return
    }

    set w .gommeasuregeometry
    catch {destroy $w}
    toplevel $w 
    wm title $w "Measure geometry"
    wm iconname $w "Measure geometry"
    
    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
	-command "destroy $w"
    button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
	-command "lulDoMeasureGeometry $w"
    button $w.buttons.help    -text Help     -font "$gomControlFont" \
	-command \
	"htmlShowHelp $gomHelpFile(measuregeom)"
    pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1
    
    label  $w.label -text "Measure geometry:" -font {Arial 12 bold}
    pack   $w.label -side top -anchor w
    
    frame  $w.atoms
    pack   $w.atoms -side left -anchor w

    frame  $w.geometry -borderwidth 2 -relief ridge
    pack   $w.geometry -side right -anchor e -fill x

    for {set set 1} {$set <= 4} {incr set} {
	# Atoms
	frame  $w.atoms.set${set} -borderwidth 2 -relief ridge
	pack   $w.atoms.set${set} -side top
	label  $w.atoms.set${set}.label -text "Atom #${set}:"
	pack   $w.atoms.set${set}.label -side top -anchor w
	lulCreateAtomInputEntries $w.atoms.set${set} 7 "" "" "" "" 0

	# Geometry
	set setm1 [expr $set - 1]
	set setm2 [expr $set - 2]
	set setm3 [expr $set - 3]
	if { $setm3 >= 1 } {
	    # Torsion
	    label  $w.geometry.tlabel${setm3} -text "Torsion #${setm3}-#${setm2}-#${setm1}-#${set}:"
	    pack   $w.geometry.tlabel${setm3} -side top -anchor w
	    entry  $w.geometry.torsion${setm3}${setm2}${setm1}${set} -width 10
	    pack   $w.geometry.torsion${setm3}${setm2}${setm1}${set} -side top -anchor w -padx 12
	}
	if { $setm2 >= 1 } {
	    # Angle
	    label  $w.geometry.alabel${setm2} -text "Angle #${setm2}-#${setm1}-#${set}:"
	    pack   $w.geometry.alabel${setm2} -side top -anchor w
	    entry  $w.geometry.angle${setm2}${setm1}${set} -width 10
	    pack   $w.geometry.angle${setm2}${setm1}${set} -side top -anchor w -padx 12
	}
	if { $setm1 >= 1 } {
	    # Distance
	    label  $w.geometry.dlabel${setm1} -text "Distance #${setm1}-#${set}:"
	    pack   $w.geometry.dlabel${setm1} -side top -anchor w 
	    entry  $w.geometry.distance${setm1}${set} -width 10
	    pack   $w.geometry.distance${setm1}${set} -side top -anchor w -padx 12
	}
    }
}

#############################################################
# PROC
proc lulDoMeasureGeometry { w } {

    # Delete all texts.
    for {set set 2} {$set <= 4} {incr set} {
	set setm1 [expr $set - 1]
	set setm2 [expr $set - 2]
	set setm3 [expr $set - 3]
	# Torsion
	if { $setm3 >= 1 } {$w.geometry.torsion${setm3}${setm2}${setm1}${set} delete 0 end}
	# Angle
	if { $setm2 >= 1 } {$w.geometry.angle${setm2}${setm1}${set} delete 0 end}
	# Distance
	if { $setm1 >= 1 } {$w.geometry.distance${setm1}${set} delete 0 end}
    }

    # Measure new values.
    for {set set 1} {$set <= 4} {incr set} {
	# Atoms
	set Seg($set) [string trim [$w.atoms.set${set}.segment.input get]]
	set Res($set) [string trim [$w.atoms.set${set}.residue.input get]]
	set Atm($set) [string trim [$w.atoms.set${set}.atom.input    get]]
	if { "  " == "$Seg($set) $Res($set) $Atm($set)" } return
	set AtomCount($set) [lulCalculateIndexRangeListLength [
		find atoms join $Seg($set) $Res($set) $Atm($set) all]]

	# Geometry
	set setm1 [expr $set - 1]
	set setm2 [expr $set - 2]
	set setm3 [expr $set - 3]
	if { $setm3 >= 1 } {
	    # Torsion
	    regsub {\.0*$} [expr \
		$AtomCount($setm3).0 * \
		$AtomCount($setm2).0 * \
		$AtomCount($setm1).0 * \
		$AtomCount($set).0] {} count
	    if { "$count.0" < 100 || [tk_messageBox \
		-parent $w -type yesno -default no -icon question \
		-title "Complex calculation" \
		-message [join [list \
			"In all $count torsion angles will be calculated." \
			"This is actually quite a much." \
			"Do you really want to do this?"] "\n"]] == "yes" } {
		# Calculate
		$w.geometry.torsion${setm3}${setm2}${setm1}${set} insert 0 \
		    [calculate torsion \
			$Seg($setm3) $Res($setm3) $Atm($setm3) \
			$Seg($setm2) $Res($setm2) $Atm($setm2) \
			$Seg($setm1) $Res($setm1) $Atm($setm1) \
			$Seg($set)   $Res($set)   $Atm($set)]
	    }
	}
	if { $setm2 >= 1 } {
	    # Angle
	    regsub {\.0*$} [expr \
		$AtomCount($setm2).0 * \
		$AtomCount($setm1).0 * \
		$AtomCount($set).0] {} count
	    if { "$count.0" < 100 || [tk_messageBox \
		-parent $w -type yesno -default no -icon question \
		-title "Complex calculation" \
		-message [join [list \
			"In all $count angles will be calculated." \
			"This is actually quite a much." \
			"Do you really want to do this?"] "\n"]] == "yes" } {
		# Calculate
		$w.geometry.angle${setm2}${setm1}${set} insert 0 \
		    [calculate angle \
			$Seg($setm2) $Res($setm2) $Atm($setm2) \
			$Seg($setm1) $Res($setm1) $Atm($setm1) \
			$Seg($set)   $Res($set)   $Atm($set)]
	    }
	}
	if { $setm1 >= 1 } {
	    # Distance
	    regsub {\.0*$} [expr \
		$AtomCount($setm1).0 * \
		$AtomCount($set).0] {} count
	    if { "$count.0" < 100 || [tk_messageBox \
		-parent $w -type yesno -default no -icon question \
		-title "Complex calculation" \
		-message [join [list \
			"In all $count distances will be calculated." \
			"This is actually quite a much." \
			"Do you really want to do this?"] "\n"]] == "yes" } {
		# Calculate
		$w.geometry.distance${setm1}${set} insert 0 \
		    [calculate distance \
			$Seg($setm1) $Res($setm1) $Atm($setm1) \
			$Seg($set)   $Res($set)   $Atm($set)]
	    }
	}
    }
}

###################################################################
# PROC
proc lulPushToCommandStack {Command} {

     global gomCommandStack
     global gomCommandStackDeepMax
     global gomCommandStackPeek
     global gomCommandStackFill
     global gomCommandStackNext
     global gomCommandStackPrev
     global gomCommandStackCounter

     if {$Command == ""} return

     set gomCommandStack($gomCommandStackFill) $Command
     set gomCommandStackPeek $gomCommandStackFill
     set gomCommandStackFill [expr [incr gomCommandStackFill] % \
                                    $gomCommandStackDeepMax]

     if {$gomCommandStackCounter < $gomCommandStackDeepMax} {

         set gomCommandStackNext 0
         set gomCommandStackPrev [expr $gomCommandStackPeek - 1]
    } else {
         set gomCommandStackNext $gomCommandStackFill
         set gomCommandStackPrev [expr $gomCommandStackPeek - 1]
         if {$gomCommandStackPeek < 0} {
             set gomCommandStackPrev [expr $gomCommandStackDeepMax - 1]
         }
    }

    incr gomCommandStackCounter
} 
#############################################################
# PROC
proc lulGetNextFromStack {w} {

     global gomCommandStack
     global gomCommandStackDeepMax
     global gomCommandStackPeek
     global gomCommandStackFill
     global gomCommandStackNext
     global gomCommandStackPrev
     global gomCommandStackCounter

     if {!$gomCommandStackCounter} {
         bell
         return
         }

     $w delete 0 end
     $w insert 0 $gomCommandStack($gomCommandStackPeek)

     
     if {$gomCommandStackCounter < $gomCommandStackDeepMax} {
         incr gomCommandStackPeek
         if {$gomCommandStackPeek == $gomCommandStackFill} {
              set gomCommandStackPeek 0
              }
    } else {
         incr gomCommandStackPeek
         if {$gomCommandStackPeek == $gomCommandStackDeepMax} {
              set gomCommandStackPeek 0
         }
    }
}

#############################################################
# PROC
proc lulGetPrevFromStack {w} {

     global gomCommandStack
     global gomCommandStackDeepMax
     global gomCommandStackPeek
     global gomCommandStackNext
     global gomCommandStackPrev
     global gomCommandStackCounter
     global gomCommandStackFill

     if {!$gomCommandStackCounter} {
         bell
         return
         }

     $w delete 0 end
     $w insert 0 $gomCommandStack($gomCommandStackPeek)
     
     if {$gomCommandStackCounter < $gomCommandStackDeepMax} {
         set gomCommandStackPeek [expr $gomCommandStackPeek - 1]
         if {$gomCommandStackPeek < 0} {
              set gomCommandStackPeek [expr $gomCommandStackFill - 1]
              }
    } else {
         set gomCommandStackPeek [expr $gomCommandStackPeek - 1]
         if {$gomCommandStackPeek < 0} {
              set gomCommandStackPeek [expr $gomCommandStackDeepMax - 1]
         }
    }
}

##################### Plot colour scale ############################
# PROC
proc lulPlotColourScale {} {

     global gomControlFont
     global gomHelpFile
     global gomPlotCScale

set w .gomplotcolourscale
catch {destroy $w}
toplevel $w 
wm title $w "Colour Scale"
wm iconname $w "Colour Scale"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulDoPlotColourScale $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(plotcolourscale)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Plot colour scale:" -font {Arial 12 bold}
pack   $w.label -side top -anchor w -pady 3

frame  $w.frame1
pack   $w.frame1 -side top -pady 4
label  $w.frame1.min -text "Min value (blue):" -width 20
pack   $w.frame1.min -side left  -padx 4
entry  $w.frame1.minv -width 10
pack   $w.frame1.minv -side left

$w.frame1.minv insert 0 [show cscale mini]

frame  $w.frame2
pack   $w.frame2 -side top -pady 4
label  $w.frame2.max -text "Max value (red):"  -width 20
pack   $w.frame2.max -side left  -padx 4
entry  $w.frame2.maxv -width 10
pack   $w.frame2.maxv -side left

$w.frame2.maxv insert 0 [show cscale maxim]

frame  $w.frame3
pack   $w.frame3 -side top -pady 4
label  $w.frame3.bins -text "Number of bins:"  -width 20
pack   $w.frame3.bins -side left  -padx 4
entry  $w.frame3.binsv -width 10
pack   $w.frame3.binsv -side left

$w.frame3.binsv insert 0 [show cscale bins]

frame  $w.frame4 -borderwidth 2 -relief ridge
label  $w.frame4.label -text "Display state:"
pack   $w.frame4.label -side top
radiobutton $w.frame4.on  -text "On"  -command \
       {plot cscale [.gomplotcolourscale.frame3.binsv get]  \
                    [.gomplotcolourscale.frame1.minv  get]  \
                    [.gomplotcolourscale.frame2.maxv  get]  \
                    ;lulInvalidateDisplay} \
                    -variable gomPlotCScale -value 1

radiobutton $w.frame4.off -text "Off" -command \
       {plot -cscale;lulInvalidateDisplay} -variable gomPlotCScale -value 0

pack $w.frame4 -side top -anchor w 
pack $w.frame4.on $w.frame4.off -side top -anchor w

if {[show cscale state]} {
  $w.frame4.on select
} else {
  $w.frame4.off select
}

}
###############################################################
# PROC
proc lulDoPlotColourScale { w } {

     global gomPlotCScale

if {$gomPlotCScale} {
  plot cscale [$w.frame3.binsv get]  \
              [$w.frame1.minv  get]  \
              [$w.frame2.maxv  get] ;lulInvalidateDisplay
} else {
  plot -cscale
  }
}

##################### run xvibs ############################
# PROC
proc lulRunXvibs {} {

     global gomHelpFile
     global gomControlFont



set w .gomrunxvibs
catch {destroy $w}
toplevel $w 
wm title $w "Run Xvibs"
wm iconname $w "Run Xvibs"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulRunXvibsCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
       {htmlShowHelp $gomHelpFile(runxvibs)}
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Input file name:" -width 18
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoRunXvibsInput $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left -padx 2

frame  $w.frame1
pack   $w.frame1 -side top -anchor w

label  $w.frame1.label     -text "Vibration(s):"  -width 19
entry  $w.frame1.vibra     -width 10
label  $w.frame1.text      -text "\{all \| 1..(3*N-6)\}"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.vibra $w.frame1.text -side left 

frame  $w.frame3
pack   $w.frame3 -side top -anchor w

text   $w.frame3.text -relief sunken -bd 2 \
      -yscrollcommand "$w.frame3.scroll set" -setgrid 1 \
	  -height 30 -width 80

scrollbar $w.frame3.scroll -command "$w.frame3.text yview"
pack $w.frame3.scroll -side right -fill y
pack $w.frame3.text   -expand yes -fill both
$w.frame3.text delete 0.0 end

}

##################### lulDoRunXvibsInput ############################
#
# PROC
proc lulDoRunXvibsInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Gaussian/Gamess output"		{.out .log}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .log -title "Input Xvibs File"]]

    if [string compare $file ""] {
     $w.frame0.filename delete 0 end
	 $w.frame0.filename insert 0 "$file"

# change to current directory
     lulChangeDirectory "$file"

    }
}

########################################################################
# PROC
proc lulRunXvibsCommand { w } {

     global env
	 global gomEnv
	 global tcl_platform

set InputFile  [string trim [$w.frame0.filename get]]
set Vibra      [string trim [$w.frame1.vibra get]]
if {$Vibra == ""} {
    set Vibra "all"
    $w.frame1.vibra insert 0 "all"
    }

if {$tcl_platform(platform) == "windows"} {
  set Program    "java.exe -jar [file join $gomEnv(GOM_BIN) xvibs.jar] $InputFile $Vibra"
} else {
  set Program    "java     -jar [file join $gomEnv(GOM_BIN) xvibs.jar] $InputFile $Vibra"
}

puts "'$Program'"

set Command $Program

puts $Command

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

    $w.frame3.text insert $i.0 "Job is done!\n"

catch {close $f}

}


##################### demo ############################
# PROC
proc lulShowDemo {} {

     global gomHelpFile
     global gomControlFont
	 global gomInterruptDemo
     global env
	 global gomEnv


set w .gomdemo
catch {destroy $w}
toplevel $w 
wm title $w "Demo"
wm iconname $w "Demo"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.stop    -text Stop  -font "$gomControlFont" \
        -command {set gomInterruptDemo 1} \
		-state disabled
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"

pack   $w.buttons.stop $w.buttons.dismiss  -side left -expand 1

label  $w.label -text "Demos:" -font {Arial 12 bold}
pack   $w.label -side top -anchor w

frame  $w.frame -borderwidth 2 -relief raised
pack   $w.frame -side top

button $w.frame.button1 -text "Simple molecule display" -width 30 \
        -command "lulPlayDemo $w demo1.tcl {Simple molecule display}"
pack   $w.frame.button1 -side top -pady 3

button $w.frame.button2 -text "Molecular dynamics" -width 30 \
        -command "lulPlayDemo $w demo2.tcl {Molecular dynamics}"
pack   $w.frame.button2 -side top -pady 3

button $w.frame.button3 -text "Orbital and density" -width 30 \
        -command "lulPlayDemo $w demo3.tcl {Orbital and density}"
pack   $w.frame.button3 -side top -pady 3

button $w.frame.button4 -text "Object display" -width 30 \
        -command "lulPlayDemo $w demo4.tcl {Object display}"
pack   $w.frame.button4 -side top -pady 3

button $w.frame.button5 -text "Manipulations" -width 30 \
        -command "lulPlayDemo $w demo5.tcl {Manipulations}"
pack   $w.frame.button5 -side top -pady 3

button $w.frame.button6 -text "Clip plane" -width 30 \
        -command "lulPlayDemo $w demo6.tcl {Clip plane}"
pack   $w.frame.button6 -side top -pady 3
}

############################################################################
# PROC
proc lulPlayDemo { w file title } {
    global gomInterruptDemo
    global gomEnv
    
    set gomInterruptDemo 0

    $w.buttons.stop configure -state normal
    catch {source [file join $gomEnv(GOM_DEMO) $file]}
    $w.buttons.stop configure -state disable

    if {$gomInterruptDemo} {
	lulMessageDialog "'$title'\ndemo is stopped."
    } else {
	lulMessageDialog "'$title'\ndemo has finished."
    }
}

############################################################################
# PROC
proc lulPauseDemo { time } {
    global gomInterruptDemo

    update
    if {$gomInterruptDemo} {error "Demo interrupted by user"}

    pause $time

    update
    if {$gomInterruptDemo} {error "Demo interrupted by user"}
}

############################################################################
# PROC
proc lulGetVisibleForegroundColour { bg } {
    scan [lulColourHex2Float $bg] "%f %f %f" red green blue

    if { 0.299*$red + 0.587*$green + 0.114 > 0.5 } {
	return "#000000"
    } else {
	return "#ffffff"
    }
}

############################################################################
# PROC
proc lulChooseButtonColour { button var text } {

    set colour [tk_chooseColor -title $text \
	-parent $button \
	-initialcolor [expr $$var]]
		 
    if {"" != $colour} {
	set $var $colour
	$button configure \
	    -fg [lulGetVisibleForegroundColour $colour] \
	    -bg $colour
    }
}

############################################################################
# PROC
proc lulTimeSeriesManipulation {} {

     global gomHelpFile
     global gomTrajFileName
     global gomControlFont
     global gomTimeSeriesMan
     global gomTimeSeriesType
     global gomControlFont

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

     if {![show monitor distance timeseries] &&
         ![show monitor angle    timeseries] &&
         ![show monitor torsion  timeseries]} {
          lulErrorDialog {ERROR: no distance, angle and torsion time series defined!}
	      return
     }

set w .gomtimeseries
catch {destroy $w}
toplevel $w 
wm title $w "Manipulate Timeseries"
wm iconname $w "Manipulate Timeseries"

#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulDoManTimeSeries  $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command {htmlShowHelp $gomHelpFile(mantimeseries)}

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1
#

label  $w.label -text "Manipulate timeseries:" -font "$gomControlFont"
pack   $w.label -side top -anchor w

frame  $w.top 
pack   $w.top -side top -anchor w

# LEFT
frame  $w.top.left -borderwidth 2 -relief ridge
pack   $w.top.left -side left -anchor w

#define     DAVERAGE      1
#define     SQUARE        2
#define     COS           3
#define     COS2          4
#define     SQRT          5
#define     DINITIAL      6
#define     COPYNR        7
#define     ADDNR         8
#define     LOG           9
#define     EXP          10
#define     POWERREAL    11
#define     MULTREAL     12
#define     DIVIDEREAL   13
#define     SHIFTREAL    14
#define     DMIN         15
#define     ABS          16
#define     DIVFIRST     17
#define     DIVMAXIMUM   18
#define     PSPECTRUM    19
#define     ZERO         20

label  $w.top.left.label -text "Action:"
pack   $w.top.left.label -side top -anchor w

# default value
set         gomTimeSeriesMan  ""

radiobutton $w.top.left.daverage  -text "DAVERAGE"  -value "daverage" \
             -variable gomTimeSeriesMan
pack        $w.top.left.daverage  -side top -anchor w 
radiobutton $w.top.left.square    -text "SQUARE"    -value "square" \
             -variable gomTimeSeriesMan              
pack        $w.top.left.square    -side top -anchor w 
radiobutton $w.top.left.cos       -text "COS"       -value "cos"  \
             -variable gomTimeSeriesMan
pack        $w.top.left.cos       -side top -anchor w 
radiobutton $w.top.left.cos2      -text "COS2"      -value "cos2" \
             -variable gomTimeSeriesMan
pack        $w.top.left.cos2      -side top -anchor w 
radiobutton $w.top.left.sqrt      -text "SQRT"      -value "sqrt" \
             -variable gomTimeSeriesMan
pack        $w.top.left.sqrt      -side top -anchor w 
radiobutton $w.top.left.dinitial  -text "DINITIAL"  -value "dinitial" \
             -variable gomTimeSeriesMan
pack        $w.top.left.dinitial  -side top -anchor w 
#
frame       $w.top.left.fcopy
pack        $w.top.left.fcopy -side top -anchor w -fill x
radiobutton $w.top.left.fcopy.copy      -text "COPY"      -value "copy" \
             -variable gomTimeSeriesMan
pack        $w.top.left.fcopy.copy      -side left -anchor w
entry       $w.top.left.fcopy.number    -width 6
pack        $w.top.left.fcopy.number    -side right -anchor e
label       $w.top.left.fcopy.label     -text "Time series nr: "
pack        $w.top.left.fcopy.label     -side right -anchor e
#
frame       $w.top.left.fadd
pack        $w.top.left.fadd -side top -anchor w -fill x
radiobutton $w.top.left.fadd.add       -text "ADD"       -value "add" \
             -variable gomTimeSeriesMan
pack        $w.top.left.fadd.add       -side left -anchor w 
entry       $w.top.left.fadd.number    -width 6
pack        $w.top.left.fadd.number    -side right -anchor e
label       $w.top.left.fadd.label     -text "Time series nr: "
pack        $w.top.left.fadd.label     -side right -anchor e
#
radiobutton $w.top.left.log       -text "LOG"       -value "log" \
             -variable gomTimeSeriesMan
pack        $w.top.left.log       -side top -anchor w 
radiobutton $w.top.left.exp       -text "EXP"       -value "exp" \
             -variable gomTimeSeriesMan
pack        $w.top.left.exp       -side top -anchor w 
#
frame       $w.top.left.fpowerreal
pack        $w.top.left.fpowerreal -side top -anchor w -fill x
radiobutton $w.top.left.fpowerreal.powerreal -text "POWERREAL" -value "powerreal" \
             -variable gomTimeSeriesMan
pack        $w.top.left.fpowerreal.powerreal -side left -anchor w 
entry       $w.top.left.fpowerreal.number    -width 6
pack        $w.top.left.fpowerreal.number    -side right -anchor e
label       $w.top.left.fpowerreal.label     -text "Exponent: "
pack        $w.top.left.fpowerreal.label     -side right -anchor e
#
frame       $w.top.left.fmultreal
pack        $w.top.left.fmultreal -side top -anchor w -fill x
radiobutton $w.top.left.fmultreal.multreal  -text "MULTREAL"  -value "multreal" \
             -variable gomTimeSeriesMan
pack        $w.top.left.fmultreal.multreal  -side left -anchor w 
entry       $w.top.left.fmultreal.number    -width 6
pack        $w.top.left.fmultreal.number    -side right -anchor e
label       $w.top.left.fmultreal.label     -text "Constant: "
pack        $w.top.left.fmultreal.label     -side right -anchor e
#
frame       $w.top.left.fdividereal
pack        $w.top.left.fdividereal -side top -anchor w -fill x
radiobutton $w.top.left.fdividereal.dividereal  -text "DIVIDEREAL"  -value "dividereal" \
             -variable gomTimeSeriesMan
pack        $w.top.left.fdividereal.dividereal  -side left -anchor w 
entry       $w.top.left.fdividereal.number      -width 6
pack        $w.top.left.fdividereal.number      -side right -anchor e
label       $w.top.left.fdividereal.label       -text "Constant: "
pack        $w.top.left.fdividereal.label       -side right -anchor e
#
frame       $w.top.left.fshiftreal 
pack        $w.top.left.fshiftreal -side top -anchor w -fill x
radiobutton $w.top.left.fshiftreal.shiftreal -text "SHIFTREAL"  -value "shiftreal" \
             -variable gomTimeSeriesMan
pack        $w.top.left.fshiftreal.shiftreal -side left -anchor w 
entry       $w.top.left.fshiftreal.number    -width 6
pack        $w.top.left.fshiftreal.number    -side right -anchor e
label       $w.top.left.fshiftreal.label     -text "Constant: "
pack        $w.top.left.fshiftreal.label     -side right -anchor e
#
radiobutton $w.top.left.dmin      -text "DMIN"  -value "dmin" \
             -variable gomTimeSeriesMan
pack        $w.top.left.dmin      -side top -anchor w 
radiobutton $w.top.left.abs       -text "ABS"  -value "abs" \
             -variable gomTimeSeriesMan
pack        $w.top.left.abs       -side top -anchor w 
radiobutton $w.top.left.divfirst  -text "DIVFIRST"  -value "divfirst" \
             -variable gomTimeSeriesMan
pack        $w.top.left.divfirst  -side top -anchor w 
radiobutton $w.top.left.divmaximum  -text "DAVERAGE"  -value "divmaximum" \
             -variable gomTimeSeriesMan
pack        $w.top.left.divmaximum  -side top -anchor w 
radiobutton $w.top.left.pspectrum  -text "PSPECTRUM"  -value "pspectrum" \
             -variable gomTimeSeriesMan
pack        $w.top.left.pspectrum  -side top -anchor w 
radiobutton $w.top.left.zero       -text "ZERO"  -value "zero" \
             -variable gomTimeSeriesMan
pack        $w.top.left.zero       -side top -anchor w 

# RIGHT

set          gomTimeSeriesType   ""

frame  $w.top.right -borderwidth 2 -relief ridge
pack   $w.top.right -side left -anchor e

label  $w.top.right.label -text "Type of time series:    "
pack   $w.top.right.label -side top -anchor w

if {[show monitor distance timeseries]} {
  radiobutton $w.top.right.distance  -text "Distance"  -value "distance" \
               -variable gomTimeSeriesType
  pack        $w.top.right.distance  -side top -anchor w 
}

if {[show monitor angle    timeseries]} {
  radiobutton $w.top.right.angle     -text "Angle"     -value "angle" \
               -variable gomTimeSeriesType
  pack        $w.top.right.angle     -side top -anchor w 
}

if {[show monitor torsion  timeseries]} {
  radiobutton $w.top.right.torsion   -text "Torsion"   -value "torsion" \
               -variable gomTimeSeriesType
  pack        $w.top.right.torsion   -side top -anchor w 
}

frame  $w.which -borderwidth 2 -relief ridge
pack   $w.which -side top -anchor w -pady 5

label  $w.which.label -text "Destination time series nr: "
pack   $w.which.label -side left -anchor w
entry  $w.which.series -width 10
pack   $w.which.series -side left -anchor w

}
#############################################################
# PROC
proc lulDoManTimeSeries { w } {

     global gomTimeSeriesMan
     global gomTimeSeriesType

  if {$gomTimeSeriesMan == ""} {
      lulErrorDialog {ERROR: Time Series manipulation action not defined!}
      return
  }

  if {$gomTimeSeriesType == ""} {
      lulErrorDialog {ERROR: Time Series type not defined!}
      return
  }

  if {[$w.which.series get] == ""} {
      lulErrorDialog {ERROR: Time Series number not defined!}
      return
  }

  set Command "manipulate timeser $gomTimeSeriesType $gomTimeSeriesMan [$w.which.series get]"

# COPY
  if {$gomTimeSeriesMan == "copy"} {

  set Which [$w.top.left.fcopy.number get]

  if {$Which == ""} {
      lulErrorDialog {ERROR: Time Series number to copy from is not defined!}
      return
  }

     set Command "$Command $Which"
# ADD
  } elseif {$gomTimeSeriesMan == "add"} {

  set Which [$w.top.left.fadd.number get]

  if {$Which == ""} {
      lulErrorDialog {ERROR: Time Series number to add is not defined!}
      return
  }

     set Command "$Command $Which"
# POWERREAL
  } elseif {$gomTimeSeriesMan == "powerreal"} {

  set Which [$w.top.left.fpowerreal.number get]

  if {$Which == ""} {
      lulErrorDialog {ERROR: exponent is not defined!}
      return
  }

     set Command "$Command $Which"
# MULTREAL
  } elseif {$gomTimeSeriesMan == "multreal"} {

  set Which [$w.top.left.fmultreal.number get]

  if {$Which == ""} {
      lulErrorDialog {ERROR: Value to be multiplied is not defined!}
      return
  } 

     set Command "$Command $Which"
# DIVIDEREAL
  } elseif {$gomTimeSeriesMan == "dividereal"} {

  set Which [$w.top.left.fdividereal.number get]

  if {$Which == ""} {
      lulErrorDialog {ERROR: Value to divide with is not defined!}
      return
  }

     set Command "$Command $Which"
#SHIFTREAL
  } elseif {$gomTimeSeriesMan == "shiftreal"} {

  set Which [$w.top.left.fshiftreal.number get]

  if {$Which == ""} {
      lulErrorDialog {ERROR: Value to shift with is not defined!}
      return
  }

     set Command "$Command $Which"

  }

# Just do it...
  eval $Command

}

################## Plot Root Mean Square Deviation ######################
#
# define RMSD plot
#
proc lulDefineRMSDplot {} {

     global gomDefaultRMSDmagnifier
     global gomControlFont
     global gomHelpFile
     global gomControlFont
     global gomPlotRMSDaction

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomplotrmsd
catch {destroy $w}
toplevel $w 
wm title $w "Plot RMSD"
wm iconname $w "Plot RMSD"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulDoRMSDplot $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(plotrmsd)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Plot RMSD:" -font {Arial 12 bold}
pack   $w.label -side top -anchor w -pady 3

frame  $w.frame8 -borderwidth 2 -relief ridge
button $w.frame8.delete -text "Delete RMSD plot" \
       -command "plot -sphere;lulInvalidateDisplay"
pack   $w.frame8 -side bottom -anchor w -padx 10 -pady 5
pack   $w.frame8.delete -side left 

frame  $w.frame1 
label  $w.frame1.segment -text "Segment:" -width 10
entry  $w.frame1.segment_input -width 30 
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.segment $w.frame1.segment_input -side left

if { [$w.frame1.segment_input get] == "" } {$w.frame1.segment_input insert 0 "*"}

frame  $w.frame2 
label  $w.frame2.residue -text "Residue:" -width 10
entry  $w.frame2.residue_input -width 30  
pack   $w.frame2 -side top -anchor w
pack   $w.frame2.residue $w.frame2.residue_input -side left 

if { [$w.frame2.residue_input get] == ""} {$w.frame2.residue_input insert 0 "*"}

frame  $w.frame3 
label  $w.frame3.atom -text "Atom:" -width 10
entry  $w.frame3.atom_input -width 30  
pack   $w.frame3 -side top -anchor w
pack   $w.frame3.atom $w.frame3.atom_input -side left

if { [$w.frame3.atom_input get] == "" } {$w.frame3.atom_input insert 0 "*"}

frame  $w.frame4 -borderwidth 2 -relief ridge
label       $w.frame4.label -text "Apply fit:"
radiobutton $w.frame4.on  -text "On"  \
            -variable gomPlotRMSDaction -value fit
radiobutton $w.frame4.off -text "Off" \
            -variable gomPlotRMSDaction -value nofit

pack $w.frame4 -side left 
pack $w.frame4.label $w.frame4.on $w.frame4.off -side top -anchor w
$w.frame4.on select

frame  $w.frame6 -borderwidth 2 -relief ridge
label  $w.frame6.label -text "Amplifier factor: "
entry  $w.frame6.factor      -width 10

pack   $w.frame6 -side bottom -anchor e -pady 3
pack   $w.frame6.label $w.frame6.factor -side left

$w.frame6.factor insert 0 $gomDefaultRMSDmagnifier

}

########################################################################
proc lulDoRMSDplot { w } {

     global gomPlotRMSDaction

     set Segment [string trim [$w.frame1.segment_input get]]
	 if {$Segment == ""} {set Segment "*"}
	 set Residue [string trim [$w.frame2.residue_input get]]
	 if {$Residue == ""} {set Residue "*"}
     set Atom    [string trim [$w.frame3.atom_input    get]]
	 if {$Atom    == ""} {set Atom    "*"}
     set Factor  [string trim [$w.frame6.factor        get]]
	 if {$Factor  == ""} {set Factor  "1.0"}

     eval "p_rmsd \{$Segment\} \{$Residue\} \{$Atom\} $Factor $gomPlotRMSDaction"

     lulInvalidateDisplay
}

##################### display props ############################
#
# display properties
#
# PROC
proc lulChangeDisplayProperties {} {

     global gomHelpFile
     global gomControlFont

     package require gom::gui::Widgets
     namespace import -force ::gom::gui::Widgets::ColorButton

# return if no molecular systems defined
#     if {[show trajectory frames] < 1} {
#          lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
#	      return
#	 }

set w .gomchangedisplayprops
catch {destroy $w}
toplevel $w 
wm title $w "Display properties"
wm iconname $w "Display properties"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss    -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply      -font "$gomControlFont" \
        -command { } -state disabled
button $w.buttons.help    -text Help       -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(displayprops)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Display Properties:"
pack   $w.label -side top -anchor w

############################################################
frame  $w.bgcolor -borderwidth 2 -relief ridge
pack   $w.bgcolor -side top -fill x -pady 4
label  $w.bgcolor.label -text "Background colour (click to change): "
pack   $w.bgcolor.label -side left -padx 3
set gom::gui::wdata($w.bgcolor.apply) [eval lulColourFloat2Hex [show bgcolour]]
ColorButton $w.bgcolor.apply \
    -title "Choose background colour" \
    -variable gom::gui::wdata($w.bgcolor.apply) \
    -command "
	define bgcolour \$gom::gui::wdata($w.bgcolor.apply)
	lulInvalidateDisplay
    "
pack   $w.bgcolor.apply -side left -anchor center -padx 3

############################################################
frame  $w.bondtype -borderwidth 2 -relief ridge
pack   $w.bondtype -side top -fill x -pady 4
label  $w.bondtype.label -text "Bond display type: "
pack   $w.bondtype.label -side left -padx 3
radiobutton $w.bondtype.half   -text "Half"   -value 1 \
              -command {define bondstyle half;lulInvalidateDisplay} \
              -variable gomBondDisplayType
radiobutton $w.bondtype.smooth -text "Smooth" -value 0 \
              -command {define bondstyle smooth;lulInvalidateDisplay} \
              -variable gomBondDisplayType
pack        $w.bondtype.half $w.bondtype.smooth -side left

if {[show bondstyle] == "half"} {
    $w.bondtype.half select
} else {
    $w.bondtype.smooth select
}
############################################################
frame  $w.colortype -borderwidth 2 -relief ridge
pack   $w.colortype -side top -fill x -pady 4
label  $w.colortype.label -text "Color display type: "
pack   $w.colortype.label -side left -padx 3
radiobutton $w.colortype.color   -text "Color"   -value 1 \
              -command {define colortype color;lulInvalidateDisplay} \
              -variable gomColorDisplayType
radiobutton $w.colortype.gray -text "Grayscale" -value 0 \
              -command {define colortype gray;lulInvalidateDisplay} \
              -variable gomColorDisplayType
pack        $w.colortype.color $w.colortype.gray -side left

if {[show colortype] == "color"} {
    $w.colortype.color select
} else {
    $w.colortype.gray select
}
#############################################################

frame  $w.colormapping -borderwidth 2 -relief ridge
pack   $w.colormapping -side top -fill x -pady 4
label  $w.colormapping.label -text "Color mapping type: "
pack   $w.colormapping.label -side left -padx 3
radiobutton $w.colormapping.texture   -text "Texture"   -value 1 \
              -command {define colormapping texture;lulInvalidateDisplay} \
              -variable gomColorMappingType
radiobutton $w.colormapping.rainbow  -text "Rainbow"   -value 0 \
              -command {define colormapping rainbow;lulInvalidateDisplay} \
              -variable gomColorMappingType
pack        $w.colormapping.texture $w.colormapping.rainbow -side left

#
# return value for 'show colormapping':
# texture
# rainbow
#
if {[show colormapping] == "texture"} {
    $w.colormapping.texture select
} elseif {[show colormapping] == "rainbow"} {
    $w.colormapping.rainbow  select
} else {
  gomError "wrong type of color mapping type defined"
}
#############################################################

frame  $w.drawbuffer -borderwidth 2 -relief ridge
pack   $w.drawbuffer -side top -fill x -pady 4
label  $w.drawbuffer.label -text "Drawing buffer: "
pack   $w.drawbuffer.label -side left -padx 3
radiobutton $w.drawbuffer.front   -text "Front"   -value 1 \
              -command {define drawbuffer front} -variable gomActiveBuffer
radiobutton $w.drawbuffer.back -text "Back"       -value 0 \
              -command {define drawbuffer back}  -variable gomActiveBuffer
pack        $w.drawbuffer.front $w.drawbuffer.back -side left

if {[show drawbuffer] == "front"} {
    $w.drawbuffer.front select
} else {
    $w.drawbuffer.back  select
}
#############################################################

frame  $w.frame1 -borderwidth 2 -relief ridge
pack   $w.frame1 -side top -fill x

label  $w.frame1.label0 -text "Edit cutplane distance to near plane:"
pack   $w.frame1.label0 -side top -anchor w

frame  $w.frame1.nearplane 
pack   $w.frame1.nearplane -side top -pady 3
button $w.frame1.nearplane.dec -text "<-\]" -command "lulDecreaseNearplane $w"
pack   $w.frame1.nearplane.dec -side left -padx 4 -pady 3
entry  $w.frame1.nearplane.value -width 15
pack   $w.frame1.nearplane.value -side left -padx 4 -pady 3
button $w.frame1.nearplane.inc -text "\[+>" -command "lulIncreaseNearplane $w"
pack   $w.frame1.nearplane.inc -side left -padx 4 -pady 3

set Value [show nearplane value]
$w.frame1.nearplane.value insert end $Value

frame  $w.frame1.step 
pack   $w.frame1.step -side top -anchor w
label  $w.frame1.step.label  -text "Step size:" -width 13
pack   $w.frame1.step.label  -side left 
entry  $w.frame1.step.value  -width 10
pack   $w.frame1.step.value  -side left

set StepValue [show nearplane step]
$w.frame1.step.value insert end $StepValue

frame  $w.frame1.update 
pack   $w.frame1.update -side top -anchor w -pady 4 -fill x
checkbutton $w.frame1.update.status -text "Instant update"         \
            -onvalue 1 -offvalue 0 -variable gomNearPlaneUpdate
pack   $w.frame1.update.status -side left
button $w.frame1.update.apply -text "Apply" -command "lulDoEditDisplay $w"
pack   $w.frame1.update.apply -side right -padx 4

$w.frame1.update.status deselect

#############################################################

frame  $w.projection -borderwidth 2 -relief ridge
pack   $w.projection -side top -fill x -pady 4
frame  $w.projection.top
pack   $w.projection.top -side top -anchor w

set DisplayOn     {.gomchangedisplayprops.projection.bottom.value configure  -state normal ;\
                    .gomchangedisplayprops.projection.bottom.button configure -state normal}
set DisplayOff    {.gomchangedisplayprops.projection.bottom.value configure  -state disabled;\
                    .gomchangedisplayprops.projection.bottom.button configure -state disabled}

label  $w.projection.top.label -text "Projection type: "
pack   $w.projection.top.label -side left -padx 3
radiobutton $w.projection.top.perspective   -text "Perspective"   -value 1 \
              -command "define projection perspective;lulInvalidateDisplay;$DisplayOn"  \
              -variable gomProjectionType
radiobutton $w.projection.top.orthographic  -text "Orthographic"   -value 0 \
              -command "define projection orthographic;lulInvalidateDisplay;$DisplayOff" \
              -variable gomProjectionType
pack        $w.projection.top.perspective $w.projection.top.orthographic -side left

frame  $w.projection.bottom
pack   $w.projection.bottom -side top -anchor w
label  $w.projection.bottom.label -text "Viewing angle (degrees): "
pack   $w.projection.bottom.label -side left -padx 3
entry  $w.projection.bottom.value  -width 10
pack   $w.projection.bottom.value  -side left
button $w.projection.bottom.button -text "Apply" -command \
       {define viewangle [.gomchangedisplayprops.projection.bottom.value get];lulInvalidateDisplay}
pack   $w.projection.bottom.button -side left -padx 3

set VAngle [show viewangle]
$w.projection.bottom.value delete 0 end
$w.projection.bottom.value insert 0 $VAngle

#
# return value for 'show projection':
# 0: perspective
# 1: orthographic
#
if {[show projection] == 0} {
    $w.projection.top.perspective select
    $w.projection.bottom.value configure  -state normal
	$w.projection.bottom.button configure -state normal
} elseif {[show projection] == 1} {
    $w.projection.top.orthographic  select
    $w.projection.bottom.value  configure -state disabled
	$w.projection.bottom.button configure -state disabled
} else {
  gomError "wrong type of Projection Type defined"
}

############################################################
frame  $w.redisplaytype -borderwidth 2 -relief ridge
pack   $w.redisplaytype -side top -fill x -pady 4
label  $w.redisplaytype.label -text "Redisplay type: "
pack   $w.redisplaytype.label -side left -padx 3
radiobutton $w.redisplaytype.fast   -text "Display stick (fast)"   -value 0 \
              -variable gomReDisplayTypeSwitch \
              -command {define redisplay fast;lulInvalidateDisplay}
radiobutton $w.redisplaytype.slow -text "Display all (slow)" -value 1 \
              -variable gomReDisplayTypeSwitch \
              -command {define redisplay slow;lulInvalidateDisplay}
pack        $w.redisplaytype.fast $w.redisplaytype.slow -side left

if {[show redisplay] == "fast"} {
    $w.redisplaytype.fast select
} else {
    $w.redisplaytype.slow select
}
#############################################################

frame  $w.scaling -borderwidth 2 -relief ridge
pack   $w.scaling -side top -fill x -pady 4
label  $w.scaling.label -text "Scaling type: "
pack   $w.scaling.label -side left -padx 3
radiobutton $w.scaling.global   -text "Global"   -value 1 \
              -variable gomScaleDisplayTypeSwitch \
              -command {define scaling global;lulInvalidateDisplay}
radiobutton $w.scaling.individual  -text "Local"   -value 0 \
              -variable gomScaleDisplayTypeSwitch \
              -command {define scaling local;lulInvalidateDisplay}
pack        $w.scaling.global $w.scaling.individual -side left

#
if {[show scaling type] == "global"} {
    $w.scaling.global select
} elseif {[show scaling type] == "local"} {
    $w.scaling.individual  select
} else {
  gomError "wrong type of scaling type defined"
}

############################################################
frame  $w.stickwidth -borderwidth 2 -relief ridge
pack   $w.stickwidth -side top -fill x -pady 4
label  $w.stickwidth.label -text "Stick display width: "
pack   $w.stickwidth.label -side left -padx 3
entry  $w.stickwidth.entry -width 5
pack   $w.stickwidth.entry -side left -padx 3
$w.stickwidth.entry insert 0 [show mlinewidth]
button $w.stickwidth.apply -text "Apply" \
        -command {define mlinewidth [.gomchangedisplayprops.stickwidth.entry get];lulInvalidateDisplay}
pack   $w.stickwidth.apply -side left -anchor e -padx 3

############################################################
frame  $w.wupdate -borderwidth 2 -relief ridge
pack   $w.wupdate -side top -fill x -pady 4
label  $w.wupdate.label -text "Window update: "
pack   $w.wupdate.label -side left -padx 3
radiobutton $w.wupdate.auto   -text "Automatic"  -value 1 -variable lulWiup \
              -command {define window update automatic;lulInvalidateDisplay}
radiobutton $w.wupdate.manu   -text "Manual"     -value 0 -variable lulWiup \
              -command {define window update manual;lulInvalidateDisplay}
pack        $w.wupdate.auto $w.wupdate.manu -side left

if {[show window update] == "automatic"} {
    $w.wupdate.auto select
} elseif {[show window update] == "manual"} {
    $w.wupdate.manu select
} else {
  gomError "wrong windowing update returned (automatic/manual)"
}
############################################################
frame  $w.wtype -borderwidth 2 -relief ridge
pack   $w.wtype -side top -fill x -pady 4
label  $w.wtype.label -text "Windowing type: "
pack   $w.wtype.label -side left -padx 3
radiobutton $w.wtype.single   -text "Single"   -value 0 -variable lulWity \
              -command {define window single;lulInvalidateDisplay}
radiobutton $w.wtype.multi -text "Multi"       -value 1 -variable lulWity \
              -command {define window multi;lulInvalidateDisplay}
pack        $w.wtype.single $w.wtype.multi -side left

if {[show window style] == "single"} {
    $w.wtype.single select
} elseif {[show window style] == "multi"} {
    $w.wtype.multi select
} else {
  gomError "wrong windowing type returned (single/multi)"
}
#############################################################

}

##################### display list props ############################
#
# display properties
#
# PROC
proc lulChangeDisplayListProperties {} {

    global gomDisplayListState
    global gomHelpFile
    global gomControlFont

    set w .gomchangedisplaylistprops
    catch {destroy $w}
    toplevel $w
    wm title $w "Display list properties"
    wm iconname $w "Display list properties"

    foreach name [array names gomDisplayListState] {
	unset gomDisplayListState($name)
    }

    # Buttons
    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
	-command "destroy $w"
    button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
	-command {} -state disabled
    button $w.buttons.help    -text Help    -font "$gomControlFont" \
	-command \
	"htmlShowHelp $gomHelpFile(displayprops)"
    pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

    label  $w.label -text "Display List Properties:"
    pack   $w.label -side top -anchor w

    # Main config
    frame  $w.state -borderwidth 2 -relief ridge
    pack   $w.state -side top -fill x -pady 4
    label  $w.state.label -text "Display list state:"
    pack   $w.state.label -side left -padx 3
    set gomDisplayListState(main) [show displaylists]
    radiobutton $w.state.on  -text "On"  -value 1 \
	-command {define displaylists on} \
	-variable gomDisplayListState(main)
    radiobutton $w.state.off -text "Off" -value 0 \
	-command {define displaylists off} \
	-variable gomDisplayListState(main)
    pack   $w.state.off $w.state.on -side right

    # Default status
    frame  $w.default -borderwidth 2 -relief ridge
    pack   $w.default -side top -fill x -pady 4
    label  $w.default.label -text "Default for new type of objects:"
    pack   $w.default.label -side left -padx 3
    set gomDisplayListState(default) [show displaylists default]
    radiobutton $w.default.yes  -text "Yes"  -value 1 \
	-command {define displaylists default on} \
	-variable gomDisplayListState(default)
    radiobutton $w.default.no -text "No" -value 0 \
	-command {define displaylists default off} \
	-variable gomDisplayListState(default)
    pack   $w.default.no $w.default.yes -side right

    frame  $w.types -borderwidth 2 -relief ridge
    pack   $w.types -side top -fill x -pady 4

    # Type specific status
    set count 0
    foreach type [show displaylists types] {
	incr count
	regsub -all {[[:punct:]]} $type { } typename

	frame  $w.types.frame$count
	pack   $w.types.frame$count -side top -fill x
	label  $w.types.frame$count.label -text "Use display lists for $typename:"
	pack   $w.types.frame$count.label -side left -padx 3
	set gomDisplayListState($type) [show displaylists $type]

	radiobutton $w.types.frame$count.yes  -text "Yes"  -value 1 \
	    -command "define displaylists $type on" \
	    -variable gomDisplayListState($type)
	radiobutton $w.types.frame$count.no -text "No" -value 0 \
	    -command "define displaylists $type off" \
	    -variable gomDisplayListState($type)
	pack   $w.types.frame$count.no $w.types.frame$count.yes -side right
    }
}

#############################################################
# Immediate display ?
# PROC
proc lulInvalidateDisplay {} {

     global gomInstantDisplay

# check for instant update
     if {$gomInstantDisplay > 0} {display}
}

##################### import atom charges ############################

# ICON8
proc lulImportAtomChargesICON8 { w } {

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"ICON8 output file"		{*}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .out -title "Read Atom Charges ICON8"]]

    if {$file == ""} return 

    lulImportChargeFromICON8output "$file"

# change to current directory
    lulChangeDirectory "$file"

}

# MOPAC6
proc lulImportAtomChargesMOPAC6 { w } {

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"MOPAC6 output file"		{*}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .out -title "Read Atom Charges MOPAC6"]]

    if {$file == ""} return 

    lulImportChargeFromMOPAC6output "$file"

# change to current directory
    lulChangeDirectory "$file"

}

# GAUSSIAN
proc lulImportAtomChargesGAUSSIAN { w  which} {

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"GAUSSIAN output file"		{.log .out .LOG .OUT}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .out -title "Read Atom Charges GAUSSIAN"]]

    if {$file == ""} return 

    lulImportChargesFromGAUSSIANoutput "$file" $which

# change to current directory
    lulChangeDirectory "$file"

}

# USER
##########################################################################
# PROC
#
# The atoms have to be in the same order as the atoms were
# defined when they were read into gOpenMol.
#
# Generic text file input by Jim Kress 04November02
# Text file must have natoms lines.  Each line consists
# of the atom number (integer) and its associated charge (float)
#
proc lulImportChargeFromUSERoutput {file} {

    if {$file == ""} {
      gomError "USER output file name is missing"
      return
    }

    puts "Importing USER atom charge file '$file'"

    set f [open "$file" r]
    set i 0
    while {![eof $f]} {
    gets $f Text
    set Text [string trim $Text]
    if {($Text == "") } break
    incr i
      if {[llength $Text] > 1} {
         set Charge($i) [lindex $Text 1]
      } else {
         gomPrint "Expecting two values (atom symbol and partial charge), now only one '$Text'! Using that value"
         set Charge($i) [lindex $Text 0]
      }
    }
    set NumC $i
    set MolS [show molstructures]
    gomPrint "Found '$NumC' atom partial charges"

    for   {set i 1} {$i <= $MolS} {incr i}               {
# if number of atoms do not match jump over
     if {[show numatoms $i] != $NumC} continue
# if structure not selected jump over
     if {![show status $i]} continue

     gomPrint "Will now put '$NumC' atom partial charges into structure
'$i'"

     for  {set j 1} {$j <=  [show numatoms $i]} {incr j} {

        define atom charge  $Charge($j) $j $i
        unset Charge($j)

     }
    }

}


########################## end of import charges #########################


# USER
proc lulImportAtomChargesUSER { w } {

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"USER output file"		{*}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .out -title "Read USER Atom Charges"]]

    if {$file == ""} return 

    lulImportChargeFromUSERoutput "$file"

# change to current directory
    lulChangeDirectory "$file"

}

##########################################################################
# PROC
# Look for the string "NET CHG." in the output file.
#
# Observe that no checking for the atoms are made!!!!
# The atoms have to be in the same order as the atoms were
# defined when they were read into gOpenMol.
#
proc lulImportChargeFromICON8output {file} {

    if {$file == ""} {
      gomError "ICON8 output file name is missing"
      return
    } 

    puts "Importing ICON8 atom charge file '$file'"

    set f [open $file r]

    set i 0

    while {![eof $f]} { 
      gets $f Text
      update idletasks
# look for "NET CHG."
      if {[string match "*NET CHG.*" $Text]} {
# two text lines first
          gets $f Text
          gets $f Text
          while {![eof $f]} { 
           gets $f Text
           set Text [string trim $Text]
           if {$Text == ""} break
           incr i
           set Charge($i) [lindex $Text 2]
          }

      }
    }

    catch {close $f}

    if {$i < 1} {
      gomError "no atom partial charges found in ICON8 output file"
      return
    }

    set NumC $i
    set MolS [show molstructures]
    gomPrint "Found '$NumC' atom partial charges"

    for   {set i 1} {$i <= $MolS} {incr i}               {
# if number of atoms do not match jump over
     if {[show numatoms $i] != $NumC} continue
# if structure not selected jump over
     if {![show status $i]} continue

     gomPrint "Will now put '$NumC' atom partial charges into structure '$i'"

     for  {set j 1} {$j <=  [show numatoms $i]} {incr j} {

        define atom charge $Charge($j) $j $i
        unset Charge($j)

     }
    }

}

# MOPAC6
##########################################################################
# PROC
# Look for the string "NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS"
# in the output file.
#
# Observe that no checking for the atoms are made!!!!
# The atoms have to be in the same order as the atoms were
# defined when they were read into gOpenMol.
#
proc lulImportChargeFromMOPAC6output {file} {

    if {$file == ""} {
      gomError "MOPAC6 output file name is missing"
      return
    } 

    puts "Importing MOPAC6 atom charge file '$file'"

    set f [open $file r]

    set i 0

    while {![eof $f]} { 
      gets $f Text
      update idletasks
# look for "NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS"
      if {[string match "*NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS*" $Text]} {
# two text lines first
          gets $f Text
          gets $f Text
          while {![eof $f]} { 
           gets $f Text
           set Text [string trim $Text]
           if {($Text == "") || ([string match "*DIPOLE*" $Text])} break
           incr i
           set Charge($i) [lindex $Text 2]
          }

      }
    }

    catch {close $f}

    if {$i < 1} {
      gomError "no atom partial charges found in MOPAC6 output file"
      return
    }

    set NumC $i
    set MolS [show molstructures]
    gomPrint "Found '$NumC' atom partial charges"

    for   {set i 1} {$i <= $MolS} {incr i}               {
# if number of atoms do not match jump over
     if {[show numatoms $i] != $NumC} continue
# if structure not selected jump over
     if {![show status $i]} continue

     gomPrint "Will now put '$NumC' atom partial charges into structure '$i'"

     for  {set j 1} {$j <=  [show numatoms $i]} {incr j} {

        define atom charge $Charge($j) $j $i
        unset Charge($j)

     }
    }

}

###
# GAUSSIANXX Natural Orbital/Merz-Kollmann/Chelpg Charges
##########################################################################
# PROC
# which == 1 (Natural Orbital)
# Look for the string "Summary of Natural Population Analysis:"
# which == 2 (Merz-Kollmann)
# Look for the string "Fitting point charges to eletrostatic potential"
# which == 3 (Chelpg)
# Look for the string "Fitting point charges to eletrostatic potential"
# in the output file.
#
# Observe that no checking for the atoms are made!!!!
# The atoms have to be in the same order as the atoms were
# defined when they were read into gOpenMol.
#
proc lulImportChargesFromGAUSSIANoutput {file which} {

    if {$file == ""} {
      gomError "GAUSSIAN output file name is missing"
      return
    } 

    puts "Importing GAUSSIAN atom charge file '$file'"

    set f [open "$file" r]

    set i 0

  if {$which == 1} {
    gomPrint "Reading charges for a gaussian Natural Orbital calculation"
    set i 0
    while {![eof $f]} { 
      gets $f Text
      update idletasks
# look for "Summary of Natural Population Analysis:"
      if {[string match "*Summary of Natural Population Analysis:*" $Text]} {
# two text lines first
          gets $f Text
          gets $f Text
          gets $f Text
          gets $f Text
          gets $f Text
          while {![eof $f]} { 
           gets $f Text
           set Text [string trim $Text]
           if {($Text == "") || ([string match "*==========*" $Text])} break
           incr i
           set Charge($i) [lindex $Text 2]
          }
      }
    }
  } elseif {$which == 2 || $which == 3} {

    if {$which == 2} {
     gomPrint "Reading charges for a Gaussian Merz-Kollmann calculation"
    } else {
     gomPrint "Reading charges for a Gaussian Chelpg calculation"
    }

    set i 0
    while {![eof $f]} { 
      gets $f Text
      update idletasks
# look for "Fitting point charges to eletrostatic potential"
      if {[string match "*Fitting point charges to eletrostatic potential*" $Text]} {
# text lines first
          gets $f Text
          gets $f Text
          gets $f Text
          while {![eof $f]} { 
           gets $f Text
           set Text [string trim $Text]
           if {($Text == "") || ([string match "*---------------------*" $Text])} break
           incr i
           set Charge($i) [lindex $Text 2]
          }
      }
    }
  } else {
    gomError "Unknown Gaussian output file charge option"
  }

    catch {close $f}

    if {$i < 1} {
      gomError "no atom partial charges found in GAUSSIAN output file"
      return
    }

    set NumC $i
    set MolS [show molstructures]
    gomPrint "Found '$NumC' atom partial charges"

    for   {set i 1} {$i <= $MolS} {incr i}               {
# if number of atoms do not match jump over
     if {[show numatoms $i] != $NumC} continue
# if structure not selected jump over
     if {![show status $i]} continue

     gomPrint "Will now put '$NumC' atom partial charges into structure '$i'"

     for  {set j 1} {$j <=  [show numatoms $i]} {incr j} {

        define atom charge $Charge($j) $j $i
        unset Charge($j)

     }
    }

}

# USER
##########################################################################
# PROC
#
# The atoms have to be in the same order as the atoms were
# defined when they were read into gOpenMol.
#
proc lulImportChargeFromUSERoutputDummy {file} {

    if {$file == ""} {
      gomError "USER output file name is missing"
      return
    } 

    puts "Importing USER atom charge file '$file'"

    gomError "USER code is missing!"
}


########################## end of import charges #########################

##################### calculate connectivity ############################
#
# calculate connectivity
#
# PROC
proc lulCalculateConnectivity {} {

     global gomHelpFile
     global gomControlFont

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomcalcconn
catch {destroy $w}
toplevel $w 
wm title $w "Calculate atom connectivity"
wm iconname $w "Calc atom conn"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulDoCalculateAtomConnection $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(calcatomconn)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Calculate Atom Connectivity:" -font "$gomControlFont"
pack   $w.label -side top -anchor w

frame  $w.atoms -borderwidth 2 -relief ridge
pack   $w.atoms -side top -anchor w
label  $w.atoms.label -text "Atom(s):"
pack   $w.atoms.label -side top -anchor w
lulCreateAtomInputEntries $w.atoms 15 "*" "*" "*" "" 1

frame       $w.options -borderwidth 2 -relief ridge
pack        $w.options -side top -anchor w -fill x -pady 5

label       $w.options.label  -text "Atom search window: "
pack        $w.options.label  -side left -anchor w -padx 3
entry       $w.options.window -width 10
pack        $w.options.window -side left -anchor w -padx 3
button      $w.options.apply  -text "Apply" \
             -command "lulDoAtomSearchWindow $w.options.window"
pack        $w.options.apply  -side left -anchor w

button      $w.reset  -text "Reset connectivity" \
             -command "reset connectivity;lulInvalidateDisplay"
pack        $w.reset -side bottom

$w.options.window insert 0 [show atom window]

}

###############################################################
# PROC
proc lulDoCalculateAtomConnection { w } {

   set Seg1 [string trim [$w.atoms.segment.input get]]
   set Res1 [string trim [$w.atoms.residue.input get]]
   set Atm1 [string trim [$w.atoms.atom.input get]]


   calculate connect $Seg1 $Res1 $Atm1

   lulInvalidateDisplay
}

###############################################################
# PROC
proc lulDoAtomSearchWindow { w } {

  set Window [string trim [$w get]]

  if {$Window == ""} {
     gomError "window value has to be > 0\n Will make it = 40"
     $w insert 0 "40"
     return
  } 

  define atom window $Window

  lulInvalidateDisplay
}
#################### end of calculate connection #################

####################################################################
# 
# edit light properties ...
#
proc lulEditLightProperties { } {

     global gomHelpFile
     global gomStructure
     global gomSwitch
     global gomControlFont
	 global gomScale
	 global gomInstantDisplay

set w .gomeditlightproperties
catch {destroy $w}
toplevel $w 
wm title $w "Edit light properties"
wm iconname $w "Editlight properties"
wm geometry $w 
#280x450

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss   -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply     -font "$gomControlFont" \
        -command "lulApplyEditLightProperties $w"
button $w.buttons.help    -text Help      -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(editlightprop)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

frame        $w.lightdiff  -borderwidth 2 -relief ridge
label        $w.lightdiff.label -text "Light diffuse:"  -font {Arial 12 bold} 
pack         $w.lightdiff       -side top -anchor w
pack         $w.lightdiff.label -side top -anchor w

#red
entry        $w.lightdiff.entryred -width 10
scale        $w.lightdiff.scalered -from 0.0 -to 1.0 -length 240 -variable gomLightDiffScaleRed \
             -orient horizontal -label "Light diffuse RED" -tickinterval 0.2               \
			 -showvalue true -fg red -digits 3 -resolution 0.01 \
			 -command "lulShowLightDiffuseValueRed $w.lightdiff.entryred"
pack         $w.lightdiff.entryred -side top
pack         $w.lightdiff.scalered -side top -fill x

bind         $w.lightdiff.entryred <Return> "lulApplyEditLightProperties $w"

#green
entry        $w.lightdiff.entrygreen -width 10
scale        $w.lightdiff.scalegreen -from 0.0 -to 1.0 -length 240 -variable gomLightDiffScaleGreen \
             -orient horizontal -label "Light diffuse GREEN" -tickinterval 0.2              \
			 -showvalue true -fg green -digits 3 -resolution 0.01 \
			 -command "lulShowLightDiffuseValueGreen $w.lightdiff.entrygreen"
pack         $w.lightdiff.entrygreen -side top
pack         $w.lightdiff.scalegreen -side top -fill x

bind         $w.lightdiff.entrygreen <Return> "lulApplyEditLightProperties $w"

#blue
entry        $w.lightdiff.entryblue -width 10
scale        $w.lightdiff.scaleblue -from 0.0 -to 1.0 -length 240 -variable gomLightDiffScaleBlue \
             -orient horizontal -label "Light diffuse BLUE" -tickinterval 0.2               \
			 -showvalue true -fg blue -digits 3 -resolution 0.01 \
			 -command "lulShowLightDiffuseValueBlue $w.lightdiff.entryblue"
pack         $w.lightdiff.entryblue -side top
pack         $w.lightdiff.scaleblue -side top -fill x

bind         $w.lightdiff.entryblue <Return> "lulApplyEditLightProperties $w"

set          LightDiff [show light diffuse red]

             $w.lightdiff.scalered set $LightDiff
             $w.lightdiff.entryred delete 0 end
             $w.lightdiff.entryred insert 0 $LightDiff

set          LightDiff [show light diffuse green]

             $w.lightdiff.scalegreen set $LightDiff
             $w.lightdiff.entrygreen delete 0 end
             $w.lightdiff.entrygreen insert 0 $LightDiff

set          LightDiff [show light diffuse blue]

             $w.lightdiff.scaleblue set $LightDiff
             $w.lightdiff.entryblue delete 0 end
             $w.lightdiff.entryblue insert 0 $LightDiff

##################################################
# light position
frame        $w.lightpos -borderwidth 2 -relief raised
pack         $w.lightpos -side top -anchor w
label        $w.lightpos.label -text "Light position:"
pack         $w.lightpos.label -side top -anchor w

frame        $w.lightpos.frame1
pack         $w.lightpos.frame1 -side top -anchor w
# first row
radiobutton  $w.lightpos.frame1.nw -text "NW" -width 2 -value 1 -variable lulJunkLP \
              -command {define light position nw;lulInvalidateDisplay}
pack         $w.lightpos.frame1.nw -side left -anchor w
radiobutton  $w.lightpos.frame1.n  -text "N "  -width 2 -value 2 -variable lulJunkLP \
              -command {define light position n;lulInvalidateDisplay}
pack         $w.lightpos.frame1.n  -side left -anchor w
radiobutton  $w.lightpos.frame1.ne -text "NE" -width 2 -value 3 -variable lulJunkLP \
              -command {define light position ne;lulInvalidateDisplay}
pack         $w.lightpos.frame1.ne -side left -anchor w

frame        $w.lightpos.frame2
pack         $w.lightpos.frame2 -side top -anchor w
#second row
radiobutton  $w.lightpos.frame2.w  -text "W "  -width 2 -value 4 -variable lulJunkLP \
              -command {define light position w;lulInvalidateDisplay}
pack         $w.lightpos.frame2.w  -side left -anchor w
radiobutton  $w.lightpos.frame2.c  -text "C "  -width 2 -value 5 -variable lulJunkLP \
              -command {define light position c;lulInvalidateDisplay}
pack         $w.lightpos.frame2.c  -side left -anchor w
radiobutton  $w.lightpos.frame2.e  -text "E "  -width 2 -value 6 -variable lulJunkLP \
              -command {define light position e;lulInvalidateDisplay}
pack         $w.lightpos.frame2.e  -side left -anchor w

frame        $w.lightpos.frame3
pack         $w.lightpos.frame3 -side top -anchor w
#third row
radiobutton  $w.lightpos.frame3.sw -text "SW" -width 2 -value 7 -variable lulJunkLP \
              -command {define light position sw;lulInvalidateDisplay}
pack         $w.lightpos.frame3.sw -side left -anchor w
radiobutton  $w.lightpos.frame3.s  -text "S "  -width 2 -value 8 -variable lulJunkLP \
              -command {define light position s;lulInvalidateDisplay}
pack         $w.lightpos.frame3.s  -side left -anchor w
radiobutton  $w.lightpos.frame3.se -text "SE" -width 2 -value 9 -variable lulJunkLP \
              -command {define light position se;lulInvalidateDisplay}
pack         $w.lightpos.frame3.se -side left -anchor w

set LightPos [show light position]

 if       {$LightPos == "nw"}  {
   $w.lightpos.frame1.nw select
 } elseif {$LightPos == "n"}   {
   $w.lightpos.frame1.n select
 } elseif {$LightPos == "ne"}  {
   $w.lightpos.frame1.ne select
 } elseif {$LightPos == "w"}   {
   $w.lightpos.frame2.w select
 } elseif {$LightPos == "c"}   {
   $w.lightpos.frame2.c select
 } elseif {$LightPos == "e"}   {
   $w.lightpos.frame2.e select
 } elseif {$LightPos == "se"}  {
   $w.lightpos.frame3.se select
 } elseif {$LightPos == "s"}   {
   $w.lightpos.frame3.s select
 } elseif {$LightPos == "se"}  {
   $w.lightpos.frame3.se select
 } else {
   gomError "wrong light position returned"
 }

##################################################
frame        $w.control -borderwidth 2 -relief raised
label        $w.control.label -text "Continuous display: "
pack         $w.control -side top -anchor w -padx 4 -pady 4
pack         $w.control.label -side left -anchor w

radiobutton  $w.control.on  -text "On"  -value 1 -variable lulJunkID \
              -command {set gomInstantDisplay 1} 
radiobutton  $w.control.off -text "Off" -value 0 -variable lulJunkID \
              -command {set gomInstantDisplay 0} 
pack         $w.control.on  -side left -anchor w
pack         $w.control.off -side left -anchor w

# pick radiobutton as default
  if {$gomInstantDisplay} {
      $w.control.on  select
  } else {
      $w.control.off select
  }

}

#############################################################
# PROC
proc lulShowLightDiffuseValueRed { w value } {

     global gomLightDiffScaleRed
     global gomInstantDisplay

     set LightDiff $value

     $w delete 0 end
     $w insert 0 [format "%.4f" $LightDiff]

# check for instant update
     if {$gomInstantDisplay > 0} {
	 
     set LightDiffRed   [$w   get]

	 if {$LightDiffRed != ""} {
	     eval define light diffuse red $LightDiffRed
     }
	 
	 display}

}

#############################################################
# PROC
proc lulShowLightDiffuseValueGreen { w value } {

     global gomLightDiffScaleGreen
     global gomInstantDisplay

     set LightDiff $value

     $w delete 0 end
     $w insert 0 [format "%.4f" $LightDiff]

# check for instant update
     if {$gomInstantDisplay > 0} {

     set LightDiffGreen [$w get]

	 if {$LightDiffGreen != ""} {
	     eval define light diffuse green $LightDiffGreen
     }

     set LightDiffGreen  [$w get]

	 if {$LightDiffGreen != ""} {
	     eval define light diffuse green $LightDiffGreen
     }
	 
	 display}

}

#############################################################
# PROC
proc lulShowLightDiffuseValueBlue { w value } {

     global gomLightDiffScaleBlue
     global gomInstantDisplay

     set LightDiff $value

     $w delete 0 end
     $w insert 0 [format "%.4f" $LightDiff]

# check for instant update
     if {$gomInstantDisplay > 0} {

     set LightDiffBlue  [$w get]

	 if {$LightDiffBlue != ""} {
	     eval define light diffuse blue $LightDiffBlue
     }
	 
	 display}

}

#############################################################
# PROC
proc lulApplyEditLightProperties { w } {

     set LightDiffRed   [$w.lightdiff.entryred   get]

	 if {$LightDiffRed != ""} {
	     eval define light diffuse red $LightDiffRed
         set  LightDiff $LightDiffRed
         $w.lightdiff.scalered set $LightDiff
     }

     set LightDiffGreen [$w.lightdiff.entrygreen get]

	 if {$LightDiffGreen != ""} {
	     eval define light diffuse green $LightDiffGreen
         set  LightDiff $LightDiffGreen
         $w.lightdiff.scalegreen set $LightDiff
     }

     set LightDiffBlue  [$w.lightdiff.entryblue  get]

	 if {$LightDiffBlue != ""} {
	     eval define light diffuse blue $LightDiffBlue
         set  LightDiff $LightDiffBlue
         $w.lightdiff.scaleblue set $LightDiff
     }

     lulInvalidateDisplay
}
##################### copy timeseries to disk ############################
#
# copy timeseries to disk
#
# PROC
proc lulExportTimeSeriesFile {} {

     global gomHelpFile
     global gomControlFont


set w .gomcopytimeseries
catch {destroy $w}
toplevel $w 
wm title $w "Save timeseries"
wm iconname $w "Save timeseries"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss   -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply     -font "$gomControlFont" \
        -command "lulDoCopyTimeSeries2Disk $w"
button $w.buttons.help    -text Help      -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(copytimeseries2disk)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Save TimeSeries:"
pack   $w.label -side top -anchor w

frame  $w.frame1 -borderwidth 2 -relief ridge -width 30
pack   $w.frame1 -side top -fill x

set DistList [show monitor distance timeseries]

if {$DistList} {
  radiobutton $w.frame1.distlist -text "Distance list #:" \
               -value 0 -variable gomSaveTimeList
  pack        $w.frame1.distlist -side left -anchor w 
  entry       $w.frame1.value -width 5 
  pack        $w.frame1.value -side left -anchor w
  label       $w.frame1.flabel -text "File name: "
  pack        $w.frame1.flabel -side left -anchor w 
  entry       $w.frame1.fname -width 30
  pack        $w.frame1.fname -side left -anchor w 
  button      $w.frame1.fbrowse -text "Browse" -command "lulBrowseTimeSeriesFiles $w 1"
  pack        $w.frame1.fbrowse -side left -anchor w -padx 4
} else {
  label       $w.frame1.label -text "No distance time series are defined"
  pack        $w.frame1.label -side top -anchor w
  }

frame  $w.frame2 -borderwidth 2 -relief ridge
pack   $w.frame2 -side top -fill x

set AngList [show monitor angle timeseries]

if {$AngList} {
  radiobutton $w.frame2.anglist -text "Angle list #:     " \
               -value 1 -variable gomSaveTimeList
  pack        $w.frame2.anglist -side left -anchor w 
  entry       $w.frame2.value -width 5
  pack        $w.frame2.value -side left -anchor w
  label       $w.frame2.flabel -text "File name: "
  pack        $w.frame2.flabel -side left -anchor w 
  entry       $w.frame2.fname -width 30
  pack        $w.frame2.fname -side left -anchor w 
  button      $w.frame2.fbrowse -text "Browse" -command "lulBrowseTimeSeriesFiles $w 2"
  pack        $w.frame2.fbrowse -side left -anchor w -padx 4
} else {
  label       $w.frame2.label -text "No angle time series are defined"
  pack        $w.frame2.label -side top -anchor w
  }

frame  $w.frame3 -borderwidth 2 -relief ridge
pack   $w.frame3 -side top -fill x

set TorsList [show monitor torsion timeseries]

if {$TorsList} {
  radiobutton $w.frame3.torslist -text "Torsion list #:  " \
               -value 2 -variable gomSaveTimeList
  pack        $w.frame3.torslist -side left -anchor w 
  entry       $w.frame3.value -width 5
  pack        $w.frame3.value -side left -anchor w
  label       $w.frame3.flabel -text "File name: "
  pack        $w.frame3.flabel -side left -anchor w 
  entry       $w.frame3.fname -width 30
  pack        $w.frame3.fname -side left -anchor w 
  button      $w.frame3.fbrowse -text "Browse" -command "lulBrowseTimeSeriesFiles $w 3"
  pack        $w.frame3.fbrowse -side left -anchor w -padx 4 
} else {
  label       $w.frame3.label -text "No torsion time series are defined"
  pack        $w.frame3.label -side top -anchor w
  }

  if {!$DistList && !$AngList && !$TorsList} {
     $w.buttons.apply configure -state disabled
  }
}
#############################################################
# PROC
proc lulDoCopyTimeSeries2Disk { w } {

     global gomSaveTimeList

  switch $gomSaveTimeList {
      0 {set FileName [string trim [$w.frame1.fname get]];
         if {$FileName == ""} {return};
         export distance list [$w.frame1.value get] "$FileName"}
      1 {set FileName [string trim [$w.frame2.fname get]];
         if {$FileName == ""} {return};
         export angle    list [$w.frame2.value get] "$FileName"}
      2 {set FileName [string trim [$w.frame3.fname get]];
         if {$FileName == ""} {return};
         export torsion  list [$w.frame3.value get] "$FileName"}
	  default {lulErrorDialog "Please choose by pressing a radiobutton and giving a time series number"}
  }
}
############################################################
# PROC
proc lulBrowseTimeSeriesFiles { w which} {
    global gomTimeSeriesFileName


# return if no molecular systems defined
     if {[show molstructures] < 1} {
          gomError {ERROR: no structure available. Read a structure first!}
	      return
	 }

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
        -title "Save Time Series File"]]

    set file [string trim $file]

    if {$file != ""} {
      switch $which {
        1 {$w.frame1.fname delete 0 end;$w.frame1.fname insert 0 "$file"}
        2 {$w.frame2.fname delete 0 end;$w.frame2.fname insert 0 "$file"}
        3 {$w.frame3.fname delete 0 end;$w.frame3.fname insert 0 "$file"}
      }
# change to current directory
    lulChangeDirectory "$file"
    }
}

##################### display atom labels ############################
#
# define atom widget
#
proc lulDisplayAtomLabels {} {

     global gomHelpFile
     global gomControlFont
     global gomLabelSelectionDisplayStyle

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomatomlabels
catch {destroy $w}
toplevel $w 
wm title $w "Display Atom Labels"
wm iconname $w "Display Atom Labels"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -command "destroy $w" \
        -font "$gomControlFont"
button $w.buttons.apply   -text Apply -command "lulDoDisplayAtomLabels" \
        -font "$gomControlFont"
button $w.buttons.help    -text Help -command \
       "htmlShowHelp $gomHelpFile(atomlabeldisplay)" \
        -font "$gomControlFont"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Display atom label(s):" -font "$gomControlFont"
pack   $w.label -side top -anchor w

frame   $w.frame1
pack    $w.frame1 -side top

frame   $w.frame1.selstate -borderwidth 2 -relief raised
pack    $w.frame1.selstate -side right -anchor w -padx 5 -pady 3 -fill x

label   $w.frame1.selstate.label -text "Display state:"
pack    $w.frame1.selstate.label -side top  -anchor w

radiobutton $w.frame1.selstate.on  -text "On"  -value 1   \
             -variable gomLabelSelectionState         \
             -command  {}
radiobutton $w.frame1.selstate.off -text "Off" -value 0   \
             -variable gomLabelSelectionState         \
             -command  {}

pack        $w.frame1.selstate.on $w.frame1.selstate.off 

            $w.frame1.selstate.on select

lulCreateAtomInputEntries $w.frame1 30 "*" "*" "*" "gomLabelSelection" 1

frame   $w.frame4 -borderwidth 2 -relief raised
pack    $w.frame4 -side top -anchor w -padx 5 -pady 10 -fill x

label   $w.frame4.label -text "Font:" 
pack    $w.frame4.label -pady 5 -side top  -anchor w

eval tk_optionMenu $w.frame4.options gomAtomLabelFont  \
BITMAP_8_BY_13 \
BITMAP_9_BY_15 \
BITMAP_TIMES_ROMAN_10 \
BITMAP_TIMES_ROMAN_24 \
BITMAP_HELVETICA_10 \
BITMAP_HELVETICA_12 \
BITMAP_HELVETICA_18

pack   $w.frame4.options -side left -anchor w -pady 5

frame   $w.frame5 -borderwidth 2 -relief ridge
pack    $w.frame5 -pady 4 -side top -anchor w

label   $w.frame5.label -text "Display type: "
pack    $w.frame5.label -side left -anchor w

radiobutton $w.frame5.residue -text "Residue" -value residue   \
             -variable gomLabelSelectionDisplayStyle           \
             -command  {}
radiobutton $w.frame5.atom    -text "Atom"    -value atom   \
             -variable gomLabelSelectionDisplayStyle        \
             -command  {}
radiobutton $w.frame5.full    -text "Full"    -value full   \
             -variable gomLabelSelectionDisplayStyle        \
             -command  {}

pack        $w.frame5.residue $w.frame5.atom $w.frame5.full -side left -anchor w -padx 4

            $w.frame5.atom select

}

###########################################################################
# react to atom labels
#
# PROC
proc lulDoDisplayAtomLabels { } {

     global gomLabelSelectionSegment gomLabelSelectionResidue gomLabelSelectionAtom
     global gomLabelSelectionState
     global gomLabelSelectionDisplayStyle

     if {$gomLabelSelectionState} {
       atom label  $gomLabelSelectionSegment $gomLabelSelectionResidue $gomLabelSelectionAtom
       atom label  $gomLabelSelectionDisplayStyle
     } else {
       atom -label $gomLabelSelectionSegment $gomLabelSelectionResidue $gomLabelSelectionAtom
     }
     lulInvalidateDisplay
}
#########################################################################

##################### center system ############################
#
# center system
#
proc lulCenterSystem {} {

     global gomHelpFile
     global gomControlFont

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomcentersystem
catch {destroy $w}
toplevel $w 
wm title $w "Center system"
wm iconname $w "Center system"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -command "destroy $w" \
        -font "$gomControlFont"
button $w.buttons.apply   -text Apply -command "lulDoCenterSystem $w" \
        -font "$gomControlFont"
button $w.buttons.help    -text Help -command \
       "htmlShowHelp $gomHelpFile(centersystem)" \
        -font "$gomControlFont"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Center system:" -font "$gomControlFont"
pack   $w.label -side top -anchor w

lulCreateAtomInputEntries $w 30 "*" "*" "*" "gomDisplay" 1

frame $w.frame4 -borderwidth 2 -relief ridge
pack $w.frame4 -side top
label       $w.frame4.label -text "Center action: "
radiobutton $w.frame4.system   -text "System"   -width 8 -value 1 -variable gomCenterValue
radiobutton $w.frame4.rotation -text "Rotation" -width 8 -value 0 -variable gomCenterValue
pack $w.frame4.label $w.frame4.system $w.frame4.rotation -padx 3 -side left
 
$w.frame4.system select

}

###########################################################################
# react to center system
#
# PROC
proc lulDoCenterSystem { w } {

    global gomCenterValue

    set Segment [$w.segment.input get]
    set Residue [$w.residue.input get]
    set Atom    [$w.atom.input get]

    if {$gomCenterValue} {
	center system   $Segment $Residue $Atom
    } else {
	center rotation $Segment $Residue $Atom
    }

    lulInvalidateDisplay
}
#########################################################################

################### text widget with text input #############
#
# Text widget
#
# PROC
proc lulDisplayEditTextFile InputFile {

     global gomControlFont

set w .gomtextedit
catch {destroy $w}
toplevel    $w
wm title    $w "Display/Edit text"
wm iconname $w "text"

frame  $w.buttons -borderwidth 2 -relief raised -bd 2
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.save -text Save -font "$gomControlFont" \
        -command "lulSaveDisplayEditTextFile $w $InputFile;lulMessageDialog {File written to disk!}"
pack   $w.buttons.save -side left -expand 1
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
pack   $w.buttons.dismiss -side left -expand 1

frame  $w.frame1 -borderwidth 2 -relief raised -bd 2
pack   $w.frame1 -side top -fill x -pady 4
label  $w.frame1.label -text "Output file written: $InputFile" -font $gomControlFont
pack   $w.frame1.label -side top -anchor w

text $w.text -relief sunken -bd 2 -yscrollcommand "$w.scroll set" -setgrid 1 \
	-height 30
scrollbar $w.scroll -command "$w.text yview"
pack $w.scroll -side right -fill y
pack $w.text -expand yes -fill both
$w.text delete 0.0 end

set File [open $InputFile r]
    while {![eof $File]} {
    $w.text insert end [read $File 1000]
    }

close $File

$w.text mark set insert 0.0
}

#####################################################################################
# PROC
proc lulSaveDisplayEditTextFile { w file } {

gomPrint "saving file: $file"

set File [open $file w]
puts $File [$w.text get 0.0 end]
close $File
}


##################### export input ############################
#
# define export input
#
proc lulCreateExportProgramInputFile {} {

     global gomControlFont
     global gomHelpDir
     global gomHelpFile
     global gomFocusStructure
     global gomExportInputType

set w .gomexportcreateinput
catch {destroy $w}
toplevel $w 
wm title $w "Export input"
wm iconname $w "Export input"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulApplyCreateExportProgramInputFile $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(exportinput)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Export Input:" -font $gomControlFont
pack   $w.label -side top -anchor w

frame        $w.frame1 -borderwidth 2 -relief ridge 
label        $w.frame1.label -text "File formats: "

radiobutton  $w.frame1.gamess -text "GAMESS"   -variable gomExportInputType -value gamess \
              -command "lulRetrieveGAMESSInputWidget $w.frame5.landscape"
radiobutton  $w.frame1.icon8  -text "ICON8"    -variable gomExportInputType -value icon8 \
              -command "lulRetrieveICON8InputWidget $w.frame5.landscape"
radiobutton  $w.frame1.mopac  -text "Mopac"    -variable gomExportInputType -value mopac \
              -command "lulRetrieveMOPACInputWidget $w.frame5.landscape"
radiobutton  $w.frame1.probesurf -text "Probesurf" -variable gomExportInputType -value probesurf \
              -command "lulRetrieveProbesurfInputWidget $w.frame5.landscape"
radiobutton  $w.frame1.user -text "USER"       -variable gomExportInputType -value user \
              -command "lulRetrieveUSERInputWidget $w.frame5.landscape"

#$w.frame1.probesurf select
               
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.gamess $w.frame1.icon8 $w.frame1.mopac $w.frame1.probesurf $w.frame1.user -side left

frame  $w.frame2
label  $w.frame2.label     -text "File name: "
entry  $w.frame2.filename  -width 50
button $w.frame2.browse    -text "Browse..." -command "lulBrowseExportInputFile $w"

pack   $w.frame2 -side top -pady 6
pack   $w.frame2.label $w.frame2.filename -side left
pack   $w.frame2.browse -padx 3 -side left
bind   $w.frame2.filename <Return> "eval lulApplyCreateExportInputFile $w"

frame       $w.frame4 -borderwidth 2 -relief ridge
pack        $w.frame4 -side top -anchor w -pady 4
label       $w.frame4.label -text "Choose structure #:"
pack        $w.frame4.label -side left -anchor w

set List ""
for {set i 1} { $i <= [show molstruct]} {incr i} {
 lappend List "$i"
}

set gomFocusStructure 1

eval tk_optionMenu $w.frame4.struct gomFocusStructure $List
pack $w.frame4.struct -side left -anchor w -pady 3

# specific details for the various export options

frame       $w.frame5 
pack        $w.frame5 -side top -anchor w -pady 4

# probesurf
if {$gomExportInputType == "gamess"} {
    lulRetrieveGAMESSInputWidget $w.frame5.landscape
} elseif {$gomExportInputType == "icon8"} {
    lulRetrieveICON8InputWidget  $w.frame5.landscape
} elseif {$gomExportInputType == "mopac"} {
    lulRetrieveMOPACInputWidget  $w.frame5.landscape
} elseif {$gomExportInputType == "probesurf"} {
    lulRetrieveProbesurfInputWidget $w.frame5.landscape
} elseif {$gomExportInputType == "user"} {
    lulRetrieveUSERInputWidget $w.frame5.landscape
} else {
    gomError "wrong program option '$gomExportInputType'"
}
 
}

#####################################################################################
# PROC
proc lulRetrieveGAMESSInputWidget { w } {

catch {destroy $w}
update idletasks

}
#####################################################################################
# PROC
proc lulRetrieveICON8InputWidget { w } {

catch {destroy $w}
update idletasks

}
#####################################################################################
# PROC
proc lulRetrieveMOPACInputWidget { w } {

catch {destroy $w}
update idletasks

}
#####################################################################################
# PROC
proc lulRetrieveUSERInputWidget { w } {

catch {destroy $w}
update idletasks

}

#####################################################################################
# PROC
proc lulRetrieveProbesurfInputWidget { w } {

catch {destroy $w}
update idletasks

frame       $w 
pack        $w -side top

set Size     [show sizeofsystem]

	set xp   "60"
	set yp   "60"
	set zp   "60"

	set minx "-$Size"
	set maxx "$Size"

	set miny "-$Size"
	set maxy "$Size"

	set minz "-$Size"
	set maxz "$Size"

if {[show contour defined]} {

  if {![lulYesNoMessageDialog {There is at least one contour file defined\nDo you want to use the same header data}]} {

    set RawData "[show contour header [show contour name 1]]"
	set xp   [lindex $RawData 0]
	set yp   [lindex $RawData 1]
	set zp   [lindex $RawData 2]

	set minx [lindex $RawData 3]
	set maxx [lindex $RawData 4]

	set miny [lindex $RawData 5]
	set maxy [lindex $RawData 6]

	set minz [lindex $RawData 7]
	set maxz [lindex $RawData 8]

  }

}

frame  $w.frame1
pack   $w.frame1 -side top -anchor w -pady 3

label  $w.frame1.label1 -text "X min: " -width 8 
entry  $w.frame1.xmin   -width 15
$w.frame1.xmin insert 0 $minx

label  $w.frame1.label2 -text "X max: " -width 8
entry  $w.frame1.xmax   -width 15
$w.frame1.xmax insert 0 $maxx

pack   $w.frame1.label1 $w.frame1.xmin  -side left
pack   $w.frame1.label2 $w.frame1.xmax  -side left

frame  $w.frame2
pack   $w.frame2 -side top -anchor w

label  $w.frame2.label1 -text "Y min: " -width 8
entry  $w.frame2.ymin   -width 15
$w.frame2.ymin insert 0 $miny

label  $w.frame2.label2 -text "Y max: " -width 8
entry  $w.frame2.ymax   -width 15
$w.frame2.ymax insert 0 $maxy

pack   $w.frame2.label1 $w.frame2.ymin  -side left
pack   $w.frame2.label2 $w.frame2.ymax  -side left

frame  $w.frame3
pack   $w.frame3 -side top -anchor w -pady 3

label  $w.frame3.label1 -text "Z min: " -width 8
entry  $w.frame3.zmin   -width 15
$w.frame3.zmin insert 0 $minz

label  $w.frame3.label2 -text "Z max: " -width 8
entry  $w.frame3.zmax   -width 15
$w.frame3.zmax insert 0 $maxz

pack   $w.frame3.label1 $w.frame3.zmin  -side left
pack   $w.frame3.label2 $w.frame3.zmax  -side left

frame  $w.frame4
pack   $w.frame4 -side top -anchor w
label  $w.frame4.label1 -text "Points x-direction:" -width 20
entry  $w.frame4.xpts   -width 10
$w.frame4.xpts insert 0 $xp
pack   $w.frame4.label1 $w.frame4.xpts -side left -padx 4

frame  $w.frame5
pack   $w.frame5 -side top -anchor w
label  $w.frame5.label2 -text "Points y-direction:" -width 20
entry  $w.frame5.ypts   -width 10
$w.frame5.ypts insert 0 $yp
pack   $w.frame5.label2 $w.frame5.ypts -side left -padx 4


frame  $w.frame6
pack   $w.frame6 -side top -anchor w
label  $w.frame6.label3 -text "Points z-direction:" -width 20
entry  $w.frame6.zpts   -width 10
$w.frame6.zpts insert 0 $zp
pack   $w.frame6.label3 $w.frame6.zpts -side left -padx 4

frame  $w.frame7
label  $w.frame7.label     -text "Probe radius (A): " -width 20
entry  $w.frame7.proberad  -width 10
$w.frame7.proberad insert 0 "2.0"
pack   $w.frame7 -side top -anchor w -pady 4
pack   $w.frame7.label $w.frame7.proberad -side left
 
}

#####################################################################################
# PROC
proc lulBrowseExportInputFile { w } {


# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"INP files"		{.inp}		TEXT}
	{"All files"		*}
    }


    set file [string trim [tk_getSaveFile -filetypes $types  \
         -defaultextension .inp -title "Save Input File"]]

    if [string compare $file ""] {
      $w.frame2.filename delete 0 end
      $w.frame2.filename insert 0 "$file"
# change to current directory
      lulChangeDirectory "$file"
    }
}

#####################################################################################
# PROC
proc lulApplyCreateExportProgramInputFile { w } {

    global gomFocusStructure
    global gomExportInputType

    set file [$w.frame2.filename get]

    if {$file != ""} {

     if {$gomExportInputType == "probesurf"} {
       set ExtraInput [lulGetDataFromProbesurfWidget $w.frame5.landscape]
       eval export input $gomExportInputType $gomFocusStructure \"$file\" $ExtraInput
     } else {
       export input $gomExportInputType $gomFocusStructure "$file"
     }
     lulDisplayEditTextFile "$file"
    }
}

#####################################################################################
# PROC
proc lulGetDataFromProbesurfWidget { w } {

set Size     [show sizeofsystem]

              set Xmin     [$w.frame1.xmin get]
              if {[string trim $Xmin] == ""} {
			      set Xmin "-$Size"
				  }
              set Xmax     [$w.frame1.xmax get]
              if {[string trim $Xmax] == ""} {
			      set Xmax $Size
				  }

              set Ymin     [$w.frame2.ymin get]
              if {[string trim $Ymin] == ""} {
			      set Ymin "-$Size"
				  }

              set Ymax     [$w.frame2.ymax get]
              if {[string trim $Ymax] == ""} {
			      set Ymax $Size
				  }

              set Zmin     [$w.frame3.zmin get]
              if {[string trim $Zmin] == ""} {
			      set Zmin "-$Size"
				  }
              set Zmax     [$w.frame3.zmax get]
              if {[string trim $Zmax] == ""} {
			      set Zmax $Size
				  }

              set Xpts     [$w.frame4.xpts get]
              if {[string trim $Xpts] == ""} {
			      set Xpts 60
				  }
              set Ypts     [$w.frame5.ypts get]
              if {[string trim $Ypts] == ""} {
			      set Ypts 60
				  }
              set Zpts     [$w.frame6.zpts get]
              if {[string trim $Zpts] == ""} {
			      set Zpts 60
				  }

              set ProbeRad [$w.frame7.proberad get]
              if {[string trim $ProbeRad] == ""} {
			      set ProbeRad 2.0
				  }

    return "$Xmin $Xmax $Ymin $Ymax $Zmin $Zmax $Xpts $Ypts $Zpts $ProbeRad"

}

################### text widget with text input #############
#
# Text widget
#
# PROC
proc lulDisplayEditRunTCLScript { } {

     global gomHelpDir
     global gomHelpFile
     global gomControlFont

set w .gomtcledit
catch {destroy $w}
toplevel    $w
wm title    $w "Display/Edit Tcl/Tk script"
wm iconname $w "text"

frame  $w.buttons -borderwidth 2 -relief ridge -bd 2
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.save -text "Run Tcl/Tk script" -font "$gomControlFont" \
        -command "lulRunTclTkScript $w"
pack   $w.buttons.save -side left -expand 1
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "lulSaveTclScriptFromWidget $w;destroy $w"
pack   $w.buttons.dismiss -side left -expand 1
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(runtclscript)"
pack   $w.buttons.help -side left -expand 1

frame  $w.frame0 -borderwidth 2 -relief raised -bd 2
pack   $w.frame0 -side top -anchor w -fill x

label  $w.frame0.label     -text "Input Tcl/Tk script name: " -font $gomControlFont
entry  $w.frame0.filename  -width 45
button $w.frame0.browse    -text "Browse..." -command "lulImportTclFile $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename $w.frame0.browse -side left -padx 5

bind   $w.frame0.filename <Return> "lulImportTclTkScript2Widget $w"

#frame  $w.frame1 -borderwidth 2 -relief raised -bd 2
#pack   $w.frame1 -side top -fill x -pady 4
#label  $w.frame1.label -text "Output file written: $InputFile" -font $gomControlFont
#pack   $w.frame1.label -side top -anchor w

text $w.text -relief sunken -bd 2 -yscrollcommand "$w.scroll set" -setgrid 1 \
	-height 30 -width 80
scrollbar $w.scroll -command "$w.text yview"
pack $w.scroll -side right -fill y
pack $w.text -expand yes -fill both

$w.text delete 0.0 end

$w.text mark set insert 0.0
}
#####################################################################################
# PROC
proc lulRunTclTkScript { w } {

      $w.text configure -bg red
      $w.text configure -state disabled
      if [catch {uplevel #0 [$w.text get 0.0 end]} RetVal] {
         gomError $RetVal
      }
      $w.text configure -bg white
      $w.text configure -state normal

}
#####################################################################################
# PROC
proc lulImportTclTkScript2Widget { w } {

    set file [string trim [$w.frame0.filename get]]

    if [string compare $file ""] {

# check to see that the file exists...
    if {![file exists $file]} {
       gomError {file: $file does not exist}
       return
    } 

     $w.text delete 0.0 end

     set File_p [open "$file" r]
     while {![eof $File_p]} {
         $w.text insert end [read $File_p 1000]
     }

     $w.text mark set insert 0.0

     close $File_p
    }
}

#####################################################################################
# PROC
proc lulSaveTclScriptFromWidget { w } {

    set script [string trim [$w.text get 0.0 end]]
# save if there is something
    if {$script != ""} {
        set InputText "Do you want to save the Tcl/Tk script?"
        after idle {.gomdialogtclwrite.msg configure -wraplength 4i}
        set i [tk_dialog .gomdialogtclwrite "gOpenMol Question" $InputText \
        question 0 "No" "Yes" ]

        switch $i {
          0 {return}
        }

    
    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Tcl files"		{.tcl .TCL}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
         -defaultextension .tcl -title "Write Tcl/Tk Script"]]

    if {$file != ""} {

      if [file isdirectory $file] {
         gomError "the file is a directory!"
         return
      }

      set f [open "$file" w]
      puts $f $script
      close $f
    }
    }
}

##################### import and execute tcl file ############################
# PROC 
proc lulImportTclFile { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Tcl files"		{.tcl .TCL}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .tcl -title "Read Tcl File"]]


    if [string compare $file ""] {

# change to current directory
     lulChangeDirectory "$file"

     $w.frame0.filename delete 0 end
     $w.frame0.filename insert 0 "$file"

     $w.text delete 0.0 end

     set File_p [open "$file" r]
     while {![eof $File_p]} {
         $w.text insert end [read $File_p 1000]
     }

     $w.text mark set insert 0.0

     close $File_p

    }
}

####################################################################
# 
# Stereo display control
#
proc lulDisplayStereoPairControl { } {

     global gomHelpDir
     global gomHelpFile
     global gomStructure
     global gomControlFont
	 global gomStereoPairStateValue

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          gomError "Read in a structure first"
	      return
	 }

set w .gomstereopair
catch {destroy $w}
toplevel $w 
wm title $w "Stereo pair control"
wm iconname $w "Stereo pair"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyStereoPair $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(plotstereopair)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1 -padx 4

set SState     [show spair state]
set SDistance  [show spair distance]
set SAngle     [show spair angle]

label  $w.label -text "Stereo pair display control:" -font {Arial 12 bold}
pack   $w.label -side top -anchor w -pady 4

frame $w.frame0 -borderwidth 2 -relief ridge
pack $w.frame0 -side top

label $w.frame0.label    -text "Distance between structures:" -width 30 
entry $w.frame0.distance -width 12
pack  $w.frame0.label $w.frame0.distance -side left

$w.frame0.distance delete 0 end
$w.frame0.distance insert 0 $SDistance

frame $w.frame1 -borderwidth 2 -relief ridge
pack $w.frame1 -side top

label $w.frame1.label    -text "Tilt angle (degrees)+-:" -width 30
entry $w.frame1.angle    -width 12
pack  $w.frame1.label $w.frame1.angle -side left

$w.frame1.angle delete 0 end
$w.frame1.angle insert 0 $SAngle

frame $w.frame2 -borderwidth 2 -relief ridge
pack $w.frame2 -side top
label       $w.frame2.label -text "Display state: "
radiobutton $w.frame2.on  -text "On"  -width 4 -value 1 -variable gomStereoPairStateValue
radiobutton $w.frame2.off -text "Off" -width 4 -value 0 -variable gomStereoPairStateValue
pack $w.frame2.label $w.frame2.on $w.frame2.off -padx 3 -side left

if {$SState == "on"} {
  $w.frame2.on  select
  $w.frame2.off deselect
} elseif {$SState == "off"} {
  $w.frame2.on  deselect
  $w.frame2.off select
} else {
  gomError "unknown state '$SState' returned from 'show spair state' command"
}

}

####################################################################
# PROC
proc lulApplyStereoPair { w } {

     global gomStereoPairStateValue

     set Distance [string trim [$w.frame0.distance get]]
     define spair distance $Distance

     set Angle    [string trim [$w.frame1.angle    get]]
     define spair angle $Angle

     if {$gomStereoPairStateValue} {
	  define spair on
	 } else {
      define spair off
	 }

     lulInvalidateDisplay
}

####################################################################
# 
# PROC
proc lulHandleClickedAtomDown { Structure AtomNr } {
    catch {lulAtomTree::pick_node $Structure $AtomNr}
}

####################################################################
# 
# PROC
proc lulHandleClickedAtomRelease { Structure AtomNr } {
}

proc lulHandlePickAtomEventStart { } {
     global lulGlobalEnteringWidgetValue

     bind all <Enter> {set lulGlobalEnteringWidgetValue %W}
}

####################################################################
# 
# PROC
proc lulHandlePickAtomEventEnd { StructureList AtomList Seg Res Atm } {

     global lulGlobalEnteringWidgetValue
     global gomAtomPickingTargetWidgets

     update

     bind all <Enter> {}

     if {[llength $StructureList] < 1} return

     set Structure [join $StructureList ,]
     set AtomNr    [join $AtomList      ,]

# Command line
    if {$lulGlobalEnteringWidgetValue == ".commandframe.commandline"} {
       set String [string trim [.commandframe.commandline get]]
       set String "$String $Seg $Res $Atm"

      .commandframe.commandline delete 0 end
      .commandframe.commandline insert 0 $String

       set lulGlobalEnteringWidgetValue ""
       lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomimpose.structure1.frame1.segment_input" ||
              $lulGlobalEnteringWidgetValue == ".gomimpose.structure1.frame2.residue_input" ||
              $lulGlobalEnteringWidgetValue == ".gomimpose.structure1.frame3.atom_input"} {

              lulHandleAtomSimpose1Picking $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomimpose.structure2.frame1.segment_input" ||
              $lulGlobalEnteringWidgetValue == ".gomimpose.structure2.frame2.residue_input" ||
              $lulGlobalEnteringWidgetValue == ".gomimpose.structure2.frame3.atom_input"} {

              lulHandleAtomSimpose2Picking $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomtraceatoms.frame1.segment_input" ||
              $lulGlobalEnteringWidgetValue == ".gomtraceatoms.frame2.residue_input" ||
              $lulGlobalEnteringWidgetValue == ".gomtraceatoms.frame3.atom_input"} {

              lulHandleAtomTracePicking $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomselectldpatoms.frame1.segment" ||
              $lulGlobalEnteringWidgetValue == ".gomselectldpatoms.frame2.residue" ||
              $lulGlobalEnteringWidgetValue == ".gomselectldpatoms.frame3.atom"} {

              lulHandleAtomSelectLDPPicking $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomselectclusteratoms.frame1.segment" ||
              $lulGlobalEnteringWidgetValue == ".gomselectclusteratoms.frame2.residue" ||
              $lulGlobalEnteringWidgetValue == ".gomselectclusteratoms.frame3.atom"} {

              lulHandleAtomClusterPicking $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomselectrdfatoms.frame1.segment" ||
              $lulGlobalEnteringWidgetValue == ".gomselectrdfatoms.frame2.residue" ||
              $lulGlobalEnteringWidgetValue == ".gomselectrdfatoms.frame3.atom"} {

              lulHandleAtomRDF1Picking $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomselectrdfatoms.frame11.segment" ||
              $lulGlobalEnteringWidgetValue == ".gomselectrdfatoms.frame12.residue" ||
              $lulGlobalEnteringWidgetValue == ".gomselectrdfatoms.frame13.atom"} {

              lulHandleAtomRDF2Picking $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gommeansquaredisplacement.frame1.segment" ||
              $lulGlobalEnteringWidgetValue == ".gommeansquaredisplacement.frame2.residue" ||
              $lulGlobalEnteringWidgetValue == ".gommeansquaredisplacement.frame3.atom"} {

              lulHandleAtomMSDPicking $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomplotrmsd.frame1.segment_input" ||
              $lulGlobalEnteringWidgetValue == ".gomplotrmsd.frame2.residue_input" ||
              $lulGlobalEnteringWidgetValue == ".gomplotrmsd.frame3.atom_input"} {

              lulHandleAtomRMSDPicking $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomcontourplanexyz.frame1.coord.p1xyz"} {

              lulHandleAtomPlaneXYZPicking 1 $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomcontourplanexyz.frame1.coord.p2xyz"} {

              lulHandleAtomPlaneXYZPicking 2 $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomcontourplanexyz.frame1.coord.p3xyz"} {

              lulHandleAtomPlaneXYZPicking 3 $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomcontourclipplanexyz.frame1.coord.p1xyz"} {

              lulHandleAtomClipPlaneXYZPicking 1 $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomcontourclipplanexyz.frame1.coord.p2xyz"} {

              lulHandleAtomClipPlaneXYZPicking 2 $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } elseif {$lulGlobalEnteringWidgetValue == ".gomcontourclipplanexyz.frame1.coord.p3xyz"} {

              lulHandleAtomClipPlaneXYZPicking 3 $Structure $AtomNr

              set lulGlobalEnteringWidgetValue ""
              lulInvalidateDisplay

    } else {
        set widget $lulGlobalEnteringWidgetValue
	while { "" != $widget } {
            if { [info exist gomAtomPickingTargetWidgets($widget)] } {

		set SegmentWidget $widget.[lindex $gomAtomPickingTargetWidgets($widget) 0]
		set ResidueWidget $widget.[lindex $gomAtomPickingTargetWidgets($widget) 1]
		set AtomWidget    $widget.[lindex $gomAtomPickingTargetWidgets($widget) 2]
		set AllowMultipleAtoms [lindex $gomAtomPickingTargetWidgets($widget) 3]

                if { $lulGlobalEnteringWidgetValue != $SegmentWidget &&
                     $lulGlobalEnteringWidgetValue != $ResidueWidget &&
                     $lulGlobalEnteringWidgetValue != $AtomWidget } return

                if {$AllowMultipleAtoms} {
                    set SegmentString ""
                    set ResidueString ""
                    set AtomString    ""
                    catch {set SegmentString [string trim [$SegmentWidget get]]}
                    catch {set ResidueString [string trim [$ResidueWidget get]]}
                    catch {set AtomString    [string trim [$AtomWidget    get]]}
                    # Segment
                    if { "" == $SegmentString || "*" == $SegmentString} {
                        set SegmentString $Seg
                    } else {
                        set SegmentString [lulJoinAtomList "$SegmentString,$Seg"]
                    }
                    # Residue
                    if { "" == $ResidueString || "*" == $ResidueString} {
                        set ResidueString $Res
                    } else {
                        set ResidueString [lulJoinAtomList "$ResidueString,$Res"]
                    }
                    # Atom
                    if { "" == $AtomString || "*" == $AtomString} {
                        set AtomString $Atm
                    } else {
                        set AtomString [lulJoinAtomList "$AtomString,$Atm"]
                    }
                } else {
                    set IStruct [lindex $StructureList 0]
                    set IAtom   [lindex $AtomList      0]
                    set SegmentString [show atom segment $IAtom $IStruct]
                    set ResidueString [show atom residue $IAtom $IStruct]
                    set AtomString    $IAtom
                }
                # segment
                catch {
                    $SegmentWidget delete 0 end
                    $SegmentWidget insert 0 $SegmentString
                }
                # residue
                catch {
                    $ResidueWidget delete 0 end
                    $ResidueWidget insert 0 $ResidueString
                }
                # atom
                catch {
                    $AtomWidget    delete 0 end
                    $AtomWidget    insert 0 $AtomString
                }

                set lulGlobalEnteringWidgetValue ""

		break
            }

            # Drop the last part.
            set widget [join [lrange [split $widget .] 0 end-1] .]
        }
    }
}

###1
############################################################################
# PROC
proc lulHandleAtomSimpose1Picking { Structure AtomNr } {

# segment
       set StringSegment [string trim [.gomimpose.structure1.frame1.segment_input get]]
       set StringLenSegment [llength $StringSegment]
       if {$StringSegment == ""} {
          .gomimpose.structure1.frame1.segment_input insert 0 "*"
       }
# residue
       set StringResidue [string trim [.gomimpose.structure1.frame2.residue_input get]]
       set StringLenResidue [llength $StringResidue]
       if {$StringResidue == ""} {
          .gomimpose.structure1.frame2.residue_input insert 0 "*"
       }
# atom
       set StringAtom [string trim [.gomimpose.structure1.frame3.atom_input get]]
       set StringLenAtom [llength $StringAtom]
       if {$StringAtom == "" || $StringAtom == "*"} {
          .gomimpose.structure1.frame3.atom_input delete 0 end 
          .gomimpose.structure1.frame3.atom_input insert 0 $AtomNr 
       } else {
          .gomimpose.structure1.frame3.atom_input delete 0 end 
          .gomimpose.structure1.frame3.atom_input insert 0 "$StringAtom,$AtomNr" 
       }

}
###2
############################################################################
# PROC
proc lulHandleAtomSimpose2Picking { Structure AtomNr } {

# segment
       set StringSegment [string trim [.gomimpose.structure2.frame1.segment_input get]]
       set StringLenSegment [llength $StringSegment]
       if {$StringSegment == ""} {
          .gomimpose.structure2.frame1.segment_input insert 0 "*"
       }
# residue
       set StringResidue [string trim [.gomimpose.structure2.frame2.residue_input get]]
       set StringLenResidue [llength $StringResidue]
       if {$StringResidue == ""} {
          .gomimpose.structure2.frame2.residue_input insert 0 "*"
       }
# atom
       set StringAtom [string trim [.gomimpose.structure2.frame3.atom_input get]]
       set StringLenAtom [llength $StringAtom]
       if {$StringAtom == "" || $StringAtom == "*"} {
          .gomimpose.structure2.frame3.atom_input delete 0 end 
          .gomimpose.structure2.frame3.atom_input insert 0 $AtomNr 
       } else {
          .gomimpose.structure2.frame3.atom_input delete 0 end 
          .gomimpose.structure2.frame3.atom_input insert 0 "$StringAtom,$AtomNr" 
       }

}

############################################################################
# PROC
proc lulHandleAtomTracePicking { Structure AtomNr } {

# segment
       set StringSegment [string trim [.gomtraceatoms.frame1.segment_input get]]
       set StringLenSegment [llength $StringSegment]
       if {$StringSegment == ""} {
          .gomtraceatoms.frame1.segment_input insert 0 "*"
       }
# residue
       set StringResidue [string trim [.gomtraceatoms.frame2.residue_input get]]
       set StringLenResidue [llength $StringResidue]
       if {$StringResidue == ""} {
          .gomtraceatoms.frame2.residue_input insert 0 "*"
       }
# atom
       set StringAtom [string trim [.gomtraceatoms.frame3.atom_input get]]
       set StringLenAtom [llength $StringAtom]
       if {$StringAtom == "" || $StringAtom == "*"} {
          .gomtraceatoms.frame3.atom_input delete 0 end 
          .gomtraceatoms.frame3.atom_input insert 0 $AtomNr 
       } else {
          .gomtraceatoms.frame3.atom_input delete 0 end 
          .gomtraceatoms.frame3.atom_input insert 0 "$StringAtom,$AtomNr" 
       }

}

############################################################################
# PROC
proc lulHandleAtomSelectLDPPicking { Structure AtomNr } {

# segment
       set StringSegment [string trim [.gomselectldpatoms.frame1.segment get]]
       set StringLenSegment [llength $StringSegment]
       if {$StringSegment == ""} {
          .gomselectldpatoms.frame1.segment insert 0 "*"
       }
# residue
       set StringResidue [string trim [.gomselectldpatoms.frame2.residue get]]
       set StringLenResidue [llength $StringResidue]
       if {$StringResidue == ""} {
          .gomselectldpatoms.frame2.residue insert 0 "*"
       }
# atom
       set StringAtom [string trim [.gomselectldpatoms.frame3.atom get]]
       set StringLenAtom [llength $StringAtom]
       if {$StringAtom == "" || $StringAtom == "*"} {
          .gomselectldpatoms.frame3.atom delete 0 end 
          .gomselectldpatoms.frame3.atom insert 0 $AtomNr 
       } else {
          .gomselectldpatoms.frame3.atom delete 0 end 
          .gomselectldpatoms.frame3.atom insert 0 "$StringAtom,$AtomNr" 
       }

}

############################################################################
# PROC
proc lulHandleAtomClusterPicking { Structure AtomNr } {

# segment
       set StringSegment [string trim [.gomselectclusteratoms.frame1.segment get]]
       set StringLenSegment [llength $StringSegment]
       if {$StringSegment == ""} {
          .gomselectclusteratoms.frame1.segment insert 0 "*"
       }
# residue
       set StringResidue [string trim [.gomselectclusteratoms.frame2.residue get]]
       set StringLenResidue [llength $StringResidue]
       if {$StringResidue == ""} {
          .gomselectclusteratoms.frame2.residue insert 0 "*"
       }
# atom
       set StringAtom [string trim [.gomselectclusteratoms.frame3.atom get]]
       set StringLenAtom [llength $StringAtom]
       if {$StringAtom == "" || $StringAtom == "*"} {
          .gomselectclusteratoms.frame3.atom delete 0 end 
          .gomselectclusteratoms.frame3.atom insert 0 $AtomNr 
       } else {
          .gomselectclusteratoms.frame3.atom delete 0 end 
          .gomselectclusteratoms.frame3.atom insert 0 "$StringAtom,$AtomNr" 
       }
}
#1
############################################################################
# PROC
proc lulHandleAtomRDF1Picking { Structure AtomNr } {

# segment
       set StringSegment [string trim [.gomselectrdfatoms.frame1.segment get]]
       set StringLenSegment [llength $StringSegment]
       if {$StringSegment == ""} {
          .gomselectrdfatoms.frame1.segment insert 0 "*"
       }
# residue
       set StringResidue [string trim [.gomselectrdfatoms.frame2.residue get]]
       set StringLenResidue [llength $StringResidue]
       if {$StringResidue == ""} {
          .gomselectrdfatoms.frame2.residue insert 0 "*"
       }
# atom
       set StringAtom [string trim [.gomselectrdfatoms.frame3.atom get]]
       set StringLenAtom [llength $StringAtom]
       if {$StringAtom == "" || $StringAtom == "*"} {
          .gomselectrdfatoms.frame3.atom delete 0 end 
          .gomselectrdfatoms.frame3.atom insert 0 $AtomNr 
       } else {
          .gomselectrdfatoms.frame3.atom delete 0 end 
          .gomselectrdfatoms.frame3.atom insert 0 "$StringAtom,$AtomNr" 
       }
}
#2
############################################################################
# PROC
proc lulHandleAtomRDF2Picking { Structure AtomNr } {

# segment
       set StringSegment [string trim [.gomselectrdfatoms.frame11.segment get]]
       set StringLenSegment [llength $StringSegment]
       if {$StringSegment == ""} {
          .gomselectrdfatoms.frame11.segment insert 0 "*"
       }
# residue
       set StringResidue [string trim [.gomselectrdfatoms.frame12.residue get]]
       set StringLenResidue [llength $StringResidue]
       if {$StringResidue == ""} {
          .gomselectrdfatoms.frame12.residue insert 0 "*"
       }
# atom
       set StringAtom [string trim [.gomselectrdfatoms.frame13.atom get]]
       set StringLenAtom [llength $StringAtom]
       if {$StringAtom == "" || $StringAtom == "*"} {
          .gomselectrdfatoms.frame13.atom delete 0 end 
          .gomselectrdfatoms.frame13.atom insert 0 $AtomNr 
       } else {
          .gomselectrdfatoms.frame13.atom delete 0 end 
          .gomselectrdfatoms.frame13.atom insert 0 "$StringAtom,$AtomNr" 
       }
}

############################################################################
# PROC
proc lulHandleAtomMSDPicking { Structure AtomNr } {

# segment
       set StringSegment [string trim [.gommeansquaredisplacement.frame1.segment get]]
       set StringLenSegment [llength $StringSegment]
       if {$StringSegment == ""} {
          .gommeansquaredisplacement.frame1.segment insert 0 "*"
       }
# residue
       set StringResidue [string trim [.gommeansquaredisplacement.frame2.residue get]]
       set StringLenResidue [llength $StringResidue]
       if {$StringResidue == ""} {
          .gommeansquaredisplacement.frame2.residue insert 0 "*"
       }
# atom
       set StringAtom [string trim [.gommeansquaredisplacement.frame3.atom get]]
       set StringLenAtom [llength $StringAtom]
       if {$StringAtom == "" || $StringAtom == "*"} {
          .gommeansquaredisplacement.frame3.atom delete 0 end 
          .gommeansquaredisplacement.frame3.atom insert 0 $AtomNr 
       } else {
          .gommeansquaredisplacement.frame3.atom delete 0 end 
          .gommeansquaredisplacement.frame3.atom insert 0 "$StringAtom,$AtomNr" 
       }
}

############################################################################
# PROC
proc lulHandleAtomRMSDPicking { Structure AtomNr } {

# segment
       set StringSegment [string trim [.gomplotrmsd.frame1.segment_input get]]
       set StringLenSegment [llength $StringSegment]
       if {$StringSegment == ""} {
          .gomplotrmsd.frame1.segment_input insert 0 "*"
       }
# residue
       set StringResidue [string trim [.gomplotrmsd.frame2.residue_input get]]
       set StringLenResidue [llength $StringResidue]
       if {$StringResidue == ""} {
          .gomplotrmsd.frame2.residue_input insert 0 "*"
       }
# atom
       set StringAtom [string trim [.gomplotrmsd.frame3.atom_input get]]
       set StringLenAtom [llength $StringAtom]
       if {$StringAtom == "" || $StringAtom == "*"} {
          .gomplotrmsd.frame3.atom_input delete 0 end 
          .gomplotrmsd.frame3.atom_input insert 0 $AtomNr 
       } else {
          .gomplotrmsd.frame3.atom_input delete 0 end 
          .gomplotrmsd.frame3.atom_input insert 0 "$StringAtom,$AtomNr" 
       }

}

############################################################################
# PROC
proc lulHandleAtomPlaneXYZPicking { Which Structure AtomNr } {

       switch $Which {

         1 {
             set Coords [show atom coord $AtomNr $Structure]
             .gomcontourplanexyz.frame1.coord.p1xyz delete 0 end
             .gomcontourplanexyz.frame1.coord.p1xyz insert 0 $Coords
           }
         2 {
             set Coords [show atom coord $AtomNr $Structure]
             .gomcontourplanexyz.frame1.coord.p2xyz delete 0 end
             .gomcontourplanexyz.frame1.coord.p2xyz insert 0 $Coords
           }
         3 {
             set Coords [show atom coord $AtomNr $Structure]
             .gomcontourplanexyz.frame1.coord.p3xyz delete 0 end
             .gomcontourplanexyz.frame1.coord.p3xyz insert 0 $Coords
           }
       }
}

############################################################################
# PROC
proc lulHandleAtomClipPlaneXYZPicking { Which Structure AtomNr } {

       switch $Which {

         1 {
             set Coords [show atom coord $AtomNr $Structure]
             .gomcontourclipplanexyz.frame1.coord.p1xyz delete 0 end
             .gomcontourclipplanexyz.frame1.coord.p1xyz insert 0 $Coords
           }
         2 {
             set Coords [show atom coord $AtomNr $Structure]
             .gomcontourclipplanexyz.frame1.coord.p2xyz delete 0 end
             .gomcontourclipplanexyz.frame1.coord.p2xyz insert 0 $Coords
           }
         3 {
             set Coords [show atom coord $AtomNr $Structure]
             .gomcontourclipplanexyz.frame1.coord.p3xyz delete 0 end
             .gomcontourclipplanexyz.frame1.coord.p3xyz insert 0 $Coords
           }
       }
}

########################## end of handle picking tcl file #########################
####################################################################
# 
# plot axis ...
#
proc lulPlotAxis { } {

     global gomHelpDir
     global gomHelpFile
     global gomStructure
     global gomSwitch
     global gomControlFont

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomplotaxis
catch {destroy $w}
toplevel $w 
wm title $w "Plot Axis on Atoms"
wm iconname $w "Plot Axis"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyPlotAxis $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(plotaxis)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1 -padx 3

set  NumStruct [show molstructures]

frame $w.frame1              -borderwidth 2 -relief ridge
label $w.frame1.segmentlabel -text "Segment (x):" -width 10
entry $w.frame1.segment      -width 15
pack  $w.frame1 -side top    -anchor w
pack  $w.frame1.segmentlabel $w.frame1.segment -side left
# default
$w.frame1.segment insert 0 "*"

frame $w.frame2              -borderwidth 2 -relief ridge
label $w.frame2.residuelabel -text "Residue (y):" -width 10
entry $w.frame2.residue      -width 15
pack  $w.frame2 -side top    -anchor w
pack  $w.frame2.residuelabel $w.frame2.residue -side left
# default
$w.frame2.residue insert 0 "*"

frame $w.frame3              -borderwidth 2 -relief ridge
label $w.frame3.atomlabel    -text "Atom (z):" -width 10
entry $w.frame3.atom         -width 15
pack  $w.frame3 -side top    -anchor w
pack  $w.frame3.atomlabel    $w.frame3.atom    -side left
# default
$w.frame3.atom insert 0 "*"

set Segment [$w.frame1.segment get]
set Residue [$w.frame2.residue get]
set Atom    [$w.frame3.atom    get]

frame       $w.frame4 -borderwidth 2 -relief ridge
button      $w.frame4.button -text "Delete" -command {plot -axis; lulInvalidateDisplay}

pack        $w.frame4 -side top -anchor w -pady 4
pack        $w.frame4.button -side top

}

####################################################################
# PROC
proc lulApplyPlotAxis { w } {

     global gomStructure
     global gomSwitch

     set Segment [string trim [$w.frame1.segment get]]
     set Residue [string trim [$w.frame2.residue get]]
     set Atom    [string trim [$w.frame3.atom    get]]

     eval "plot axis \{$Segment\} \{$Residue\} \{$Atom\}"

     lulInvalidateDisplay
}
# 
# Quad Stereo
#
proc lulQuadStereoControl { } {
 
    global gomHelpDir
    global gomHelpFile
    global gomStructure
    global gomControlFont
    global gomQuadStereoStateValue
    global gomControlFont8    
   global gomControlFont12    
 
# return if no molecular systems defined
     if {[show molstructures] < 1} {
          gomError "Read in a structure first"
	      return
	 }
 
set w .gomquadstereo
catch {destroy $w}
toplevel $w 
wm title $w "Quad Stereo control"
wm iconname $w "Quad Stereo"
 
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyQuadStereo $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(quadstereo)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1 -padx 4
pack   $w.buttons.dismiss $w.buttons.apply -side left -expand 1 -padx 4

set SState     [show quadstereo state]
set SAngle     [show quadstereo angle]

label  $w.label -text "Quad Stereo control:" -font "$gomControlFont12"
pack   $w.label -side top -anchor w -pady 4
 
frame $w.frame0 -borderwidth 2 -relief ridge
pack $w.frame0 -side top
 
frame $w.frame1 -borderwidth 2 -relief ridge
pack $w.frame1 -side top
 
label $w.frame1.label    -text "Tilt angle (degrees):" -width 30
entry $w.frame1.angle    -width 12
pack  $w.frame1.label $w.frame1.angle -side left
 
$w.frame1.angle delete 0 end
$w.frame1.angle insert 0 $SAngle
 
frame $w.frame2 -borderwidth 2 -relief ridge
pack $w.frame2 -side top
label       $w.frame2.label -text "Display state: "
radiobutton $w.frame2.on  -text "On"  -width 4 -value 1 -variable gomQuadStereoStateValue
radiobutton $w.frame2.off -text "Off" -width 4 -value 0 -variable gomQuadStereoStateValue
pack $w.frame2.label $w.frame2.on $w.frame2.off -padx 3 -side left
 
if {$SState == "on"} {
  $w.frame2.on  select
  $w.frame2.off deselect
} elseif {$SState == "off"} {
  $w.frame2.on  deselect
  $w.frame2.off select
} else {
  gomError "unknown state '$SState' returned from 'show quadstereo state' command"
}
 
}
 
####################################################################
# PROC
proc lulApplyQuadStereo { w } {
 
    global gomQuadStereoStateValue
 
     set Angle    [string trim [$w.frame1.angle    get]]
     define quadstereo angle $Angle
 
     if {$gomQuadStereoStateValue} {
	  define quadstereo on
	 } else {
      define quadstereo off
	 }
 
     lulInvalidateDisplay
}

#
############################################################################
# PROC
proc lulContourClipPlaneControl {} {

     global gomControlFont
     global gomHelpDir
     global gomHelpFile
     global gomContourID
     global gomClipPlaneX
     global gomClipPlaneY
     global gomClipPlaneZ
     global gomContourClipDirectionX
     global gomContourClipDirectionY
     global gomContourClipDirectionZ

# return if no contours defined
    if {[show contour defined]   < 1} {
        lulErrorDialog {ERROR: no contour available. Read a contour first!}
        return
	}

set w .gomcontourclipplane
catch {destroy $w}
toplevel $w 
wm title $w "Contour Clip Plane Control"
wm iconname $w "Contour Clip Plane Control"
# if this window is closed the ".gomanimatecontclipplane" window
# has to be closed as well!
wm protocol .gomcontourclipplane WM_DELETE_WINDOW  {catch {destroy .gomcontourclipplane; destroy .gomanimatecontclipplane}}

#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w;catch {destroy .gomanimatecontplane}"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyContourClipPlaneCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(contourclipplane)"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

#
    scan [show sizeofsystem] "%f" SizeSystem
    set SizeSystem [expr $SizeSystem / 2.0]
    set Minx  "-$SizeSystem"
    set Maxx  "$SizeSystem" 
    set Miny  "-$SizeSystem" 
    set Maxy  "$SizeSystem" 
    set Minz  "-$SizeSystem" 
    set Maxz  "$SizeSystem" 

    set ClipValueX  0.0
    set ClipValueY  0.0
    set ClipValueZ  0.0

    set gomClipPlaneX 0
    set gomClipPlaneY 0
    set gomClipPlaneZ 0

    frame $w.frame0
	pack  $w.frame0 -side top -anchor w

    frame $w.frame0.left -borderwidth 2 -relief ridge
    label $w.frame0.left.label0 -text "Box\[cube\] dimensions: "
    label $w.frame0.left.label1 -text "Min x: [format %6.3e $Minx], \
            max x: [format %6.3e $Maxx]"
    label $w.frame0.left.label2 -text "Min y: [format %6.3e $Miny], \
            max y: [format %6.3e $Maxy]"
    label $w.frame0.left.label3 -text "Min z: [format %6.3e $Minz], \
            max z: [format %6.3e $Maxz]"
   	pack  $w.frame0.left -side left       -anchor w
	pack  $w.frame0.left.label0 -side top -anchor w
	pack  $w.frame0.left.label1 -side top
	pack  $w.frame0.left.label2 -side top
	pack  $w.frame0.left.label3 -side top
	 
    frame $w.box -borderwidth 2 -relief raised

	pack $w.box   -side top   -anchor n -fill both

    frame       $w.box.frame1 
	label       $w.box.frame1.labelplane   -text "Plane" -width 10
	label       $w.box.frame1.direction    -text "Direction" -width 12
	label       $w.box.frame1.labelcoord   -text "  \[X/Y/Z\]-coord"

	pack        $w.box.frame1              -side top    -anchor w
	pack        $w.box.frame1.labelplane   -side left   -anchor w
	pack        $w.box.frame1.direction    -side left   -anchor w
   	pack        $w.box.frame1.labelcoord   -side left   -anchor w

# 1
	frame       $w.box.frame2
	checkbutton $w.box.frame2.yz -text "YZ-plane" -width 10 \
	            -variable gomClipPlaneX \
                -command "lulHandleClipCheckButton x $w;lulInvalidateDisplay"
    frame        $w.box.frame2.display -borderwidth 2 -relief ridge 
    radiobutton  $w.box.frame2.display.pl    -text "<+>"       \
                -variable gomContourClipDirectionX  \
	            -value "+" -command "lulHandleClipCheckButton x $w;lulInvalidateDisplay"
    radiobutton  $w.box.frame2.display.mi    -text "<->"      \
                -variable gomContourClipDirectionX  \
	            -value "-" -command "lulHandleClipCheckButton x $w;lulInvalidateDisplay"  

	entry       $w.box.frame2.x          -width 10
	pack        $w.box.frame2            -side top
	pack        $w.box.frame2.yz         -side left
    pack        $w.box.frame2.display    -side left -anchor n -padx 4
    pack        $w.box.frame2.display.pl -side top
    pack        $w.box.frame2.display.mi -side top
	pack        $w.box.frame2.x          -side left

    $w.box.frame2.display.pl  select
	$w.box.frame2.display.mi deselect

	$w.box.frame2.x delete 0 end
	$w.box.frame2.x inser  0 [format "%.5f" $ClipValueX]

# 2
	frame       $w.box.frame3
	checkbutton $w.box.frame3.xz -text "XZ-plane" -width 10 \
	            -variable gomClipPlaneY \
                -command "lulHandleClipCheckButton y $w;lulInvalidateDisplay"
    frame        $w.box.frame3.display -borderwidth 2 -relief ridge
    radiobutton  $w.box.frame3.display.pl    -text "<+>"       \
                -variable gomContourClipDirectionY  \
	            -value "+" -command "lulHandleClipCheckButton y $w;lulInvalidateDisplay"
    radiobutton  $w.box.frame3.display.mi    -text "<->"      \
                -variable gomContourClipDirectionY  \
	            -value "-" -command "lulHandleClipCheckButton y $w;lulInvalidateDisplay"
                     
	entry       $w.box.frame3.y  -width 10
	pack        $w.box.frame3            -side top
	pack        $w.box.frame3.xz         -side left
    pack        $w.box.frame3.display    -side left -anchor n -padx 4
    pack        $w.box.frame3.display.pl -side top
    pack        $w.box.frame3.display.mi -side top
	pack        $w.box.frame3.y          -side left

    $w.box.frame3.display.pl  select
	$w.box.frame3.display.mi  deselect

	$w.box.frame3.y delete 0 end
	$w.box.frame3.y inser  0 [format "%.5f" $ClipValueY]

# 3
	frame       $w.box.frame4
	checkbutton $w.box.frame4.xy -text "XY-plane" -width 10 \
	            -variable gomClipPlaneZ \
                -command "lulHandleClipCheckButton z $w;lulInvalidateDisplay"
    frame        $w.box.frame4.display -borderwidth 2 -relief ridge
    radiobutton  $w.box.frame4.display.pl    -text "<+>"       \
                -variable gomContourClipDirectionZ  \
	            -value "+" -command "lulHandleClipCheckButton z $w;lulInvalidateDisplay"
    radiobutton  $w.box.frame4.display.mi    -text "<->"      \
                -variable gomContourClipDirectionZ  \
	            -value "-" -command "lulHandleClipCheckButton z $w;lulInvalidateDisplay"

	entry       $w.box.frame4.z         -width 10
	pack        $w.box.frame4            -side top
	pack        $w.box.frame4.xy         -side left
    pack        $w.box.frame4.display    -side left -anchor n -padx 4
    pack        $w.box.frame4.display.pl -side top
    pack        $w.box.frame4.display.mi -side top
	pack        $w.box.frame4.z          -side left

    $w.box.frame4.display.pl  select
	$w.box.frame4.display.mi  deselect

	$w.box.frame4.z delete 0 end
	$w.box.frame4.z inser  0 [format "%.5f" $ClipValueZ]

# arbitrary plane(s)

	frame       $w.box.frame5           -borderwidth 2 -relief raised
	pack        $w.box.frame5           -side top
#1
    frame       $w.box.frame5.xyz1
	pack        $w.box.frame5.xyz1      -side top

    label       $w.box.frame5.xyz1.label1  -text "Define XYZ plane #1:"
    button      $w.box.frame5.xyz1.button1 -text "Click here" \
                  -command "lulDefineContourClipPlaneXYZ 1 $w"
    pack        $w.box.frame5.xyz1.label1  -side left -anchor w
    pack        $w.box.frame5.xyz1.button1 -side left -anchor w -padx 5
#2
    frame       $w.box.frame5.xyz2
	pack        $w.box.frame5.xyz2      -side top

    label       $w.box.frame5.xyz2.label1  -text "Define XYZ plane #2:"
    button      $w.box.frame5.xyz2.button1 -text "Click here" \
                  -command "lulDefineContourClipPlaneXYZ 2 $w"
    pack        $w.box.frame5.xyz2.label1  -side left -anchor w
    pack        $w.box.frame5.xyz2.button1 -side left -anchor w -padx 5
#3
    frame       $w.box.frame5.xyz3
	pack        $w.box.frame5.xyz3      -side top

    label       $w.box.frame5.xyz3.label1  -text "Define XYZ plane #3:"
    button      $w.box.frame5.xyz3.button1 -text "Click here" \
                  -command "lulDefineContourClipPlaneXYZ 3 $w"
    pack        $w.box.frame5.xyz3.label1  -side left -anchor w
    pack        $w.box.frame5.xyz3.button1 -side left -anchor w -padx 5
    
}

############################################################################
# PROC
proc lulDefineContourClipPlaneXYZ { Which w } {

   global gomClipXYZ1
   global gomClipXYZ2
   global gomClipXYZ3
   global gomControlFont
   global gomContourColour
   global gomHelpDir
   global gomHelpFile

set w .gomcontourclipplanexyz
catch {destroy $w}
toplevel $w 
wm title $w "Contour XYZ Clip Plane Control"
wm iconname $w "Contour XYZ Clip Plane Control"
#
#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w;catch {destroy .gomanimatecontplane}"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyClipplaneXYZ $Which $w;lulInvalidateDisplay"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(contourxyzclipplane)"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Define Contour Clip Plane #$Which:" -font $gomControlFont
pack   $w.label  -side top -anchor w -pady 4

frame  $w.inplabel
pack   $w.inplabel -side top
label  $w.inplabel.text  -text " (ON/OFF)                    (x, y, z)"
pack   $w.inplabel.text  -side left

frame  $w.frame1 -borderwidth 2 -relief ridge
pack   $w.frame1 -side top -anchor w

	checkbutton $w.frame1.xyz -text "XYZ$Which-plane" -width 10 \
	            -variable gomClipXYZ$Which
	pack        $w.frame1.xyz                  -side left

    frame       $w.frame1.coord
	pack        $w.frame1.coord                -side left

	entry       $w.frame1.coord.p1xyz          -width 25
	entry       $w.frame1.coord.p2xyz          -width 25
	entry       $w.frame1.coord.p3xyz          -width 25
	pack        $w.frame1.coord.p1xyz          -side top
	pack        $w.frame1.coord.p2xyz          -side top
	pack        $w.frame1.coord.p3xyz          -side top

}

############################################################################
# PROC
proc lulApplyClipplaneXYZ { Which w } {

   global gomClipXYZ1
   global gomClipXYZ2
   global gomClipXYZ3

          set XYZ1 [string trim [$w.frame1.coord.p1xyz get]]
          if {[string index $XYZ1 0] == "\["} {
              eval set XYZ1 $XYZ1
          }
          set XYZ2 [string trim [$w.frame1.coord.p2xyz get]]
          if {[string index $XYZ2 0] == "\["} {
              eval set XYZ2 $XYZ2
          }
          set XYZ3 [string trim [$w.frame1.coord.p3xyz get]]
          if {[string index $XYZ3 0] == "\["} {
              eval set XYZ3 $XYZ3
          }

   switch $Which {

     1 { 
        if {$gomClipXYZ1} {
          eval contour clipplane xyz1 $XYZ1 $XYZ2 $XYZ3
          contour clipplane xyz1 on
        } else {
          contour clipplane xyz1 off
        }
     }
     2 { 
        if {$gomClipXYZ2} {
          eval contour clipplane xyz2 $XYZ1 $XYZ2 $XYZ3
          contour clipplane xyz2 on
        } else {
          contour clipplane xyz2 off
        }
     }
     3 { 
        if {$gomClipXYZ3} {
          eval contour clipplane xyz3 $XYZ1 $XYZ2 $XYZ3
          contour clipplane xyz3 on
        } else {
          contour clipplane xyz3 off
        }
     }
   }

}

############################################################################
# PROC
proc lulHandleClipCheckButton {Axis w} {

     global gomClipPlaneX
     global gomClipPlaneY
     global gomClipPlaneZ
     global gomContourClipDirectionX
     global gomContourClipDirectionY
     global gomContourClipDirectionZ

if {$Axis == "x"} {

  if {$gomClipPlaneX} {
   set ClipValueX [string trim [$w.box.frame2.x get]]
   if {$gomContourClipDirectionX == "+"} {
       contour clip x+ $ClipValueX;contour clip x on
   } elseif {$gomContourClipDirectionX == "-"} {
       contour clip x- $ClipValueX;contour clip x on
   } else {
       gomError "wrong direction option (allowed +/-)"
       return
   }
  } else {
   contour clip x off
  }
 
} elseif {$Axis == "y"} {

  if {$gomClipPlaneY} {
   set ClipValueY [string trim [$w.box.frame3.y get]]
   if {$gomContourClipDirectionY == "+"} {
       contour clip y+ $ClipValueY;contour clip y on
   } elseif {$gomContourClipDirectionY == "-"} {
       contour clip y- $ClipValueY;contour clip y on
   } else {
       gomError "wrong direction option (allowed +/-)"
       return
   }
  } else {
   contour clip y off
  }

} elseif {$Axis == "z"} {

  if {$gomClipPlaneZ} {
   set ClipValueZ [string trim [$w.box.frame4.z get]]
   if {$gomContourClipDirectionZ == "+"} {
       contour clip z+ $ClipValueZ;contour clip z on
   } elseif {$gomContourClipDirectionZ == "-"} {
       contour clip z- $ClipValueZ;contour clip z on
   } else {
       gomError "wrong direction option (allowed +/-)"
       return
   }
  } else {
   contour clip z off
  }

} else {

    gomError "unknown axis direction '$Axis'"
    return
}


}

############################################################################
proc lulApplyContourClipPlaneCommand { w } {

     global gomClipPlaneX
     global gomClipPlaneY
     global gomClipPlaneZ

# return if no contours defined
    if {[show contour defined]   < 1} {
        return
	}


  if {$gomClipPlaneX} {
   set ClipValueX [string trim [$w.box.frame2.x get]]
   contour clip x+ $ClipValueX;contour clip x on
  } else {
   contour clip x off
  }
 
  if {$gomClipPlaneY} {
   set ClipValueY [string trim [$w.box.frame3.y get]]
   contour clip y+ $ClipValueY;contour clip y on
  } else {
   contour clip y off
  }

  if {$gomClipPlaneZ} {
   set ClipValueZ [string trim [$w.box.frame4.z get]]
   contour clip z+ $ClipValueZ;contour clip z on
  } else {
   contour clip z off
  }

  lulInvalidateDisplay
}
##################### calculate hydrogen bonds ############################
#
# calculate hydrogen bonds
#
# PROC
proc lulCalculateHydrogenBonds {} {

     global gomHelpDir
     global gomHelpFile
     global gomControlFont
     global gomHBonding
     global gomHbondColour
     global gomHBondingHydro

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

set w .gomhbonds
catch {destroy $w}
toplevel $w 
wm title $w "Calculate hydrogen bonds"
wm iconname $w "Calc hyd bonds"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulDoCalculateHydrogenBonds $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(hydrogenbonds)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Calculate Hydrogen Bonds:" -font "$gomControlFont"
pack   $w.label -side top -anchor w

frame  $w.atoms
pack   $w.atoms -side top -anchor w

# Atoms
labelframe  $w.atoms.atoms -borderwidth 2 -relief ridge -text "Atoms (donor or acceptor)"
pack   $w.atoms.atoms -side left -anchor n
lulCreateAtomInputEntries $w.atoms.atoms 16 "*" "*" "*" "" 1

# Subset
labelframe  $w.atoms.subset -borderwidth 2 -relief ridge -text "Search subset"
pack   $w.atoms.subset -side left -anchor n
lulCreateAtomInputEntries $w.atoms.subset 16 "*" "*" "*" "gomHbondSubsetText" 1
frame  $w.atoms.subset.buttons
button $w.atoms.subset.buttons.around -text "Around" \
	-command "lulDoCalculateHydrogenSubsetAround $w"
button $w.atoms.subset.buttons.apply -text "Apply" \
	-command "edit hbsubset \
		\[$w.atoms.subset.segment.input get\] \
		\[$w.atoms.subset.residue.input get\] \
		\[$w.atoms.subset.atom.input get\]"
button $w.atoms.subset.buttons.reset -text "Reset" \
	-command {edit -hbsubset}
pack   $w.atoms.subset.buttons -side top -anchor e
pack   $w.atoms.subset.buttons.around \
	$w.atoms.subset.buttons.apply \
	$w.atoms.subset.buttons.reset \
	-side left -padx 2 -pady 2

frame  $w.left
frame  $w.right
pack   $w.left $w.right -side left -anchor n

# Colour
frame  $w.left.colour -borderwidth 2 -relief ridge
pack   $w.left.colour -side top -anchor w -pady 2 -padx 2
label  $w.left.colour.text -text "Choose colour:   "
pack   $w.left.colour.text -side left  -anchor w
button  $w.left.colour.button -text "Click to change colour" \
	-command "::lulChooseButtonColour $w.left.colour.button \
	    ::gomHbondColour {Choose hbond colour}" \
	-fg [::lulGetVisibleForegroundColour $gomHbondColour] \
	-bg $gomHbondColour
pack   $w.left.colour.button -side left

frame $w.left.hydro -borderwidth 2 -relief ridge
pack  $w.left.hydro -side top -anchor w -pady 2 -padx 2
label $w.left.hydro.text -text "Hydrogen atoms in structure:   "
pack  $w.left.hydro.text -side left  -anchor w
radiobutton $w.left.hydro.yes  -text "Yes"  -value hydro   -variable gomHBondingHydro  
radiobutton $w.left.hydro.no   -text "No"   -value nohydro -variable gomHBondingHydro
pack        $w.left.hydro.yes $w.left.hydro.no -side left -anchor w

$w.left.hydro.yes select

# Remove & criteria
button $w.right.remove -text "Remove hbonds" \
	-command {calculate -hbond;lulInvalidateDisplay}
button $w.right.criteria -text "Change criteria" \
	-command {lulChangeHydrogenBondCriteria}
pack   $w.right.remove $w.right.criteria -side top \
    -anchor w -fill x -padx 2 -pady 4
}
###############################################################
# PROC
proc lulChangeHydrogenBondCriteria { } {
    global gomHelpFile
    global gomControlFont
    global gomHydrogenBondingParams
    global gomHydrogenBondingParamsInGUI

    set w .gomhbondcriteria
    catch {destroy $w}
    toplevel $w
    wm title $w "Change hydrogen bond criteria"
    wm iconname $w "Change hydrogen bond criteria"
    
    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
	-command "destroy $w"
    button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
	-command {array set gomHydrogenBondingParams [
		array get gomHydrogenBondingParamsInGUI]}
    button $w.buttons.help    -text Help     -font "$gomControlFont" \
	-command \
	"htmlShowHelp $gomHelpFile(hbondcriteria)"
    pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

    if {[llength [array names gomHydrogenBondingParamsInGUI]] == 0} {
	array set gomHydrogenBondingParamsInGUI {
		d-a 3.90
		d-a_nohydrogens 3.50
		h-a 2.50
		d-h 1.00
		d-h-a 90
		d-a-aa 90
		h-a-aa 90
		aromatic 20}
    }
    array set gomHydrogenBondingParamsInGUI [array get gomHydrogenBondingParams]

# Distance criteria
    frame $w.distance
    pack  $w.distance -side left -anchor n -padx 4

    labelframe $w.distance.hyd -text "Maximum distances with hydrogens"
    pack       $w.distance.hyd -side top -anchor w

    foreach {entry label} {
	d-a {Donor - acceptor:}
	h-a {Hydrogen - acceptor:}} {
	frame $w.distance.hyd.$entry
	pack  $w.distance.hyd.$entry -side top -anchor w
	label $w.distance.hyd.$entry.label -text $label -width 20 -anchor w
	entry $w.distance.hyd.$entry.value \
	    -textvariable gomHydrogenBondingParamsInGUI($entry) -width 10
	pack  $w.distance.hyd.$entry.label $w.distance.hyd.$entry.value -side left
    }

    labelframe $w.distance.nohyd -text "Maximum and calculatory distances"
    pack       $w.distance.nohyd -side top -anchor w

    foreach {entry label} {
	d-a_nohydrogens {Donor - acceptor:}
	d-h             {Donor - hydrogen:}} {
	frame $w.distance.nohyd.$entry
	pack  $w.distance.nohyd.$entry -side top -anchor w
	label $w.distance.nohyd.$entry.label -text $label -width 20 -anchor w
	entry $w.distance.nohyd.$entry.value \
	    -textvariable gomHydrogenBondingParamsInGUI($entry) -width 10
	pack  $w.distance.nohyd.$entry.label $w.distance.nohyd.$entry.value -side left
    }

# Angle criteria
    frame $w.angle
    pack  $w.angle -side left -anchor n -padx 4

    labelframe $w.angle.min -text "Maximum angles"
    pack       $w.angle.min -side top -anchor w

    foreach {entry label} {
	d-h-a  {Donor - hydrogen - acceptor:}
	d-a-aa {Donor - acceptor - 4th atom:}
	h-a-aa {Hydrogen - acceptor - 4th atom:}} {
	frame $w.angle.min.$entry
	pack  $w.angle.min.$entry -side top -anchor w
	label $w.angle.min.$entry.label -text $label -width 25 -anchor w
	entry $w.angle.min.$entry.value \
	    -textvariable gomHydrogenBondingParamsInGUI($entry) -width 10
	pack  $w.angle.min.$entry.label $w.angle.min.$entry.value -side left
    }

    label $w.angle.min.info1 -text "4th atom is any atom bonded to aceptor."
    label $w.angle.min.info2 -text "Usually it is carbon atom."
    pack  $w.angle.min.info1 $w.angle.min.info2 -side top -anchor w

    labelframe $w.angle.max -text "Maximum angles:"
    pack       $w.angle.max -side top -anchor w
    label $w.angle.max.label -text "Aromatic angle:" -width 25 -anchor w
    entry $w.angle.max.value \
	    -textvariable gomHydrogenBondingParamsInGUI(aromatic) -width 10
    pack  $w.angle.max.label $w.angle.max.value -side left

#   Maximum distances:
#     d-a              Distance between donor and acceptor
#     d-a_nohydrogens  Distance between donor and acceptor
#                      (used then "nohy$drogens" is applied)
#     h-a              Distance between hydrogen and acceptor
#     d-h              Distance between donor and hydrogen
#   Minimum angles:
#     d-h-a            donor-hydrogen-acceptor
#     d-a-aa           donor-acceptor-&lt;atom bonded to acceptor&gt;
#     h-a-aa           hydrogen-acceptor-&lt;atom bonded to acceptor&gt;
#   Maximum angles:
#     aromatic         Angle between aromatic plane and vector donor-hydrogen
}
###############################################################
# PROC
proc lulDoCalculateHydrogenBonds { w } {

    global gomHbondColour
    global gomHBondingHydro
    
    foreach level {segment residue atom} {
	set ${level} [string trim [$w.atoms.atoms.${level}.input get]]
	if {[expr \$${level}] == ""} {
	    set ${level} "*"
	}
    }

    calculate hbond \
	$segment    $residue    $atom \
	[lulColourHex2Float $gomHbondColour] $gomHBondingHydro

    lulInvalidateDisplay
}
##################### plot vector file ############################
# PROC
proc lulDoCalculateHydrogenSubsetAround { w } {
    # Get current atoms
    set Segment [$w.atoms.atoms.segment.input get]
    set Residue [$w.atoms.atoms.residue.input get]
    set Atom    [$w.atoms.atoms.atom.input    get]

    # Get Current subset
    set SubsetSegment [$w.atoms.subset.segment.input get]
    set SubsetResidue [$w.atoms.subset.residue.input get]
    set SubsetAtom    [$w.atoms.subset.atom.input    get]
    if { "" == $SubsetAtom || "*" == $SubsetAtom } {
	set SubsetAtom {*N*,*O*}
    }

    # Find atom around current atoms
    set NewSubsetSegment [find segment around \
	$SubsetSegment $SubsetResidue $SubsetAtom 5.0 \
	$Segment $Residue $Atom all]
    set NewSubsetResidue [find residue around \
	$SubsetSegment $SubsetResidue $SubsetAtom 5.0 \
	$Segment $Residue $Atom all]

    # Create new subset
    $w.atoms.subset.segment.input delete 0 end
    $w.atoms.subset.segment.input insert 0 [find segment join \
	$Segment $Residue $Atom \
	$NewSubsetSegment $NewSubsetResidue $SubsetAtom \
	all]
    $w.atoms.subset.residue.input delete 0 end
    $w.atoms.subset.residue.input insert 0 [find residue join \
	$Segment $Residue $Atom \
	$NewSubsetSegment $NewSubsetResidue $SubsetAtom \
	all]
    $w.atoms.subset.atom.input    delete 0 end
    $w.atoms.subset.atom.input    insert 0 $SubsetAtom
}

##################### plot vector file ############################
# PROC
proc lulPlotVectorFile {} {

     global gomHelpDir
     global gomHelpFile
     global gomControlFont
     global gomVectorDefaultType
     global gomCurrentVectorDisplayStatus


set w .gomplotvectorfile
catch {destroy $w}
toplevel $w 
wm title $w "Plot Vector File"
wm iconname $w "Plot"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
        "htmlShowHelp $gomHelpFile(plotvectorfile)"
pack   $w.buttons.dismiss $w.buttons.help -side left -expand 1

label  $w.label  -text "Plot vector file:" -font $gomControlFont
pack   $w.label  -side top -anchor w

frame  $w.frame0
pack   $w.frame0 -side top -anchor w

label  $w.frame0.label     -text "Vector input text file name:" -width 25
entry  $w.frame0.filename  -width 40
button $w.frame0.browse    -text "Browse..." -command "lulDoPlotVectorFileInput $w;lulUpdateVectorDisplay $w"
pack   $w.frame0 -side top -anchor w
pack   $w.frame0.label $w.frame0.filename -side left -anchor w
pack   $w.frame0.browse -side left -anchor w -padx 4

frame       $w.frame1
pack        $w.frame1 -side top -anchor w
label       $w.frame1.label -text "Vector type:"
radiobutton $w.frame1.line  -text "Line"  -command "unset gomVectorDefaultType;lulInvalidateDisplay"  \
             -value 0 -variable lulJunkID1
radiobutton $w.frame1.solid -text "Solid" -command "lulSetSolidVectorDefaultType $w" \
             -value 1 -variable lulJunkID1
label       $w.frame1.lrad  -text "Radius: "
entry       $w.frame1.rad   -width 15
pack        $w.frame1 -side top -pady 4 -anchor w 
pack        $w.frame1.label $w.frame1.line $w.frame1.solid -side left -anchor w
pack        $w.frame1.lrad -side left -anchor w -padx 4
pack        $w.frame1.rad  -side left -anchor w -padx 4

$w.frame1.rad insert 0 "0.1"

if {[info exist gomVectorDefaultType]} {
  if {[lindex $gomVectorDefaultType 0]} {
    $w.frame1.solid select
  } else { 
    $w.frame1.line  select
  }
} else {
    $w.frame1.line  select
}

labelframe  $w.frame2 -text "Parameter control" -relief ridge -padx 2 -pady 2
pack        $w.frame2 -side top -anchor w
label       $w.frame2.label  -text "Display range: (min) "
entry       $w.frame2.min    -width 12
label       $w.frame2.label1 -text " (max) "
entry       $w.frame2.max    -width 12 
label       $w.frame2.label2 -text " Scaling: "
entry       $w.frame2.scale  -width 12
button      $w.frame2.button -text "Apply" -command "lulSetPlotVectorRange $w"
pack        $w.frame2 -side top -pady 4 -anchor w 
pack        $w.frame2.label $w.frame2.min $w.frame2.label1 $w.frame2.max \
             $w.frame2.label2 $w.frame2.scale -side left -anchor w
pack        $w.frame2.button -side left -anchor w -padx 4

set RValue [show vector]

$w.frame2.min insert 0 [lindex RValue 2]
$w.frame2.max insert 0 [lindex RValue 3]

labelframe        $w.frame3 -text "Display state" -borderwidth 2 -relief ridge -padx 2 -pady 2
#label        $w.frame3.label -text "Display state: "
pack         $w.frame3 -side top -anchor w -padx 4 -pady 6
#pack         $w.frame3.label -side left -anchor w

radiobutton  $w.frame3.on  -text "On"  -value 1 -variable gomCurrentVectorDisplayStatus \
              -command {plot vector on;lulInvalidateDisplay} 
radiobutton  $w.frame3.off -text "Off" -value 0 -variable gomCurrentVectorDisplayStatus \
              -command {plot vector off;lulInvalidateDisplay} 
pack         $w.frame3.on  -side left -anchor w
pack         $w.frame3.off -side left -anchor w

if {$gomCurrentVectorDisplayStatus} {
   $w.frame3.on  select
} else {
   $w.frame3.off select
}

frame        $w.frame4 -borderwidth 2 -relief raised
label        $w.frame4.label -text "Delete all vectors and space: "
button       $w.frame4.button -text "Click here" -command "plot -vector;destroy $w;lulInvalidateDisplay"
pack         $w.frame4 -side top -anchor w -padx 4 -pady 6
pack         $w.frame4.label $w.frame4.button -side left -anchor w -padx 4 -pady 4

}

############################################################################
# PROC
proc lulUpdateVectorDisplay { w } {

set RValue [show vector]
puts $RValue
$w.frame2.scale insert 0 [lindex $RValue 0]
$w.frame2.min   insert 0 [lindex $RValue 2]
$w.frame2.max   insert 0 [lindex $RValue 3]

}
############################################################################
# PROC
proc lulSetSolidVectorDefaultType { w } {

     global gomVectorDefaultType

set gomVectorDefaultType "1 [string trim [$w.frame1.rad get]]"
lulInvalidateDisplay

}

############################################################################
# PROC
proc lulSetPlotVectorRange { w } {

set min [string trim [$w.frame2.min get]]
set max [string trim [$w.frame2.max get]]

if {$min == ""} {
   gomError "min value missing"
   return
}

plot vector range $min $max

set pscale [string trim [$w.frame2.scale get]]

if {$pscale == ""} {
   gomError "scale value missing"
   return
}

plot vector scale $pscale

lulInvalidateDisplay

}

##################### lulDoPlotVectorFileInput ############################
#
# PROC
proc lulDoPlotVectorFileInput { w } {

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set types {
	{"Plot file"		{.txt .TXT}		TEXT}
	{"All files"		*}
    }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension .plt -title "Vector file"]]

    if [string compare $file ""] {

        $w.frame0.filename delete 0 end
	    $w.frame0.filename insert 0 $file

        import vector flatfile $file
# change to current directory
        lulChangeDirectory "$file"

    }
}

##################### lulRunPlotVectorFileCommand ############################
#
# PROC
proc lulRunPlotVectorFileCommand { w } {

     global env
     global gomOutFileAction
	 global gomEnv

set InputFile1  [string trim [$w.frame0.filename get]]
set InputFile2  [string trim [$w.frame4.filename get]]
set OutputFile [string trim [$w.frame1.filename get]]

if {$InputFile1 == "" || $InputFile2 == "" || $OutputFile == ""} {

   gomError "either input1: '$InputFile1', input2: '$InputFile2'  or output: '$OutputFile' file name is missing"
   return
}

# check to see that output is not a directory
  if {[file isdirectory $OutputFile]} {
	    lulErrorDialog "ERROR - output is a directory"
		return
  }

  if {$gomOutFileAction == 0} {
     set Program    "contman -i$InputFile1-$InputFile2 -o$OutputFile"
  } elseif {$gomOutFileAction == 1} {
     set Program    "contman -i$InputFile1+$InputFile2 -o$OutputFile"
  } else {
     gomError "wrong action type: '$gomOutFileAction'. Allowed values are 0 or 1"
     return
  }

set Command [file join $gomEnv(GOM_BIN) $Program]

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    $w.frame3.text insert $i.0 "$Text\n"
    update idletasks
    puts $Text
    incr i
  }

  catch {close $f}

}
#############################################################
# PROC
proc lulCorrectDistanceListWidget { } {

set w .gommonitordistance

if {[winfo exists $w]} {

set DefDist [string trim [show monitor distance list]]

if {[llength $DefDist] < 2} return

set DistValues [lindex $DefDist 0]

# delete old list first
     for {set i 1} {$i <= $DistValues} {incr i} {
       catch {destroy $w.distmon.text.frame$i}
     }
# now build new list
     for {set i 1} {$i <= $DistValues} {incr i} {

        set Type  [string trim [show monitor distance type $i]]
        set Color [show monitor distance colour $i]
        scan $Color "%f %f %f" Red Green Blue

        frame  $w.distmon.text.frame$i
		entry  $w.distmon.text.frame$i.entry$i -width 34

        set AtomIndex   [lindex $DefDist $i]
		set AtomIndex1  [lindex $AtomIndex 0]
		set AtomIndex2  [lindex $AtomIndex 1]

        set AtomInfo "\[$i\]\
		              ([show atom segmentname $AtomIndex1 1]:\
		               [show atom residuename $AtomIndex1 1]:\
					   [show atom atomname    $AtomIndex1 1]($AtomIndex1)) - \
					  ([show atom segmentname $AtomIndex2 1]:\
					   [show atom residuename $AtomIndex2 1]:\
					   [show atom atomname    $AtomIndex2 1]($AtomIndex2)) "
        $w.distmon.text.frame$i.entry$i delete 0 end
		$w.distmon.text.frame$i.entry$i insert 0 $AtomInfo

        menubutton  $w.distmon.text.frame$i.type$i  -text "Line type"   \
                    -menu $w.distmon.text.frame$i.type$i.menu                       \
					-cursor hand2                                       \
                    -borderwidth 2 -relief raised -padx 2
        menu        $w.distmon.text.frame$i.type$i.menu
        set m       $w.distmon.text.frame$i.type$i.menu
        $m add radiobutton -label "* * * *"   -variable gomDistMonLineType -value 1 \
           -command "monitor distance type 1 $i;lulInvalidateDisplay"
        $m add radiobutton -label "- - - -"   -variable gomDistMonLineType -value 2 \
           -command "monitor distance type 2 $i;lulInvalidateDisplay"
        $m add radiobutton -label "* - - *"   -variable gomDistMonLineType -value 3 \
           -command "monitor distance type 3 $i;lulInvalidateDisplay"
        $m add radiobutton -label "-------"   -variable gomDistMonLineType -value 4 \
           -command "monitor distance type 4 $i;lulInvalidateDisplay"
        $m invoke [show monitor distance type $i]

		button $w.distmon.text.frame$i.color$i  -text "Colour"    \
		       -cursor hand2                                      \
			   -bg [lulColourFloat2Hex $Red $Green $Blue] -padx 2 \
			   -command "lulChangeMonitorDistanceLineColor $i $w.distmon.text.frame$i.color$i;lulInvalidateDisplay"  

        button $w.distmon.text.frame$i.edit$i -text "Edit"        \
               -cursor hand2                                      \
               -command "lulMonitorEditDistance $i"

		pack   $w.distmon.text.frame$i.entry$i -side left -padx 3
		pack   $w.distmon.text.frame$i.type$i  -side left
		pack   $w.distmon.text.frame$i.color$i -side left
        pack   $w.distmon.text.frame$i.edit$i  -side left
        $w.distmon.text window create $i.0 -window $w.distmon.text.frame$i
    }

}
}
#############################################################
# PROC
proc lulMonitorEditDistance { which } {

     global gomHelpDir
     global gomHelpFile
     global gomControlFont
     global gomAtomPickingTargetWidgets

set w .gommonitoreditdistance
catch {destroy $w}
toplevel $w 
wm title $w "Monitor Edit Distance"
wm iconname $w "Monitor Edit Distance"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss    -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply      -font "$gomControlFont" \
        -command "lulDoMonitorEditDistance $w $which"
button $w.buttons.help    -text Help       -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(monitordistance)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Edit monitor distance:"
pack   $w.label -side top -anchor w

frame  $w.atoms
pack   $w.atoms -side top -anchor w

set DefList   [string trim [show monitor distance list]]
set DefList1  [lindex [lindex $DefList $which ] 0] 
set DefList2  [lindex [lindex $DefList $which ] 1] 

# set #1,#2
foreach set {1 2} {
    frame  $w.atoms.set$set -borderwidth 2 -relief ridge
    pack   $w.atoms.set$set -side left
    label  $w.atoms.set$set.label -text "Atom #$set:" -width 10
    pack   $w.atoms.set$set.label -side left 
    entry  $w.atoms.set$set.atom_input -width 15 
    pack   $w.atoms.set$set.atom_input -side left

    $w.atoms.set$set.atom_input delete 0 end
    $w.atoms.set$set.atom_input insert 0 [expr \$DefList$set]

    set gomAtomPickingTargetWidgets($w.atoms.set$set) {{} {} atom_input 0}
}

frame       $w.options 
pack        $w.options -side top -anchor w -fill x

frame       $w.options.left
pack        $w.options.left   -side left -anchor w

#frame       $w.options.right 
#pack        $w.options.right  -side right -anchor e

frame       $w.options.left.frame8 -borderwidth 2 -relief ridge
button      $w.options.left.frame8.delete -text "Delete" \
            -command "monitor delete distance $which;destroy .gommonitoreditdistance;lulInvalidateDisplay"

pack        $w.options.left.frame8 -side bottom -anchor w
pack        $w.options.left.frame8.delete -side left

}

#############################################################
# PROC
proc lulDoMonitorEditDistance { w which } {

   set Atom1 [string trim [$w.atoms.set1.atom_input get]]
   set Atom2 [string trim [$w.atoms.set2.atom_input get]]

   if {$Atom1 != "" && $Atom2 != "" } {
     monitor edit distance $which $Atom1 $Atom2 -1 -1.0 -1.0 -1.0
   }

   lulInvalidateDisplay
}

#############################################################
# PROC
proc lulCorrectAngleListWidget { } {

set w .gommonitorangle

if {[winfo exists $w]} {

set DefAng [string trim [show monitor angle list]]

if {[llength $DefAng] < 2} return

set AngValues [lindex $DefAng 0]

# delete old list first
     for {set i 1} {$i <= $AngValues} {incr i} {
       catch {destroy $w.angmon.text.frame$i}
     }

     for {set i 1} {$i <= $AngValues} {incr i} {

        set Color [show monitor angle colour $i]
        scan $Color "%f %f %f" Red Green Blue

        frame  $w.angmon.text.frame$i
		entry  $w.angmon.text.frame$i.entry$i -width 34

        set AtomIndex   [lindex $DefAng $i]
		set AtomIndex1  [lindex $AtomIndex 0]
		set AtomIndex2  [lindex $AtomIndex 1]
		set AtomIndex3  [lindex $AtomIndex 2]

        set AtomInfo "\[$i\]\
		              ([show atom segmentname $AtomIndex1 1]:\
		               [show atom residuename $AtomIndex1 1]:\
					   [show atom atomname    $AtomIndex1 1]($AtomIndex1)) - \
					  ([show atom segmentname $AtomIndex2 1]:\
					   [show atom residuename $AtomIndex2 1]:\
					   [show atom atomname    $AtomIndex2 1]($AtomIndex2)) - \
					  ([show atom segmentname $AtomIndex3 1]:\
					   [show atom residuename $AtomIndex3 1]:\
					   [show atom atomname    $AtomIndex3 1]($AtomIndex3))"
        $w.angmon.text.frame$i.entry$i delete 0 end
		$w.angmon.text.frame$i.entry$i insert 0 $AtomInfo

        menubutton  $w.angmon.text.frame$i.type$i  -text "Line type"   \
                    -menu $w.angmon.text.frame$i.type$i.menu                       \
					-cursor hand2                                       \
                    -borderwidth 2 -relief raised -padx 2
        menu        $w.angmon.text.frame$i.type$i.menu
        set m       $w.angmon.text.frame$i.type$i.menu
        $m add radiobutton -label "* * * *"   -variable gomAngMonLineType -value 1 \
           -command "monitor angle type 1 $i;lulInvalidateDisplay"
        $m add radiobutton -label "- - - -"   -variable gomAngMonLineType -value 2 \
           -command "monitor angle type 2 $i;lulInvalidateDisplay"
        $m add radiobutton -label "* - - *"   -variable gomAngMonLineType -value 3 \
           -command "monitor angle type 3 $i;lulInvalidateDisplay"
        $m add radiobutton -label "-------"   -variable gomAngMonLineType -value 4 \
           -command "monitor angle type 4 $i;lulInvalidateDisplay"
        $m invoke [show monitor angle type $i]

		button $w.angmon.text.frame$i.color$i  -text "Colour"    \
		       -cursor hand2                                      \
			   -bg [lulColourFloat2Hex $Red $Green $Blue] -padx 2 \
			   -command "lulChangeMonitorAngleLineColor $i $w.angmon.text.frame$i.color$i;lulInvalidateDisplay"  

        button $w.angmon.text.frame$i.edit$i -text "Edit"        \
               -cursor hand2                                      \
               -command "lulMonitorEditAngle $i"
		pack   $w.angmon.text.frame$i.entry$i -side left
		pack   $w.angmon.text.frame$i.type$i  -side left -padx 3
		pack   $w.angmon.text.frame$i.color$i -side left
        pack   $w.angmon.text.frame$i.edit$i  -side left
        $w.angmon.text window create $i.0 -window $w.angmon.text.frame$i
     }

}

}
#############################################################
# PROC
proc lulMonitorEditAngle { which } {

     global gomHelpDir
     global gomHelpFile
     global gomControlFont

set w .gommonitoreditangle
catch {destroy $w}
toplevel $w 
wm title $w "Monitor Edit Angle"
wm iconname $w "Monitor Edit Angle"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss    -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply      -font "$gomControlFont" \
        -command "lulDoMonitorEditAngle $w $which"
button $w.buttons.help    -text Help       -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(monitorangle)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Edit monitor angle:"
pack   $w.label -side top -anchor w

frame  $w.atoms
pack   $w.atoms -side top -anchor w

set DefList   [string trim [show monitor angle list]]
set DefList1  [lindex [lindex $DefList $which ] 0] 
set DefList2  [lindex [lindex $DefList $which ] 1] 
set DefList3  [lindex [lindex $DefList $which ] 2] 

# set #1,#2,#3
foreach set {1 2 3} {
    frame  $w.atoms.set$set -borderwidth 2 -relief ridge
    pack   $w.atoms.set$set -side left
    label  $w.atoms.set$set.label -text "Atom #$set:" -width 10
    pack   $w.atoms.set$set.label -side left 
    entry  $w.atoms.set$set.atom_input -width 15 
    pack   $w.atoms.set$set.atom_input -side left

    $w.atoms.set$set.atom_input delete 0 end
    $w.atoms.set$set.atom_input insert 0 [expr \$DefList$set]

    set gomAtomPickingTargetWidgets($w.atoms.set$set) {{} {} atom_input 0}
}

frame       $w.options 
pack        $w.options -side top -anchor w -fill x

frame       $w.options.left
pack        $w.options.left   -side left -anchor w

#frame       $w.options.right 
#pack        $w.options.right  -side right -anchor e

frame       $w.options.left.frame8 -borderwidth 2 -relief ridge
button      $w.options.left.frame8.delete -text "Delete" \
            -command "monitor delete angle $which;destroy .gommonitoreditangle;lulInvalidateDisplay"

pack        $w.options.left.frame8 -side bottom -anchor w
pack        $w.options.left.frame8.delete -side left

}

#############################################################
# PROC
proc lulDoMonitorEditAngle { w which } {

   set Atom1 [string trim [$w.atoms.set1.atom_input get]]
   set Atom2 [string trim [$w.atoms.set2.atom_input get]]
   set Atom3 [string trim [$w.atoms.set3.atom_input get]]

   if {$Atom1 != "" && $Atom2 != "" && $Atom3 != "" } {
     monitor edit angle $which $Atom1 $Atom2 $Atom3 -1 -1.0 -1.0 -1.0
   }

   lulInvalidateDisplay
}


#############################################################
# PROC
proc lulCorrectTorsionListWidget { } {

set w .gommonitortorsion

if {[winfo exists $w]} {

set DefTors [string trim [show monitor torsion list]]

if {[llength $DefTors] < 2} return

set TorsValues [lindex $DefTors 0]

# delete old list first
     for {set i 1} {$i <= $TorsValues} {incr i} {
       catch {destroy $w.angmon.text.frame$i}
     }

     for {set i 1} {$i <= $TorsValues} {incr i} {

        set Color [show monitor torsion colour $i]
        scan $Color "%f %f %f" Red Green Blue

        frame  $w.torsmon.text.frame$i
		entry  $w.torsmon.text.frame$i.entry$i -width 34

        set AtomIndex   [lindex $DefTors $i]
		set AtomIndex1  [lindex $AtomIndex 0]
		set AtomIndex2  [lindex $AtomIndex 1]
		set AtomIndex3  [lindex $AtomIndex 2]
		set AtomIndex4  [lindex $AtomIndex 3]

        set AtomInfo "\[$i\]\
		              ([show atom segmentname $AtomIndex1 1]:\
		               [show atom residuename $AtomIndex1 1]:\
					   [show atom atomname    $AtomIndex1 1]($AtomIndex1)) - \
					  ([show atom segmentname $AtomIndex2 1]:\
					   [show atom residuename $AtomIndex2 1]:\
					   [show atom atomname    $AtomIndex2 1]($AtomIndex2)) - \
					  ([show atom segmentname $AtomIndex3 1]:\
					   [show atom residuename $AtomIndex3 1]:\
					   [show atom atomname    $AtomIndex3 1]($AtomIndex3)) - \
					  ([show atom segmentname $AtomIndex4 1]:\
					   [show atom residuename $AtomIndex4 1]:\
					   [show atom atomname    $AtomIndex4 1]($AtomIndex4))"
        $w.torsmon.text.frame$i.entry$i delete 0 end
		$w.torsmon.text.frame$i.entry$i insert 0 $AtomInfo

        menubutton  $w.torsmon.text.frame$i.type$i  -text "Line type"   \
                    -menu $w.torsmon.text.frame$i.type$i.menu                       \
					-cursor hand2                                       \
                    -borderwidth 2 -relief raised -padx 2
        menu        $w.torsmon.text.frame$i.type$i.menu
        set m       $w.torsmon.text.frame$i.type$i.menu
        $m add radiobutton -label "* * * *"   -variable gomTorsMonLineType -value 1 \
           -command "monitor torsion type 1 $i;lulInvalidateDisplay"
        $m add radiobutton -label "- - - -"   -variable gomTorsMonLineType -value 2 \
           -command "monitor torsion type 2 $i;lulInvalidateDisplay"
        $m add radiobutton -label "* - - *"   -variable gomTorsMonLineType -value 3 \
           -command "monitor torsion type 3 $i;lulInvalidateDisplay"
        $m add radiobutton -label "-------"   -variable gomTorsMonLineType -value 4 \
           -command "monitor torsion type 4 $i;lulInvalidateDisplay"
        $m invoke [show monitor torsion type $i]

		button $w.torsmon.text.frame$i.color$i  -text "Colour"    \
		       -cursor hand2                                      \
			   -bg [lulColourFloat2Hex $Red $Green $Blue] -padx 2 \
			   -command "lulChangeMonitorTorsionLineColor $i $w.torsmon.text.frame$i.color$i;lulInvalidateDisplay"  

        button $w.torsmon.text.frame$i.edit$i -text "Edit"        \
               -cursor hand2                                      \
               -command "lulMonitorEditTorsion $i"
		pack   $w.torsmon.text.frame$i.entry$i -side left
		pack   $w.torsmon.text.frame$i.type$i  -side left -padx 3
		pack   $w.torsmon.text.frame$i.color$i -side left
        pack   $w.torsmon.text.frame$i.edit$i  -side left
        $w.torsmon.text window create $i.0 -window $w.torsmon.text.frame$i
     }

}

}
#############################################################
# PROC
proc lulMonitorEditTorsion { which } {

     global gomHelpDir
     global gomHelpFile
     global gomControlFont

set w .gommonitoredittorsion
catch {destroy $w}
toplevel $w 
wm title $w "Monitor Edit Torsion"
wm iconname $w "Monitor Edit Torsion"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss    -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply      -font "$gomControlFont" \
        -command "lulDoMonitorEditTorsion $w $which"
button $w.buttons.help    -text Help       -font "$gomControlFont" \
        -command \
         "htmlShowHelp $gomHelpFile(monitortorsion)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Edit monitor torsion:"
pack   $w.label -side top -anchor w

frame  $w.atoms
pack   $w.atoms -side top -anchor w

set DefList   [string trim [show monitor torsion list]]
set DefList1  [lindex [lindex $DefList $which ] 0] 
set DefList2  [lindex [lindex $DefList $which ] 1]
set DefList3  [lindex [lindex $DefList $which ] 2] 
set DefList4  [lindex [lindex $DefList $which ] 3] 

# set #1,#2,#3,#4
foreach set {1 2 3 4} {
    frame  $w.atoms.set$set -borderwidth 2 -relief ridge
    pack   $w.atoms.set$set -side left
    label  $w.atoms.set$set.label -text "Atom #$set:" -width 10
    pack   $w.atoms.set$set.label -side left 
    entry  $w.atoms.set$set.atom_input -width 15 
    pack   $w.atoms.set$set.atom_input -side left

    $w.atoms.set$set.atom_input delete 0 end
    $w.atoms.set$set.atom_input insert 0 [expr \$DefList$set]

    set gomAtomPickingTargetWidgets($w.atoms.set$set) {{} {} atom_input 0}
}

frame       $w.options 
pack        $w.options -side top -anchor w -fill x

frame       $w.options.left
pack        $w.options.left   -side left -anchor w

#frame       $w.options.right 
#pack        $w.options.right  -side right -anchor e

frame       $w.options.left.frame8 -borderwidth 2 -relief ridge
button      $w.options.left.frame8.delete -text "Delete" \
            -command "monitor delete torsion $which;destroy .gommonitoredittorsion;lulInvalidateDisplay"

pack        $w.options.left.frame8 -side bottom -anchor w
pack        $w.options.left.frame8.delete -side left

}

#############################################################
# PROC
proc lulDoMonitorEditTorsion { w which } {

   set Atom1 [string trim [$w.atoms.set1.atom_input get]]
   set Atom2 [string trim [$w.atoms.set2.atom_input get]]
   set Atom3 [string trim [$w.atoms.set3.atom_input get]]
   set Atom4 [string trim [$w.atoms.set4.atom_input get]]

   if {$Atom1 != "" && $Atom2 != "" && $Atom3 != "" && $Atom4 != ""} {
     monitor edit torsion $which $Atom1 $Atom2 $Atom3 $Atom4 -1 -1.0 -1.0 -1.0
   }

   lulInvalidateDisplay
}



#############################################################
# PROC
proc lulPrepareAtomLabel { i j option} {

   global gomAtomLabelStringDef
   global gomAtomLabelStringExt

   if {$option} {
# sprintf(text,"%s:%d:%s(%d)",GetSegName(i,j),res1[j],GetAtmName(i,j),(j+1));
     eval define gtext "$gomAtomLabelStringExt"
   } else {
# sprintf(text,"%s",GetAtmName(i,j));
     eval define gtext "$gomAtomLabelStringDef" 
   }
}
###############################################################################
#
# special yes/no dialog widget
#
proc lulYesNoMessageDialog {InputText} {

after idle {.dialogyn.msg configure -wraplength 4i}
set i [tk_dialog .dialogyn "gOpenMol Message" $InputText \
info 0 "Yes" "No" ]

}

##################### Socket server/client control ############################
#
# define socket server/client widgets
#
proc lulSocketServerControl {} {

     global gomHelpDir
     global gomHelpFile
     global gomControlFont
     global gomServerState
     global gomSocketPort
     global gomAllowConnection

set w .gomsocketservercontrol
catch {destroy $w}
toplevel $w 
wm title $w "Socket Control"
wm iconname $w "Socket Control"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulDoSocketServerControl $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(atomtype)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Socket server/client control:"
pack   $w.label -side top -anchor w

labelframe  $w.frame1 -text "Allowed host(s):" -padx 2 -pady 2
entry  $w.frame1.host_input -width 30 -textvariable gomHostName 
button $w.frame1.host_allowed -text "Apply" -command "gomUpdateAllowedConnections $w"
pack   $w.frame1 -side top -anchor w
pack   $w.frame1.host_input -side left
pack   $w.frame1.host_allowed -side left -padx 3

if { [$w.frame1.host_input get] == "" } {$w.frame1.host_input insert 0 $gomAllowConnection}

labelframe  $w.frame2 -text "Socket port:" -padx 2 -pady 2
entry  $w.frame2.socket_input -width 30 -textvariable gomSocketPort 
pack   $w.frame2 -side top -anchor w
pack   $w.frame2.socket_input -side left 

if { [$w.frame2.socket_input get] == ""} {$w.frame2.socket_input insert 0 "*"}

labelframe  $w.frame4 -borderwidth 2 -relief ridge -text "Socket state" -padx 2 -pady 2
radiobutton $w.frame4.on  -text "On"  -command "gomEnableSocketServer $w" -variable gomServerState -value 1
radiobutton $w.frame4.off -text "Off" -command gomDisableSocketServer  -variable gomServerState -value 0
pack $w.frame4 -side left 
pack $w.frame4.on $w.frame4.off -side top -anchor w
if {$gomServerState} {
  $w.frame4.on  select
} else {
  $w.frame4.off select
}
}
###############################################################################
#
proc gomUpdateAllowedConnections { w } {
   global   gomAllowConnection

set gomAllowConnection [$w.frame1.host_input get]
   
}
###############################################################################
#
proc lulDoSocketServerControl { w } {
}
###############################################################################
#
proc gomEnableSocketServer { w } {

set socket [$w.frame2.socket_input get]

lulSocketComm::Server $socket

}
###############################################################################
#
proc gomDisableSocketServer { } {

lulSocketComm::Server 0

}

#############################################################
# 3D ==> 2D
#source [file join $gomEnv(GOM_DATA) mol3dto2d.tcl]
#############################################################
# Make Movie
#source [file join $gomEnv(GOM_DATA) make_movie.tcl]
#############################################################
# VRML filter
#source [file join $gomEnv(GOM_DATA) vrml.tcl]
#
# small utility routine to plot a truncated octahedron in
# an arbitrary place in space, and another to place one at (0,0,0)
# given the "side-to-side" cell dimension.
#
#source [file join $gomEnv(GOM_DATA) plotocta.tcl]
#source [file join $gomEnv(GOM_DATA) atomtree.tcl]
#source [file join $gomEnv(GOM_DATA) trajectory.tcl]

#############################################################
#
# source automatically all files in the autotcl directory
#
#puts "Enabling all tcl files in the '[file join $gomEnv(GOM_DATA) autotcl gui]' directory"
#lappend auto_path [file join $gomEnv(GOM_DATA) autotcl gui]
#puts "Sourcing all tcl files in the '[file join $gomEnv(GOM_DATA) inittcl gui]' directory ..."
#foreach i [glob -directory [file join $gomEnv(GOM_DATA) inittcl gui] *.tcl] {
#    gomPrint "    $i"
#    catch {source "$i"} errRet
#    if {"" != $errRet} {
#	gomError "Parse error in file $i"
#    }
#}
#
#
#
# help system
#
#puts "HTML based help system..."
# source the help function
#source [file join $gomEnv(GOM_DATA) html_help.tcl]
#
