##############################################################################
#                       Copyright (c) 1994 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Enhancements 2003 - 2004 by: Eero Häkkinen
##############################################################################

namespace eval ::gom::gui::Contour {

############################################################################
# PROC
proc Control {} {
    package require gom::gui::Widgets
    namespace import -force ::gom::gui::Widgets::*

    variable displayMethod
    set w .gomcontour
    if { ! [ InitDlg $w "Contour Control"] } return

    grid columnconfigure $w 1 -weight 1
    grid rowconfigure    $w 2 -weight 1
    # import
    grid \
	[label $w.filelabel -text "File:" -anchor w] \
	[ns FileEntry $w.filename \
	     -textvariable fileName \
	     -command      "ImportContour $w \$fileName" \
	     -type         open \
	     -title        "Select a contour file" \
	     -filetypes    {
		 {{Contour files} {.plt .PLT} TEXT}
		 {{All files}     *}} \
	     -defaultextension .plt] \
	[ns button $w.fileimport \
	     -text    "Import" \
	     -command "ImportContour $w \$fileName"] \
	-sticky we
    # contours and levels
    grid \
	[labelframe $w.contours -text "Contours" -height 20] - - \
	-ipady 3 -padx 3 -sticky we
    grid \
	[labelframe $w.levels   -text "Contour levels" -height 20] - - \
	-padx 3 -sticky swen
    # display method
    grid \
	[label $w.dmethodlabel -text "Display method:" -anchor w] \
	[ns RadioGroup $w.dmethodgroup \
	     -texts    {Direct      Cached {Display lists}} \
	     -values   {direct      cached displist} \
	     -command  {SelectDisplayMethod $displayMethod} \
	     -variable displayMethod] - \
	-sticky we
    # main buttons
    grid \
	[ns MainButtonFrame $w.buttons \
	     -applycommand [list Apply $w] \
	     -helpcommand  {htmlShowHelp $gomHelpFile(contour)}] - - \
	-sticky we
    # update variables
    if { [show contour polygon method] == "save" } {
	set displayMethod cached
    } else {
	set displayMethod direct
    }
    catch {
	if { [show displaylists] && [show displaylists contour] } {
	    set displayMethod displists
	}
    }

    UpdateFileList $w
}

proc InitContourDlg { w title } {
    variable dlgs
    lappend dlgs $w
    # return if no molecular systems defined
    if { [show contour defined] < 1 } {
	lulErrorDialog {ERROR: no contours available. Import a contour first!}
	return 0
    }
    return [InitDlg $w $title]
}

namespace export InitContourDlg

############################################################################
# PROC
proc ImportContour { w FileName } {
    contour file $FileName [expr [show contour defined] + 1]
    UpdateFileList $w
}

############################################################################
# PROC
proc UpdateFileList { w } {
    regrid labelframe $w.contours
    regrid labelframe $w.levels
    set NumContours [show contour defined]
    if { $NumContours <= 0 } return
    set   f $w.contours
    grid columnconfigure $f 4 -weight 1
    # titles
    grid \
	[label $f.contourlabel -text "File name:" -anchor w] \
	[label $f.minlabel     -text "Min:"     -anchor w] \
	[label $f.maxlabel     -text "Max:"     -anchor w] \
	x \
	-sticky we
    # contours
    set NumContours [show contour defined]
    for { set i 1 } { $i <= $NumContours } { incr i } {
	grid \
	    [ns radiobutton $f.contour$i \
		 -anchor w \
		 -text "$i: [file tail [show contour filename $i]]" \
		 -variable contour -value $i \
		 -command  [list Level::UpdateList $w.levels $i]] \
	    [label $f.min$i \
		 -text [format %.3f [show contour minimum $i]] -anchor e] \
	    [label $f.max$i \
		 -text [format %.3f [show contour maximum $i]] -anchor e] \
	    [ns button $f.details$i \
		 -text "Details ..." \
		 -command [list EditDetails $i]] \
	    -padx 2 -pady 1 -sticky we
    }
    # buttons
    grid \
	[ns button $f.delete \
	     -text "Delete all contours" \
	     -command [list DeleteContours $w]] \
	[ns button $f.material \
	     -text "Edit material properties" \
	     -command EditMaterialProperties] - - \
	-padx 2 -pady 1 -sticky we
    # select the last one
    $f.contour$NumContours invoke
}

############################################################################
# PROC
proc EditDetails { Contour } {
    global gomStructureNames
    variable map
    variable glue

    set w .gomcontour_${Contour}
    if { ! [InitContourDlg $w "Contour $Contour"] } return

    set NumContours [show contour defined]
    for { set i 1 } { $i <= $NumContours } { incr i } {
	lappend maplist "$i: [file tail [show contour filename $i]]"
    }
    lset maplist [expr $Contour - 1] "$Contour: no mapping"
    set i [expr [lindex [show contour mapping $Contour] 0] - 1]
    set  map($Contour) [lindex $maplist $i]

    set NumStruct [show molstructures]
    for { set i 1 } { $i <= $NumStruct } { incr i } {
	lappend gluelist "$i: $gomStructureNames($i)"
    }
    set i [expr [show contour structure $Contour] - 1]
    set glue($Contour) [lindex $gluelist $i]

    eval grid [lb ns ComboBox $w.map \
		   -label        "Map source:" \
		   -editable     false \
		   -values       $maplist \
		   -textvariable map($Contour)] -sticky we
    eval grid [lb ns ComboBox $w.glue \
		   -label        "Glue to structure:" \
		   -editable     false \
		   -values       $gluelist \
		   -textvariable glue($Contour)] -sticky we
#    eval grid [lb ns RadioGroup $w.type \
		   -label    "Type:" \
		   -texts    {curve surface volume} \
		   -values   {curve surface volume} \
		   -variable type($Contour)] -sticky we
    grid [ns MainButtonFrame $w.buttons \
	      -applycommand [list ApplyDetails $Contour]] - -sticky we
}

############################################################################
# PROC
proc ApplyDetails { Contour } {
    variable glue
    variable map
    contour combine $Contour [lindex [split $glue($Contour) :] 0]
    contour mapping $Contour [lindex [split $map($Contour) :] 0]
    lulInvalidateDisplay
}

############################################################################
# PROC
proc DeleteContours { w } {
    contour delete
    variable dlgs
    lappend dlgs
    foreach dlg $dlgs { catch { destroy $dlg } }
    unset dlgs
    UpdateFileList $w
    lulInvalidateDisplay
}

############################################################################
proc Apply { w } {
    variable contour
    variable displayMethod

    set plotargs [list]

# Fixed 1998-07-28/LUL to return right contour mapping name
# set MappingContourName [show contour name [lindex $gomMappingIndex 0]]
# New Fix 1999-01-10/LUL. Introduced new command "contour mapping Name1 Name2
# where the values from Name2 will be colour mapped on Name1

    set MappingContourName [lindex [show contour mapping $contour] 0]

    if { $contour == $MappingContourName } {
	lappend plotargs $contour
	for { set i 1 } {  "" != [string trim $Level::value($contour,$i)] } { incr i } {
	    # <value> <colour>
	    lappend plotargs \
		[string trim $Level::value($contour,$i)] \
		$Level::colour($contour,$i)
	}
    } else { 
# apply contour mapping...
	lappend plotargs [list $contour $MappingContourName]
	for { set i 1 } {  "" != [string trim $Level::value($contour,$i)] } { incr i } {
	    # <value> <min> <max>
	    lappend plotargs \
		[lindex $Level::value($contour,$i) 0] \
		[lindex $Level::value($contour,$i) 1] \
		[lindex $Level::value($contour,$i) 2]
	}
    }
    eval contour plot $plotargs
    if { [show contour defined] } {
	SelectDisplayMethod $displayMethod
    }
    Level::UpdateList $w.levels $contour
    lulInvalidateDisplay
}

############################################################################
# PROC
proc EditMaterialProperties {} {
    variable specular
    variable shininess
    variable instantUpdateMaterial

    set w .gommaterialproperties
    if { ! [ InitDlg $w "Edit material properties"] } return

    set instantUpdateMaterial false

    foreach colour {red green blue} {
	grid \
	    [ns scale $w.${colour}scale \
		 -from         0.0 \
		 -to           1.0 \
		 -length       240 \
		 -variable     specular($colour) \
		 -orient       horizontal \
		 -label        "Material specular $colour" \
		 -tickinterval 0.2 \
		 -showvalue    false \
		 -foreground   dark$colour \
		 -digits       3 \
		 -resolution   0.01 \
		 -command      [list InstantUpdateMaterialProperty \
				    define material specular $colour]] \
	    [ns entry $w.${colour}entry \
		 -textvariable specular($colour) \
		 -width        5] \
	    -sticky w
	set specular($colour) [format %.2f [show material specular $colour]]
    }
    grid \
	[ns scale $w.shininessscale \
	     -from         0.0 \
	     -to           128.0 \
	     -length       240 \
	     -variable     shininess \
	     -orient       horizontal \
	     -label        "Material shininess" \
	     -tickinterval 32 \
	     -showvalue    false \
	     -digits       4 \
	     -resolution   0.1 \
	     -command      [list InstantUpdateMaterialProperty \
				define material shininess]] \
	[ns entry $w.shininessentry \
	     -textvariable shininess \
	     -width        5] \
	-sticky w
    set shininess [show material shininess]
    grid [ns checkbutton $w.instant \
	      -text     "Instant update" \
	      -variable instantUpdateMaterial] - -sticky w
    # main buttons
    grid \
	[ns MainButtonFrame $w.buttons \
	     -applycommand ApplyMaterialDetails \
	     -helpcommand  {htmlShowHelp $gomHelpFile(editmatprop)}] - \
	-sticky we
}

############################################################################
# PROC
proc ApplyMaterialProperties {} {
    variable specular
    variable shininess
    foreach colour { red green blue } {
	define material specular $colour $specular($colour)
    }
    define material shininess $shininess
}

############################################################################
# PROC
proc InstantUpdateMaterialProperty { args } {
    variable instantUpdateMaterial
    if { ! $instantUpdateMaterial } return
    eval $args
    display nocheck
}

############################################################################
# PROC
proc SelectDisplayMethod { method } {
    switch $method {
	direct {
	    contour method direct
	    catch {define displaylists contour off}
	}
	cached {
	    contour method save
	    catch {define displaylists contour off}
	}
	displist {
	    contour method direct
	    define displaylists on
	    catch {define displaylists contour on}
	}
    }
}

}
