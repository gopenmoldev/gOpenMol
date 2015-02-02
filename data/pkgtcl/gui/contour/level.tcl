##############################################################################
#                       Copyright (c) 1994 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Enhancements 2003 - 2004 by: Eero Häkkinen
##############################################################################

namespace eval ::gom::gui::Contour::Level {

############################################################################
# NOTE:
#   UpdateList will unset all namespace level variables.
#   Thus do not define any.
############################################################################

array set planeCoord [list]

############################################################################
# PROC
proc UpdateList { w Contour } {
    package require gom::gui::Widgets
    namespace import -force ::gom::gui::Widgets::*
    namespace import -force [namespace parent]::InitContourDlg
    namespace export InitContourDlg

    # reset variables
    catch {eval unset [info vars [namespace current]::*]} err
    # reset widgets
    catch { eval destroy [winfo children $w] }
    pack  [ScrollFrame $w.frame] -fill both -expand true
    set f [$w.frame getframe]
    # create widgets
    set Max 10
    catch { set Max [expr [show contour levels $Contour] + 10] }
    for { set i 1 } { $i <= $Max } { incr i } {
	grid \
	    [label $f.label$i -text "Value $i:" -anchor w] \
	    [ns entry $f.value$i \
		 -textvariable value($Contour,$i)] \
	    [ns ColorButton $f.colour$i \
		 -variable     colour($Contour,$i)] \
	    [ns button $f.details$i \
		 -text         "Details ..." \
		 -command [list Details::Edit $Contour $i]] \
	    -padx 2 -sticky w
	# Destroy buttons which cannot be used, anyway.
	if { $i > $Max - 10 } { destroy $f.details$i }
    }
    # update variables
    variable value
    variable colour
    for { set i [show contour levels $Contour] } { $i > 0 } { incr i -1 } {
	set val [show contour values $Contour $i]
	set value($Contour,$i)  [lindex $val 0]
	set colour($Contour,$i) [lrange $val 1 3]
    }
    # This is to simplify a for loop in Contour::Apply
    set value($Contour,[expr $Max + 1]) ""
}

namespace eval Details {
    catch {eval unset [info vars [namespace current]::*]} err
############################################################################
# PROC
proc Edit { Contour Level } {
    package require gom::gui::Widgets
    namespace import -force ::gom::gui::Widgets::*
    namespace import -force [namespace parent]::InitContourDlg
    # Check that $Contour $Level exists.
    set w .gomcontour_${Contour}_${Level}
    if { ! [InitContourDlg $w "Contour $Contour, level $Level"] } return
    if { [show contour levels $Contour] < $Level } {
	destroy $w
	lulErrorDialog "ERROR: A contour level $Level is not defined yet."
	return
    }

    set ID $Contour,$Level

    # frames
    pack [labelframe      $w.state   -text "State"  ] \
	-side top -fill x -expand true -padx 5
    pack [labelframe      $w.opacity -text "Opacity"] \
	-side top -fill x -expand true -padx 5
    pack [frame           $w.plane] \
	-side top -fill x -expand true -padx 5
    pack [ns MainButtonFrame $w.buttons \
	      -applycommand [list Apply $Contour $Level] \
	      -helpcommand  {htmlShowHelp $gomHelpFile(contour_details)}] \
	-side top -fill both -expand true

    # states
    grid columnconfigure $w.state 1 -weight 1
    eval grid [lb ns RadioGroup $w.state.type \
		   -label  "Contour type:" \
		   -texts  {Solid Mesh Line} \
		   -values {solid mesh line} \
		   -variable type($ID)] \
	-sticky we
    foreach state {display smooth cullface} {
	eval grid [lb ns RadioGroup $w.state.$state \
		       -label "[string totitle $state] state:" \
		       -texts {On Off} -values {on off} \
		       -variable ${state}($ID)] \
	    -sticky we
    }

    # opacity
    grid columnconfigure $w.opacity 0 -weight 1
    grid columnconfigure $w.opacity 2 -weight 1
    grid [ns scale $w.opacity.scale \
	      -from 0.0 -to 1.0 -length 240 \
	      -variable opacity($ID) \
	      -orient horizontal -tickinterval 0.2 \
	      -showvalue true -digits 3 -resolution 0.01] - - \
	-sticky we
    grid x [ns entry $w.opacity.entry -width 5 \
		-textvariable opacity($ID)] x \
	-sticky we

    # plane
    pack [labelframe $w.plane.lf -text "Plane"] -fill both -expand true
    grid columnconfigure $w.plane.lf 1 -weight 1
    eval grid [lb ns RadioGroup $w.plane.lf.axis \
		   -label "Plane:" \
		   -texts {YZ XZ XY} -values {x y z} \
		   -variable planeAxis($ID)] - -padx 2 -sticky we
    eval grid [lb ns entry $w.plane.lf.coord \
		   -label Coordinate: -textvariable planeCoord($ID)] -sticky we

    # update variables
    variable type
    set type($ID)    [lindex \
			  [$w.state.type cget -values] \
			  [show contour type $Contour $Level]]
    foreach state {display smooth cullface} {
	variable $state
	set ${state}($ID) [expr [show contour $state $Contour $Level] ? \
			       {"on"} : {"off"}]
    }
    variable opacity
    set opacity($ID) [show contour alpha $Contour $Level]
    variable planeAxis
    variable planeCoord
    set planeAxis($ID)  [show contour clipplane direction $Contour $Level]
    set planeCoord($ID) [show contour clipplane position  $Contour $Level]
}

############################################################################
# PROC
proc Apply { Contour Level } {
    set ID $Contour,$Level
    foreach state {display smooth cullface type} {
	variable $state
	contour $state $Contour [expr \$${state}($ID)] $Level
    }
    variable opacity
    contour alpha $Contour $opacity($ID) $Level
    set opacity($ID) [show contour alpha $Contour $Level]
    variable planeAxis
    variable planeCoord
    contour clipplane $planeAxis($ID) $Contour $planeCoord($ID) $Level
    set planeCoord($ID) [show contour clipplane position $Contour $Level]
    lulInvalidateDisplay
}

}

}
