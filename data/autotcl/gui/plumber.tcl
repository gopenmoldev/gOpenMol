##############################################################################
#                           Copyright (c) 2002 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################
#
# Make plumber list
#

namespace eval lulPlumber {

variable plumberColour           "#ffffff"
variable defaultPlumberRadius    0.5
variable defaultPlumberWidth     2.0
variable defaultPlumberThickness 0.5
variable plumberAction           0
variable AppendPlumber           0

##################### Plumber type display ############################
#
# define plumber display
#
proc DefinePlumber {} {
    package require gom::gui::Widgets
    namespace import -force ::gom::gui::Widgets::*

    global gomControlFont
    global gomHelpDir
    global gomHelpFile
    global gomControlFont

    variable secondaryStructureRecord
    variable secondaryStructureRecords
    variable defaultPlumberRadius
    variable defaultPlumberWidth
    variable defaultPlumberThickness
    variable AppendPlumber
    variable GluePlumber
    variable plumberType
    variable plumberColour

    # return if no molecular systems defined
    if {[show molstructures] < 1} {
	lulErrorDialog {ERROR: no structure available. Read a structure first!}
	return
    }

    set ns [namespace current]

    set w .gomplumber
    catch {destroy $w; destroy .gomplumberseco}
    toplevel $w
    wm title $w "Plumber"
    wm iconname $w "Plumber"
    
    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
	-command "destroy $w"
    button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
	-command "${ns}::DoPlumberType $w; lulInvalidateDisplay"
    button $w.buttons.help    -text Help    -font "$gomControlFont" \
	-command "htmlShowHelp $gomHelpFile(plumber)"
    pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help \
	-side left -expand 1
    
    label  $w.label -text "Plumber type display:" -font {Arial 12 bold}
    pack   $w.label -side top -anchor w -pady 3

    ::lulCreateAtomInputEntries $w 30 "*" "*" "*" "" 1
    
    frame  $w.frame4
    label  $w.frame4.label -text "Secondary structure: "
    button $w.frame4.button -text "Click for information" \
	-command ${ns}::ShowSecondaryStructureRecords
    pack   $w.frame4 -side top -anchor w -pady 4
    pack   $w.frame4.label $w.frame4.button -side left
    
    if {[array exists secondaryStructureRecord] ||
        [find chains "*" "*" selected] != ""} {
	$w.frame4.label  configure -state normal
	$w.frame4.button configure -state normal
    } else {
	$w.frame4.label  configure -state disabled
	$w.frame4.button configure -state disabled
    }

    frame  $w.edit -borderwidth 2 -relief ridge
    button $w.edit.list -text "List plumber(s)" \
	-command ${ns}::ShowPlumberList
    button $w.edit.delete -text "Delete plumber(s)" \
	-command "
	    plumber display off
	    plumber -atoms
	    $w.display.off select
	    lulInvalidateDisplay"
    pack   $w.edit -side bottom -anchor w 
    pack   $w.edit.list $w.edit.delete -padx 4 -pady 2 -side left
    
    frame  $w.display -borderwidth 2 -relief ridge
    label       $w.display.label -text "Display:"
    radiobutton $w.display.on  -text "On"  \
	-variable ${ns}::PlumberAction -value 1 \
	-command {plumber display on;lulInvalidateDisplay}
    radiobutton $w.display.off -text "Off" \
	-variable ${ns}::PlumberAction -value 0 \
	-command {plumber display off;lulInvalidateDisplay}
    pack   $w.display -side left 
    pack   $w.display.label $w.display.on $w.display.off -side top -anchor w
    if {[show plumber display]} {
	$w.display.on select
    } else {
	$w.display.off select
    }
    
    frame  $w.type -borderwidth 2 -relief ridge

    radiobutton $w.type.tube    -text "Tube"     \
	-variable ${ns}::plumberType -value tube \
	-command {gomPrint "Tube mode is on"}
    
    radiobutton $w.type.ribbon   -text "Ribbon"    \
	-variable ${ns}::plumberType -value ribbon \
	-command {gomPrint "Ribbon mode is on"}

    radiobutton $w.type.shelix   -text "Solid helix"    \
	-variable ${ns}::plumberType -value shelix \
	-command {gomPrint "Solid helix mode is on"}

    radiobutton $w.type.fhelix   -text "Flat helix"    \
	-variable ${ns}::plumberType -value fhelix \
	-command {gomPrint "Flat helix mode is on"}

    radiobutton $w.type.arrow   -text "Arrow"    \
	-variable ${ns}::plumberType -value arrow \
	-command {gomPrint "Arrow mode is on"}

    radiobutton $w.type.strand  -text "Strand"    \
	-variable ${ns}::plumberType -value strand \
	-command {gomPrint "Strand mode is on"}

    radiobutton $w.type.trace   -text "Trace"    \
	-variable ${ns}::plumberType -value trace \
	-command {gomPrint "Trace mode is on"}

    pack   $w.type -side left
    pack   $w.type.tube $w.type.ribbon  \
	$w.type.shelix $w.type.fhelix   \
	$w.type.arrow  $w.type.strand   \
	$w.type.trace    \
	-side top -anchor w
    $w.type.tube select

    frame  $w.radius -borderwidth 2 -relief ridge
    label  $w.radius.label -text "Radius: "
    entry  $w.radius.value -width 17
    pack   $w.radius -side top -anchor e
    pack   $w.radius.label $w.radius.value -side left
    
    $w.radius.value insert 0 $defaultPlumberRadius

    frame  $w.colour -borderwidth 2 -relief ridge
    ColorButton $w.colour.button \
	-width 23 \
	-text "Choose plumber colour" \
	-variable ${ns}::plumberColour
#    button $w.colour.button -text "Colour" -width 23 \
	-command "::lulChooseButtonColour $w.colour.button \
	    ${ns}::plumberColour {Choose plumber colour}" \
	-fg [::lulGetVisibleForegroundColour $plumberColour] \
	-bg $plumberColour
    pack   $w.colour -side top -anchor e
    pack   $w.colour.button -side left
    
    frame  $w.width -borderwidth 2 -relief ridge
    label  $w.width.label -text "Helix major axis: "
    entry  $w.width.value -width 10
    pack   $w.width -side top -anchor e
    pack   $w.width.label $w.width.value -side left
    
    $w.width.value insert 0 $defaultPlumberWidth
    
    frame  $w.thickness -borderwidth 2 -relief ridge
    label  $w.thickness.label -text "Helix minor axis: "
    entry  $w.thickness.value -width 10
    pack   $w.thickness -side top -anchor e
    pack   $w.thickness.label $w.thickness.value -side left
    
    $w.thickness.value insert 0 $defaultPlumberThickness
    
    frame  $w.glue -borderwidth 2 -relief ridge
    checkbutton $w.glue.button -text "Glue to atoms" -width 20 \
	-variable ${ns}::GluePlumber
    pack   $w.glue -side top -anchor e
    pack   $w.glue.button -side left
    
    $w.glue.button select

    frame  $w.append -borderwidth 2 -relief ridge
    checkbutton $w.append.button -text "Append to stack" -width 20 \
	-variable ${ns}::AppendPlumber
    pack   $w.append -side bottom -anchor e
    pack   $w.append.button -side left
    
    $w.append.button deselect
}

########################################################################
proc DoPlumberType { w } {
    variable plumberColour
    variable plumberType
    variable PlumberAction
    variable GluePlumber
    variable AppendPlumber

    set Segment   [string trim [$w.segment.input   get]]
    if {$Segment == ""} {set Segment "*"}
    set Residue   [string trim [$w.residue.input   get]]
    if {$Residue == ""} {set Residue "*"}
    set Atom      [string trim [$w.atom.input      get]]
    if {$Atom    == ""} {set Atom    "*"}
    set Radius    [string trim [$w.radius.value    get]]
    if {$Radius  == ""} {set Radius  "0.5"}
    set Width     [string trim [$w.width.value     get]]
    if {$Width   == ""} {set Width   "2.0"}
    set Thickness [string trim [$w.thickness.value get]]
    if {$Thickness == "" } {set Thickness "0.5"}

    if {!$AppendPlumber} {
	plumber -atoms
    }

    set Command [list plumber atoms $Segment $Residue $Atom \
		     $Radius  $plumberColour $plumberType $Width $Thickness]

    if {$GluePlumber} {set Command "$Command glue"}

    eval $Command
    if {$PlumberAction} {
	plumber display on
    } else {
	plumber display off
    }
}
##########################################################################
proc ShowSecondaryStructureRecords {} {
    global gomHelpDir
    global gomHelpFile
    global gomControlFont
    global gomControlFontFixed

    variable secondaryStructureRecord
    variable secondaryStructureRecords
    variable secStructInformation

    array unset secondaryStructureInformation

    # return if no molecular systems defined
    if {[show molstructures] < 1} {
	lulErrorDialog {ERROR: no structure available. Read a structure first!}
	return
    }

    set ns [namespace current]

    set w .gomplumberseco
    catch {destroy $w}
    toplevel $w 
    wm title $w "Secondary Structure"
    wm iconname $w "Secondary"

    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
	-command "destroy $w"
    button $w.buttons.select -text Select -font "$gomControlFont" \
	-command "${ns}::ProcessSecondaryStructureRecord $w 1"
    button $w.buttons.help    -text Help    -font "$gomControlFont" \
	-command "htmlShowHelp $gomHelpFile(plumberseco)"
    pack   $w.buttons.dismiss $w.buttons.select $w.buttons.help -side left -expand 1

    button $w.applyByType \
	-text "Create all plumbers of this kind" \
	-font "$gomControlFont" \
	-command "${ns}::ProcessSecondaryStructureRecordsByType $w"
    pack   $w.applyByType -side bottom -pady 3

    label  $w.label -text "Plumber secondary structure information:" \
	-font {Arial 12 bold}
    pack   $w.label -side top -anchor w -pady 3

    frame $w.lines -borderwidth 2 -relief ridge
    pack  $w.lines -side top -anchor n -fill both -expand 1

    text $w.lines.text \
	-xscrollcommand "$w.lines.scrollx set" \
	-yscrollcommand "$w.lines.scrolly set" \
	-width 50 -height 30
    scrollbar $w.lines.scrollx -command "$w.lines.text xview" \
	-orient horizontal
    scrollbar $w.lines.scrolly -command "$w.lines.text yview"
    pack $w.lines.scrolly -side right  -fill y
    pack $w.lines.scrollx -side bottom -fill x
    pack $w.lines.text    -side left   -fill both -expand 1

    set loop 0

    for {set i 1} {$i <= [show molstructures]} {incr i} {

	if {![show status $i]} continue

	if {[info exists secondaryStructureRecords($i)]} {
	    for {set j 1} {$j <= $secondaryStructureRecords($i) } {incr j} {
		ShowSecondaryStructureRecord $w loop \
		    "$secondaryStructureRecord($i,$j)"
	    }
	}
	ShowSecondaryStructureChains $w loop $i
    }

    if {$loop} {$w.lines.text.box1.line select}
}

##########################################################################
proc ShowSecondaryStructureChains { w loopVar iStruct } {
    variable secondaryStructureRecord
    variable secondaryStructureRecords
    upvar $loopVar loop

    set count 0

    set chains [find chains "*" "*" $iStruct]
    if {[llength $chains] < 1} return

    foreach chain $chains {
	# split list
	set ID    [lindex $chain 0]
	set chain [lindex $chain 1]
	if {$ID == ""} {set ID "*"}
	# add new trace
	ShowSecondaryStructureChain \
	    $w loop count "TRACE" $iStruct $ID $chain
	# append to chain list
	lappend ChainList($ID) $chain
    }

    # create record list
    if {[info exists secondaryStructureRecords($iStruct)]} {
	for {set j 1} {$j <= $secondaryStructureRecords($iStruct) } {incr j} {
	    set line "$secondaryStructureRecord($iStruct,$j)"
	    switch -glob $line {
		HELIX* {
		    set iStart [string trim [string range $line 21 24]]
		    set iStop  [string trim [string range $line 33 36]]
		    set ID     [string trim [string range $line 19 19]]
		    if {$ID == ""} {set ID "*"}
		    lappend RecordList($ID) $iStart-$iStop
		}
		SHEET* {
		    set iStart [string trim [string range $line 22 25]]
		    set iStop  [string trim [string range $line 33 36]]
		    set ID     [string trim [string range $line 21 21]]
		    if {$ID == ""} {set ID "*"}
		    lappend RecordList($ID) $iStart-$iStop
		}
	    }
	}
    }

    # sort by residue number (range value)
    lappend RecordList(*)
    if {[llength [array names RecordList]] == 1} {
	foreach ID [array names ChainList] {
	    set RecordList($ID) [lsort -dictionary $RecordList(*)]
	}
    } else {
	foreach ID [array names RecordList] {
	    set RecordList($ID) [lsort -dictionary [
		concat $RecordList($ID) $RecordList(*)]]
	}
    }
    unset RecordList(*)

    set count 0
    foreach ID [array names ChainList] {
	foreach residuelist $ChainList($ID) {
	    if {[info exists RecordList($ID)]} {
		set residuelist [split $residuelist ,]
		set chain ""
		while {[llength $residuelist] > 0} {
		    set iStart [lindex [split $residuelist -]   0]
		    set iStop  [lindex [split $residuelist -] end]
		    
		    set iiStart $iStart
		    set iiStop  $iStop
		    while {$iiStart <= $iiStop} {
			foreach record $RecordList($ID) {
			    scan $record "%d-%d" rStart rStop
			    if {$rStart >= $iiStart && $rStart < $iiStop } {
				# Helix or sheet overlaps the end or
				# middle part of the range.
				set iiStop  $rStart
			    }
			    if {$rStart <= $iiStart && $rStop  > $iiStart} {
				# Helix or sheet overlaps the beginning of the
				# range.
				set iiStart $rStop 
			    }
			    if {$iiStop < $iStop} break
			}
			if {$iiStart > $iStart} {
			    # Create a line.
			    ShowSecondaryStructureChain \
				$w loop count "LOOP" $iStruct $ID $chain
			    set chain ""
			}
			if {$iiStart <= $iiStop} {
			    if {"" != $chain} {set chain "$chain,"}
			    if {$iiStart < $iiStop} {
				set chain "${chain}$iiStart-$iiStop"
			    } else {
				set chain "${chain}$iiStart"
			    }
			}
			if {$iiStop < $iStop} {
			    # Create a line.
			    ShowSecondaryStructureChain \
				$w loop count "LOOP" $iStruct $ID $chain
			    set chain ""
			}
			set iiStart [expr $iiStop+1]
			set iiStop  $iStop
			set iStart  $iiStart
		    }
		    set residuelist [lrange $residuelist 1 end]
		}
		if {"" != $chain} {
		    # Create a line.
		    ShowSecondaryStructureChain \
			$w loop count "LOOP" $iStruct $ID $chain
		}
	    } else {
		ShowSecondaryStructureChain \
		    $w loop count "LOOP" $iStruct $ID $chain
	    }
	}
    }
}

##########################################################################
proc ShowSecondaryStructureChain { w loopVar countVar text iStruct ID chain} {
    upvar $loopVar loop
    upvar $countVar count

    set residuelist [split $chain ,-]
    if {[llength $residuelist] < 2} return

    set iStart [lindex $residuelist   0]
    set iStop  [lindex $residuelist end]

    set rStart "*"
    set rStop  "*"
    catch {
	set iAStart [show atomnumber "$ID" $iStart "CA" $iStruct]
	set iAStop  [show atomnumber "$ID" $iStop  "CA" $iStruct]
	set rStart  [show atom residue $iAStart $iStruct]
	set rStop   [show atom residue $iAStop  $iStruct]
    }

    incr count
#1234567890123456789012345678901234567890123456789012345678901234567890123456
#         1         2         3         4         5         6         7
#HELIX. SER .ID RES C RESN  RES C RESN T COMMENT....................... COUNT
#LOOP.. SER RES CCC RESN  RES CCC RESN   ..................RESIDUE.LIST COUNT
    ShowSecondaryStructureRecord $w loop [format \
	"%-40s%30s %5s    " \
	[format "%-6s %3s %-3s %-3s %4s  %-3s %-3s %4s" \
	    $text $count \
	    $rStart $ID $iStart \
	    $rStop  $ID $iStop] \
	$chain [::lulCalculateIndexRangeListLength $chain]]
}

##########################################################################
proc ShowSecondaryStructureRecord { w loopVar text } {
    global gomControlFontFixed

    variable secondaryStructureInformation
    variable secondaryStructureSelection

    upvar $loopVar loop

    set ns [namespace current]

    incr loop
    frame $w.lines.text.box$loop
    pack  $w.lines.text.box$loop -side top -anchor w

    set secondaryStructureInformation($loop) "$text"
    radiobutton $w.lines.text.box$loop.line \
	-text "$text" -cursor hand2 -value "$loop" \
	-font $gomControlFontFixed -width 100 -anchor w \
	-variable ${ns}::secondaryStructureSelection
    pack $w.lines.text.box$loop.line  -side top -anchor w
    $w.lines.text window create $loop.0 -window $w.lines.text.box$loop
}

##########################################################################
proc ProcessSecondaryStructureRecordsByType { w } {
    variable secondaryStructureInformation
    variable secondaryStructureSelection
    variable AppendPlumber

    set wp    .gomplumber
    set sel   $secondaryStructureSelection
    set app   $AppendPlumber
    set info  "$secondaryStructureInformation($secondaryStructureSelection)"
    set type  [lindex [split "$info"] 0]
    set first 1

    foreach i [array names secondaryStructureInformation] {
	if {[lindex [split $secondaryStructureInformation($i)] 0] == $type} {
	    if {$first} {
		set first         0
	    } else {
		set AppendPlumber 1
	    }
	    set secondaryStructureSelection  $i
	    ProcessSecondaryStructureRecord $w 0
	    DoPlumberType                    $wp
	}
    }

    set secondaryStructureSelection $sel
    set AppendPlumber               $app

    lulInvalidateDisplay
}

##########################################################################
proc ProcessSecondaryStructureRecord { w setType } {
    variable secondaryStructureInformation
    variable secondaryStructureSelection

    set info     "$secondaryStructureInformation($secondaryStructureSelection)"
    set wp       .gomplumber
    set residues ""

    switch -glob $info {
	"HELIX*" {
	    set iStart [string trim [string range $info 21 24]]
	    set iStop  [string trim [string range $info 33 36]]
	    set iD     [string trim [string range $info 19 19]]
	    set type    shelix
	}
	"SHEET*" {
	    set iStart [string trim [string range $info 22 25]]
	    set iStop  [string trim [string range $info 33 36]]
	    set iD     [string trim [string range $info 21 21]]
	    set type   arrow
	}
	"TURN*" {
	    return
	}
	"TRACE*" {
		set iD       [lindex $info 3]
		set residues [lindex $info 8]
	    set type trace
	}
	"LOOP*" {
		set iD       [lindex $info 3]
		set residues [lindex $info 8]
	    set type tube
	}
	default {
	    gomError "recognizes only records starting with 'HELIX', 'SHEET', 'TURN', 'TRACE' and 'LOOP'"
	    return
	}
    }

    if {$iD == ""} {set iD "*"}
    if {$residues == ""} {set residues $iStart-$iStop}

    $wp.segment.input delete 0 end
    $wp.segment.input insert 0 $iD
    $wp.residue.input delete 0 end
    $wp.residue.input insert 0 "$residues"
    $wp.atom.input    delete 0 end
    $wp.atom.input    insert 0 "CA"

    if {$setType} {$wp.type.$type select}
}

########################################################################
# PROC
proc ShowPlumberList {} {
    global gomHelpFile
    global gomControlFont

    set ns [namespace current]

    set PlumberSets [show plumber count]

    set w .gomplumberlist
    catch {destroy $w}
    toplevel $w 
    wm title $w "Plumber list"
    wm iconname $w "Plumber list"

    UpdatePlumberSet $w 1

    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.update -text Update  -font "$gomControlFont" \
	-command "${ns}::UpdatePlumberSet $w 0"
    button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
	-command "destroy $w"
    button $w.buttons.help    -text Help     -font "$gomControlFont" \
	-command "htmlShowHelp $gomHelpFile(plumberlist)"
    pack   $w.buttons.update $w.buttons.dismiss $w.buttons.help \
	-side left -expand 1

    frame  $w.select -borderwidth 2 -relief ridge
    pack   $w.select -side left -anchor nw -pady 1m
    button $w.select.all    -text "Select all" \
	-font "$gomControlFont" \
	-command "${ns}::SetPlumberSelection $w select"
    button $w.select.none   -text "Select none" \
	-font "$gomControlFont" \
	-command "${ns}::SetPlumberSelection $w deselect"
    button $w.select.invert -text "Invert selection" \
	-font "$gomControlFont" \
	-command "${ns}::SetPlumberSelection $w toggle"
    button $w.select.bytype -text "Select by type" \
	-font "$gomControlFont" \
	-command "${ns}::SetPlumberSelectionByType $w"
    pack   $w.select.all $w.select.none \
	$w.select.invert $w.select.bytype \
	-side top -fill x

    frame  $w.edit -borderwidth 2 -relief ridge
    pack   $w.edit -side left -anchor nw -padx 5m -pady 1m
    button $w.edit.changecolour -text "Change colour" \
	-font "$gomControlFont" \
	-command "${ns}::ChangePlumberColour $w"
    button $w.edit.invertcolour -text "Invert colour" \
	-font "$gomControlFont" \
	-command "${ns}::InvertPlumberColour $w"
    button $w.edit.unglue       -text "Unglue from atoms" \
	-font "$gomControlFont" \
	-command "${ns}::UngluePlumber $w"
    button $w.edit.destroy      -text "Destroy" \
	-font "$gomControlFont" \
	-command "${ns}::DestroyPlumber $w"
    pack   $w.edit.changecolour $w.edit.invertcolour \
	$w.edit.unglue $w.edit.destroy \
	-side top  -fill x
}

proc UpdatePlumberSet { w create } {
    global gomControlFontFixed
    variable plumberSelectionList

    set ns [namespace current]

    set PlumberSets [show plumber count]

    array unset plumberSelectionList

    if {$create} {
	frame $w.plumbers -borderwidth 2 -relief ridge
	pack  $w.plumbers -side top -fill both -expand 1

	label $w.plumbers.label -anchor w -text "Plumbers:"
	scrollbar $w.plumbers.scrollx -command "$w.plumbers.text xview" \
	    -orient horizontal
	scrollbar $w.plumbers.scrolly -command "$w.plumbers.text yview"

	pack  $w.plumbers.label -side top -anchor w
	pack  $w.plumbers.scrolly -side right  -fill y
	pack  $w.plumbers.scrollx -side bottom -fill x
    } else {
	destroy $w.plumbers.text
    }

    text $w.plumbers.text \
	-xscrollcommand "$w.plumbers.scrollx set" \
	-yscrollcommand "$w.plumbers.scrolly set" \
	-width 50 -height 30
    pack $w.plumbers.text -side left -fill both -expand 1
    

    for {set i 1} {$i <= $PlumberSets} {incr i} {
	scan [show plumber color $i] "%f %f %f" red green blue
	set bgColour    [lulColourFloat2Hex $red $green $blue]
	set bgColourInv [lulColourFloat2Hex [
	    expr 1.0-$red] [expr 1.0-$green] [expr 1.0-$blue]]

	set type [show plumber type $i]
	switch $type {
	    "tube"   {set type "Tube"}
	    "ribbon" {set type "Ribbon"}
	    "shelix" {set type "Solid helix"}
	    "fhelix" {set type "Flat helix"}
	    "arrow"  {set type "Arrow"}
	    "strand" {set type "Strand"}
	    "trace"  {set type "Trace"}
	}

	set plumberSelectionList($i) 0
	frame $w.plumbers.text.frame$i
	checkbutton $w.plumbers.text.plumber$i \
	    -text [format "%-11s %s" $type [GetPlumberDescription $i]] \
	    -font $gomControlFontFixed \
	    -width 100 -anchor w \
	    -fg [lulGetVisibleForegroundColour $bgColour] \
	    -bg $bgColour -selectcolor $bgColour \
	    -activeforeground [lulGetVisibleForegroundColour $bgColourInv] \
	    -activebackground $bgColourInv \
	    -onvalue 1 -offvalue 0 \
	    -variable ${ns}::plumberSelectionList($i)
	pack  $w.plumbers.text.plumber$i -side top
	$w.plumbers.text window create $i.0 -window $w.plumbers.text.plumber$i
    }
}

proc GetPlumberDescription { Which } {
    global gomStructureNames

    catch {
	set NStruct [show plumber structure $Which]
	set Struct  "$gomStructureNames($NStruct)"
	set atoms   [show plumber atoms $Which]

	if { $atoms != "" } {
	    set first [lindex [split $atoms ,-] 0]
	    set last  [lindex [split $atoms ,-] end]
	    set Seg   [show atom segmentname $first $NStruct]
	    set Desc  "$Struct   [format \
		"%-3s %-3s %4s  %-3s %-3s %4s" \
		[show atom residuename $first $NStruct] \
		$Seg \
		[show atom resnumber   $first $NStruct] \
		[show atom residuename $last  $NStruct] \
		$Seg \
		[show atom resnumber   $last  $NStruct]]" \
	} else {
	    set Desc "$Struct   <<< not gluead to atoms >>>"
	}
	return "$Desc"
    }
    return ""
}

proc DestroyPlumber { w } {
    variable plumberSelectionList

    set PlumberSets [show plumber count]
    set ypos        [lindex [$w.plumbers.text yview] 0]

    for {set i $PlumberSets} {$i >= 1} {set i [expr $i-1]} {
	if {$plumberSelectionList($i)} {plumber destroy $i}
    }
    if {[show plumber count] <= 0} {plumber display off}
    lulInvalidateDisplay
    UpdatePlumberSet $w 0
    $w.plumbers.text yview moveto $ypos
}

proc UngluePlumber { w } {
    variable plumberSelectionList

    set PlumberSets [show plumber count]
    set tempList    [array get plumberSelectionList]
    set ypos        [lindex [$w.plumbers.text yview] 0]

    for {set i 1} {$i <= $PlumberSets} {incr i} {
	if {$plumberSelectionList($i)} {
	    plumber structure [show plumber structure $i] $i
	}
    }
    UpdatePlumberSet $w 0
    array set plumberSelectionList $tempList
    $w.plumbers.text yview moveto $ypos
}

proc ChangePlumberColour { w } {
    variable plumberSelectionList

    set colour [tk_chooseColor \
		    -title {Choose plumber colour} \
		    -parent $w]
    if {"" == "$colour"} {return}

    set PlumberSets [show plumber count]
    set tempList    [array get plumberSelectionList]
    set ypos        [lindex [$w.plumbers.text yview] 0]

    for {set i 1} {$i <= $PlumberSets} {incr i} {
	if {$plumberSelectionList($i)} {
	    plumber colour $colour $i
	}
    }
    lulInvalidateDisplay
    UpdatePlumberSet $w 0
    array set plumberSelectionList $tempList
    $w.plumbers.text yview moveto $ypos
}

proc InvertPlumberColour { w } {
    variable plumberSelectionList
    set PlumberSets [show plumber count]
    set tempList    [array get plumberSelectionList]
    set ypos        [lindex [$w.plumbers.text yview] 0]

    for {set i 1} {$i <= $PlumberSets} {incr i} {
	if {$plumberSelectionList($i)} {
	    scan [show plumber colour $i] "%f %f %f" red green blue
	    plumber colour "[
		expr 1.0-$red] [expr 1.0-$green] [expr 1.0-$blue]" $i
	}
    }
    lulInvalidateDisplay
    UpdatePlumberSet $w 0
    array set plumberSelectionList $tempList
    $w.plumbers.text yview moveto $ypos
}

proc SetPlumberSelection { w method } {
    set PlumberSets [show plumber count]

    for {set i 1} {$i <= $PlumberSets} {incr i} {
	$w.plumbers.text.plumber$i $method
    }
}

proc SetPlumberSelectionByType { w } {
    variable plumberSelectionList

    set PlumberSets [show plumber count]
    set types ""

    for {set i 1} {$i <= $PlumberSets} {incr i} {
	if {$plumberSelectionList($i)} {
	    set type [show plumber type $i]
	    if {[lsearch -exact $types $type] < 0} {
		lappend types $type
	    }
	}
    }

    if {[llength $types] == 0} {
	lulErrorDialog {Select at least one plumber first!}
	return
    }

    for {set i 1} {$i <= $PlumberSets} {incr i} {
	set type [show plumber type $i]
    	if {[lsearch -exact $types $type] >= 0} {
	    $w.plumbers.text.plumber$i select
	}
    }
}

}
