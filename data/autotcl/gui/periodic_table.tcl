##############################################################################
#                           Copyright (c) 1997 - 2002 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Enhancements by: Eero Häkkinen
##############################################################################
#
# Periodic system
#
proc lulDefinePeriodicTable {} {
    global gomControlFont
    global gomHelpDir
    global gomHelpFile
    global gomAtom_SN
    
    set w .periodictbl
    catch {destroy $w}
    toplevel $w 
    wm title $w "Periodic Table"
    wm iconname $w "Periodic"
    
    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss \
	-command "destroy $w" -font $gomControlFont
    button $w.buttons.apply   -text Apply  \
	-state disabled -font $gomControlFont
    button $w.buttons.help    -text Help -command \
	"htmlShowHelp $gomHelpFile(periodictable)" \
	-font $gomControlFont
    pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1
	
    frame   $w.frame -borderwidth 2 -relief ridge
    pack    $w.frame -side top
    frame   $w.frame.main
    pack    $w.frame.main -side top

    set line         0
    set elementCount 0
    foreach info [list {1 0 1} {2 0 6}  {2 0 6} {2 0 16} {2 0 16} \
	{2 15 15} \
	[list 2 15 [expr [llength [array names gomAtom_SN]] - 103]]] {

	incr line
	scan $info "%d %d %d" left middle right
	set  count  [expr $left + $middle + $right]

	frame $w.frame.main.line$line
	pack  $w.frame.main.line$line -side top -fill x
	frame $w.frame.main.line$line.left
	frame $w.frame.main.line$line.right
	pack  $w.frame.main.line$line.left  -side left
	if {$middle} {
	    pack  $w.frame.main.line$line.right -side left
	} else {
	    pack  $w.frame.main.line$line.right -side right
	}

	for {set i 1} {$i <= $left} {incr i} {
	    set element  [expr $elementCount + $i]
	    set symbol   [lindex $gomAtom_SN($element) 0]
	    set bgColour [eval lulColourFloat2Hex [lindex [show atom info $symbol] 2]]
	    set fgColour [lulGetVisibleForegroundColour $bgColour]
	    set name     [string tolower $symbol]
	    button $w.frame.main.line$line.left.$name -text $symbol \
		-fg $fgColour -bg $bgColour -width 2 -font $gomControlFont \
		-command "lulDisplayPeriodicTableElement $element"
	    pack   $w.frame.main.line$line.left.$name -side left
	}
	if {$middle} {
	    set element [expr $elementCount + $i]
	    set symbol  [lindex $gomAtom_SN($element) 0]
	    set name     [string tolower $symbol]
	    button $w.frame.main.line$line.left.$name -text $symbol \
		-width 2 -state disabled -font $gomControlFont
	    pack   $w.frame.main.line$line.left.$name -side left
	}
	for {set i [expr $i + $middle]} {$i <= $count} {incr i} {
	    set element [expr $elementCount + $i]
	    set symbol  [lindex $gomAtom_SN($element) 0]
	    set bgColour [eval lulColourFloat2Hex [lindex [show atom info $symbol] 2]]
	    set fgColour [lulGetVisibleForegroundColour $bgColour]
	    set name     [string tolower $symbol]
	    button $w.frame.main.line$line.right.$name -text $symbol \
		-fg $fgColour -bg $bgColour -width 2 -font $gomControlFont \
		-command "lulDisplayPeriodicTableElement $element"
	    pack   $w.frame.main.line$line.right.$name -side left
	}

	set elementCount [expr $elementCount + $count]
    }

    foreach info {{{Lanthanides} 57 15} {{Actinides}   89 15}} {
	incr line
	scan $info "%s %d %d" title start count
	set  title [lindex $info 0]

	frame $w.frame.line$line
	pack  $w.frame.line$line -side top -fill x
	label $w.frame.line$line.label -text $title:
	pack  $w.frame.line$line.label -side left
	frame $w.frame.line$line.atoms
	pack  $w.frame.line$line.atoms -side right

	for {set i 0} {$i < $count} {incr i} {
	    set element [expr $start + $i]
	    set symbol  [lindex $gomAtom_SN($element) 0]
	    set bgColour [eval lulColourFloat2Hex [lindex [show atom info $symbol] 2]]
	    set fgColour [lulGetVisibleForegroundColour $bgColour]
	    set name     [string tolower $symbol]
	    button $w.frame.line$line.atoms.$name -text $symbol \
		-fg $fgColour -bg $bgColour -width 2 -font $gomControlFont \
		-command "lulDisplayPeriodicTableElement $element"
	    pack   $w.frame.line$line.atoms.$name -side left
	}
    }
}

proc lulDisplayPeriodicTableElement { Element } {
    global gomControlFont
    global gomHelpFile
    global gomAtom_SN
    
    set w .periodicelem$Element
    catch {destroy $w}
    toplevel $w 
    wm title    $w "Element info ([lindex $gomAtom_SN($Element) 0])"
    wm iconname $w "Element info ([lindex $gomAtom_SN($Element) 0])"
    
    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss \
	-command "destroy $w" -font $gomControlFont
    button $w.buttons.apply   -text Apply \
	-state disabled -font $gomControlFont
    button $w.buttons.help    -text Help -command \
	"htmlShowHelp $gomHelpFile(elementinfo)" \
	-font $gomControlFont
    pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

    set info [show atom info $Element]

    frame  $w.symbol
    pack   $w.symbol -side top -anchor w
    label  $w.symbol.label -text "Symbol:" -width 15 -anchor w
    label  $w.symbol.value -text [lindex $gomAtom_SN($Element) 0]
    pack   $w.symbol.label $w.symbol.value -side left

    frame  $w.name
    pack   $w.name -side top -anchor w
    label  $w.name.label -text "Name:" -width 15 -anchor w
    label  $w.name.value -text [lindex $gomAtom_SN($Element) 1]
    pack   $w.name.label $w.name.value -side left

    frame  $w.element
    pack   $w.element -side top -anchor w
    label  $w.element.label -text "Atomic number:" -width 15 -anchor w
    label  $w.element.value -text $Element
    pack   $w.element.label $w.element.value -side left

    frame  $w.basisset
    pack   $w.basisset -side top -anchor w
    label  $w.basisset.label -text "Gaussian basis set:" -width 15
    label  $w.basisset.value -text [lindex $info 3]
    pack   $w.basisset.label $w.basisset.value -side left

    frame  $w.covar
    pack   $w.covar -side top -anchor w
    label  $w.covar.label -text "Covar radius:" -width 15 -anchor w
    entry  $w.covar.value -width 10
    button $w.covar.apply -text "Apply" \
	-command "define atom covar $Element \[$w.covar.value get\]"
    pack   $w.covar.label $w.covar.value $w.covar.apply -side left
    $w.covar.value insert 0 [lindex $info 4]

    labelframe $w.gaussian -text "Gaussian data"
    pack       $w.gaussian -side top -anchor w -fill x

    set gaussian [lindex $info 5]

    if {[llength $gaussian] > 0} {
	foreach {name value} $gaussian {
	    frame  $w.gaussian.$name
	    pack   $w.gaussian.$name -side top -anchor w
	    label  $w.gaussian.$name.label -text "$name:" -width 10 -anchor w
	    label  $w.gaussian.$name.value -text $value
	    pack   $w.gaussian.$name.label $w.gaussian.$name.value -side left
	}
    } else {
	label $w.gaussian.missing -text "Not available"
	pack  $w.gaussian.missing -side left
    }
}
