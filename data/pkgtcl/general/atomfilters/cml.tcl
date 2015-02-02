##############################################################################
#                        Copyright (c) 2002 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded 2002 - 2004 by: Eero Häkkinen
##############################################################################
##########################################################################
#
# This file contains import and export filters for CML coordinate files.
#
##########################################################################

namespace eval gom::FileFilters::AtomCoordinates {

;
# see pkgtcl/atomfiles/main.tcl
array set coordTypes {
    cml {{1 1} CML {CML files} {.xml .XML}}
}

namespace eval Cml {

;
###########################################################################
# PROC
proc Write {Which FileName VisibleOnly} {
    set NAtoms [show numatoms $Which]
    if {$NAtoms < 1} {
	gomError "No atoms available in structure #$Which. Read in structure first"
	return 1
    }   

    # open file ...
    set File_p [open $FileName w]

    puts "Writing structure #$Which to a CML file '$FileName'"

    puts $File_p "<?xml version=\"1.0\" ?>"
    puts $File_p "<cml xmlns:cml=\"http://www.xml-cml.org/dtd/cml1_0_1.dtd\" id=\"c1\">"
    puts $File_p "  <molecule id=\"m$Which\">"

    # loop the atoms for atom tags ...
    puts $File_p "    <atomArray>"
    for {set i 1} {$i <= $NAtoms} {incr i} {
	# apply the display mask...
	if {$VisibleOnly && [show atom displaystate $i $Which] == 0} continue

	puts $File_p "      <atom id=\"a$i\" title=\"[
	show atom atomname $i $Which]\">"
	puts $File_p "\t<string title=\"chain\"   dictRef=\"cml:chain\"  >[
	show atom segmentname $i $Which]</string>"
	puts $File_p "\t<string title=\"resCode\" dictRef=\"cml:resType\">[
	show atom residuename $i $Which]</string>"
	puts $File_p "\t<string title=\"hetId\"   dictRef=\"cml:resId\"  >[
	show atom resnumber $i $Which]</string>"
	puts $File_p "\t<string title=\"atId\"    dictRef=\"cml:atomId\" >[
	show atom atomname $i $Which]</string>"
	puts $File_p "\t<string builtin=\"elementType\">[
	show atom type $i $Which]</string>"
	puts $File_p "\t<float builtin=\"formalCharge\">[
	show atom nuclearcharge $i $Which]</float>"
	puts $File_p "\t<coordinate3 builtin=\"xyz3\">[
	show atom coordinates $i $Which]</coordinate3>"
	puts $File_p "\t<float title=\"bFactor\" dictRef=\"cml:bFactor\">[
	show atom bvalue $i $Which]</float>"
	puts $File_p "      </atom>"
    }
    puts $File_p "    </atomArray>"

    # loop the atoms for bond tags ...
    puts $File_p "    <bondArray>"
    for {set i 1} {$i <= $NAtoms} {incr i} {
	# apply the display mask...
	if {$VisibleOnly && [show atom displaystate $i $Which] == 0} continue

	set connections [show atom connectivity $i $Which]

	foreach j $connections {
	    # Write bonds only once.
	    if {$j < $i} continue
	    # apply the display mask...
	    if {[show atom displaystate $j $Which] == 0 && $Mask} continue

	    puts $File_p "      <bond>"
	    puts $File_p "\t<string builtin=\"atomRef\">a$i</string>"
	    puts $File_p "\t<string builtin=\"atomRef\">a$j</string>"
	    puts $File_p "      </bond>"
	}
    }
    puts $File_p "    </bondArray>"
    puts $File_p "  </molecule>"
    puts $File_p "</cml>"

    # close file 
    close $File_p

    puts "Done!"
    return 0
}

namespace eval ReadUtilities {

proc GetElementName {name} {
    regsub "^.*:" $name "" name
    return $name
}

###########################################################################
# Called for every opening tags.
###########################################################################
proc ElemStart {name attlist args} {
    variable cmlInputState
    switch -regexp $name {
	"(.*:|)molecule" {
	    # Save attributes.
	    set cmlInputState(moleculeData) $attlist
	    if {$cmlInputState(atomIndex) >= 0} {
		set cmlInputState(atomIndex) 0
	    }
	    set cmlInputState(atomCount) 0
	}
	"(.*:|)(atom|atomArray)" {
	    # Save attributes.
	    set cmlInputState([GetElementName $name]Data) $attlist
	}
    }
    # Put parametres to stacks.
    lappend cmlInputState(elements) $name
    lappend cmlInputState(atts) $attlist
}

###########################################################################
# Called for data between opening and closing tags.
###########################################################################
proc CharData {data} {
    variable cmlInputState
    # Omit empty data.
    if {"" == $data} {
	return
    }
    # Handle only direct sub elements of
    # atom and atomArray tags.
    if {[llength $cmlInputState(elements)] < 2} {
	return
    }
    set lastElement [lindex $cmlInputState(elements) end-1]
    switch -regexp $lastElement {
	"(.*:|)(atom|atomArray)" {
	    set attlist [lindex $cmlInputState(atts) end]
	    foreach {var value} $attlist {
		if {[string equal "builtin" $var] ||
		    [string equal "title" $var]} {
		    # Append the attribute to the attribute list.
		    lappend cmlInputState([
			GetElementName $lastElement]Data) $value "$data"
		    break
		}
	    }
	}
    }
}

###########################################################################
# Called for closing "molecule" tag.
###########################################################################
proc EndMolecule {attlist} {
    global lulIdentifyAtoms
    variable cmlInputState

    if {$cmlInputState(atomIndex) >= 0 } {
	# We are reading atoms. Molecules are already created.
	set struct [lindex $cmlInputState(moleculeIndeces) 0]

	set cmlInputState(moleculeIndeces) [
	    lrange $cmlInputState(moleculeIndeces) 1 end]
	return
    }
    if {[lsearch $cmlInputState(elements) "molecule"] >= 0 } {
	# This isn't top level molecule.
	return
    }
    if {$cmlInputState(atomCount) <= 0 } {
	# gOpenMol doesn't accept empty structures.
	return
    }

    # Find molecule title.
    set title ""
    foreach {var value} $attlist {
	switch $var {
	    "id" {if {"" == $title} {set title $value}}
	    "title" {set title $value}
	}
    }
    if {"" == $title} {
	set title "Molecule [expr [
      llength $cmlInputState(moleculeIndeces)]+1]"
    }
    # Create new structure.
    lappend cmlInputState(moleculeIndeces) [
	define structure $title $cmlInputState(atomCount) \
	    $cmlInputState(action)]
}

###########################################################################
# Action depends on value of $cmlInputState(atomIndex):
#  -value < 0 :   increments cmlInputState(atomCount)
#  -value >= 0 :  appends atom data to cmlInputState(atoms)
###########################################################################
proc AddAtom {seg res resN elem title x y z charge occupancy} {
    variable cmlInputState
    if {"" == $x || "" == $y} {
	# Atom doesn't have coordinates.
	return
    }
    incr cmlInputState(atomCount)
    if {$cmlInputState(atomIndex) < 0 } {
	# We are just counting atoms.
	return
    }

    if {"" == $seg} {set seg "Seg1"}
    if {"" == $res} {set res "Res1"}
    if {"" == $resN} {set resN 1}
    if {"" == $elem} {set elem "000"}
    if {"" == $title} {set title $elem}
    if {"" == $x} {set x 0.0}
    if {"" == $y} {set y 0.0}
    if {"" == $z} {set z 0.0}
    if {"" == $charge} {set charge 0.0}
    if {"" == $occupancy} {set occupancy 1.0}

    incr cmlInputState(atomIndex)
    set structIndex [lindex $cmlInputState(moleculeIndeces) 0]
    set atomIndex $cmlInputState(atomIndex)

    eval [concat {define segment label} \
	      $seg $atomIndex $structIndex]    
    eval [concat {define residue label} \
	      $res $atomIndex $structIndex]    
    eval [concat {define atom resnumber} \
	      $resN $atomIndex $structIndex]
    # Specify atom by its element type.
    eval [concat {define atom label} \
	      $elem $atomIndex $structIndex]
    eval [concat {define atom coordinates} \
	      $x $y $z $atomIndex $structIndex]
    eval [concat {define atom charge} \
	      $charge $atomIndex $structIndex]
    # Fill atom information according to the atom label.
    fill atom $atomIndex $structIndex
    # Set real atom label - not the element type.
    eval [concat {define atom label} \
	      $title $atomIndex $structIndex]
}

###########################################################################
# Called for closing "atomArray" tag.
###########################################################################
proc EndAtomArray {attlist} {
    set segList [list]
    set resList [list]
    set resNList [list]
    set elemList [list]
    set titleList [list]
    set xList [list]
    set yList [list]
    set zList [list]
    set chargeList [list]
    set occupancyList [list]
    foreach {builtin data} $attlist {
	switch $builtin {
	    "x2" {if {[llength $xList] == 0} {set xList [split $data]}}
	    "y2" {if {[llength $yList] == 0} {set yList [split $data]}}
	    "x3" {set xList [split $data]}
	    "y3" {set yList [split $data]}
	    "z3" {set zList [split $data]}
	    "elementType" {set elemList [split $data]}
	    "title" {set titleList [split $data]}
	    "formalCharge" {set chargeList [split $date]}
	    "occupancy" {set occupancyList [split $data]}
	}
    }
    foreach \
	seg    $segList \
	res    $resList \
	resN   $resNList \
	elem   $elemList \
	title  $titleList \
	x      $xList \
	y      $yList \
	z      $zList \
	charge $chargeList \
	occupancy $occupancyList {
	AddAtom $seg $res $resN $elem $title $x $y $z $charge $occupancy
    }
}

###########################################################################
# Called for closing "atom" tag.
###########################################################################
proc EndAtom {attlist} {
    set seg ""
    set res ""
    set resN ""
    set elem ""
    set title ""
    set x ""
    set y ""
    set z ""
    set charge ""
    set occupancy ""
    foreach {builtin data} $attlist {
	switch $builtin {
	    # Set 2D coordinates only if 3D coordinates
	    # aren't set yet.
	    "x2" {if {"" == $x} {set x $data}}
	    "y2" {if {"" == $y} {set y $data}}
	    "x3" {set x $data}
	    "y3" {set y $data}
	    "z3" {set z $data}
	    "xy2" {
		if {"" == $x$y} {
		    scan $data "%f %f" x y
		}
	    }
	    "xyz3" {scan $data "%f %f %f" x y z}
	    "elementType" {set elem $data}
	    "title" {set title $data}
	    "resCode" {set res $data}
	    "chain" {set seg $data}
	    "hetId" {set resN $data}
	    "formalCharge" {set charge $data}
	    "occupancy" {set occupancy $data}
	}
    }
    AddAtom $seg $res $resN $elem $title \
	$x $y $z $charge $occupancy
}

###########################################################################
# Called for every closing tags.
###########################################################################
proc ElemEnd {name args} {
    variable cmlInputState
    set attlist [lindex $cmlInputState(atts) end]
    set cmlInputState(elements) [
	lrange $cmlInputState(elements) 0 end-1]
    set cmlInputState(atts) [
	lrange $cmlInputState(atts) 0 end-1]
    switch -regexp $name {
	"(.*:|)molecule" {
	    EndMolecule $attlist
	}
	"(.*:|)atomArray" {
	    EndAtomArray $cmlInputState(atomArrayData)
	}
	"(.*:|)atom" {
	    EndAtom $cmlInputState(atomData)
	}
    }
}

# namespace ReadUtilities
}

###########################################################################
# PROC
proc Read {FileName Append} {
    global errorInfo

    set ns [namespace current]::ReadUtilities

    puts "Reading CML file '$FileName':"
    puts "    preparing..."
    package require xml

    upvar "#0" ${ns}::cmlInputState cmlInputState
    array set cmlInputState [list \
	"elements"        {} \
	"atts"            {} \
	"moleculeData"    {} \
	"atomData"        {} \
	"atomArrayData"   {} \
	"moleculeIndeces" {} \
	"atomIndex"       -1 \
	"atomCount"       0 \
	"action"          $Append]

    # open file ...
    set File_p [open $FileName]

    # Create parser.
    set parser [xml::parser \
		    -elementstartcommand  ${ns}::ElemStart \
		    -characterdatacommand ${ns}::CharData \
		    -elementendcommand    ${ns}::ElemEnd \
		    -final 0]

    puts "    reading structure information..."
    # Create structures.
    $parser configure -final 0
    set lnro 1
    if { [catch {
	while {![eof $File_p]} {
	    $parser parse [read $File_p 1024]
	    incr lnro
	}
	$parser configure -final 1
	$parser parse ""
	if {[llength $cmlInputState(moleculeIndeces)] == 0} {
	    # We don't have any molecule tag but
	    # we might still have some atom tags.
	    # Create a new structure if needed.
	    ${ns}::EndMolecule {}
	}
    } errRet] } {
	close $File_p
	error "parse error occured at line $lnro while reading CML file '$FileName'\n$errRet" $errorInfo
    }
    if {[llength $cmlInputState(moleculeIndeces)] == 0} {
	close $File_p
	error "no atoms in XML file '$FileName'"
    }

    puts "    reading atom information..."
    # Create atoms.
    seek $File_p 0 start
    set cmlInputState(atomIndex) 0
    $parser reset
    $parser configure -final 0
    while {![eof $File_p]} {
	$parser parse [read $File_p 1024]
    }
    $parser configure -final 1
    $parser parse ""
    
    if {[llength $cmlInputState(moleculeIndeces)] > 0} {
	${ns}::EndMolecule {}
    }

    $parser free
    close $File_p
    puts "Done!"

    # Atoms are already identified.
    set ::lulIdentifyAtoms "no"
    return 0
}

proc PostRead { FileName Action errorCode } {
    set ::lulIdentifyAtoms "yes"
}

}; # end of namespace Cmd
}; # end of namespace gom::FileFilters::AtomCoordinates
