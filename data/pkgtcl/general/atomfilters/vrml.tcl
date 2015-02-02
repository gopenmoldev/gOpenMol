##############################################################################
#                           Copyright (c) 1999 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Enhancements 2002 - 2004 by: Eero Häkkinen
##############################################################################

namespace eval gom::FileFilters::AtomCoordinates {

;
# see pkgtcl/atomfiles/main.tcl
array set coordTypes {
    vrml {{0 1} VRML {VRML files} {.wrl .WRL}}
}

namespace eval Vrml {

;
##########################################################################
# 
# This is a filter to create VRML files out of gOpenMol
#
##########################################################################

variable lulNumAtoms
variable lulTransX
variable lulTransY
variable lulTransZ

###########################################################################
# PROC
proc Write {Which FileName VisibleOnly} {

    variable lulNumAtoms
    variable lulAtomEntryFlag
    variable lulTransX
    variable lulTransY
    variable lulTransZ

    set lulNumAtoms [show numat $Which]
    if {$lulNumAtoms < 1} {
	gomError "No atoms available in structure #$Which. Read in structure first"
	return 1
    }   

    # open file ...
    set File_p [open $FileName w]
    if {$File_p == ""} {
	catch {gomError {ERROR: can't open VRML file!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open VRML file!"
	}
	error "ERROR: can't open VRML file!"
	return 1
    }
    # Now we know the number of atoms ...

    puts "Writing structure #$Which to a VRML file '$FileName'"

    # calculate the translation array

    set TransXYZ [show transl array]
    set lulTransX [lindex $TransXYZ 0]
    set lulTransY [lindex $TransXYZ 1]
    set lulTransZ [lindex $TransXYZ 2]

    set Zshift [show sizeofsystem]
    set ViewA  [expr [show viewangle]/180. * 3.14159266]

    # write first out the file header...

    puts $File_p "#VRML V2.0 utf8"
    puts $File_p "#"
    puts $File_p "# File generated from gOpenMol (http://laaksonen.csc.fi/gopenmol/)"
    puts $File_p "# [clock format [clock seconds] -format "%A %B %d %H:%M:%S %Y"]"
    puts $File_p "# Leif Laaksonen/CSC"
    puts $File_p "#"

    # start Group and children ...
    puts $File_p "Group \{"
    puts $File_p "  children \["
    puts $File_p "  Viewpoint \{"
    puts $File_p "   position 0.0 0.0 [expr 1.5 * $Zshift]"
    puts $File_p "   fieldOfView $ViewA"
    # viewpoint ends
    puts $File_p "  \}"

    set lulAtomEntryFlag 1

    # loop the atoms for CPK entries ...

    for {set i 1} {$i <= $lulNumAtoms} {incr i} {

	# apply the display mask...
	if {$VisibleOnly && [show atom displaystate $i $Struct] == 0} continue

	if {[show atom displaystat $i $Which]} {
	    PlotStickType $File_p $i $Which
	}

	if {[show atom cpkstate $i $Which]} {
	    set Coords [show atom coord $i $Which]
	    set Xc [expr [lindex $Coords 0] - $lulTransX]
	    set Yc [expr [lindex $Coords 1] - $lulTransY]
	    set Zc [expr [lindex $Coords 2] - $lulTransZ]
	    set vDW [show atom vdw $i 1]
	    set cpkScale [show atom cpkscale $i $Which]
	    set Color [show atom colo $i $Which]
	    # start elements
	    if {$lulAtomEntryFlag} {
		puts $File_p ","
	    }

	    puts $File_p "    Transform \{"
	    puts $File_p "      translation $Xc $Yc $Zc"
	    puts $File_p "      children Shape \{"
	    puts $File_p "        appearance Appearance \{"
	    puts $File_p "          material Material \{"
	    puts $File_p "           diffuseColor $Color"
	    puts $File_p "          \}"
	    puts $File_p "        \}"
	    puts $File_p "        geometry Sphere \{"
	    puts $File_p "# vDW value: $vDW, vDW scale: $cpkScale"
	    puts $File_p "           radius [expr $cpkScale * $vDW]"
	    puts $File_p "        \}"
	    puts $File_p "      \}"
	    puts $File_p "    \}"
	}
    }

    # children ends
    puts $File_p "  \]"
    # group ends
    puts $File_p "\}"

    # close file 
    close $File_p

    puts "Done!"
    return 0

}
###########################################################################
# PROC
proc PlotStickType {File_p AtomIndex StructureIndex} {

    variable lulAtomEntryFlag
    variable lulTransX
    variable lulTransY
    variable lulTransZ

    set AList [show atom connect $AtomIndex $StructureIndex]
    set Loop [lindex $AList 0]

    set licosphere [show licosphere]
    set licocylinder [show licocyl]

    set Coords [show atom coord $AtomIndex $StructureIndex]
    set Xc1 [expr [lindex $Coords 0] - $lulTransX]
    set Yc1 [expr [lindex $Coords 1] - $lulTransY]
    set Zc1 [expr [lindex $Coords 2] - $lulTransZ]
    set Color1 [show atom colo $AtomIndex $StructureIndex]

    if {[show atom licoricestate $AtomIndex $StructureIndex]} {

	if {$lulAtomEntryFlag} {
	    puts $File_p ","
	}

	puts $File_p "    Transform \{"
	puts $File_p "      translation $Xc1 $Yc1 $Zc1"
	puts $File_p "      children Shape \{"
	puts $File_p "        appearance Appearance \{"
	puts $File_p "          material Material \{"
	puts $File_p "           diffuseColor $Color1"
	puts $File_p "          \}"
	puts $File_p "        \}"
	puts $File_p "        geometry Sphere \{"
	puts $File_p "# lico sphere radius"
	puts $File_p "           radius $licosphere"
	puts $File_p "        \}"
	puts $File_p "      \}"
	puts $File_p "    \}"

	set lulAtomEntryFlag 1

    }

    if {$Loop < 1} {
	return
    }

    for {set j 1} {$j <= $Loop} {incr j} {
	
	set AListj [lindex $AList $j]

	if {$AListj > $AtomIndex} {
	    continue
	}

	if {[show atom displaystat $AListj $StructureIndex]} {
	    set Coords [show atom coord $AListj $StructureIndex]
	    set Xc2 [expr [lindex $Coords 0] - $lulTransX]
	    set Yc2 [expr [lindex $Coords 1] - $lulTransY]
	    set Zc2 [expr [lindex $Coords 2] - $lulTransZ]
	    set Color2 [show atom colo $AListj $StructureIndex]

	    set Xc12 [expr ($Xc1 + $Xc2) / 2.0]
	    set Yc12 [expr ($Yc1 + $Yc2) / 2.0]
	    set Zc12 [expr ($Zc1 + $Zc2) / 2.0]

	    if {$lulAtomEntryFlag} {
		puts $File_p ","
	    }

	    puts $File_p " Shape \{"
	    puts $File_p "    appearance Appearance \{"
	    puts $File_p "        material Material \{"
	    puts $File_p "           emissiveColor $Color1"
	    puts $File_p "           diffuseColor  $Color1"
	    puts $File_p "        \}"
	    puts $File_p "    \}"
	    puts $File_p "    geometry IndexedLineSet \{"
	    puts $File_p "       coord Coordinate \{"
	    puts $File_p "         point \["
	    puts $File_p "              $Xc1 $Yc1 $Zc1, $Xc12 $Yc12 $Zc12"
	    puts $File_p "         \]"
	    puts $File_p "       \}"
	    puts $File_p "       coordIndex \[0 , 1\]"
	    puts $File_p "    \}"
	    puts $File_p "  \},"
	    puts $File_p " Shape \{"
	    puts $File_p "    appearance Appearance \{"
	    puts $File_p "        material Material \{"
	    puts $File_p "           emissiveColor $Color2"
	    puts $File_p "           diffuseColor  $Color2"
	    puts $File_p "        \}"
	    puts $File_p "    \}"
	    puts $File_p "    geometry IndexedLineSet \{"
	    puts $File_p "       coord Coordinate \{"
	    puts $File_p "         point \["
	    puts $File_p "              $Xc12 $Yc12 $Zc12, $Xc2 $Yc2 $Zc2"
	    puts $File_p "         \]"
	    puts $File_p "       \}"
	    puts $File_p "       coordIndex \[0 , 1\]"
	    puts $File_p "    \}"
	    puts $File_p "  \}"

	    set lulAtomEntryFlag 1

	    if {[show atom licoricestate $AListj $StructureIndex]} {

		if {$lulAtomEntryFlag} {
		    puts $File_p ","
		}

		puts $File_p "    Transform \{"
		puts $File_p "      translation $Xc2 $Yc2 $Zc2"
		puts $File_p "      children Shape \{"
		puts $File_p "        appearance Appearance \{"
		puts $File_p "          material Material \{"
		puts $File_p "           diffuseColor $Color2"
		puts $File_p "          \}"
		puts $File_p "        \}"
		puts $File_p "        geometry Sphere \{"
		puts $File_p "# lico sphere radius"
		puts $File_p "           radius $licosphere"
		puts $File_p "        \}"
		puts $File_p "      \}"
		puts $File_p "    \}"

		set lulAtomEntryFlag 1

		if {$lulAtomEntryFlag} {
		    puts $File_p ","
		}

		set Xl [expr $Xc2 - $Xc1]
		set Yl [expr $Yc2 - $Yc1]
		set Zl [expr $Zc2 - $Zc1]

		set cLength [expr sqrt($Xl * $Xl + $Yl * $Yl + $Zl * $Zl) / 2.0]

		set Angles [gomVector2Angles $Xc1 $Yc1 $Zc1 $Xc2 $Yc2 $Zc2]

		set alpha [lindex $Angles 0]
		set beta  [lindex $Angles 1]

		puts $File_p "Transform \{"
		puts $File_p "      translation $Xc1 $Yc1 $Zc1"
		puts $File_p "      rotation 0.0   0.0  1.0 $beta"
		puts $File_p " children \["
		puts $File_p " Transform \{"
		puts $File_p "   rotation 0.0   1.0  0.0 $alpha"
		puts $File_p "   children \["
		puts $File_p "    Transform \{"
		puts $File_p "      translation 0.0 0.0 [expr $cLength / 2.0]"
		puts $File_p "      rotation 1.0   0.0  0.0 1.57079632679"
		puts $File_p "      children Shape \{"
		puts $File_p "        appearance Appearance \{"
		puts $File_p "          material Material \{"
		puts $File_p "           diffuseColor $Color1"
		puts $File_p "          \}"
		puts $File_p "        \}"
		puts $File_p "        geometry Cylinder \{"
		puts $File_p "# lico sphere and cylinder radius"
		puts $File_p "           radius $licocylinder"
		puts $File_p "           height $cLength"
		puts $File_p "        \}"
		puts $File_p "      \}"
		puts $File_p "    \}"
		puts $File_p "   \]"
		puts $File_p " \}"
		puts $File_p " \]"
		puts $File_p "\},"
		puts $File_p "Transform \{"
		puts $File_p "      translation $Xc12 $Yc12 $Zc12"
		puts $File_p "      rotation 0.0   0.0  1.0 $beta"
		puts $File_p " children \["
		puts $File_p " Transform \{"
		puts $File_p "   rotation 0.0   1.0  0.0 $alpha"
		puts $File_p "   children \["
		puts $File_p "    Transform \{"
		puts $File_p "      translation 0.0 0.0 [expr $cLength / 2.0]"
		puts $File_p "      rotation 1.0   0.0  0.0 1.57079632679"
		puts $File_p "      children Shape \{"
		puts $File_p "        appearance Appearance \{"
		puts $File_p "          material Material \{"
		puts $File_p "           diffuseColor $Color2"
		puts $File_p "          \}"
		puts $File_p "        \}"
		puts $File_p "        geometry Cylinder \{"
		puts $File_p "# lico sphere and cylinder radius"
		puts $File_p "           radius $licocylinder"
		puts $File_p "           height $cLength"
		puts $File_p "        \}"
		puts $File_p "      \}"
		puts $File_p "    \}"
		puts $File_p "   \]"
		puts $File_p " \}"
		puts $File_p " \]"
		puts $File_p "\}"

		set lulAtomEntryFlag 1

	    }

	    set lulAtomEntryFlag 1

	}
    }
}

}; # end of namespace Vrml
}; # end of namespace gom::FileFilters::AtomCoordinates
