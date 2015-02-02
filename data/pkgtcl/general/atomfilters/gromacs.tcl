##############################################################################
#                           Copyright (c) 1994 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Enhancements by: Eero Häkkinen
##############################################################################        
namespace eval gom::FileFilters::AtomCoordinates {

;
# see pkgtcl/atomfiles/main.tcl
array set coordTypes {
    gromacs {{1 0} GROMACS {GROMACS files} {.gro .GRO}}
}

namespace eval Gromacs {

####################### import external coordinates #############
# PROC
# Input file name:                          FileName
# Read new structure or append to old list: Action (new / append)
proc Read {FileName Action} {

    # open file ...
    set File_p [open $FileName r]
    if {$File_p == ""} {
	catch {lulErrorDialog {ERROR: can't open GROMACS file!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open GROMACS file!"
	}
	error "ERROR: can't open GROMACS file!"
	return
    }

    # do the job ...
    puts "Reading GROMACS coordinate file '$FileName'..."

    # starting from title line
    gets $File_p InputData
    gomPrint "Title line: $InputData"

    # number of atoms
    gets $File_p InputData
    scan $InputData "%d" lulNumAtoms

    if {$lulNumAtoms < 1} {
	gomError "Number of atoms < 1"
	return
    }

    define structure "$FileName" $lulNumAtoms $Action

    set StrucNum [show molstructures]

    # initialize atom counter
    
    for {set i 1} {$i <= $lulNumAtoms} {incr i} {

	gets $File_p InputData

	set ResNum   [string trim [string range $InputData 0   4]]
	set ResName  [string trim [string range $InputData 5   9]]
	set AtomName [string trim [string range $InputData 10 14]]
	set RunNum   [string range $InputData 15 19]
	set Xc       [expr 10.0 * [string trim [string range $InputData 20 27]]]
	set Yc       [expr 10.0 * [string trim [string range $InputData 28 35]]]
	set Zc       [expr 10.0 * [string trim [string range $InputData 36 44]]]
	# put the values into the structure ...
	define atom    label    $AtomName   $i $StrucNum
	define residue label    $ResName    $i $StrucNum
	define segment label    GRMA        $i $StrucNum
	define atom coordinates $Xc $Yc $Zc $i $StrucNum
	define atom charge      0.0         $i $StrucNum
	define atom resnumber   $ResNum     $i $StrucNum 
    }
    gets $File_p InputData
    scan $InputData "%f %f %f" Xbox Ybox Zbox 

    set Xbox [expr 10.0 * $Xbox]
    set Ybox [expr 10.0 * $Ybox]
    set Zbox [expr 10.0 * $Zbox]
    define cell dimensions $Xbox $Ybox $Zbox
    gomPrint "Defining cell size: A $Xbox, B $Ybox, C $Zbox"   

    # close file and return ...
    puts Done!
    close $File_p
}
}; # end of namespace Gromacs
}; # end of namespace gom::FileFilters::AtomCoordinates
