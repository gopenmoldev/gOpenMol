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
    gromos96 {{1 0} GROMOS96 {GROMOS96 files} {.dat .DAT}}
}

namespace eval Gromos96 {

###########################################################################
# PROC
proc SolveGROMOS96Atoms {FileName} {
    set i 0
    
    # open file ...
    set File_p [open $FileName r]
    if {$File_p == ""} {
	error "ERROR: can't open GROMOS96 file!"
    }

    while {![eof $File_p]} {

	gets $File_p InputData

	if {[string match POSITION* $InputData]} {
	    while {![eof $File_p]} {
		gets $File_p InputData

		# break condition		  
		if {[string match END* $InputData]} {
		    close $File_p
		    return $i
		}
		
		if {[string index $InputData 0] == "#" ||
		    [string index $InputData 1] == "#"} continue
		incr i
	    }
	}
    }
    # close file and return ...
    puts Done!
    close $File_p

    return $i
}

####################### import external coordinates #############
# PROC
# Input file name:                          FileName
# Read new structure or append to old list: Action (new / append)
proc Read {FileName Action} {

    set lulCoordAmp [show gromos96 coordamplifier]

    # check first the number of atoms
    set lulNumAtoms [SolveGROMOS96Atoms $FileName]

    # return if no atoms defined
    if {!$lulNumAtoms} {
	error "ERROR: no atoms in GROMOS96 file"
    }
    
    # open file ...
    set File_p [open $FileName r]
    if {$File_p == ""} {
	error "ERROR: can't open GROMOS96 file!"
    }

    # do the job ...
    puts "Reading GROMOS96 coordinate file '$FileName'..."

    define structure $FileName $lulNumAtoms $Action

    set StrucNum [show molstructures]

    # initialize atom counter
    set i 1
    
    while {![eof $File_p]} {

	gets $File_p InputData

	if {[string match POSITION* $InputData]} {
	    while {![eof $File_p]} {
		gets $File_p InputData
		if {[string index $InputData 0] == "#" ||
		    [string index $InputData 1] == "#"} continue

		# break condition		  
		if {[string match END* $InputData]} {
		    break
		}

		set ResNum   [lindex $InputData 0]
		set ResName  [lindex $InputData 1]
		set AtomName [lindex $InputData 2]
		set RunNum   [lindex $InputData 3]
		set Xc       [expr $lulCoordAmp * [lindex $InputData 4]]
		set Yc       [expr $lulCoordAmp * [lindex $InputData 5]]
		set Zc       [expr $lulCoordAmp * [lindex $InputData 6]]

		# put the values into the structure ...
		define atom    label    $AtomName   $i $StrucNum
		define residue label    $ResName    $i $StrucNum
		define segment label    GROM        $i $StrucNum
		define atom coordinates $Xc $Yc $Zc $i $StrucNum
		define atom charge      0.0         $i $StrucNum
		define atom resnumber   $ResNum     $i $StrucNum
		incr i
	    }
	} elseif {[string match BOX* $InputData]} {
	    while {![eof $File_p]} {
		gets $File_p InputData

		if {[string index $InputData 0] == "#" ||
		    [string index $InputData 1] == "#"} continue

		# break condition		  
		if {[string match END* $InputData]} {
		    break
		}
		scan $InputData "%f %f %f" Xbox Ybox Zbox 
		# Apply the coordinate amplifier
		set Xbox [expr $lulCoordAmp * $Xbox]
		set Ybox [expr $lulCoordAmp * $Ybox]
		set Zbox [expr $lulCoordAmp * $Zbox]
		define cell dimensions $Xbox $Ybox $Xbox
		puts "Defining cell size: A $Xbox, B $Ybox, C $Zbox"   
	    }
	}
    }
    # close file and return ...
    puts Done!
    close $File_p
    return 0
}
}; # end of namespace Gromos96
}; # end of namespace gom::FileFilters::AtomCoordinates
