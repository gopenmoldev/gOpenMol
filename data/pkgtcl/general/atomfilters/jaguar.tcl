##############################################################################
#                           Copyright (c) 1994 - 2002 by:
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
    jaguar {{1 0} Jaguar {Jaguar files} {.out .OUT}}
}

namespace eval Jaguar {

;
#################################################################
# This filter program reads an JAGUAR output (log) file and reads
# on the atom and coordinate information.
#
# It also produces an XMOL file of all the coordinate sets
# that is found in the file. This can then be displayed
# using the trajectory tools.
#################################################################
# PROC 
proc Read {FileName Action} {
    #
    # Read JAGUAR log file 
    #
    # open file ...
    puts "$FileName $Action"
    set File_p [open $FileName r]
    if {$File_p == ""} {
	error "ERROR: can't open Jaguar log file!"
    }

    # get number of atoms in file ...
    set sets    0
    set Trigger 0
    set natoms  0

    while {![eof $File_p]} { 

	gets $File_p Text

	if {[string match "*atom            x                 y                 z*" $Text]} {

	    incr sets
	    set natoms 0
	    while {![eof $File_p]} {
		gets $File_p Text

		if {[string trim $Text] == ""} {
		    break
		}
		incr natoms
	    }

	}
    }
    # rewind
    seek $File_p 0 start

    if {$natoms < 1} {
	close $File_p
	error "can't find atoms in the log file"
    }

    gomPrint "Found '$sets' coordinate sets and '$natoms' atoms in each set"

    define structure "$FileName" $natoms $Action

    set StrucNum [show molstructures]

    set XmolOut "[file rootname $FileName].xmol"

    set File_px [open $XmolOut w]
    if {$File_px == ""} {
	error "ERROR: can't open XMOL output file!"
	return
    }

    # do the job ...
    gomPrint "Reading Jaguar output file '$FileName'..."
    gomPrint "Writing XMOL file '$XmolOut'"

    set j 0

    while {![eof $File_p]} { 

	gets $File_p Text

	if {[string match "*atom            x                 y                 z*" $Text]} {

	    incr j

	    gomPrint "########################################"
	    gomPrint "#    M O L E C U L A R  S E T # $j"
	    gomPrint "########################################"
	    gomPrint $Text

	    puts $File_px "$natoms"
	    puts $File_px "Conformation/structure #$j" 
	    
	    set i 0
	    while {![eof $File_p]} {
		gets $File_p Text

		if {[string trim $Text] == ""} {
		    break
		}
		incr i

		gomPrint $Text

		# read now the atoms ...

		set ResN   1
		set Res    "JAG"
		set Atm    [string toupper [lindex $Text 0] 0 0]
		set Seg    "JAG"
		set Xc     [lindex $Text 1]
		set Yc     [lindex $Text 2]
		set Zc     [lindex $Text 3]
		set Ch     0.0
		# put the values into the structure ...
		define atom    label    $Atm        $i $StrucNum
		define residue label    $Res        $i $StrucNum
		define segment label    $Seg        $i $StrucNum
		define atom coordinates $Xc $Yc $Zc $i $StrucNum
		define atom charge      $Ch         $i $StrucNum
		define atom resnumber   $ResN       $i $StrucNum

		# make the XMOL file
		puts $File_px "$Atm $Xc $Yc $Zc"
		set Trigger 1
	    }
	}
    }
    # close file and return ...
    close $File_p
    close $File_px

    gomPrint "########################################"
    gomPrint "Coordinate set(s) (total #'$j') also saved in file '$XmolOut'"

    if {!$Trigger} {
	gomError "Could not find the atoms in the Jaguar file"
	error "Could not find the atoms in the Jaguar file"
    }
    puts Done!
}
}; # end of namespace Jaguar
}; # end of namespace gom::FileFilters::AtomCoordinates
