##############################################################################
#                           Copyright (c) 1994 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
##############################################################################        
namespace eval gom::FileFilters::AtomCoordinates {

;
# see pkgtcl/atomfiles/main.tcl
array set coordTypes {
    adf {{1 0} ADF {ADF files} {.log .LOG .out .OUT}}
}

namespace eval Adf {

;
####################### import external coordinates #############
#
# This filter program reads an ASF output (log) file and reads
# on the atom and coordinate information.
#
# It also produces an XMOL file of all the coordinate sets
# that is found in the file. This can then be displayed
# using the trajectory tools.
#################################################################
# PROC 
proc Read {FileName Action} {
    #
    # Read ADF log file 
    #
    # open file ...
    puts "$FileName $Action"
    set File_p [open $FileName r]

    # get number of atoms in file ...
    set sets    0
    set Trigger 0

    while {![eof $File_p]} { 

	gets $File_p Text

	if {[string match "*Coordinates (Cartesian)*" $Text]} {

	    while {![eof $File_p]} {
		gets $File_p Text
		if {[string match "*-------------------*" $Text]} {
		    incr sets
		    break
		}
	    }

	    set natoms 0
	    while {![eof $File_p]} {
		gets $File_p Text

		if {[lindex $Text 1] == "XX"} {
		    continue
		} elseif {[string match "*-------------------*" $Text]} {
		    break
		}
		incr natoms
	    }

	}
    }
    # rewind
    seek $File_p 0 start

    gomPrint "Found '$sets' coordinate sets and '$natoms' atoms in each set"

    if {$natoms < 1} {
	gomError "can't find atoms in the log file"
	return
    }

    define structure "$FileName" $natoms $Action

    set StrucNum [show molstructures]

    set XmolOut "[file rootname $FileName].xmol"

    set File_px [open $XmolOut w]
    if {$File_px == ""} {
	catch {lulErrorDialog {ERROR: can't open XMOL output file!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open XMOL output file!"
	}
	error "ERROR: can't open XMOL output file!"
	return
    }

    # do the job ...
    gomPrint "Reading ADF log coordinate file '$FileName'..."
    gomPrint "Writing all found coordinate sets into a XMOL file '$XmolOut'"
    gomPrint "!!!Last coordinate set found will be the one displayed!!!"

    lulMessageDialog \
	"Reading ADF log coordinate file '$FileName'...\n\
Writing all found coordinate sets into a XMOL file '$XmolOut'\n\
!!!Last coordinate set found will be the one displayed!!!"

    set j 0

    while {![eof $File_p]} { 

	gets $File_p Text

	if {[string match "*Coordinates (Cartesian)*" $Text]} {

	    while {![eof $File_p]} {
		gets $File_p Text
		if {[string match "*-------------------*" $Text]} {
		    incr j
		    puts $File_px "$natoms"
		    puts $File_px "Conformation/structure #$j" 
		    break
		}
	    }

	    set i 0
	    while {![eof $File_p]} {
		gets $File_p Text

		if {[lindex $Text 1] == "XX"} {
		    continue
		} elseif {[string match "*-------------------*" $Text]} {
		    break
		}
		incr i

		# read now the atoms ...

		set ResN   1
		set Res    "ADF"
		set Atm    [lindex $Text 1]
		set Seg    "ADF"
		set Xc     [lindex $Text 5]
		set Yc     [lindex $Text 6]
		set Zc     [lindex $Text 7]
		set Ch     0.0

		# put the values into the structure ...
		define atom    label    $Atm        $i $StrucNum
		define residue label    $Res        $i $StrucNum
		define segment label    $Seg        $i $StrucNum
		define atom coordinates $Xc $Yc $Zc $i $StrucNum
		define atom charge      $Ch         $i $StrucNum
		define atom resnumber   $ResN       $i $StrucNum
		# make the XMOLE file
		puts $File_px "$Atm $Xc $Yc $Zc"
		set Trigger 1
	    }
	}
    }
    # close file and return ...
    close $File_p
    close $File_px
    if {!$Trigger} {
	gomError "Could not find the atoms in the ADF file"
	error "Could not find the atoms in the ADF file"
    }
    puts Done!

}
}; # end of namespace Adf
}; # end of namespace gom::FileFilters::AtomCoordinates
