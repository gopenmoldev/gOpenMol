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
    turbomole {{1 0} Turbomole {Turbomole files} {.txt .TXT}}
}

##########################################################################
# NAMESPACE Turbomole
#
namespace eval Turbomole {

variable lulNumAtoms
variable lulTurbomoleLogFileName
variable xc
variable yc
variable zc
variable atom
variable setnr 0

############################################################################
# PROC
proc Read { FileName Action } {

    variable setnr

    if {$FileName == ""} {
	gomError "file name is not supplied"
	return 1
    }

    set f [open $FileName r]

    if {$f == ""} {
	catch {lulErrorDialog {ERROR: can't open Turbomole coordinate file for reading!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open Turbomole file for reading!"
	}
	return 1
    }


    # change to current directory
    lulChangeDirectory $FileName

    # scan the coordinate file to get the atoms ...

    set setnr 0
    while {![eof $f]} { 
	gets $f Text
	update idletasks
	# look for "$coord"
	if {[string match "\$coord*" $Text]} {
	    getCoordinates $f $FileName $Action
	} elseif {[string match "\$grad*cartesian gradients*" $Text]} {
	    getCoordsGrads $f $FileName $Action
	}
        
    }

    close $f

    puts "Done!"

    return 0
    # end of proc
}
#########################################################################
# just the coordinates
proc getCoordinates {file_p FileName Action} {

    variable xc
    variable yc
    variable zc
    variable atom
    variable setnr

    set au2a  0.52917715
    set  i 0
    incr setnr
    gomPrint "Read set #: $setnr"

    while {![eof $file_p]} { 
	gets $file_p Text
	set Text [string trim $Text]
	if {($Text == "") || ([string match "\$*" $Text])} break
	incr i
	set xc($i)   [expr $au2a * [lindex $Text 0]]
	set yc($i)   [expr $au2a * [lindex $Text 1]]
	set zc($i)   [expr $au2a * [lindex $Text 2]]
	set atom($i) [string toupper [lindex $Text 3] 0 0]
    }

    set NumAtoms $i

    define structure "$FileName" $NumAtoms $Action

    set StrucNum [show molstructures]

    for {set i 1} {$i <= $NumAtoms} {incr i} {
	# put the values into the structure ...
	define atom    label    $atom($i)                $i $StrucNum
	define residue label    tbmo                     $i $StrucNum
	define segment label    tbmo                     $i $StrucNum
	define atom coordinates $xc($i) $yc($i) $zc($i)  $i $StrucNum
	define atom charge      0.0                      $i $StrucNum
	define atom resnumber   1                        $i $StrucNum

	unset atom($i)
	unset xc($i)
	unset yc($i)
	unset zc($i)
    }

}

#
#########################################################################
# coordinates and gradients
proc getCoordsGrads {file_p FileName Action} {

    variable xc
    variable yc
    variable zc
    variable atom
    variable setnr

    set au2a  0.52917715
    set FirstTime 0

    set fx [open [file rootname $FileName]_traj.xmol w]
    if {$fx == ""} {
	gomError "can't open xmol trajectory output file"
	return
    }

    gomPrint "writing coordinate sets as a trajectory into file '[file rootname $FileName]_traj.xmol'"
    set setnr 0
    
    while {![eof $file_p]} { 
	incr setnr
	gomPrint "Reading set #: $setnr"
	
	set i 0
	while {![eof $file_p]} { 
	    gets $file_p Text
	    set TextLabel [string trim $Text]
	    if {($TextLabel == "") || ([string match "\$*" $TextLabel])} {
            close $fx
            gomPrint "Number of sets read: [expr $setnr - 1]"        
            return
        }
	    # check to see if it is the label
	    if {([string match "cycle =*" $TextLabel])} {
                set TextLabelSave $Text
                gets $file_p Text
	    } 
	    # atoms
	    if {[llength $Text] == 4} {
		incr i
		set xc($i)   [expr $au2a * [lindex $Text 0]]
		set yc($i)   [expr $au2a * [lindex $Text 1]]
		set zc($i)   [expr $au2a * [lindex $Text 2]]
		set atom($i) [string toupper [lindex $Text 3] 0 0]
		# gradients
	    } else {
		set NumAtoms $i
		if {!$FirstTime} {
		    define structure "$FileName" $NumAtoms $Action

		    set FirstTime 1
		}
#
# This is a patch for Linux (25.09.2002)
#
        regsub  -all {\-\.} $Text {-0.} Text
        regsub  -all {D}    $Text {E}   Text
# end of patch
		set xg(1)   [lindex $Text 0]
		set yg(1)   [lindex $Text 1]
		set zg(1)   [lindex $Text 2]
		
		for {set j 2} {$j <= $NumAtoms} {incr j} {
		    gets $file_p Text
#
# This is a patch for Linux (25.09.2002)
#
             regsub  -all {\-\.} $Text {-0.} Text
             regsub  -all {D}    $Text {E}   Text
# end of patch
		    set xg($j)   [lindex $Text 0]
		    set yg($j)   [lindex $Text 1]
		    set zg($j)   [lindex $Text 2]
		}

		# write out gradients
		set GradFile "[file rootname $FileName]_$setnr\_grad.txt"
		set fg [open $GradFile w]
		if {$fg == ""} {
		    gomError "can't open gradient output file"
		    return
		}
		puts $fg "Dummy"
		puts $fg "Dummy"
		puts $fg "Dummy"
		puts $fg "Dummy"
		puts $fg "Dummy"
		for {set j 1} {$j <= $NumAtoms} {incr j} {
		    puts $fg [format "%f %f %f %f %f %f" $xc($j) $yc($j) $zc($j) $xg($j) $yg($j) $zg($j)]
		}
		close $fg
		# done
		set StrucNum [show molstructures]

		puts $fx $NumAtoms
		puts $fx "$TextLabelSave"

		for {set i 1} {$i <= $NumAtoms} {incr i} {
		    # put the values into the structure ...
		    define atom    label    $atom($i)                $i $StrucNum
		    define residue label    tbmo                     $i $StrucNum
		    define segment label    tbmo                     $i $StrucNum
		    define atom coordinates $xc($i) $yc($i) $zc($i)  $i $StrucNum
		    define atom charge      0.0                      $i $StrucNum
		    define atom resnumber   1                        $i $StrucNum
		    # xmol file
		    puts $fx "$atom($i) $xc($i) $yc($i) $zc($i)"

		    unset atom($i)
		    unset xc($i)
		    unset yc($i)
		    unset zc($i)
        }

                break
	    }
	}
    }
}

}; # end of namespace Turbomole
}; # end of namespace gom::FileFilters::AtomCoordinates
