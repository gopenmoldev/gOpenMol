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
    txyz      {{1 1} Tinker    {Tinker files}         {.txyz .TXYZ}}
}

namespace eval Txyz {

array unset txyzTinkerInputConnectivity

#####################################################################
# PROC
proc Read {FileName Action} {

    variable txyzTinkerInputConnectivity

    # open file ...
    set File_p [open $FileName r]
    if {$File_p == ""} {
	catch {lulErrorDialog {ERROR: can't open TinkerXYZ file!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open TinkerXYZ file!"
	}
	error "ERROR: can't open TinkerXYZ file!"
	return 1
    }

    # do the job ...
    puts "Reading Tinker XYZ coordinate file '$FileName'..."

    gets $File_p InputData

    set NumAtoms [lindex $InputData 0]

    define structure "$FileName" $NumAtoms $Action

    set StrucNum [show molstructures]
    set ConnectivityLines 0
    set Trigger 0

    for {set i 1} {$i <= $NumAtoms} {incr i} {

	gets $File_p InputData

	set AtomName [lindex $InputData 1]
	set Xc       [lindex $InputData 2]
	set Yc       [lindex $InputData 3]
	set Zc       [lindex $InputData 4]

	if {[llength $InputData] > 6} {
	    for {set ii 6} {$ii < [llength $InputData]} {incr ii} {
		set iatom [lindex $InputData $ii]
		set txyzTinkerInputConnectivity($ConnectivityLines) "$i $iatom"
		incr ConnectivityLines
	    }
	}

	# put the values into the structure ...
	define atom    label    $AtomName   $i $StrucNum
	define residue label    TINK        $i $StrucNum
	define segment label    TINK        $i $StrucNum
	define atom coordinates $Xc $Yc $Zc $i $StrucNum
	define atom charge      0.0         $i $StrucNum
	define atom resnumber   1           $i $StrucNum
	set Trigger 1
    }
    # close file and return ...
    close $File_p
    if {!$Trigger} {
	gomError "Could not find the atoms in the TINKER file"
	error "Could not find the atoms in the TINKER file"
    }
    puts Done!

    return 0
}

###################################################################
# this is for the TINKER coordinate files that has the connectivity
# included
#
# the array is set during the process when the coordinates are read in
# the connectivity can not be set at that stage because the connectivity
# will later be reset and recalculated
#
proc PostRead { FileName Action errorCode} {
    variable txyzTinkerInputConnectivity

    catch {
    for {set i 1} {$i <= [array size txyzTinkerInputConnectivity]} {incr i} {

	edit bond create * * [lindex $txyzTinkerInputConnectivity($i) 0] \
	                 * * [lindex $txyzTinkerInputConnectivity($i) 1];
    }
    }

    array unset txyzTinkerInputConnectivity
}

######################## export coordinates #####################
# PROC
proc Write {Struct FileName Mask} {

    set def_title  "Default title for: Tinker"

    set File_p [open "$FileName" w]

    if {$File_p == ""} {
	catch {lulErrorDialog {ERROR: can't open TINKER output file for writing!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open TINKER output file for writing!"
	}
	return 1
    }

    puts "Writing Tinker coordinate file '$FileName' ..."

    set TotAtoms [show numatoms $Struct]

    puts $File_p " $TotAtoms  $def_title"

    for {set i 1} {$i <= $TotAtoms } {incr i} {

	# apply the display mask...
	if {$Mask && [show atom displaystate $i $Struct] == 0} continue

	scan [show atom coordinates   $i $Struct] "%f %f %f" xc yc zc

	set Connectivity [lrange [show atom connectivity  $i $Struct] 1 end]

	puts $File_p " $i [format "%s" [show atom atomname      $i $Struct]] \
                     [format "%.5f" $xc] [format "%.5f" $yc] [format "%.5f" $zc] \
                     [show atom type          $i $Struct] \
                     $Connectivity"
    }

    close $File_p

    puts "Done!"

    return 0
}

}; # end of namespace Txyz
}; # end of namespace gom::FileFilters::AtomCoordinates
