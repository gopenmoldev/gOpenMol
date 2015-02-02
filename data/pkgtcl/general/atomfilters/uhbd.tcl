##############################################################################
#                           Copyright (c) 1994 - 2004 by:
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
    uhbd      {{1 1} UHBD      {UHBD files}           {.qcd .QCD}}
}

namespace eval Uhbd {

;
####################### import external coordinates #############
# PROC 
# Input file name:                          FileName
# Read new structure or append to old list: Action (new / append)

proc Read {FileName Action} {
    #
    # Read UHBD qcd coordinate file 
    #
    # open file ...
    set File_p [open $FileName r]
    if {$File_p == ""} {
	catch {lulErrorDialog {ERROR: can't open UHBD qcd file!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open UHBD qcd file!"
	}
	error "ERROR: can't open UHBD qcd file!"
	return 1
    }

    # do the job ...
    puts "Reading UHBD qcd coordinate file '$FileName'..."

    set NumAtoms [GetUHBDqcdAtoms $File_p]

    define structure "$FileName" $NumAtoms $Action

    set StrucNum [show molstructures]

    set i 0

    while {![eof $File_p]} {

	gets $File_p Text

	set Option [lindex $Text 0]

	if {([string trim $Text] == "") || ($Option != "ATOM")} {
	    continue}

        incr i

        set ResN   [lindex $Text 1]
        set Res    [lindex $Text 2]
        set Atm    [lindex $Text 3]
        set Seg    "UHBD"
        set Xc     [lindex $Text 4]
        set Yc     [lindex $Text 5]
        set Zc     [lindex $Text 6]
        set Ch     [lindex $Text 7]

	# put the values into the structure ...
	define atom    label    $Atm        $i $StrucNum
	define residue label    $Res        $i $StrucNum
	define segment label    $Seg        $i $StrucNum
	define atom coordinates $Xc $Yc $Zc $i $StrucNum
	define atom charge      $Ch         $i $StrucNum
	define atom resnumber   $ResN       $i $StrucNum

    }
    # close file and return ...
    puts Done!
    close $File_p
    return 0
}

########################################################################
# PROC
proc GetUHBDqcdAtoms { File_p } {

    set i 0

    # rewind it to be able to start reading from beginning...
    seek $File_p 0 start

    while {![eof $File_p]} { 

	gets $File_p Text

	set Option [lindex $Text 0]

	if {([string trim $Text] == "") || ($Option != "ATOM")} {
	    continue}

	incr i
    }

    # rewind it to be able to start reading from beginning...
    seek $File_p 0 start
    return $i
}

###########################################################################
# PROC
proc Write {Struct FileName Mask} {

    set def_title  "Default title for: Tinker"

    set File_p [open $FileName w]

    if {$File_p == ""} {
	catch {lulErrorDialog {ERROR: can't open UHBD (qcd) output file for writing!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open UHBD (qcd) output file for writing!"
	}
	return 1
    }

    puts "Writing UHBD (qcd) coordinate file '$FileName' ..."

    set TotAtoms [show numatoms $Struct]

    #   puts $File_p " $TotAtoms  $def_title"

    for {set i 1} {$i <= $TotAtoms } {incr i} {

	# apply the display mask...
	if {[show atom displaystate $i $Struct] == 0 && $Mask} continue

	scan [show atom coordinates   $i $Struct] "%f %f %f" xc yc zc
	set  charge [show atom charge $i $Struct]
	set  vdw    [show atom vdw    $i $Struct]
	set  resn   [show atom resnum $i $Struct]
	set  res    [show atom residu $i $Struct]
	set  atom   [show atom atomna $i $Struct]

	puts $File_p "ATOM $resn $res $atom $xc $yc $zc $charge $vdw"
    }

    close $File_p

    puts "Done!"

    return 0

}

}; # end of namespace Uhbd

namespace eval Qcd {

proc Read { FileName Action } {
    return [[namespace parent]::Uhbd::Read $FileName $Action]
}

proc Write { Struct FileName Mask } {
    return [[namespace parent]::Uhbd::Write $Struct $FileName $Mask]
}

}; # end of namespace Qcd

}; # end of namespace gom::FileFilters::AtomCoordinates
