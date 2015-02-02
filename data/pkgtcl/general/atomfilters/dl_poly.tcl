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
    dl_poly {{1 0} DL_Poly {DL_Poly CONFIG files} {.config .CONFIG}}
}

namespace eval Dl_poly {

;
#################################################################
# PROC 
proc Read {FileName Action} {
    #
    # Read DL_Poly CONFIG file 
    #
    # open file ...
    puts "$FileName $Action"
    set File_p [open $FileName r]

    # title
    gets $File_p Text
    gomPrint "Title: $Text"
    # control
    gets $File_p Text
    set levcfg [lindex $Text 0]
    set imcon  [lindex $Text 1]

    if {$imcon > 0} {
	gets $File_p Text
	gets $File_p Text
	gets $File_p Text
    }

    set Trigger 0
    set i       0
    while {![eof $File_p]} {

	# read now the atoms ...
	gets $File_p Text
	set Text [string trim $Text]
	if {$Text == ""} {
	    continue
	}
	incr i
	set Atm($i)    [lindex $Text 0]
	set Charge($i) "0.0"
	gets $File_p Text
	set Xc($i)     [lindex $Text 0]
	set Yc($i)     [lindex $Text 1]
	set Zc($i)     [lindex $Text 2]
	if {$levcfg > 0} {
	    gets $File_p Text
	}
	if {$levcfg > 1} {
	    gets $File_p Text
	}
    }

    set natms $i

    define structure "$FileName" $natms $Action

    set StrucNum [show molstructures]

    for {set i 1} {$i <= $natms} {incr i} {
	# put the values into the structure ...
	define atom    label    $Atm($i)                $i $StrucNum
	unset Atm($i)
	define residue label    "DL_P"                  $i $StrucNum
	define segment label    "DL_P"                  $i $StrucNum
	define atom coordinates $Xc($i) $Yc($i) $Zc($i) $i $StrucNum
	unset Xc($i)
	unset Yc($i)
	unset Zc($i)
	define atom charge      $Charge($i)             $i $StrucNum
	unset Charge($i)
	define atom resnumber   "1"                     $i $StrucNum
	set Trigger 1
    }

    # close file and return ...
    close $File_p
    if {!$Trigger} {
        gomError "Could not find the atoms in the DL_Poly file"
        error "Could not find the atoms in the DL_Poly file"
    }
    puts Done!
}
}; # end of namespace Dl_poly
}; # end of namespace gom::FileFilters::AtomCoordinates
