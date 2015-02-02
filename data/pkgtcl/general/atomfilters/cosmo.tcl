##############################################################################
#                           Copyright (c) 1994 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#
# 1) Place the script in the data/autotcl directory.
# 2) When you start gOpenMol you will find a new
#    input coordinate file. You will see it after
#    File => Inport => Coors... Remember to uncheck
#    the "Select file format by extension"!
# 3) Doing the previous gives you a list of supported
#    coordinate input filters. Click the "Cosmo".
# 4) Click browse and file your file.
# 5) Read in the coordinates.
#     * If the file contains the field "$segment_information"
#       it reads in the surface otherwise it only reads
#       the coordinates in the "$coord_rad" record.
#     * By default it reads in the charge and plots spheres with 
#       a radius of 0.1 A and colour codes the sphere from blue to red.
#
# 6) If you want to have spheres of a different radius than 0.1
#    you have to give the command "set lulCosmoRadius your.value"
#    BEFORE you read in the file! 
# 7) If you want to display the potential you have to give the command
#    "set lulCosmoType potential" (the default type is "charge") BEFORE
#    you read in the file.
#
#
#
# If you want to remove the spheres the command is "plot -sphere".
# The procedure adds always spheres to the display... until you
# remove them all.
#
##############################################################################

namespace eval gom::FileFilters::AtomCoordinates {

;
# see pkgtcl/atomfiles/main.tcl
array set coordTypes {
    cosmo {{1 0} Cosmo {Cosmo files} {*}}
}

##########################################################################
# NAMESPACE gom::FileFilters::AtomCoordinates::Cosmo
#
namespace eval Cosmo {

;
variable lulNumAtoms
variable lulCosmoLogFileName
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
	catch {lulErrorDialog {ERROR: can't open Cosmo coordinate file for reading!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open Cosmo file for reading!"
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
    if {[string match "\$coord_rad*" $Text]} {
      getCoordinates   $f $FileName $Action
    } elseif {[string match "\$segment_information*" $Text]} {
# check for graphics
      if {[show graphics]} {
	    if {![lulYesNoMessageDialog "do you want to read the surface?"]} {
         getCoordsSurface $f $FileName $Action
        }
	  } else {
# no graphics so read everything
         getCoordsSurface $f $FileName $Action
	  } 
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
#
# label
#
	gets $file_p Text

    while {![eof $file_p]} { 
	gets $file_p Text
	set Text [string trim $Text]
	if {($Text == "") || ([string match "\$*" $Text])} break
	incr i
	set xc($i)   [expr $au2a * [lindex $Text 1]]
	set yc($i)   [expr $au2a * [lindex $Text 2]]
	set zc($i)   [expr $au2a * [lindex $Text 3]]
	set atom($i) [string toupper [lindex $Text 4] 0 0]
    }

    set NumAtoms $i

    define structure "$FileName" $NumAtoms $Action

    set StrucNum [show molstructures]

    for {set i 1} {$i <= $NumAtoms} {incr i} {
	# put the values into the structure ...
	define atom    label    $atom($i)                $i $StrucNum
	define residue label    cosm                     $i $StrucNum
	define segment label    cosm                     $i $StrucNum
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
proc getCoordsSurface {file_p FileName Action} {

    variable xc
    variable yc
    variable zc
    variable atom
    variable setnr

	global lulCosmoRadius
    global lulCosmoType

    set au2a  0.52917715
    set FirstTime 0

# set defaults
    if {![info exists lulCosmoRadius]} {
	  set lulCosmoRadius 0.1
	}
    if {![info exists lulCosmoType]} {
	  set lulCosmoType "charge"
	}

	set i 0
    while {![eof $file_p]} { 

	    gets $file_p Text
	    set TextLabel [string trim $Text]
        if {[string match "#*" $TextLabel]} {
		   continue
		}
	    if {($TextLabel == "") || ([string match "\$*" $TextLabel])} {
           break
		}

        incr i
        set atom($i)     [lindex $Text 1]
		set xc($i)       [expr $au2a * [lindex $Text 2]]
		set yc($i)       [expr $au2a * [lindex $Text 3]]
		set zc($i)       [expr $au2a * [lindex $Text 4]]

        if {$lulCosmoType == "charge"} {
           set charge($i)   [lindex $Text 5]
		} elseif {$lulCosmoType == "potential"} {
           set charge($i)   [lindex $Text 8]
		} else {
		  gomError "unknow option: $lulCosmoType"
		  return
		}

	}

    set NumAtoms $i

    set MinVal  1.0E+30
	set MaxVal -1.0E+30

    for {set i 1} {$i <= $NumAtoms} {incr i} {
	    if {$charge($i) > $MaxVal} {
		    set MaxVal $charge($i)
		}
		if {$charge($i) < $MinVal} {
		    set MinVal $charge($i)
		}
    }

    for {set i 1} {$i <= $NumAtoms} {incr i} {

	   set value [expr ($charge($i) - $MinVal) / ($MaxVal - $MinVal)]

       plot sphere $xc($i) $yc($i) $zc($i) $lulCosmoRadius "[show rainbow $value]" append

       unset atom($i)
       unset xc($i)
       unset yc($i)
	   unset zc($i)
	   unset charge($i)
    }

}

}; # end of namespace Cosmo
}; # end of namespace gom::FileFilters::AtomCoordinates
