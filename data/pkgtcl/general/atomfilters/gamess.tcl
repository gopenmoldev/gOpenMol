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
    gamess {{1 0} GAMESS {GAMESS files} {.out .OUT .dat .DAT .irc .IRC .log .LOG}}
}

namespace eval Gamess {

;
#################################################################
#
#  First version of a GAMESS reader by Leif Laaksonen/CSC
#
#  Extensive rewrite and additions 2001 by:
#
#  Dr. Germund Höjer
#  Departamento de Física y Química Teórica
#  Facultad de Química
#  UNAM
#  México D.F.
#  germund@eros.pquim.unam.mx

###### import coordinates from an external gamess file ##########
#
# The process determines file type, .out, .dat, .irc or .log, and
# chooses the appropriate filter to do the job
#
#################################################################
# PROC
# Read GAMESS coordinate file
#
# Input file name:                          FileName
# Read new structure or append to old list: Action (new / append)

proc Read {FileName Action} {

    gomPrint "Reading GAMESS coordinate file '$FileName'..."

    # determine file type, .out, .dat, .irc or .log, and choose appropriate filter.
    set ext [file extension $FileName]
#    gomPrint "$ext"

    if {[string equal -nocase $ext ".DAT"]} {
	  ReadDATCoordinates $FileName $Action
      return
    }

    if {[string equal -nocase $ext ".IRC"]} {
	  ReadIRCCoordinates $FileName $Action
      return
    }

    if {[string equal -nocase $ext ".LOG"]} {
	  ReadLOGCoordinates $FileName $Action
      return
    }

    if {[string equal -nocase $ext ".OUT"]} {
	  ReadLOGCoordinates $FileName $Action
      return
    }

    gomError "unknown file extension '$ext'"
    error "unknown file extension '$ext'"

}

##### import coordinates from an external gamess dat file #######
#
# This filter program reads a GAMESS DAT output file (.dat).
#
# It reads the final coordinate set in the file produced with the
# MOLPLT option set in the $CONTRL group.
#
#################################################################
# PROC
# Input file name:                          FileName
# Read new structure or append to old list: Action (new  / append)

proc ReadDATCoordinates {FileName Action} {

    global gomAtom_SN
    #
    # Read GAMESS DAT file with MOLPLT coordinates
    #
    # open file ...
    set File_p [open $FileName r]
    if {$File_p == ""} {
	catch {lulErrorDialog {ERROR: can't open GAMESS DAT file!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open GAMESS DAT file!"
	}
	error "ERROR: can't open GAMESS DAT file!"
	return 1
    }

    # do the job ...
    puts "Reading GAMESS DAT file with MOLPLT coordinates '$FileName'..."

    set Trigger  0

    while {![eof $File_p]} {

	gets $File_p Text

	if {[string match "*START OF -MOLPLT- INPUT FILE*" $Text]} {

	    set Trigger 1
	    gets $File_p Text

	    # some essentials ...
	    scan $Text "%*s %d %*s %d %*s %d" natoms nkinds nbonds

	    define structure "$FileName" $natoms $Action

	    set StrucNum [show molstructures]

	    # title ...
	    gets $File_p Text
	    gomPrint "Title: $Text"

	    # some dummy ...
	    for {set i 0} {$i < $nkinds} {incr i} {
		gets $File_p Text
	    }

	    # read now the atoms ...
	    for {set i 1} {$i <= $natoms} {incr i} {

		gets $File_p Text

		set ResN   1
		set Res    "GMS"
		set Seg    "GMS"
		set Atm    [string trim [lindex $Text 0]]

		# assure proper notation for chemical symbols, for example Fe and not FE or fe ...
		if {[string length $Atm] == 2} {
		    set char1 [string toupper [string index $Atm 0]]
		    set char2 [string tolower [string index $Atm 1]]
		    set Atm [concat $char1$char2]
		}

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
	    }
	    break
	}
    }

    close $File_p
    if {!$Trigger} {
        gomError "Could not find the MOLPLT information in file"
        error "Could not find the MOLPLT information in file"
    }
    # close file and return ...
    puts Done!
    return 0
}

##### import coordinates from an external gamess irc file #######
#
# This filter program reads a GAMESS IRC output file (.irc).
#
# It also produces an XMOL file of all the coordinate sets
# that is found in the file. This can then be displayed
# using the trajectory tools.
#
# G Höjer, 2001-05-21
#
#################################################################
# PROC
# Read GAMESS IRC coordinate file
#
# Input file name:                          FileName
# Read new structure or append to old list: Action (new / append)

proc ReadIRCCoordinates {FileName Action} {

    global gomAtom_SN

    # open input and output files ...

    set File_p [open $FileName r]
    if {$File_p == ""} {
	catch {lulErrorDialog {ERROR: can't open GAMESS IRC file!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open GAMESS IRC file!"
	}
	error "ERROR: can't open GAMESS IRC file!"
	return
    }

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

    # find first coordinate set and show its structure ...
    set Coord_sets  0
    set Trigger     0

    while {![eof $File_p]} {
	gets $File_p Text

	if {[string match "***** BEGIN IRC INFORMATION PACKET *****" $Text]} {

	    # title ...
	    gets $File_p Text
	    set Title $Text
	    gomPrint "Title: $Title"

	    # line with text "CARTESIAN COORDINATES (BOHR)"
	    gets $File_p Text
	    set i 0

	    while {![eof $File_p]} {
		gets $File_p Text
		if {[lindex $Text 0] == "MASS-WEIGHTED"} {
		    break
		}
		gomPrint $Text
		incr i

		set AtmNum  [expr int([lindex $Text 1])]
		set Atm($i) [lindex $gomAtom_SN($AtmNum) 0]
		set Xc($i)  [expr [lindex $Text 2] * 0.52917715]
		set Yc($i)  [expr [lindex $Text 3] * 0.52917715]
		set Zc($i)  [expr [lindex $Text 4] * 0.52917715]
		set Ch($i)  0.0
	    }

	    set Natoms $i

	    define structure "$FileName" $Natoms $Action

	    set StrucNum [show molstructures]

	    # add set to XMOLE file and
	    # put the values into the structure ...
	    puts $File_px $Natoms
	    puts $File_px $Title
	    set ResN    1
	    set Res     "GMS"
	    set Seg     "GMS"
	    for {set i 1} {$i <= $Natoms} {incr i} {
		puts $File_px "$Atm($i) $Xc($i) $Yc($i) $Zc($i)"
		define atom    label    $Atm($i)                $i $StrucNum
		define residue label    $Res                    $i $StrucNum
		define segment label    $Seg                    $i $StrucNum
		define atom coordinates $Xc($i) $Yc($i) $Zc($i) $i $StrucNum
		define atom charge      $Ch($i)                 $i $StrucNum
		define atom resnumber   $ResN                   $i $StrucNum
	    }

	    # skip mass-weighted gradient
	    for {set i 1} {$i <= $Natoms} {incr i} {
		gets $File_p Text
	    }
	    incr Coord_sets
	    set Trigger 1
	}
	if {$Trigger == 1} {
	    break
	}
    }

    # find next coordinate set(s) ...
    while {![eof $File_p]} {
	gets $File_p Text

	if {[string match "***** BEGIN IRC INFORMATION PACKET *****" $Text]} {

	    # title ...
	    gets $File_p Text
	    set Title $Text
	    gomPrint "Title: $Title"

	    # line with text "CARTESIAN COORDINATES (BOHR)"
	    gets $File_p Text

	    for {set i 1} {$i <= $Natoms} {incr i} {
		gets $File_p Text
		gomPrint $Text
		set AtmNum  [expr int([lindex $Text 1])]
		set Atm($i) [lindex $gomAtom_SN($AtmNum) 0]
		set Xc($i)  [expr [lindex $Text 2] * 0.52917715]
		set Yc($i)  [expr [lindex $Text 3] * 0.52917715]
		set Zc($i)  [expr [lindex $Text 4] * 0.52917715]
	    }

	    # add set to XMOLE file ...
	    puts $File_px $Natoms
	    puts $File_px $Title
	    for {set i 1} {$i <= $Natoms} {incr i} {
		puts $File_px "$Atm($i) $Xc($i) $Yc($i) $Zc($i)"
	    }

	    # skip mass-weighted gradient
	    gets $File_p Text
	    for {set i 1} {$i <= $Natoms} {incr i} {
		gets $File_p Text
	    }
	    incr Coord_sets
	}
    }

    gomPrint "Found $Coord_sets coordinate sets and $Natoms atoms in each set"

    # close file and return ...
    close $File_p
    close $File_px
    puts Done!
}

##### import coordinates from an external gamess log file #######
#
# This filter program reads a GAMESS LOG output file (.log).
#
# It also produces an XMOL file of all the coordinate sets
# that is found in the file. This can then be displayed
# using the trajectory tools.
#
# G Höjer, 2001-05-21
#
#################################################################
# PROC
# Read GAMESS LOG file
#
# Input file name:                          FileName
# Read new structure or append to old list: Action (new / append)

proc ReadLOGCoordinates {FileName Action} {

    global gomAtom_SN

    # open input and output files ...

    set File_p [open $FileName r]
    if {$File_p == ""} {
	catch {lulErrorDialog {ERROR: can't open GAMESS LOG file!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open GAMESS LOG file!"
	}
	error "ERROR: can't open GAMESS LOG file!"
	return
    }

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

    # get title ...
    set Trigger 0

    while {![eof $File_p]} {
	gets $File_p Text
	set Text [string trim $Text]
	if {[string match "RUN TITLE" $Text]} {
	    gets $File_p Text
	    gets $File_p Text
	    set Title $Text
	    gomPrint "Title: $Title"
	    set Trigger 1
	}
	if {$Trigger == 1} {
	    break
	}
    }
    puts " Title: $Title"


    # find original coordinate set ...
    set Coord_sets  0
    set Trigger     0

    while {![eof $File_p]} {
	gets $File_p Text
	set Text [string trim $Text]
	if {[string match "ATOM      ATOMIC                      COORDINATES (BOHR)" $Text]} {
	    gets $File_p Text
	    set i 0

	    while {![eof $File_p]} {
		gets $File_p Text
		gomPrint $Text
		set Text [string trim $Text]
		if {$Text == ""} {
		    break
		}
		incr i

		set AtmNum  [expr int([lindex $Text 1])]
		set Atm($i) [lindex $gomAtom_SN($AtmNum) 0]
		set Xc($i)  [expr [lindex $Text 2] * .52917715]
		set Yc($i)  [expr [lindex $Text 3] * .52917715]
		set Zc($i)  [expr [lindex $Text 4] * .52917715]
		set Ch($i)  0.0
	    }

	    set Natoms $i

	    define structure "$FileName" $Natoms $Action

	    set StrucNum [show molstructures]

	    # add set to XMOLE file and
	    # put the values into the structure ...
	    puts $File_px $Natoms
	    puts $File_px $Title
	    set ResN    1
	    set Res     "GMS"
	    set Seg     "GMS"
	    for {set i 1} {$i <= $Natoms} {incr i} {
		puts $File_px "$Atm($i) $Xc($i) $Yc($i) $Zc($i)"
		define atom    label    $Atm($i)                $i $StrucNum
		define residue label    $Res                    $i $StrucNum
		define segment label    $Seg                    $i $StrucNum
		define atom coordinates $Xc($i) $Yc($i) $Zc($i) $i $StrucNum
		define atom charge      $Ch($i)                 $i $StrucNum
		define atom resnumber   $ResN                   $i $StrucNum
	    }
	    incr Coord_sets
	    set Trigger 1
	}
	if {$Trigger == 1} {
	    break
	}
    }

    # check for more coordinate sets
    set More_sets   0
    set Trigger     0
    set CtrlOpt "?CONTRL OPTIONS"



    while {![eof $File_p]} {
	gets $File_p Text
	set Text1 [string trim $Text]
	if {[string match $CtrlOpt $Text1]} {
	    gets $File_p Text
	    gets $File_p Text
	    set RunOpt [lindex $Text 1]
	    if {($RunOpt == "RUNTYP=OPTIMIZE") || ($RunOpt == "RUNTYP=SADPOINT")} {
		set More_sets 1
		gomPrint "Option: $RunOpt"
	    }
	    set Trigger 1
	}
	if {$Trigger == 1} {
	    break
	}
    }

    if {$More_sets == 1} {

	# find next coordinate set

	while {![eof $File_p]} {
	    gets $File_p Text
	    if {[lindex $Text 0] == "1NSERCH="} {
		gomPrint $Text
		set Trigger 0
		while {![eof $File_p]} {
		    gets $File_p Text
		    if {$Text == " COORDINATES OF ALL ATOMS ARE (ANGS)"} {
			gets $File_p Text
			gets $File_p Text
			for {set i 1} {$i <= $Natoms} {incr i} {
			    gets $File_p Text
			    gomPrint $Text
			    set AtmNum  [expr int([lindex $Text 1])]
			    set Atm($i) [lindex $gomAtom_SN($AtmNum) 0]
			    set Xc($i)  [lindex $Text 2]
			    set Yc($i)  [lindex $Text 3]
			    set Zc($i)  [lindex $Text 4]
			}
			set Trigger 1
		    }
		    if {$Trigger == 1} {
			break
		    }
		}

		set Trigger 0

		while {![eof $File_p]} {
		    gets $File_p Text
		    if {[string first " FINAL ENERGY IS" $Text 0]
                    || [string first " FINAL MCSCF ENERGY IS" $Text 0] != -1} {
			gomPrint "Title: $Text"

			# add set to XMOLE file ...
			puts $File_px $Natoms
			puts $File_px $Text
			for {set i 1} {$i <= $Natoms} {incr i} {
			    puts $File_px "$Atm($i) $Xc($i) $Yc($i) $Zc($i)"
			}
			incr Coord_sets
			set Trigger 1
		    }
		    if {$Trigger == 1} {
			break
		    }
		}
	    }
	}
    }

    gomPrint "Found $Coord_sets coordinate sets and $Natoms atoms in each set"

    # close file and return ...
    close $File_p
    close $File_px
    puts Done!
}


####################### import external coordinates #############
# PROC
# Input file name:                          FileName
# Read new structure or append to old list: Action (= 0 new  != 0 append)

proc lulReadGAMESSCoordinatesLUL {FileName Action} {

    global gomAtom_SN
    #
    # Read GAMESS MOLPLT coordinate file
    #
    # open file ...
    set File_p [open $FileName r]
    if {$File_p == ""} {
	catch {lulErrorDialog {ERROR: can't open GAMESS MOLPLT file!}} errRet
	if {$errRet != ""} {
	    puts "ERROR: can't open GAMESS MOLPLT file!"
	}
	error "ERROR: can't open GAMESS MOLPLT file!"
	return
    }

    # do the job ...
    puts "Reading GAMESS MOLPLT coordinate file '$FileName'..."

    set Trigger  0
    set OnlyOnce 0
    set MolSets  1

    while {![eof $File_p]} {

	gets $File_p Text

	if {[string match "*START OF -MOLPLT- INPUT FILE*" $Text]} {

	    set Trigger 1
	    gets $File_p Text

	    # some essentials ...
	    scan $Text "%*s %d %*s %d %*s %d" natoms nkinds nbonds

	    if {!$Action} {
		define structure "$FileName" $natoms new
	    } else {
		define structure "$FileName" $natoms append
	    }

	    set StrucNum [show molstructures]

	    # title ...
	    gets $File_p Text
	    gomPrint "Title: $Text"

	    # some dummy ...
	    for {set i 0} {$i < $nkinds} {incr i} {
		gets $File_p Text
	    }

	    # read now the atoms ...
	    for {set i 1} {$i <= $natoms} {incr i} {

		gets $File_p Text

		set ResN   1
		set Res    "GMS"
		set Atm    [lindex $Text 0]
		set Seg    "GMS"
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
	    }
	    break
	}
    }

    close $File_p
    if {!$Trigger} {
        gomError "Could not find the MOLPLT information in file"
        error "Could not find the MOLPLT information in file"
    }
    # close file and return ...
    puts Done!
}
}; # end of namespace Gamess
}; # end of namespace gom::FileFilters::AtomCoordinates
