##############################################################################
#                           Copyright (c) 1994 - 2002 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Enhancements by: Eero Häkkinen
##############################################################################
namespace eval lulFile {

variable exportCoordStructure

#<<< Options
set appendStructure 0

set selectedCoordType(import) "charmm"
set selectedCoordType(export) "charmm"

set selectCoordFileFormatByExtension(import) 1
set selectCoordFileFormatByExtension(export) 0
#>>>

#<<< Import procedures

###################################################################
# PROC
proc ImportCoordFile {} {
    global gomHelpDir
    global gomHelpFile
    global gomControlFont
    global lulCalculateAtomConnectivity

    set ns [namespace current]
    
    set w .gomimportcoordinates
    catch {destroy $w}
    toplevel $w
    wm title $w "Import coordinates"
    wm iconname $w "Import coordinates"

    # Buttons
    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
	-command "destroy $w"
    button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
	-command "${ns}::DoImportCoordinates $w"
    button $w.buttons.help    -text Help    -font "$gomControlFont" \
	-command \
	"htmlShowHelp $gomHelpFile(importcoordinates)"
    pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1
    
    frame $w.file  
    pack  $w.file -side left  -fill both


    # File formats
    CreateCoordFileFormatFrame        import $w.formats
    SelectCoordFileFormatByExtension  import $w.formats
    
    # File input
    labelframe  $w.file.input     -text  "Input file name" \
	-borderwidth 2 -relief ridge -padx 2 -pady 2
    entry       $w.file.input.filename -width 40
    button      $w.file.input.browse   -text "Browse..."    \
	-command "${ns}::BrowseImportCoordFile $w.file.input.filename" \
	-width 10

    button      $w.file.peek      -text "Peek file" \
	-width 10 -command "${ns}::PeekCoordFile $w.file.input.filename"

    checkbutton $w.file.append    -text "Append structure" \
	-variable ${ns}::appendStructure
    checkbutton $w.file.calcconn  -text "Calculate connectivity" \
	-variable lulCalculateAtomConnectivity \
	-onvalue "yes" -offvalue "no" 
    checkbutton $w.file.typebyext -text "Select file format by extension" \
	-variable ${ns}::selectCoordFileFormatByExtension(import) \
	-onvalue 1 -offvalue 0 \
	-command  "${ns}::SelectCoordFileFormatByExtension import $w.formats"

    # bind the <return>
    bind $w.file.input.filename <Return> "${ns}::DoImportCoordinates $w"
    
    pack   $w.file.input          -side top
    pack   $w.file.input.filename -side left
    pack   $w.file.input.browse   -padx 2 -side left
    
    pack   $w.file.peek           -padx 2 -pady 4 -side right -anchor n
    pack   $w.file.append         -padx 2 -pady 4 -side top   -anchor w
    pack   $w.file.calcconn       -padx 2 -pady 4 -side top   -anchor w
    pack   $w.file.typebyext      -padx 2 -pady 4 -side top   -anchor w
    
    if {$lulCalculateAtomConnectivity == "no"}  {
	$w.file.calcconn deselect
    } else {
	$w.file.calcconn select
    }
}

###################################################################
# PROC
proc BrowseImportCoordFile { w } {
    set fileTypes [GetCoordFileFormatList import]
    
    set file [tk_getOpenFile -filetypes $fileTypes \
		  -parent $w -title "Read Coordinate File"]

    if {"" != $file} {
	$w delete 0 end
	$w insert 0 $file
    }
}

##################################################################
# PROC
proc PeekCoordFile { w } {
    set fileName [$w get]

    if {"" != $fileName} {
	::lulDisplayTextFile $fileName
    }
}

###################################################################
# PROC
proc DoImportCoordinates { w } {
    global gomURLFileName
    variable appendStructure

    set fileName [$w.file.input.filename get]

    if {"" == $fileName} return

    set coordType [GetCoordFileFormat import $fileName]
    if {"" == $coordType} {
	gomError "No import filter found for '$fileName'"
	return
    }

    if {$appendStructure} {
	set append "append"
    } else {
	set append "new"
    }

    set errorValue [import coord $coordType $fileName $append]

    if {0 == $errorValue} {
	# was it an URL?
	if {[string match -nocase "http://*" $fileName]} {
	    set fileName $gomURLFileName
	    if {$gomURLFileName == ""} {
		lulInvalidateDisplay
		return
	    }
	}
	# change to current directory
	lulChangeDirectory "$fileName"

	lulInvalidateDisplay
    }
}

    #>>>

#<<< Export procedures

##################### export coordinate file ############################
# PROC
proc ExportCoordFile {} {
    global gomHelpDir
    global gomHelpFile
    global gomControlFont

    variable exportCoordStructure

    set ns [namespace current]

    # return if no molecular systems defined
    if {[show molstructures] < 1} {
	lulErrorDialog {ERROR: no structure available. Read a structure first!}
	return
    }

    set w .gomexportcoordinates
    catch {destroy $w}
    toplevel $w 
    wm title $w "Export coordinates"
    wm iconname $w "Export coordinates"

    # Buttons
    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
    button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "${ns}::DoExportCoordinates $w"
    button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
	"htmlShowHelp $gomHelpFile(exportcoordinates)"
    pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

    frame $w.file
    pack  $w.file -side left -fill both

    # File formats
    CreateCoordFileFormatFrame        export $w.formats
    SelectCoordFileFormatByExtension  export $w.formats

    # File export
    labelframe $w.file.input -text "Export to file" \
	-borderwidth 2 -relief ridge -padx 2 -pady 2
    entry      $w.file.input.filename -width 40
    button     $w.file.input.browse   -text "Browse..." \
	-command "${ns}::BrowseExportCoordFile $w.file.input.filename"

    # bind the <return>
    bind $w.file.input.filename <Return> "${ns}::DoExportCoordinates $w"

    set NumStruct [show molstructures]

    labelframe $w.file.structures -text "Choose structure" \
	-borderwidth 2 -relief ridge -padx 2 -pady 2
    text       $w.file.structures.text \
	-yscrollcommand "$w.file.structures.scrollbar set" \
	-width 20 -height 5 -bg gray
    scrollbar  $w.file.structures.scrollbar \
	-command "$w.file.structures.text yview"

    for {set i 1} {$i <= $NumStruct} {incr i} {
	frame       $w.file.structures.text.frame$i -borderwidth 2 -relief ridge
	radiobutton $w.file.structures.text.frame$i.check -text "Structure nr: $i" \
	    -variable ${ns}::exportCoordStructure -value $i \
	    -cursor hand2 -bg gray
	pack        $w.file.structures.text.frame$i.check -side left
	$w.file.structures.text window create $i.0 \
	    -window $w.file.structures.text.frame$i
    }
    $w.file.structures.text.frame1.check select

    checkbutton $w.file.typebyext -text "Select file format by extension" \
	-variable ${ns}::selectCoordFileFormatByExtension(export) \
	-onvalue 1 -offvalue 0 \
	-command "${ns}::SelectCoordFileFormatByExtension export $w.formats"

    pack $w.file.input           -side top
    pack $w.file.input.filename  -side left
    pack $w.file.input.browse    -side left -padx 2
    pack $w.file.structures      -side top -anchor w
    pack $w.file.structures.text -side left
    pack $w.file.typebyext       -side top -anchor w
}

###################################################################
# PROC
proc BrowseExportCoordFile { w } {
    variable selectedCoordType
    variable selectCoordFileFormatByExtension

    set fileTypes [GetCoordFileFormatList export]
    if {$selectCoordFileFormatByExtension(export)} {
	set defext ""
    } else {
	set patternlist [lindex \
	    $gom::FileFilters::AtomCoordinates::coordTypes($selectedCoordType(export)) \
	    $gom::FileFilters::AtomCoordinates::coordTypeListIndeces(patterns)]
	set defext      [lindex $patternlist 0]
    }	

    set file [tk_getSaveFile -filetypes $fileTypes \
		  -defaultextension $defext \
		  -parent $w -title "Export Coordinate File"]

    if {"" != $file} {
	$w delete 0 end
	$w insert 0 $file
    }
}

#################################################################
# PROC
proc DoExportCoordinates { w } {
    variable exportCoordStructure

    set fileName [$w.file.input.filename get]

    if {"" == $fileName} return

    set coordType [GetCoordFileFormat export $fileName]
    if {"" == $coordType} {
	gomError "No export filter found for '$fileName'"
	return
    }
 
    set errorCode [export coord $exportCoordStructure $coordType "$fileName"]

    if {0 == $errorCode} {
	# change to current directory
	lulChangeDirectory "$fileName"
	lulMessageDialog "Coordinate file exported!"
    }
}

#>>>

#<<< Common procedures
###################################################################
# PROC
proc CreateCoordFileFormatFrame { which w } {
    set ns [namespace current]

    labelframe $w -text "Coordinate formats:" \
	-borderwidth 2 -relief ridge -padx 2 -pady 2

    set types ""

    foreach type [array names gom::FileFilters::AtomCoordinates::coordTypes] {
	# Check if importable or exportable
	set context [lindex \
			 $gom::FileFilters::AtomCoordinates::coordTypes($type) \
			 $gom::FileFilters::AtomCoordinates::coordTypeListIndeces(context)]
	if {![lindex $context \
		  $gom::FileFilters::AtomCoordinates::coordTypeListIndeces(context_$which)]} {
	    continue
	}

	lappend types [list $type [lindex \
	    $gom::FileFilters::AtomCoordinates::coordTypes($type) \
	    $gom::FileFilters::AtomCoordinates::coordTypeListIndeces(name)]]
    }

    set typeCount 0
    foreach type [lsort -dictionary -index 1 $types] {

	set column col[expr $typeCount/10+1]
	if {$typeCount % 10 == 0} {
	    frame $w.$column
	    pack  $w.$column -side left -anchor n
	}

	set name [lindex $type 1]
	set type [lindex $type 0]

	radiobutton $w.$column.$type -text $name \
	    -variable ${ns}::selectedCoordType($which) -value $type
	pack        $w.$column.$type -side top -anchor w

	incr typeCount
    }
}

###################################################################
# PROC
proc SelectCoordFileFormatByExtension { which w } {
    variable selectCoordFileFormatByExtension

    if {$selectCoordFileFormatByExtension($which)} {
	pack forget $w
    } else {
	pack $w -side right -padx 2
    }
}

###################################################################
# PROC
proc GetCoordFileFormat { which fileName } {
    variable selectedCoordType
    variable selectCoordFileFormatByExtension

    if {$selectCoordFileFormatByExtension($which)} {
	return [gom::GetAtomCoordinateFileFormatByExtension $which $fileName]
    } else {
	return $selectedCoordType($which)
    }
}

###################################################################
# PROC
proc GetCoordFileFormatList { which } {
    variable selectedCoordType
    variable selectCoordFileFormatByExtension

    if {$selectCoordFileFormatByExtension($which)} {
	set fileTypes ""
	set patterns  ""
	set omitted   ""

	foreach type [lsort [array names gom::FileFilters::AtomCoordinates::coordTypes]] {
	    # Check if importable or exportable
	    set context [lindex \
			     $gom::FileFilters::AtomCoordinates::coordTypes($type) \
			     $gom::FileFilters::AtomCoordinates::coordTypeListIndeces(context)]
	    if {![lindex $context \
		      $gom::FileFilters::AtomCoordinates::coordTypeListIndeces(context_$which)]} {
		continue
	    }

	    set description [lindex \
		$gom::FileFilters::AtomCoordinates::coordTypes($type) \
		$gom::FileFilters::AtomCoordinates::coordTypeListIndeces(patternname)]
	    set patternlist [lindex \
		$gom::FileFilters::AtomCoordinates::coordTypes($type) \
		$gom::FileFilters::AtomCoordinates::coordTypeListIndeces(patterns)]

	    foreach pattern $patternlist {
		# Accept only suffixes
		if {[string first * $pattern] >= 0} continue

		set index [lsearch -exact $patterns $pattern]
		if {$index >= 0} {
		    # suffix is amphibious
		    set patterns  [lreplace $patterns  $index $index]
		    set fileTypes [lreplace $fileTypes $index $index]
		    lappend omitted $pattern
		} elseif {[lsearch -exact $omitted $pattern] < 0} {
		    lappend fileTypes [list $description $pattern TEXT]
		    lappend patterns  $pattern
		}
	    }
	}

	set fileTypes [linsert $fileTypes 0 \
			   [list \
				{Unambiguous molecule files} \
				$patterns \
				TEXT]]
    } else {
	set description [lindex \
	    $gom::FileFilters::AtomCoordinates::coordTypes($selectedCoordType($which)) \
	    $gom::FileFilters::AtomCoordinates::coordTypeListIndeces(patternname)]
	set patternlist [lindex \
	    $gom::FileFilters::AtomCoordinates::coordTypes($selectedCoordType($which)) \
	    $gom::FileFilters::AtomCoordinates::coordTypeListIndeces(patterns)]

	if {"" == $patternlist} {
	    set fileTypes ""
	} else {
	    set fileTypes [list \
			       [list $description $patternlist TEXT] \
			       {{All files} {*}}]
	}
    }

    return $fileTypes
}
#>>>

}
