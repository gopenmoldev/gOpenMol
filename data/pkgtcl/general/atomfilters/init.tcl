##############################################################################
#                        Copyright (c) 2002 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded 2002 - 2004 by: Eero Häkkinen
##############################################################################

package provide gom::FileFilters::AtomCoordinates 1.0
package require gom::PackageUtilities

namespace eval gom {
namespace eval FileFilters::AtomCoordinates {

;
# File format declarations
set coordTypeListIndeces(context)        0
set coordTypeListIndeces(context_import) 0
set coordTypeListIndeces(context_export) 1
set coordTypeListIndeces(name)           1
set coordTypeListIndeces(patternname)    2
set coordTypeListIndeces(patterns)       3

# Format is following:
# {
#   {
#     can_be_imported
#     can_be_exported
#   }
#   filter_name
#   file_type_name
#   file_pattern_list
# }

# It's totally unclear what this Amber file format really is
# so it has been disabled but can be included again by removing
# the comment character. (LUL 16.09.2002)
#set coordTypes(amber)     {{1 0} Amber     {Amber files}     {*}}
#
# 01.04.2003
# Don't swim against the stream!
# The .xyz extension points to XMOL and the .gxyz is the specific gOpenMol xyz file format
#
# Filters for the following formats are implemented in C source.
array set coordTypes {
    charmm    {{1 1} CHARMM    {CHARMM files}    {.crd .CRD}}
    frame     {{1 0} Frame     {}                {}}
    gaussian  {{1 0} GAUSSIAN  {GAUSSIAN files}  {.fch .FCH .fchk .FCHK}}
    hyperchem {{1 0} HyperChem {HyperChem files} {.hin .HIN}}
    insight   {{1 0} Insight   {Insight files}   {.car .CAR}}
    mol2      {{1 0} MOL2      {MOL2 files}      {.mol2 .MOL2}}
    mopac     {{1 0} MOPAC     {MOPAC files}     {.gpt .GPT}}
    mumod     {{1 0} MUMOD     {MUMOD files}     {*}}
    openmol   {{1 1} OpenMol   {OpenMol files}   {cli*.data CLI*.DATA}}
    pdb       {{1 1} PDB       {PDB files}       {.pdb .ent .PDB .ENT}}
    xmol      {{1 0} Xmol      {XMOL files}      {.xyz .XYZ .xmol .XMOL}}
    xyz       {{0 1} XYZ       {XMOL files}      {.xyz .XYZ .xmol .XMOL}}
    gxyz      {{1 0} GXYZ      {GXYZ files}      {.gxyz .GXYZ}}
    yasp      {{1 0} YASP      {YASP files}      {*}}
}

###################################################################
# PROC
proc GetAtomCoordinateFileFormatByExtension { which fileName } {
    variable coordTypeListIndeces
    variable coordTypes

    set suffix [file extension $fileName]
    if { "" == $suffix } {
	return ""
    }

    foreach {type spec} [array get coordTypes] {
	# Check if importable or exportable
	if { ! [lindex $spec \
		    $coordTypeListIndeces(context) \
		    $coordTypeListIndeces(context_$which)] } continue

	foreach pattern [lindex $spec $coordTypeListIndeces(patterns)] {
	    # Accept only suffixes
	    if { [string first * $pattern] >= 0 } continue

	    if { [string equal $suffix $pattern] } {
		return $type
	    }
	}
    }

    return ""
}

##################################################################
# PROC
proc PreReadAtomCoordinates { type fileName append } {
    set ns [namespace current]

    set function ${ns}::[string totitle $type]::PreRead
    set function [info procs $function]
    if { "" != $function } {
	return [$function $fileName $append]
    }
    return 0
}

##################################################################
# PROC
proc ReadAtomCoordinates { type fileName append } {
    set ns [namespace current]

    set function ${ns}::[string totitle $type]::Read
    set function [info procs $function]
    if { "" != $function } {
       return [$function $fileName $append]
    }
    error "No import filter found for coordinate type '$type'"
}

##################################################################
# PROC
proc PostReadAtomCoordinates { type fileName append errorExit } {
    set ns [namespace current]

    set function ${ns}::[string totitle $type]::PostRead
    set function [info procs $function]
    if { "" != $function } {
	return [$function $fileName $append $errorExit]
    }
    return 0
}

##################################################################
# PROC
proc WriteAtomCoordinates { struct type fileName visibleOnly } {
    set ns [namespace current]

    set function ${ns}::[string totitle $type]::Write
    set function [info procs $function]
    if { "" != $function } {
	return [$function $struct $fileName $visibleOnly]
    }
    error "No export filter found for coordinate type '$type'"
}

namespace export {[A-Z]*}

}; # end of namespace FileFilters::AtomCoordinates


namespace import FileFilters::AtomCoordinates::*

}; # end of namespace gom

gom::SourcePackageFiles general/atomfilters
