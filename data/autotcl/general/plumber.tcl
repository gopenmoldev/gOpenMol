##############################################################################
#                           Copyright (c) 2001 - 2002 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################        
namespace eval lulPlumber {

variable secondaryStructureRecord
variable secondaryStructureRecords

array unset secondaryStructureRecord
array unset secondaryStructureRecords

#################################################################
# PROC
proc SecondaryStructureSaver {Struc Position Input} {

    variable secondaryStructureRecord
    variable secondaryStructureRecords

    set Input [string replace "$Input" \
		   [expr [string length "$Input"] - 1] end ""]
    set secondaryStructureRecord($Struc,$Position) "$Input"
    set secondaryStructureRecords($Struc)          $Position

    gomPrint "Record added ($Struc,$Position): $secondaryStructureRecord($Struc,$Position)"
}

}
