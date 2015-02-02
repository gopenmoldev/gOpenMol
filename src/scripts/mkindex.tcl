# /bin/sh
##############################################################################
#                           Copyright (c) 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################
# The next line is part of a Tcl comment but is a Posix shell command. \
exec tclsh "$0" "$@"
if { [info exists env(GOM_DATA)] } {
    set datadir $env(GOM_DATA)
} elseif { [info exists env(GOM_ROOT)] } {
    set datadir [file join $env(GOM_ROOT) data]
} else {
    set datadir [file join [file dirname $argv0] .. .. data]
}
foreach dir [glob -types d -directory $datadir/pkgtcl * */* */*/* */*/*/*] {
    puts stderr $dir
    file delete [file join $dir pkgIndex.tcl]
    if { [llength [glob -nocomplain -directory $dir *.tcl]] } {
	pkg_mkIndex -verbose $dir *.tcl
#	auto_mkindex $dir *.tcl
    }
}
