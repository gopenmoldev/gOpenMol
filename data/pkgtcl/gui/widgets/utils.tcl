##############################################################################
#                           Copyright (c) 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

namespace eval ::gom::gui::Widgets {

############################################################################
# SYNOPSIS
#   lb <cmd> ?<subcmd> ...? <path> ?<options>? -label <label> ?<options>?
# FUNCTION
#   Evaluates following two commands
#       label <path>label -text <label> -anchor w
#       <cmd> ?<subcmd> ...? <path> ?<options>? ?<options>?
# RETURN VALUE
#   [list <path>label <path>]
############################################################################
proc lb { args } {
    set i     [lsearch  $args -label]
    set label [lindex   $args [expr $i + 1]]
    set args  [lreplace $args $i [expr $i + 1]]
    set path  [lsearch -glob -inline $args .*]
    label ${path}label -text $label -anchor w
    uplevel $args
    return [list ${path}label $path]
}

############################################################################
# SYNOPSIS
#   ns <cmd> ?<subcmd> ...? ?<path>? ?<options>?
# FUNCTION
#   Replace values of -*command* and -*variable* options with values
#   referring to a namespace of a caller.
#   Return a value returned by <cmd>.
############################################################################
# ns <cmd> ?<subcmd> ...? ?<path>? ?<options>?
proc ns { args } {
    foreach i [lsearch -all -glob $args -*command*] {
	set j [expr $i + 1]
	lset args $j [uplevel namespace code [lrange $args $j $j]]
    }
    foreach i [lsearch -all -glob $args -*variable*] {
	set j [expr $i + 1]
	regexp {^(.*?)(|\([^()]*\))$} [lindex $args $j] dummy var1 var2
	# Allow global variables.
	set var [uplevel namespace which -variable $var1]
	if { "" == $var } {
	    set var [uplevel namespace current]::$var1
	}
	lset args $j $var$var2
    }
    return [uplevel $args]
}

proc regrid { widget path } {
    set gridopts [grid info $path]
    set opts     [list]
    foreach spec [$path configure] {
	if { [llength $spec] == 5 } {
	    lappend opts [lindex $spec 0] [lindex $spec 4]
	}
    }
    destroy $path
    eval [list $widget $path] $opts
    eval [list grid    $path] $gridopts
}

namespace export lb ns regrid

}
