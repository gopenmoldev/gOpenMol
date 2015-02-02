package provide gom::general 1.0
package require gom::FileFilters::AtomCoordinates

if { "" == [info commands lassign] } {
    # lassing is part of Tcl 8.5
    proc lassign { list args } {
	foreach val $list var $args {
	    upvar $var v
	    set v $val
	}
	if { [llength $list] <= [llength $args] } {
	    return ""
	} else {
	    return [lrange $list [llength $args] end]
	}
    }
}
