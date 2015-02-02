##############################################################################
#                            Copyright (c) 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded 2004 by: Eero Häkkinen
##############################################################################

package provide gom::gui::Plugins 1.0

namespace eval gom::gui {

namespace eval Plugins {

proc UpdatePluginsMenu {} {
    global gomPluginsMenu
    global gomPluginsMenuItems

    $gomPluginsMenu delete 0 end

    set items {}

    foreach {plugin list} [array get gomPluginsMenuItems] {
	set label   [lindex $list 0]
	set command [lindex $list 1]
	if { "" == $command } {
	    lappend items [list cascade -label $label \
			       -menu    $gomPluginsMenu.$plugin]
	} else {
	    lappend items [list command -label $label \
			       -command $command]
	}
    }
    foreach item [lsort -dictionary -index 2 $items] {
	eval [list $gomPluginsMenu add] $item
    }
}

namespace export {[A-Z]*}

}

namespace import Plugins::*

}

proc lulUpdatePluginsMenu {} {
    gom::gui::UpdatePluginsMenu
}
