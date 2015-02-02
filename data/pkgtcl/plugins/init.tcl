##############################################################################
#                            Copyright (c) 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded 2004 by: Eero Häkkinen
##############################################################################

package provide gom::Plugins 1.0

namespace eval gom {

namespace eval Plugins {

# Called by an initialization script (in GUI mode) or by user.
proc _InitPackage {} {
    package require gom::PackageUtilities
    set g  1
    set $g [show graphics]
    if { $g } {
	package require gom::gui::Plugins
    }
    foreach file [gom::ListPackageFiles plugins] {
	# Check if a plugin is compiled.
	if { [catch {_FindPlugin [file rootname [file tail $file]]} err] } {
	    puts $err
	    continue
	}
	# Init a plugin.
	namespace eval :: source [list $file]
    }
}

# Called by $gomEnv(GOM_DATA)/inittcl/plugins/$plugin.tcl.
proc RegisterPlugin { plugin what } {
    if { ! [show graphics] } return

    set ::gomPluginsMenuItems($plugin) \
	[list \
	     "Enable $what" \
	     [namespace code [list _EnableRegisteredPlugin $plugin $what]]]
    gom::gui::UpdatePluginsMenu
}

# Called by a menu command or by RegisterPlugin.
proc EnablePlugin { plugin } {
    # Load a plugin.
    LoadPlugin $plugin
    # Call a procedure ::<Plugin>Plugin::Enable once if it exists.
    set p [info procs ::[string totitle $plugin]Plugin::Enable]
    if { "" != $p } {
	$p
	rename $p ""
    }
    gomPrint "Plugin '$plugin' is enabled!"
}

proc LoadPlugin { plugin } {
    load [_FindPlugin $plugin]
}

proc _FindPlugin { plugin } {
    global gomEnv

    set libext [info sharedlibextension]
    foreach lib [list ${plugin}${libext} lib${plugin}${libext}] {
	set file [file join $gomEnv(GOM_PLUGINS) $lib]
	if { [file exists $file] } { return $file }
    }
    error "ERROR: No plugin '$plugin' in a directory '$gomEnv(GOM_PLUGINS)'."
}

proc _EnableRegisteredPlugin { plugin what } {
    global gomPluginsMenuItems
    EnablePlugin $plugin
    unset gomPluginsMenuItems($plugin)
    ::lulUpdatePluginsMenu
    ::lulMessageDialog "[string totitle $what] is enabled!"
}

namespace export {[A-Z]*}

}

namespace import Plugins::*

}

gom::Plugins::_InitPackage
