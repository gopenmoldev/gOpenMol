##############################################################################
#                           Copyright (c) 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

package require BWidget

namespace eval gom::gui::Widgets {
############################################################################
# SYNOPSIS
#   MainButtonFrame <pathName> ?<options>?
# OPTIONS
#     -dismisscommand DismissCommand
#     -applycommand   ApplyCommand
#     -acceptcommand  AcceptCommand
#     -helpcommand    HelpCommand
#         Commands to assosiate to buttons.
#         Default DismissCommand is to destroy a toplevel window.
#         Create key bindings <Escape>, <Alt-A>, <Return> and <F1>.
# COMMANGS
#   <pathName> cget <option>
#   <pathName> configure ?<option> ?<value> <option> <value> ...??
############################################################################
namespace eval MainButtonFrame {
    variable Class [string range [namespace current] 2 end]
    Widget::define  $Class buttonfrm Button
    foreach i {apply accept dismiss help} {
	Widget::bwinclude $Class Button .f.$i \
	    prefix  [list $i -command]
    }
    proc create { path args } {
	global   gomControlFont
	set toplevel .[join [lrange [split $path .] 1 end-1] .]
	# init
	variable Class
	Widget::init $Class $path \
	    [linsert $args 0 -dismisscommand [list destroy $toplevel]]
	# frame
	frame $path
	# sub frame
	frame $path.f -borderwidth 2 -relief raised
	pack  $path.f -fill x -expand true -pady 10
	# sub widgets
	foreach \
	    w     {dismiss apply  accept help} \
	    text  {Dismiss Apply  Accept Help} \
	    index {-1      0      1      -1  } \
	    key   {Escape  {}     {}     F1  } \
	{
	    set opts [Widget::subcget $path .f.$w]
	    if { "" == $opts } continue
	    eval [list Button $path.f.$w -text $text \
		      -underline $index -font $gomControlFont] $opts
	    if { "" != $key } {
		bind $toplevel <$key> [list $path.f.$w invoke]
	    }
	    pack $path.f.$w -side left -expand true -padx 10
	}
	return [Widget::create $Class $path]
    }
    proc cget      { path option } { return [Widget::cget      $path $option] }
    proc configure { path args   } { return [Widget::configure $path $args  ] }
}
namespace export MainButtonFrame
}
