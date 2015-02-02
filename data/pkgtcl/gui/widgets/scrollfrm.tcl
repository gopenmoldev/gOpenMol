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
#   ScrollFrame <pathName> ?<options>?
# OPTIONS
#   OPTIONS from ScrolledWindow
#   OPTIONS from ScrollableFrame
# COMMANGS
#   <pathName> getframe
#   <pathName> cget <option>
#   <pathName> configure ?<option> ?<value> <option> <value> ...??
############################################################################
proc ScrollFrame {} {}; # for pkg_mkIndex
namespace eval ScrollFrame {
    variable Class [string range [namespace current] 2 end]
    Widget::define    $Class scrollfrm ScrolledWindow ScrollableFrame
    Widget::bwinclude $Class ScrolledWindow  .w
    Widget::bwinclude $Class ScrollableFrame .f

    proc create { path args } {
	# init
	variable Class
	Widget::init $Class $path $args
	# frame
	frame $path -class ScrollFrame
	# sub widgets
	eval [list ScrolledWindow  $path.w] [Widget::subcget $path .w  ]
	eval [list ScrollableFrame $path.f] [Widget::subcget $path .f  ]
	# pack
	$path.w setwidget $path.f
	pack $path.w -fill both -expand true
	# create a megawidget
	return [Widget::create $Class $path]
    }

    proc getframe  { path        } { return [$path.f getframe] }
    proc cget      { path option } { return [Widget::cget      $path $option] }
    proc configure { path args   } { return [Widget::configure $path $args  ] }
}
namespace export ScrollFrame
}
