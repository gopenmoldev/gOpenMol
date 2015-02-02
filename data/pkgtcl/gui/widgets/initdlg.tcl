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

proc InitOptDlg { w title } {
    if { [ winfo exists $w ] } {
	raise $w
	focus $w
	return 0
    } else {
	toplevel $w
	wm title    $w $title
	wm iconname $w $title
	return 1
    }
}    

proc InitDlg { w title } {
    # return if no molecular systems defined
    if { [show molstructures] < 1 } {
	lulErrorDialog {ERROR: no structure available. Read a structure first!}
	return 0
    }
    return [InitOptDlg $w $title]
}

namespace export InitDlg InitOptDlg

}
