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
#   FileEntry <pathName> ?<options>?
# OPTIONS
#   OPTIONS from Entry
#   OPTIONS from Button
#     -buttontext                  <Text> 
#         See -text. Defaults to "Browse ...".
#     -buttontextvariable          <Variable>
#         See -textvariable
#     -buttonbg, -buttonbackground <BackGround>
#         See -background
#   OPTIONS from tk_getOpenFile
#   OPTIONS from tk_getSaveFile
#   OPTIONS from tk_chooseDirectory
#   WIDGET-SPECIFIC OPTIONS
#     -type <Type>
#         <Type> may be "open", "save" or "directory". Defaults to "open".
# COMMANGS
#   <pathName> cget <option>
#   <pathName> configure ?<option> ?<value> <option> <value> ...??
############################################################################
namespace eval FileEntry {
    variable Class [string range [namespace current] 2 end]
    Widget::define    $Class fileentry Entry Button
    Widget::bwinclude $Class Button .browse \
	remove {-command -height -width} \
	prefix {button -text -textvariable -bg -background}
    Widget::bwinclude $Class Entry  .file \
	remove {-height -width}
    array set dlgCmd {
	open      tk_getOpenFile
	save      tk_getSaveFile
	directory tk_chooseDirectory
    }
    set dlgOpts {
	{-defaultextension String "" 0}
	{-filetypes        String "" 0}
	{-initialdir       String "" 0}
	{-initialfile      String "" 0}
	{-multiple         String "" 0}
	{-message          String "" 0}
	{-title            String "" 0}
    }
    Widget::declare   $Class [concat $dlgOpts {
	{-height           Int 0 0 {%d >= 0}}
	{-width            Int 0 0 {%d >= 0}}
	{-type             Enum   open 0 {open save directory}}
    }]

    proc create { path args } {
	# init
	variable Class
	Widget::init $Class $path \
	    [linsert $args 0 \
		 -buttontext "Browse ..." \
		 -underline  0]
	# frame
	frame $path \
	    -height [Widget::cget $path -height] \
	    -width  [Widget::cget $path -width ] \
	    -class  FileEntry
	# sub widgets
	set file [Widget::cget $path -initialdir \
		     ][Widget::cget $path -initialfile]
	eval [list Entry  $path.file  ] [Widget::subcget $path .file  ] \
	    [list -text $file]
	eval [list Button $path.browse] [Widget::subcget $path .browse] \
	    [list -command [namespace code [list $path  _browse]]]
	# pack
	pack $path.file   -side left -padx 2 -fill x -expand true
	pack $path.browse -side left -padx 2
	# create a megawidget
	return [Widget::create $Class $path]
    }
    proc _browse { path } {
	variable dlgCmd
	variable dlgOpts
	# Collect dialog options.
	# Not all dialog supports all options thus drop options with
	# empty values.
	set opts(-parent) $path
	foreach spec $dlgOpts {
	    set opt [lindex $spec 0]
	    set val [$path cget $opt]
	    if { "" != $val } { set opts($opt) $val }
	}
	# Set initial file/directory if possible.
	set file [$path cget -text]
	set type [$path cget -type]
	if { [file isfile $file] && $type != "directory" } {
	    set opts(-initialfile) $file
	} elseif { [file isdirectory $file] } {
	    set opts(-initialdir)  $file
	}
	# Show a dialog.
	set file [eval [list $dlgCmd($type)] [array get opts]]
	if { "" == $file } return
	cd [file dirname $file]
	$path configure -text $file
    }
    proc cget      { path option } { return [Widget::cget      $path $option] }
    proc configure { path args   } { return [Widget::configure $path $args  ] }
}
namespace export FileEntry
}
