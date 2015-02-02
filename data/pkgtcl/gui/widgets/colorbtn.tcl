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

namespace eval SelectColor {
    # BWidget uses lazy sourcing.Force SelectColor::dialog to be
    # sourced (path is missing thus dialog is not shown).
    catch { SelectColor::dialog }
    # Extend SelectColor::dialog.
    rename dialog _orig_dialog
    proc dialog { path args } {
	variable _use_alt_dlg false
	set opts(-parent) .[join [lrange [split $path .] 1 end-1] .]
	foreach { opt val } $args {
	    switch -- $opt {
		-color  { set opts(-initialcolor) $val }
		-parent -
		-title  { set opts($opt) $val }
	    }
	}
	if { "unix" == $::tcl_platform(platform) } {
	    set cmd "set _use_alt_dlg true; destroy $path"
	    after idle [list $path add \
			    -text "Alternative dialog" \
			    -command [namespace code $cmd]]
	    set colour [eval [list _orig_dialog $path] $args]
	    if { ! $_use_alt_dlg } { return $colour }
	}
	return [eval tk_chooseColor [array get opts]]
    }
}

namespace eval gom::gui::Widgets {
############################################################################
# SYNOPSIS
#   ColorButton <pathName> ?<options>?
# OPTIONS
#   OPTIONS from Button
#     -text <Text>
#         A text to be displayed in a button. Defaults to "Colour ...".
#   WIDGET-SPECIFIC OPTIONS
#     -command <Command>
#         A command to be executed after user has selected a new colour.
#     -colour, -color <Colour>
#         A colour value. Default to "white".
#     -variable, <Variable>
#         A variable to which a colour value is stored
#     -title <Title>
#         A title to be display in a colour selection dialog.
#         Defaults to "Choose a colour".
# COMMANGS
#   <pathName> cget <option>
#   <pathName> configure ?<option> ?<value> <option> <value> ...??
############################################################################
namespace eval ColorButton {
    variable Class [string range [namespace current] 2 end]
    Widget::define    $Class colorbtn Button
    Widget::bwinclude $Class Button .b \
	remove {
	    -command
	    -bg   -background -fg   -foreground
	    -activebackground -activeforeground
	}
    # NOTE: Be careful with synonyms.
    # E.g. Widget::hasChanged does not handle synonyms.
    # User of this widget is allowed to use synonyms but inside this
    # namespace we use only a form "colour".
    Widget::declare   $Class {
	{-command  String  "" 0}
	{-colour   String  "" 0}
	{-color    Synonym -colour}
	{-variable String  "" 0}
	{-title    String  "Choose a colour" 0}
    }
    Widget::syncoptions $Class Button .b {-colour -background}

    proc create { path args } {
	# init
	variable Class
	Widget::init $Class $path \
	    [linsert $args 0 -text "Colour ..."]
	# frame
	frame $path -class ColorButton
	# sub widgets
	eval [list Button $path.b \
		  -command [namespace code [list _selectColour $path]]] \
	    [Widget::subcget $path .b]
	# create a megawidget.
	Widget::create $Class $path
	# here we need a megawidget
	$path _setTrace  [$path cget -variable]
	if { "" == [$path cget -colour] } {
	    # get a colour from the variable
	    $path _trace
	    if { "" == [$path cget -colour] } {
		# use default colour
		$path _setColour "white"
	    }
	} else {
	    $path _setColour [$path cget -colour]
	}
	# pack
	pack $path.b -fill both -expand true
	# return $path as normal
	return $path
    }
    proc _trace { path } {
	catch {$path configure -colour [expr \$::[$path cget -variable]]}
    }
    proc _setTrace { path variable } {
	if { "" == $variable } return
	namespace eval :: \
	    [list trace add variable $variable write "$path _trace; #"]
    }
    proc _selectColour { path } {
	# show popup menu.
	set colour \
	    [SelectColor $path.dlg \
		 -color     [$path cget -colour] \
		 -parent    $path \
		 -placement below \
		 -title     [$path cget -title] \
		 -type      popup]
	if { "" != $colour } {
	    # change colour and execute a command.
	    $path configure -colour $colour
	    namespace eval :: [$path cget -command]
	}
    }
    proc _setColour { path bg } {
	if { [llength $bg] == 3 } {
	    lassign $bg red green blue
	    if { [string is int -strict $red  ] &&
		 [string is int -strict $green] &&
		 [string is int -strict $blue ] } {
		# bg is {Ired Igreen Iblue}; Ired, Igreen, Iblue in 0...255
		set bg [format "#%02x%02x%02x" $red $green $blue]
	    } else {
		# bg is {Fred Fgreen Fblue}; Fred, Fgreen, Fblue in 0.0...1.0
		set bg [format "#%04x%04x%04x" \
			    [expr round($red   * 0xffff)] \
			    [expr round($green * 0xffff)] \
			    [expr round($blue  * 0xffff)]]
	    }
	    $path configure -colour $bg
	    return
	}
	lassign [winfo rgb $path $bg] red green blue
	# red, green and blue will be in range of 0...0xffff.
	if { 0.299 * $red + 0.587 * $green + 0.114 > 0x8000 } {
	    set fg black
	} else {
	    set fg white
	}
	# update a button
	$path.b configure \
	    -foreground $fg -activeforeground $fg \
	    -background $bg -activebackground $bg
	# update a variable
	set variable [$path cget -variable]
	if { "" != $variable } { set ::$variable $bg }
    }
    proc configure { path args } {
	# a base configuration
	set res [Widget::configure $path $args]
	# Configure more only if necessary.
	# This is a) to avoid unnessary work
	#         b) to avoid race condition with tracing
	if { [Widget::hasChanged $path -variable variable] } {
	    $path _setTrace $variable
	}
	if { [Widget::hasChanged $path -colour colour] } {
	    $path _setColour $colour
	}
	return $res
    }
    proc cget      { path option } { return [Widget::cget      $path $option] }
}
namespace export ColorButton
}
