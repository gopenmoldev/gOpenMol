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
#   RadioGroup <pathName> ?<options>?
# OPTIONS
#   OPTIONS from frame
#     -background <BackGround>
#     -height     <Height>
#     -width      <Width>
#   OPTIONS from labelframe
#     -labelanchor <Anchor>
#     -label       <Label>
#         Use labelframe widget. See -text.
#   WIDGET-SPECIFIC OPTIONS
#     -orient   <Orient>
#         May be "horizontal" (default) or "vertical".
#     -rows     <Rows>
#         A maximum number of rows. Default is no maximum.
#     -columns  <Columns>
#         A maximum number of columns. Default is no maximum.
#     -texts    {<Text>    <Text>    <Text>    ...}
#     -commands {<Command> <Command> <Command< ...}
#     -values   {<Value>   <Value>   <Value>   ...}
#         Lists containing option values for radiobutton widgets.
#     -variable <Variable>
#         A target variable for Values. Defaults to the last part of
#         <pathName>.
# COMMANGS
#   <pathName> cget <option>
#   <pathName> configure ?<option> ?<value> <option> <value> ...??
#   <pathName> add            ?<option> <value> ...?
#   <pathName> insert <index> ?<option> <value> ...?
#   <pathName> delete <index>
#   <pathName> index  <index>
#   <pathName> itemcget <index> <option>
#   <pathName> itemconfigure <index> ?<option> ?<value> <option> <value> ...??
#   <pathName> invoke   <index>
#   <pathName> select   <index>
#   <pathName> deselect <index>
#   <pathName> flash    <index>
############################################################################
namespace eval RadioGroup {
    variable Class [string range [namespace current] 2 end]
    Widget::define  $Class koe
    set frameOpts {
	{-background  TkResource "" 0 frame}
	{-bg          Synonym    -background}
	{-height      TkResource "" 0 frame}
	{-width       TkResource "" 0 frame}
    }
    set labelframeOpts {
	{-labelanchor TkResource "" 0 labelframe}
    }
    set radiobuttonOpts {
	{-commands    String     "" 0}
	{-command     String     "" 0}
	{-texts       String     "" 0}
	{-values      String     "" 0}
	{-variable    String     "" 0}
    }
    Widget::declare $Class \
	[concat $frameOpts $labelframeOpts $radiobuttonOpts {
	    {-orient      Enum       horizontal 1 {horizontal vertical}}
	    {-columns     Int        0  1 "%d >= 0"}
	    {-rows        Int        0  1 "%d >= 0"}
	    {-label       String     "" 0}
	}]

    # NOTE: An idea is to put radiobuttons to a upper-left corner of
    # a expandable frame.
    # To implement this, we create a sub frame into an upper-left
    # bounding box of 2x2 grid in a main frame. The upper-left
    # bounding box is kept fixed. The three other bounding boxes are
    # configured to fill the rest of the main frame.
    proc create { path args } {
	# init
	variable Class
	Widget::init $Class $path \
	    [linsert $args 0 -variable [lindex [split $path .] end]]
	# frame
	eval [list frame $path -class RadioGroup] [_getOpts $path frame]
	# sub frame
	set label [Widget::cget $path -label]
	if { "" == $label } {
	    frame $path.lf
	} else {
	    eval [list labelframe $path.lf -text $label] \
		[_getOpts $path labelframe]
	}
	# sub widgets
	Widget::getVariable $path incrCount
	set incrCount 0
	set variable [Widget::cget $path -variable]
	set command0 [Widget::cget $path -command ]
	set i 0
	foreach \
	    text    [Widget::cget $path -texts   ] \
	    value   [Widget::cget $path -values  ] \
	    command [Widget::cget $path -commands] \
	{
	    if { "" == $command } { set command $command0 }
	    _add $path \
		[list \
		     -text     $text \
		     -command  $command \
		     -variable $variable \
		     -value    $value]
	}
	# pack
	_repack $path
	grid $path.lf
	grid columnconfigure $path 1 -weight 1
	grid rowconfigure    $path 1 -weight 1
	# create a megawidget
	return [Widget::create $Class $path]
    }
    proc _getOpts { path class } {
	# Get options but drop options with empty values.
	variable Class
	upvar \#0 ${Class}::${class}Opts classOpts
	set opts [list]
	foreach spec $classOpts {
	    set opt [lindex $spec 0]
	    set val [Widget::cget $path $opt]
	    if { "" != $val } { lappend opts $opt $val }
	}
	return $opts
    }
    proc _add { path options } {
	# buttons is a list of radiobuttons.
	# Because user is allowed to delete and insert radio buttons,
	# $buttons is not sorted.
	Widget::getVariable $path buttons
	Widget::getVariable $path incrCount; # A total count of widget created
	set button $path.lf.r$incrCount
	eval [list radiobutton $button] $options
	lappend buttons $button
	incr    incrCount
    }
    proc _repack { path } {
	Widget::getVariable $path buttons
	catch {eval grid forget [grid slaves $path.lf]}
	# Extract layout options.
	foreach opt {columns rows} {
	    set lim($opt) [Widget::cget $path -$opt]
	}
	switch [Widget::cget $path -orient] {
	    horizontal { array set dim {1 column 2 row} }
	    vertical   { array set dim {1 row 2 column} }
	}
	array set pos {column 0 row 0}
	foreach button $buttons {
	    grid $button -column $pos(column) -row $pos(row) -sticky nw
	    # Move to a next column/row.
	    incr pos($dim(1))
	    if { [expr $pos($dim(1))] == $lim($dim(1)s) } {
		# A row/column if full. Move to the next row/column.
		set  pos($dim(1)) 0
		incr pos($dim(2))
	    }
	}
    }
    proc cget      { path option } { return [Widget::cget      $path $option] }
    proc configure { path args   } { return [Widget::configure $path $args  ] }
    proc add { path args } {
	$path _add $args
	$path _repack
    }
    proc delete { path index } {
	Widget::getVariable $path buttons
	set  button  [lrange   $buttons $index $index]
	set  buttons [lreplace $buttons $index $index]
	$path _repack
	eval destroy $button
    }
    proc index { path index } {
	Widget::getVariable $path buttons
	return [lindex $buttons $index]
    }
    proc insert { path index args } {
	$path _add $args
	# Re-organize radio buttons.
	Widget::getVariable $path buttons
	set buttons \
	    [linsert [lrange $buttons 0 end-1] $index [lindex $buttons end-1]]
	$path _repack
    }
    # Create simple sub commands.
    # Use format to make it look nicer.
    foreach i {invoke select deselect flash} {
	proc $i { path index } [format {
	    [$path index $index] %s
	} $i]
    }
    foreach i {cget configure} {
	proc item$i { path index args } [format {
	    return [eval [list [$path index $index] %s] $args]
	} $i]
    }
}
namespace export RadioGroup
}
