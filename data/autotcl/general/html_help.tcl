
# This code is based on:
#
# Here is a sample html viewer to demonstrate the library usage
# Copyright (c) 1995 by Sun Microsystems
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#
# This REQUIRES Tk4.0 -- make sure "wish" on the next line is a 4.0 version
# The next line is a TK comment, but a shell command \
#  exec wish4.0 -f "$0" "$@" & exit 0

#if {$tk_version < 4.0 || [regexp {b[123]} $tk_patchLevel] } {
#	puts stderr "This library requires TK4.0, this is only $tk_version, \
#			patchlevel $tk_patchLevel"
#	exit 1
#}

if {[catch {array get env *}]} {
	puts stderr "This library requires tcl7.4, this version is too old!"
	exit 1
}

#
# source first the library
#
source [file join $gomEnv(GOM_DATA) html_library.tcl]

# construct a simple user interface

######################################################################
# PROC
proc htmlSetup {} {

     global htmlUrl
     global gomControlFont

        set w .html
        catch {destroy $w}
        toplevel $w 
        wm title $w "HTML help"
        wm iconname $w "HTML"

	frame      $w.frame
	menubutton $w.frame.menu   -relief raised -bd 2   -text options... \
                    -menu $w.frame.menu.m
	entry      $w.frame.entry  -textvariable htmlUrl  -width 35
	label      $w.frame.file   -text file:
	label      $w.frame.status -textvariable htmlRunning -width 6 -relief ridge \
			    -bd 2 -padx 9 -pady 3
    button     $w.frame.back    -text "<=" -font {Arial 14 bold} \
                -command {htmlPreviousURL}
    button     $w.frame.forward -text "=>" -font {Arial 14 bold} \
                -command {htmlNextURL}

	label      $w.msg -textvariable htmlMessage -font {Arial 10 bold}
	scrollbar  $w.scrollbary  -command "$w.text yview"  -orient v
	scrollbar  $w.scrollbarx  -command "$w.text xview"  -orient h
	text       $w.text    -yscrollcommand "$w.scrollbary set" \
                -xscrollcommand "$w.scrollbarx set"           \
                -padx 3 -pady 3 -takefocus 0                  \
                -width 90 -height 30 -wrap none

	pack       $w.frame $w.msg -side top
    pack       $w.frame.back    -side left 
    pack       $w.frame.forward -padx 4 -side left 
	pack       $w.frame.file $w.frame.entry -side left
	pack       $w.frame.status -padx 4 -side left
	pack       $w.frame.menu -side left

    frame      $w.buttons -borderwidth 2 -relief raised -bd 2
    pack       $w.buttons -side bottom -fill x -pady 2m
    button     $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
                -command "destroy $w"
    pack       $w.buttons.dismiss -side left -expand 1
	pack       $w.scrollbarx -side bottom -expand 0 -fill x

	pack       $w.scrollbary -side right  -expand 0 -fill y
	pack       $w.text -side top -fill both -expand 1

	# set up some sample keyboard bindings for the text widget
	bind $w.frame.entry <Return> {htmlRender $htmlUrl}
	bind $w <End>   "$w.text yview end"
	bind $w <Home>  "$w.text yview 0.0"
	bind $w <Next>  "$w.text yview scroll  1 page"
	bind $w <Prior> "$w.text yview scroll -1 page"

	menu $w.frame.menu.m
	$w.frame.menu.m add command -label "option menu"
	$w.frame.menu.m add separator
	$w.frame.menu.m add command -label "font size" -foreground red
	$w.frame.menu.m add radiobutton -label small -value  2   -variable htmlSize \
		-command {HMset_state .html.text -size $htmlSize; htmlRender $htmlUrl}
	$w.frame.menu.m add radiobutton -label medium -value 4   -variable htmlSize \
		-command {HMset_state .html.text -size $htmlSize; htmlRender $htmlUrl}
	$w.frame.menu.m add radiobutton -label large -value  6   -variable htmlSize \
		-command {HMset_state .html.text -size $htmlSize; htmlRender $htmlUrl}
	$w.frame.menu.m add separator
	$w.frame.menu.m add command -label "indent level" -foreground red
	$w.frame.menu.m add radiobutton -label small  -value 0.6  -variable htmlIndent \
		-command {HMset_indent .html.text $htmlIndent}
	$w.frame.menu.m add radiobutton -label medium -value 1.2  -variable htmlIndent \
		-command {HMset_indent .html.text $htmlIndent}
	$w.frame.menu.m add radiobutton -label large  -value 2.4  -variable htmlIndent \
		-command {HMset_indent .html.text $htmlIndent}

}

# Go render a page.  We have to make sure we don't render one page while
# still rendering the previous one.  If we get here from a recursive 
# invocation of the event loop, cancel whatever we were rendering when
# we were called.
# If we have a fragment name, try to go there.

proc htmlRender {file} {
	global HM.text htmlUrl
	global htmlRunning htmlMessage
    global htmlStack
    global htmlStackDeep
    global htmlStackCurrent

	set fragment ""
	regexp {([^#]*)#(.+)} $file dummy file fragment
	if {$file == "" && $fragment != ""} {
		HMgoto .html.text $fragment
		return
	}
	HMreset_win .html.text
	set htmlRunning busy
	set htmlMessage "Displaying $file"
	update idletasks
	if {$fragment != ""} {
		HMgoto .html.text $fragment
	}
	set htmlUrl $file
	HMparse_html [get_html $file] {HMrender .html.text}
	set htmlRunning ready
	HMset_state .html.text -stop 1	;# stop rendering previous page if busy
	set htmlMessage ""
#
# Update stack
#
#       set  htmlStack($htmlStackDeep) $htmlUrl
#       set  htmlStackCurrent $htmlStackDeep
#       incr htmlStackDeep

}

# given a file name, return its html, or invent some html if the file can't
# be opened.

proc get_html {file} {
	global htmlHome

	if {[catch {set fd [open $file]} msg]} {
		return "
			<title>Bad file $file</title>
			<h1>Error reading $file</h1><p>
			$msg<hr>
			<a href=$htmlHome>Go htmlHome</a>
		"
	}
	set result [read $fd]
	close $fd
	return $result
}

# Override the library link-callback routine for the sample app.
# It only handles the simple cases.

proc HMlink_callback {win href} {
    global htmlUrl
    global htmlStack
    global htmlStackDeep
    global htmlStackCurrent

	if {[string match #* $href]} {
		htmlRender $href
		return
	}
	if {[string match /* $href]} {
		set htmlUrl $href
	} else {
		set htmlUrl [file dirname $htmlUrl]/$href
	}
        update
        htmlRender $htmlUrl
        set  htmlStack($htmlStackDeep) $htmlUrl
        set  htmlStackCurrent $htmlStackDeep
        incr htmlStackDeep
}

# Supply an image callback function
# Read in an image if we don't already have one
# callback to library for display

proc HMset_image {win handle src} {
	global htmlUrl htmlMessage
	if {[string match /* $src]} {
		set image $src
	} else {
		set image [file dirname $htmlUrl]/$src
	}
	set htmlMessage "fetching image $image"
	update
	if {[string first " $image " " [image names] "] >= 0} {
		HMgot_image $handle $image
	} else {
		set type photo
		if {[file extension $image] == ".bmp"} {set type bitmap}
		catch {image create $type $image -file $image} image
		HMgot_image $handle $image
	}
}

# Handle base tags.  This breaks if more than 1 base tag is in the document

proc HMtag_base {win param text} {
	global htmlUrl
	upvar #0 HM$win var
	HMextract_param $param href htmlUrl
}

# downloading fonts can take a long time.  We'll override the default
# font-setting routine to permit better user feedback on fonts.  We'll
# keep our own list of installed fonts on the side, to guess when delays
# are likely

proc HMset_font {win tag font} {
	global htmlMessage Fonts
	if {![info exists Fonts($font)]} {
		set Fonts($font) 1
		.html.msg configure -fg blue
		set htmlMessage "downloading font $font"
		update
	}
	.html.msg configure -fg black
	set htmlMessage ""
	catch {$win tag configure $tag -font $font} htmlMessage
}

# Lets invent a new HTML tag, just for fun.
# Change the color of the text. Use html tags of the form:
# <color value=blue> ... </color>
# We can invent a new tag for the display stack.  If it starts with "T"
# it will automatically get mapped directly to a text widget tag.

proc HMtag_color {win param text} {
	upvar #0 HM$win var
	set value bad_color
	HMextract_param $param value
	$win tag configure $value -foreground $value
	HMstack $win "" "Tcolor $value"
}

proc HMtag_/color {win param text} {
	upvar #0 HM$win var
	HMstack $win / "Tcolor {}"
}

# Add a font size manipulation primitive, so we can use this sample program
# for on-line presentations.  sizes prefixed with + or - are relative.
#  <font size=[+-]3>  ..... </font>.  Note that this is not the same as
# Netscape's <font> tag.

proc HMtag_font {win param text} {
	upvar #0 HM$win var
	set size 0; set sign ""
	HMextract_param $param size
	regexp {([+-])? *([0-9]+)} $size dummy sign size
	if {$sign != ""} {
		set size [expr [lindex $var(size) end] $sign $size]
	}
	HMstack $win {} "size $size"
}

# commented by LUL
# This version is closer to what Netscape does

#proc HMtag_font {win param text} {
#	upvar #0 HM$win var
#	set size 0; set sign ""
#	HMextract_param $param size
#	regexp {([+-])? *([0-9]+)} $size dummy sign size
#	if {$sign != ""} {
#		set size [expr [lindex $var(size) end] $sign  $size*2]
#		HMstack $win {} "size $size"
#	} else {
#		HMstack $win {} "size [expr 10 + 2 * $size]"
#	}
#}

proc HMtag_/font {win param text} {
	upvar #0 HM$win var
	HMstack $win / "size {}"
}

proc htmlShowHelp { htmlFile } {

     global htmlSize
     global htmlIndent
     global htmlHome
     global htmlUrl
     global htmlRunning
     global htmlMessage
     global gomHelpDir
     global htmlStack
     global htmlStackDeep
     global htmlStackCurrent

# set initial values
     set htmlStackDeep 0               ;# set stack depth to 0 (reset)
     set htmlStackCurrent 0            ;# current stack pointer
     set htmlSize 2                    ;# font size adjustment
     set htmlIndent 1.2                ;# tab spacing (cm)
     set htmlHome gopenmol.html        ;# htmlHome document
     set htmlUrl $htmlHome             ;# current file
     set htmlRunning busy              ;# page status
     set htmlMessage ""                ;# message line

# make the interface and render the htmlHome page
     catch htmlSetup		;# the catch lets us re-source this file
     HMinit_win   .html.text
     HMset_state  .html.text -size $htmlSize
     HMset_indent .html.text $htmlIndent

     if {$htmlFile == ""} {
         set  htmlStack($htmlStackDeep) [file join $gomHelpDir $htmlHome]
         set  htmlStackCurrent $htmlStackDeep
         incr htmlStackDeep
         htmlRender   [file join $gomHelpDir $htmlHome]
     } else {
         set  htmlStack($htmlStackDeep) [file join $gomHelpDir $htmlFile]
         set  htmlStackCurrent $htmlStackDeep
         incr htmlStackDeep
         htmlRender   [file join $gomHelpDir $htmlFile]
     }
}

proc htmlPreviousURL {} {

    global htmlStack
    global htmlStackDeep
    global htmlStackCurrent

set place [expr $htmlStackCurrent - 1]

if {$place < 0} {
    return
    }

#puts "URL: $htmlStack($place)"

set htmlStackCurrent $place

	htmlRender $htmlStack($place)
}

proc htmlNextURL {} {

    global htmlStack
    global htmlStackDeep
    global htmlStackCurrent

set place [expr $htmlStackCurrent + 1]

if {$place >= $htmlStackDeep} {
    return
    }

#puts "URL: $htmlStack($place)"

set htmlStackCurrent $place

	htmlRender $htmlStack($place)

}
