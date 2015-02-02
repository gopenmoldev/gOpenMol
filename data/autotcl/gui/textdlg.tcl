################### text widget with text input #############
#
# Text widget
#
proc lulCreateTextWidget {w Title} {
    global gomControlFont

    catch {destroy $w}
    toplevel $w 
    wm title $w $Title
    wm iconname $w "text"

    frame $w.buttons -borderwidth 2 -relief raised -bd 2
    pack $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
	-command "destroy $w"
    pack $w.buttons.dismiss -side left -expand 1

    text $w.text -relief sunken -bd 2 -yscrollcommand "$w.scrolly set" \
	-xscrollcommand "$w.scrollx set" -setgrid 1                  \
	-height 30 -font {Courier 10 normal} -wrap none
    scrollbar $w.scrolly -command "$w.text yview"
    scrollbar $w.scrollx -orient horizontal -command "$w.text xview"
    pack $w.scrolly -side right -fill y
    pack $w.text -expand yes -fill both
    pack $w.scrollx -side bottom -fill x
    
    return $w.text
}
################### text widget with text input #############
#
# Text widget
#
proc lulDisplayTextText InputText {
    set w [lulCreateTextWidget .gomtext "Display text"]
    
    $w delete 0.0 end
    $w insert 0.0 $InputText
    $w mark set insert 0.0
}
################### text widget with text input #############
#
# Text widget
#
# PROC
proc lulDisplayTextFile InputFile {
    set w [lulCreateTextWidget .gomtext "Display text"]

    $w delete 0.0 end

    set File [open $InputFile r]
    while {![eof $File]} {
	$w insert end [read $File 1000]
    }

    $w mark set insert 0.0

    close $File
}

set gomSkipErrorMsg ""
###############################################################################
#
# dialog widget
#
proc lulErrorDialog InputText {
    global gomSkipErrorMsg
    set    gomSkipErrorMsg $InputText

    after idle {.errdialog.msg configure -wraplength 5i}
    set i [tk_dialog .errdialog "gOpenMol ERROR" $InputText error 0 OK Debug]
    switch $i {
	1 {set gomSkipErrorMsg ""}
    }
}

###############################################################################
#
# dialog widget
#
proc lulMessageDialog InputText {
    after idle {.dialog1.msg configure -wraplength 4i}
    tk_dialog .dialog1 "gOpenMol Message" $InputText info 0 OK
}

###############################################################################
#
# Replacer for the default Tcl background error handler.
#
proc lulErrorStackDialog InputText {
    global errorInfo

    set code  ok
    set stack $errorInfo

    after idle {.errstackdialog.msg configure -wraplength 5i}
    set i [tk_dialog .errstackdialog "gOpenMol ERROR" $InputText error 0 \
	       OK "Stack trace"]
    switch $i {
	1 {
	    # Stack trace
	    set top .gomerrorstack
	    set w [lulCreateTextWidget $top "gOpenMol Error - Stack trace"]
	    $w configure -wrap word
	    $w delete 0.0 end
	    $w insert 0.0 $stack
	    $w mark set insert 0.0
	    grab          $top
	    tkwait window $top
	    grab release  $top
	}
    }
}

###############################################################################
#
# Replacer for the default Tcl background error handler.
#
proc bgerror InputText {
    global errorInfo
    global gomSkipErrorMsg

    set code  ok
    set stack $errorInfo

    # Check if this error message is already shown by lulErrorDialog.
    if {"x$gomSkipErrorMsg" != "x$InputText"} {
	after idle {.bgerrdialog.msg configure -wraplength 5i}
	set i [tk_dialog .bgerrdialog "gOpenMol ERROR" $InputText error 0 \
		   OK "Stack trace" Skip]
	switch $i {
	    1 {
		# Stack trace
		set top .gomerrorstack
		set w [lulCreateTextWidget $top "gOpenMol Error - Stack trace"]
		$w configure -wrap word
		$w delete 0.0 end
		$w insert 0.0 $stack
		$w mark set insert 0.0
		grab          $top
		tkwait window $top
		grab release  $top
	    }
	    2 {
		# Skip pending background error messages.
		set code break
	    }
	}
    }

    set gomSkipErrorMsg ""
    return -code $code
}
