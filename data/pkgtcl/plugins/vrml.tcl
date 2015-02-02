namespace eval VrmlPlugin {

variable static          0
variable staticstyle     0
variable trajectorystyle 0

proc makevrmlt {fn delay bs} {
	set NumFrames [show trajectory frames]
	if {$NumFrames < 1} {
           lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
           return
	}
	VRMLmovie $NumFrames $bs
	for {set i 1} {$i<=$NumFrames} {incr i } {
		import coord frame $i
		VRMLframe $i 
	}
	VRMLwrite $fn $delay
}

proc makevrmls {fn bs} {
	VRMLStill $fn $bs
}

proc dovrml { w } {
    variable static
    variable staticstyle
    variable trajectorystyle
    
    set filename [$w.file.filename get]
    if {"" == $filename} return

    if {$static} {
	makevrmls $filename $staticstyle
    } else {
	set delay [$w.trajectory.delay.value get]
	if {$delay<100} {
	    set delay 100
	}
	makevrmlt $filename $delay $trajectorystyle
    }
}	

proc browseFileName { w } {
    set filename [tk_getSaveFile -defaultextension "wrl" \
	-filetypes [list \
	    {{VRML files} {.wrl .WRL} TEXT} \
	    {{All files}  {*}         TEXT}]]

    if { "" != $filename } {
	$w.file.filename insert 0 $filename
    }
}

proc makevrml {} {
    global   gomControlFont

    variable static

    gom::LoadPlugin "vrml"

    set w       .vrmlrc
    catch       {destroy $w}
    toplevel    $w
    wm title    $w "Export VRML animation"
    wm iconname $w "Export VRML animation"

    set ns [namespace current]

    # Create buttons.
    frame       $w.buttons -borderwidth 2 -relief raised
    pack        $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
    button $w.buttons.apply   -text Apply  -font "$gomControlFont" \
        -command "${ns}::dovrml $w"
    pack        $w.buttons.dismiss $w.buttons.apply -side left -expand 1

    # Create a type frame.
    frame       $w.type -borderwidth 2 -relief ridge
    label       $w.type.label      -text "Type:"
    radiobutton $w.type.static     -text "Static Image" \
	-variable ${ns}::static -value 1 \
	-command "pack $w.static -side top -anchor w;     pack forget $w.trajectory"
    radiobutton $w.type.trajectory -text "Trajectory" \
	-variable ${ns}::static -value 0 \
	-command "pack $w.trajectory -side top -anchor w; pack forget $w.static"

    pack        $w.type -side top -anchor w -fill x
    pack        $w.type.label $w.type.static $w.type.trajectory -side left

    labelframe  $w.file -text "Filename" -borderwidth 2 -relief ridge -padx 2 -pady 2
    entry       $w.file.filename -width 60
    button      $w.file.browse   -text "Browse..." -command "${ns}::browseFileName $w"

    pack        $w.file -side top -anchor w -pady 2m
    pack        $w.file.filename $w.file.browse -side left -padx 2

    # Create a static frame.
    frame       $w.static
    frame       $w.static.style -borderwidth 2 -relief ridge
    label       $w.static.style.label -text "Display style:"
    radiobutton $w.static.style.licor -text "Licorice" \
	-variable ${ns}::staticstyle -value 0
    radiobutton $w.static.style.stick -text "Sticks" \
	-variable ${ns}::staticstyle -value 1
    pack        $w.static.style -side top -anchor w
    pack        $w.static.style.label $w.static.style.licor $w.static.style.stick -side left

    # Create a trajectory frame.
    frame       $w.trajectory
    frame       $w.trajectory.style -borderwidth 2 -relief ridge
    label       $w.trajectory.style.label -text "Display style:"
    radiobutton $w.trajectory.style.cpk   -text "CPK Only" \
	-variable ${ns}::trajectorystyle -value 0
    radiobutton $w.trajectory.style.ball  -text "Ball and Stick" \
	-variable ${ns}::trajectorystyle -value 1
    radiobutton $w.trajectory.style.licor -text "Licorice" \
	-variable ${ns}::trajectorystyle -value 2
    radiobutton $w.trajectory.style.stick -text "Sticks" \
	-variable ${ns}::trajectorystyle -value 3
    frame       $w.trajectory.delay -borderwidth 2 -relief ridge
    label       $w.trajectory.delay.label -text "Frame Delay (ms):"
    entry       $w.trajectory.delay.value -width 15
    pack        $w.trajectory.style -side top -anchor w
    pack        $w.trajectory.style.label $w.trajectory.style.cpk \
			$w.trajectory.style.ball $w.trajectory.style.licor \
			$w.trajectory.style.stick -side left
    pack        $w.trajectory.delay -side top -anchor w
    pack        $w.trajectory.delay.label $w.trajectory.delay.value -side left

    # Display the right frame.
    if {$static} {
	$w.type.static     invoke
    } else {
	$w.type.trajectory invoke
    }
}

if {[show graphics]} {
# Correction by Dan Macks, change  set ::gomPluginsMenuItems(vrml) "{VRML...} command [namespace current]::makevrml" to
# 20-10-2003
    set ::gomPluginsMenuItems(vrml) "{VRML...} [namespace current]::makevrml"
    gom::gui::UpdatePluginsMenu
}

}
