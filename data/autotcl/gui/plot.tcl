
#######################################################################
# PROC
proc lulPlot { MinX MaxX MinY MaxY } {

set w .plot
catch {destroy $w}
toplevel $w
wm title $w "Plot"
wm iconname $w "Plot"
#positionWindow $w
set c $w.c

frame $w.buttons
pack $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -command "destroy $w"
pack $w.buttons.dismiss -side left -expand 1

set width  600
set height 500

canvas $c -relief raised -width $width -height $height
pack $w.c -side top -fill x

set plotFont -*-Helvetica-Medium-R-Normal--*-180-*-*-*-*-*-*

# window pixels
set Wx [expr ($width  - 150)/10]
set Wy [expr ($height - 100)/10]

set Ox  100
set Oy [expr $height - 50]

$c create line  $Ox                   $Oy \
               [expr $Ox + 10 * $Wx]  $Oy  -width 2

$c create line  $Ox  $Oy \
                $Ox  [expr $Oy - 10 * $Wy] -width 2

$c create text 225 20 -text "A Simple Plot" -font $plotFont -fill brown

set Dx [expr ($MaxX - $MinX)/10.]
set Dy [expr ($MaxY - $MinY)/10.]

# first x-axis
for {set i 0} {$i <= 10} {incr i} {
    set x [expr {$Ox + ($i * $Wx)}]
    $c create line $x $Oy $x [expr $Oy - 5] -width 2
    $c create text $x $Oy -text [expr "$MinX + $i * $Dx"] -anchor n -font $plotFont
}
# then y-axis
for {set i 0} {$i <= 10} {incr i} {
    set y [expr {$Oy - ($i * $Wy)}]
    $c create line $Ox  $y [expr $Ox + 5] $y -width 2
    $c create text $Ox  $y -text [expr $MinX + $i * $Dx] -anchor e -font $plotFont
}

foreach point {{12 56} {20 94} {33 98} {32 120} {61 180}
	{75 160} {98 223}} {
    set x [expr {100 + (3*[lindex $point 0])}]
    set y [expr {250 - (4*[lindex $point 1])/5}]
    set item [$c create oval [expr $x-6] [expr $y-6] \
	    [expr $x+6] [expr $y+6] -width 1 -outline black \
	    -fill SkyBlue2]
    $c addtag point withtag $item
}

$c bind point <Any-Enter> "$c itemconfig current -fill red"
$c bind point <Any-Leave> "$c itemconfig current -fill SkyBlue2"
$c bind point <1> "plotDown $c %x %y"
$c bind point <ButtonRelease-1> "$c dtag selected"
bind $c <B1-Motion> "plotMove $c %x %y"

set plot(lastX) 0
set plot(lastY) 0
}
# plotDown --
# This procedure is invoked when the mouse is pressed over one of the
# data points.  It sets up state to allow the point to be dragged.
#
# Arguments:
# w -		The canvas window.
# x, y -	The coordinates of the mouse press.

proc plotDown {w x y} {
    global plot
    $w dtag selected
    $w addtag selected withtag current
    $w raise current
    set plot(lastX) $x
    set plot(lastY) $y
}

# plotMove --
# This procedure is invoked during mouse motion events.  It drags the
# current item.
#
# Arguments:
# w -		The canvas window.
# x, y -	The coordinates of the mouse.

proc plotMove {w x y} {
    global plot
    $w move selected [expr $x-$plot(lastX)] [expr $y-$plot(lastY)]
    set plot(lastX) $x
    set plot(lastY) $y
}

proc lulWorld2Pixel {Maxx Minx Maxy Miny Xc Yc } {

     set wx [expr ($Xc - $Minx)/($Maxx - $Minx)]
     set wy [expr ($Yc - $Miny)/($Maxy - $Miny)]

}