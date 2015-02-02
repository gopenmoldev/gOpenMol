###########################################################################
# demo3
#
# Copyright Leif Laaksonen 2001
#
#
    global env
#
# PROCS
proc lulResetText {} {
 plot -text
}

plot -text
plot text green 0.1 0.9 "Density display (no smoothing)"
import coord karpl [file join $gomEnv(GOM_DEMO) orbital.crd]

atom lico
gscale display 3.0
define licocyli 0.1
define licosphe 0.12
rotate display 0.0 -45. 0.0

lulPauseDemo 5

#
# electron density
#
plot -text
plot text green 0.1 0.9 "Total density display"
contour file [file join $gomEnv(GOM_DEMO) density.plt] toka
contour plot toka 0.1 violet
contour smooth toka on

display

lulPauseDemo 5

#
# transparency
#
plot -text
plot text green 0.1 0.9 "How to display covered different isosurfaces?"

display

lulPauseDemo 5

#
# Clip plane
#
plot -text
plot text green 0.1 0.9 "Use a clip plane cutting the positive z-axis"

contour plot toka 0.1 violet 0.2 red 0.3 blue 0.4 green 0.05 yellow 0.5 cyan

contour clipplane z+ 0.0
contour clipplane z  on
display

lulPauseDemo 5

#
# moving clip plane
#

    set zvalue  0.0

for {set i 0} {$i <= 10} {incr i} {
    plot -text
    plot text green 0.1 0.9 "Clip plane at z = $zvalue"
    contour clipplane z+ $zvalue
    display sleep 1000
    set  zvalue [expr $zvalue + 0.1]
}

    set  zvalue [expr $zvalue - 0.1]

for {set i 0} {$i <= 10} {incr i} {
    plot -text
    plot text green 0.1 0.9 "Clip plane at z = $zvalue"
    contour clipplane z+ $zvalue
    display sleep 1000
    set  zvalue [expr $zvalue - 0.1]
}

lulPauseDemo 5

plot -text
plot text green 0.1 0.9 "Use a clip plane cutting the negative x-axis and positive z-axis"

    contour clipplane x- 0.0
    contour clipplane x  on
    set xvalue  0.0

for {set i 0} {$i <= 10} {incr i} {
    plot -text
    plot text green 0.1 0.9 "Clip plane at x = $xvalue (z = 0.0)"
    contour clipplane x- $xvalue
    display sleep 1000
    set  xvalue [expr $xvalue - 0.1]
}

    set  xvalue [expr $xvalue + 0.1]

for {set i 0} {$i <= 10} {incr i} {
    plot -text
    plot text green 0.1 0.9 "Clip plane at x = $xvalue (z = 0.0)"
    contour clipplane x- $xvalue
    display sleep 1000
    set  xvalue [expr $xvalue + 0.1]
}

#
# show contour profile
#

plot -text
plot text green 0.1 0.9 "Isosurface profile using a clip and cut plane"

contour clipplane x off
contour clipplane y off
contour clipplane z off

contour clipplane z+ 0.0
contour clipplane z on

plot cutplane z toka -0.1 0.0 1.0

display

lulPauseDemo 5

contour -file
plot -text

#
# end of show
#
