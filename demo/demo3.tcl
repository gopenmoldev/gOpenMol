###########################################################################
# demo3
#
# Copyright Leif Laaksonen 1997
#
#
    global env
#
# PROCS
proc lulResetText {} {
 plot -text
}

plot -text
plot text green 0.1 0.9 "Orbital display (no smoothing)"
import coord karpl [file join $gomEnv(GOM_DEMO) orbital.crd]

atom lico
gscale display 3.0
define licocyli 0.1
define licosphe 0.12
rotate display 0.0 -45. 0.0
contour file [file join $gomEnv(GOM_DEMO) orbital.plt] eka
contour plot eka -0.15 blue 0.15 red
display

lulPauseDemo 5

plot -text
plot text green 0.1 0.9 "Orbital display (smoothing)"
contour smooth eka on
display

lulPauseDemo 5

#
# electron density
#
plot -text
plot text green 0.1 0.9 "Orbital + total density display"
contour file [file join $gomEnv(GOM_DEMO) density.plt] toka
contour plot toka 0.1 violet
contour smooth toka on

display

lulPauseDemo 5

#
# transparency
#
plot -text
plot text green 0.1 0.9 "Transparent density (60%)"
contour alpha toka 0.6

display

lulPauseDemo 5

#
# cut plane through the density
#
plot -text
plot text green 0.1 0.9 "Total density as a cut plane"
plot cscale 10 0.0 1.0

contour display toka off

plot cutplane z toka 0.0 0.0 1.0

display

lulPauseDemo 5

#
# total density as two cut planes
#

plot -text
plot text green 0.1 0.9 "Total density with two cut planes"

contour display eka off

plot cutplane x toka 0.0 0.0 1.0

display

lulPauseDemo 5

#
# total density as two cut planes and move one
#

plot -text
plot text green 0.1 0.9 "Total density with two cut planes"

contour display eka off

    set xvalue  0.0

for {set i 0} {$i <= 10} {incr i} {
    plot cutplane x toka $xvalue 0.0 1.0
    display sleep 1000
    set  xvalue [expr $xvalue - 0.1]
}

for {set i 0} {$i <= 18} {incr i} {
    plot cutplane x toka $xvalue 0.0 1.0
    display sleep 1000
    set  xvalue [expr $xvalue + 0.1]
}
    lulPauseDemo 5

define cutplane damping -1.0
plot cutplane 3dx
set xvalue 0.0

for {set i 0} {$i <= 10} {incr i} {
    plot cutplane x toka $xvalue 0.0 1.0
    display sleep 1000
    set  xvalue [expr $xvalue - 0.1]
}

for {set i 0} {$i <= 18} {incr i} {
    plot cutplane x toka $xvalue 0.0 1.0
    display sleep 1000
    set  xvalue [expr $xvalue + 0.1]
}

plot cutplane 2dx
contour -file
plot -cscale
plot -text

plot -text
plot text green 0.1 0.9 "The electron density coloured according to the orbital values (mapping)"
import coord karpl [file join $gomEnv(GOM_DEMO) density.crd]
contour file [file join $gomEnv(GOM_DEMO) density.plt] eka
contour file [file join $gomEnv(GOM_DEMO) orbital.plt] toka
contour plot {eka toka} 0.1 -0.15 0.15
gscale display 3.0
rotate display 0.0 -45. 0.0
display

lulPauseDemo 5

contour -file
plot -text

#
# end of show
#
