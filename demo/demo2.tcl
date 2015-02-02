############################################################################
# demo2
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
plot text green 0.1 0.9 "Molecular dynamics"
import coord karpl [file join $gomEnv(GOM_DEMO) opt4pti.crd]
trajectory file charmm [file join $gomEnv(GOM_DEMO) dyn4pti.dcd]
display

lulPauseDemo 5

#
# show frames
#
#lulResetText
#plot text green 0.1 0.9 "Show frames"

for {set i 1} { $i <= [show traject frames] } {incr i 2} {
 import coord frame $i
# lulResetText
# plot text green 0.1 0.9 "Show frame # $i"
display sleep 100
}

#
# show frames as a distance plot
#
lulResetText
plot text green 0.1 0.9 "Show frames as a series of LDP's"
lulPauseDemo 5

plot ldp atoms * * CA
plot ldp on

for {set i 1} { $i <= [show traject frames] } {incr i 2} {
 import coord frame $i
 lulResetText
 plot text green 0.1 0.9 "Show frame # $i"
 display sleep 100
}

plot ldp off

lulPauseDemo 5

#
# show frames (distance)
#
lulResetText
plot text green 0.1 0.9 "Show frames with a distance"
monitor distance * * 653 * * 677
monitor display distance on
gscale display 1.5
lulPauseDemo 3

for {set i 1} { $i <= [show traject frames] } {incr i 2} {
 import coord frame $i
 lulResetText
 plot text green 0.1 0.9 "Show frame # $i, Dist: [calculate distan * * 653 * * 677]"
 display sleep 100
}

monitor -distance

#
# show frames (angle)
#
lulResetText
plot text green 0.1 0.9 "Show frames with an angle"
monitor angle * * 653 * * 677 * * 474
monitor display angle on
lulPauseDemo 3

for {set i 1} { $i <= [show traject frames] } {incr i 2} {
 import coord frame $i
 lulResetText
 plot text green 0.1 0.9 "Show frame # $i, Ang: [calculate angle * * 653 * * 677 * * 474]"
 display sleep 100
}

monitor -angle

#
# show frames (torsion)
#
lulResetText
plot text green 0.1 0.9 "Show frames with a torsion angle "
monitor torsion * * 653 * * 677 * * 474 * * 569
monitor display torsion on
lulPauseDemo 3

for {set i 1} { $i <= [show traject frames] } {incr i 2} {
 import coord frame $i
 lulResetText
 plot text green 0.1 0.9 "Show frame # $i, Tors: [calculate torsion * * 653 * * 677 * * 474 * * 569]"
 display sleep 100
}

monitor -torsion
lulResetText

atom -displ * * *
atom  displ * TIP3 *
atom  cpk   * TIP3 *

gscale display 0.8

for {set i 1} { $i <= [show traject frames] } {incr i 10} {
 import coord frame $i
 lulResetText
 plot text green 0.1 0.9 "Show frame # $i (backbone as a tube)"
 plumber atoms * * CA 0.4 blue tube
 plumber display on
 display sleep 100
 plumber -atoms
}

lulResetText
atom -cpk

#
# show molecular movement as ellipsoids
#
atom -displ
plumber atoms * * CA 0.2 blue tube
plumber display on
lulResetText
plot text green 0.1 0.9 "Show C(alpha) movement as ellipsoids (* 5.0)"
p_rmsd * * CA 5.0 fit
display
lulPauseDemo 6
plumber -atoms

#
# show force at minimum
#

plot -text
plot text green 0.1 0.9 "Force acting on the atoms (CA) at min. struct."
import coord karpl [file join $gomEnv(GOM_DEMO) opt4pti.crd]
trajectory file charmm [file join $gomEnv(GOM_DEMO) dyn4pti.dcd]
import vector charmm [file join $gomEnv(GOM_DEMO) min_force.crd]

plot vector atoms * * CA 0.3 15.0
plot vector on
#atom colo vector
atom -displ
atom displ * * CA,N,C
atom lico  * * *
gscale displ 1.3
define bondstyle smooth

display

lulPauseDemo 5

#
# show force at a dyn frame
#

plot -text
plot text green 0.1 0.9 "Force acting on the atoms (CA) for a frame"
import coord karpl [file join $gomEnv(GOM_DEMO) opt4pti.crd]
trajectory file charmm [file join $gomEnv(GOM_DEMO) dyn4pti.dcd]
import vector charmm [file join $gomEnv(GOM_DEMO) f100_force.crd]
import coord frame 100

plot vector atoms * * CA 0.3 0.1
plot vector on
atom colo vector
atom -displ
atom displ * * CA,N,C
atom lico  * * *
gscale displ 1.3
plot cscale 10 0.0 110.

display

lulPauseDemo 5

plot -vector
plot -cscale
plot -text
atom -lico
atom displ
define bondstyle half

#
# end of show
#
