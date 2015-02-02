##########################################################################
# demo1
#
# Copyright Leif Laaksonen 1997
#
    global env
#
# PROCS
proc lulResetText {} {
 plot -text
}

plot -text
plot text green 0.1 0.9 "Simple stick display"
import coord karpl [file join $gomEnv(GOM_DEMO) oxaz_min.crd]
display

lulPauseDemo 5

#
# rotate
#
lulResetText
plot text green 0.1 0.9 "Rotate molecule"

for {set i 0} { $i <= 360 } {incr i 5} {
 rotate display 0.0 5.0 0.0
 display sleep 100
}

lulPauseDemo 5

#
# licorice display
#
lulResetText
plot text green 0.1 0.9 "Licorice display"

atom lico
display

lulPauseDemo 5
atom -lico

#
# CPK display
#
lulResetText
plot text green 0.1 0.9 "CPK display"

atom cpk
display

lulPauseDemo 5

#
# CPK run ...
#
lulResetText
plot text green 0.1 0.9 "CPK display with a running atom"
atom -cpk
for {set i 1} { $i <= [show numatoms] } {incr i} {
 atom cpk * * $i
 display sleep 500
}

lulPauseDemo 5

#
# CPK scaling ...
#
lulResetText
plot text green 0.1 0.9 "CPK atom scaling"

for {set i 1} { $i <= [show numatoms] } {incr i} {
 atom scale cpk 0.5 * * $i
 display sleep 500
}
