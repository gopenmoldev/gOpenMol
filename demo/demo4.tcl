#
# This demo will demonstrate the usage of the various graphics objects
# available in gOpenMol'
#
# Leif Laaksonen 1999
#

    global env

# PROCS
proc lulResetText {} {
 plot -text
}

plot -cylinder
plot -sphere
plot -arrow
plot -text
plot text green 0.1 0.9 "Display of a cylinder..."
import coord karpl [file join $gomEnv(GOM_DEMO) oxaz_min.crd]
atom -display
plot cylinder 0.0 -4.0 0.0 0.0 4.0 0.0 1. red
display

lulPauseDemo 5

plot cylinder 0.0 -4.0 0.0 0.0 4.0 0.0 1. {red green}
display

lulPauseDemo 5

plot cylinder 0.0 -4.0 0.0 0.0 4.0 0.0 {1. 0.5} {red green}
display

lulPauseDemo 5

plot -cylinder
plot -text
plot text green 0.1 0.9 "Display of an arrow..."
plot arrow 0.0 4.0 0.0 0.0 -4.0 0.0 1. red
display

lulPauseDemo 5

plot -arrow
plot -cylinder
plot -text
plot text green 0.1 0.9 "Display of a sphere..."

plot sphere 0.0 -4.0 0.0 3. red
display 

lulPauseDemo 5

plot sphere 0.0 -4.0 0.0 3. red
display 

lulPauseDemo 5

set scalex 1.0
for {set i 1} {$i <= 50} {incr i} {

plot -arrow
plot arrow 0.0 4.0 0.0 0.0 [expr 3.0 * $scalex - 4.0] 0.0 0.5 blue
plot arrow 0.0 -12.0 0.0 0.0 [expr -3.0 * $scalex - 4.0] 0.0 0.5 blue append
plot sphere 0.0 -4.0 0.0 3. red 1. $scalex 1.0
set scalex [expr $scalex - 0.01]
display sleep 1000

}
lulPauseDemo 5

for {set i 1} {$i <= 50} {incr i} {

plot -arrow
plot arrow 0.0 4.0 0.0 0.0 [expr 3.0 * $scalex - 4.0] 0.0 0.5 blue
plot arrow 0.0 -12.0 0.0 0.0 [expr -3.0 * $scalex - 4.0] 0.0 0.5 blue append
plot sphere 0.0 -4.0 0.0 3. red 1. $scalex 1.0
set scalex [expr $scalex + 0.01]
display sleep 1000

}

lulPauseDemo 5

plot -cylinder
plot -sphere
plot -arrow
plot -text

atom display
