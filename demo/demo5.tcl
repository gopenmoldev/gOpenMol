#
# manipulations
#
# Leif Laaksonen 1999
#

    global env

# PROCS
proc lulResetText {} {
 plot -text
}

plot -text
plot text green 0.1 0.9 "Translate an atom..."
import coord karpl [file join $gomEnv(GOM_DEMO) oxaz_min.crd]
plot -axis
atom -lico
atom -cpk
plot -cylinder
plot -sphere
plot -arrow
display

select atoms * * N11
set Coords [show atom coord 11 1]
set xc [lindex $Coords 0]
set yc [lindex $Coords 1]
set zc [lindex $Coords 2]

lulPauseDemo 5

for {set i 0} {$i < 20} {incr i} {
translate selection 0.1 0.0 0.0
display sleep 100
}

lulPauseDemo 5

for {set i 0} {$i < 20} {incr i} {
translate selection -0.1 0.0 0.0
display sleep 100
}

plot -text
plot text green 0.1 0.9 "Rotate a set of atoms..."
select atoms * * 1-6

lulPauseDemo 5

for {set i 0} {$i < 30} {incr i} {
rotate selection 0.0 0.0 1.0
display sleep 100
}

lulPauseDemo 5

for {set i 0} {$i < 30} {incr i} {
rotate selection 0.0 0.0 -1.0
display sleep 100
}

lulPauseDemo 5

plot -axis
atom -lico
atom -cpk
plot -cylinder
plot -sphere
plot -arrow
plot -text
plot text green 0.1 0.9 "Rotate systems ..."
import coord karpl [file join $gomEnv(GOM_DEMO) water.crd]
atom lico
plot axis * * 1,4
monitor distance * * 1 * * 4
monitor display distance on
display

lulPauseDemo 5

# first water molecule down ...
select atoms * * 1-3
for {set i 1} {$i < 50} {incr i} {
  translate selection 0.0 -0.1 0.0 
  display sleep 100
}

# second water molecule up ...
select atoms * * 4-6
for {set i 1} {$i < 50} {incr i} {
  translate selection 0.0 0.1 0.0 
  display sleep 100
}

# first water molecule right ...
select atoms * * 1-3
for {set i 1} {$i < 50} {incr i} {
  translate selection 0.12 0.0 0.0 
  display sleep 100
}

# second water molecule up ...
select atoms * * 4-6
for {set i 1} {$i < 50} {incr i} {
  translate selection -0.12 0.0 0.0 
  display sleep 100
}

# rotate system
for {set i 1} {$i < 90} {incr i} {
  rotate display 0.0 1.0 0.0
  display sleep 100
}

# translate system
for {set i 1} {$i < 40} {incr i} {
  translate selection 0.0 0.0 0.1
  display sleep 100
}

# rotate water molecule
for {set i 1} {$i < 180} {incr i} {
  rotate select 0.0 1.0 0.0
  display sleep 100
}

for {set i 1} {$i < 20} {incr i} {
  gscale display 1.05
  display sleep 100
}
plot -axis
atom -lico
atom -cpk
plot -cylinder
plot -sphere
plot -arrow
plot -text
monitor -distance
