##########################################################################
# simple.tcl
#
#NOTE:  comment lines start with a # character
# simple demo tcl script 
#
    global env
#
# PROCS   This is here just to show how you can define procedures (functions)
#                    to avoid having to repeat sets of commands.  In this case the procedure
#                    only does one thing  -- erase any text that is on the display
proc lulResetText {} {
 plot -text
}

#here is a somewhat more complicated procedure.  It rotates the display in 5 degree increments 
#through a total angle of "RotAngle" degrees, pausing for 10 msec at each step
#note that in the procedure, the parameter is referenced with a preceeding $
#note also that my proc name is prefaced with sla (my initials)  The important point is to make
#sure that your proc name is unique, otherwise you might inadvertently interfere with some system procedure

proc slaRotateDisplay {RotAngle} {
	for {set i 0} { $i <= $RotAngle } {incr i 5} {
		rotate display 0.0 5.0 0.0
		display sleep 10
	 }
}



#  I like the orthographic, rather than perspective view
define projection orthographic

#now import a trajectory file -- the display command just displays the structure, not the animation
trajectory file xmol h2o001.xmol
display



# rotate the display orientation to make the symmetry axis vertical
rotate display 90.0 0.0 0.0
rotate display 0.0 90.0 0.0
display
pause 1

#white backgroud
define bgcolour white
#improve display quality for publication
define cylinderquality 100
define spherequality 100

#make atoms black
atom  color * * * Black

# CPK display
atom cpk
atom scale cpk 0.1 * * *

#make matte finish
define material specular red 0.0 
define material specular blue 0.0 
define material specular green 0.0

     set NumAtoms   [show numatoms 1]

       import coord frame 5

       for {set i 1} {$i <= $NumAtoms} {incr i} {

          set Valuesc [show atom coord $i 1]
          scan $Valuesc "%f %f %f" Xc Yc Zc
          plot arrow $Xc $Yc $Zc [expr $Xc + $Xv] [expr $Yc + $Yv] [expr $Zc + $Zv] 0.2 \
               [show atom colo $i 1] append
       }
          display
          pause 2
          plot -arrow
