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

# import coordinates for water from a gaussian formatted checkpoint file
import coord gaus h2o.fch
# you have to issue a display command to see the last action
display
#wait a bit
pause 1

# rotate the display orientation to make the symmetry axis vertical
rotate display 90.0 0.0 0.0
display
pause 1

# CPK display
atom cpk
display
# call procedure to rotate -- note that you just give the parameter value without parentheses
slaRotateDisplay  450

#scaled cpk
atom scale cpk 0.1 * * *
display 
# now call our procedure to rotate the display
slaRotateDisplay  360

gscale   display   4.0
display
pause 1


#licorice display
atom  licorice  
display
pause 2

#define sphere and cylinder diameters to get a ball and stick display
define   licosphere   0.15
define   licocylinder 0.05
display
# now call our procedure to rotate the display
slaRotateDisplay  360

#load up a contour file
contour file mo.plt mo5
#plot contours -- -0.2 is red, +0.2 is blue
contour plot mo5 -0.2 red 0.2 blue
display
pause 2
#now make the contours 50% transparent
contour alpha mo5 .6
display

# now call our procedure to rotate the display
slaRotateDisplay  85

#now convert the contours to mesh type
contour   type        mo5 mesh
display
# now call our procedure to rotate the display
slaRotateDisplay  80

#turn off contours
contour display mo5 off

#now import a trajectory file -- the display command just displays the structure, not the animation
trajectory file xmol h2o001.xmol
display

pause 1

#plot some text on the display
plot text green 0.1 0.9 "Show frames with a distance"
#set up to monitor and display distance between atoms 1 and 3
monitor distance * * 1 * * 3
monitor display distance on


#now we step through the frames of the trajectory
#for this first trajectory, we are going to explicitly display the frame number so
#the next line turns off the (default) automatic frame number display
define  trajectory       fid           off
#the first for loop is just to repeat the animation 3 times
#the 2nd for loop actually steps through the trajectory
for {set j 1} { $j <= 3 } {incr j 1} {
	for {set i 1} { $i <= [show traject frames] } {incr i 1} {
		 import coord frame $i
# the next line clears the text
		lulResetText
#this line demonstrates plotting text and using a variable in it
		plot text green 0.1 0.9 "Show frame # $i"
		display sleep 80
	}
}
#turn off distance monitor
monitor display distance off
pause 1

#for this 2nd trajectory, we are going let gOpenMol do the automatic frame number display
define  trajectory       fid           on

#another way to clear the text in the display
plot -text
plot text green 0.1 0.9 "Show frames with an angle"
pause 1
#set up to display the angle between atoms 1, 2 and 3
monitor angle * * 1 * * 2 * * 3
monitor display angle on
for {set j 1} { $j <= 3 } {incr j 1} {
	for {set i 1} { $i <= [show traject frames] } {incr i 1} {
		 import coord frame $i
		display sleep 80
	}
}
