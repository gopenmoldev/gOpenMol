#
# Demo program to spawn the contour program
# to generate the surface, read it into gOpenMol
# display the surface.
#
# This is done for all frames
#
proc lulDemoReaction { value } {

    global gomEnv

# return if no molecular systems defined
     set NumFrames [show trajectory frames]
    if {$NumFrames < 1} {
	gomError {ERROR: no trajectory is defined. Define a trajectory first!}
            return
    }

     define display off

    for {set i 1} { $i <= [show numfr]} {incr i} {

     import coord frame  $i

     set RootName "[file join $gomEnv(GOM_TEMP) frame.inp]"
     export input probesurf 1 $RootName 
     run probesurf $RootName
     contour file [file join $gomEnv(GOM_TEMP) probesurf.plt] eka
     contour plot eka $value red [expr $value - 0.2] green
     contour clip z+ 2.5
     contour clip z on
#     calc conn * * C
     display
#
     contour -file
#
     update idle
    }

}
