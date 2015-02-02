##########################################################################
# 
# 2000-07-31
#
##########################################################################
##########################################################################
# NAMESPACE Trajectory
#
namespace eval lulTrajectory {

###########################################################################
# PROC
proc ReadVelocities {} {

# return if no molecular systems defined
     if {[show molstructures] < 1} {
          lulErrorDialog {ERROR: no structure available. Read a structure first!}
	      return
	 }

     trajectory retrieve veloci on
     set NumAtoms   [show numatoms 1]
     set NumFrames  [show numframes]

      for {set j 1} {$j <= $NumFrames} {incr j} {

       import coord frame $j

       for {set i 1} {$i <= $NumAtoms} {incr i} {

          set Valuesv [show atom velocity $i]
          scan $Valuesv "%f %f %f" Xv Yv Zv
          set Valuesc [show atom coord $i 1]
          scan $Valuesc "%f %f %f" Xc Yc Zc
          plot arrow $Xc $Yc $Zc [expr $Xc + $Xv] [expr $Yc + $Yv] [expr $Zc + $Zv] 0.2 \
               [show atom colo $i 1] append
       }
          display
          pause 2
          plot -arrow
      }
}
}
