############################################################################
# This is a very simple solution to the following problem:
#
# What I would like to do is read in the structure for the 
# first frame (equilibrium geometry), then load up frame 5 (one of the 
# extrema of the vibrational trajectory) and the draw arrows from 
# the center of the each atom in the equilibrium geometry to the 
# position that the atom would have in the extremum frame.
#
# Copyright Leif Laaksonen/CSC 2001 Version 0.1
#
# Input text string can be:
# 1) frame number (numbering goes from 1 ...) where it is frame number + 1 in the file
# 2) frame number, arrow radius
# 3) frame number, arrow radius and colour
#
# Example
#                  lulVectorDemo   1          # plot frame nr 1
#                  lulVectorDemo  {1 0.1}     # plot frame number 1 vector radius 0.1
#                  lulVectorDemo  {1 0.1 red} # plot frame number 1 vector radius 0.1 in red
#
# If the colour is not given the default atom colour is used for the
# vectors
#
# This procedure deletes all previously defined arrows!
#
proc lulVectorDemo { input } {

     global gomControlFont
     global gomPictureType
     global gomPictureExt
     global gomHelpDir
     global gomHelpFile
     global env

# return if no molecular systems defined
     set NumFrames [show trajectory frames]
	 if {$NumFrames < 1} {
            gomError {ERROR: no trajectory is defined. Define a trajectory first!}
            return
	 }

# now read in the desired frame
         set frame [expr [lindex $input 0] + 1]
         set rad   [lindex $input 1]
         set color [lindex $input 2]

         if {$rad == ""}    {set rad   0.1}

         if {$frame < 1 || $frame > $NumFrames} {
            gomError "frame index ('$frame') out of allowed ramge (1 and $NumFrames)"
            return
         }
# read first frame #frame and save the coordinates
         import coord frame $frame
# now loop the atoms
         for {set i 1} { $i <= [show numatoms 1]} { incr i} {
            set xc($i) [lindex [show atom coordin $i 1] 0]
            set yc($i) [lindex [show atom coordin $i 1] 1]
            set zc($i) [lindex [show atom coordin $i 1] 2]
         }
         import coord frame 1
# calculate vector data
         plot -arrow
         for {set i 1} { $i <= [show numatoms 1]} { incr i} {
            set xcc [lindex [show atom coordin $i 1] 0]
            set ycc [lindex [show atom coordin $i 1] 1]
            set zcc [lindex [show atom coordin $i 1] 2]
# arrow
            if {$color == ""} {
                set colort [show atom color $i 1]
            } else {
                set colort "$color"
            }
            plot arrow $xcc $ycc $zcc $xc($i) $yc($i) $zc($i) $rad "$colort"  append 
            unset xc($i)
            unset yc($i)
            unset zc($i)
         }
# throw it on the screen
         display
   
}
