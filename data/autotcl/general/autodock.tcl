##########################################################################
# 
# This is a first approximation to create an analysis and display
# tool for AutoDock trajectories.
# Currently it is not helpful for more than just the display of
# the conformations. 
#
# Copyright Leif Laaksonen CSC 1998
# Leif.Laaksonen@csc.fi
#
# This is the alpha version of the new AutoDock Crucher that
# can also handle AutoDock 3 output files.
#
# I have no personal knowledge of the information in an output
# file so I may interpreted something wrong. Please let me know
# and I will correct it.
#
# 1999-12-07
#
# 2001-03-02: A small fix for changed Autodock v. 3.05 output format 
#
##########################################################################
#
# AutoDock wrapper
#
#################################################################
# AutoDock command
proc autodock args {

    global env

if {[lindex $args 0] == ""} {
	    gomError "AutoDock parameters missing"
		return
}

set Action [lindex $args 0]

# ligand
  if {[string match liga* $Action]} {

# check that there is an input file given...

    set InputFile [lindex $args 1]

    if { $InputFile == ""} {
          gomError "input file name (ligand) missing!"
	      return
    }

    import coordi pdb $InputFile

# change to current directory
    lulChangeDirectory $InputFile

# move information to AutoDock namespace
    set lulAutoDock::lulAutoDockLigandFileName $InputFile

    return
  } elseif {[string match syst* $Action]} {
# system
# check that there is an input file given...

    set InputFile [lindex $args 1]

    if { $InputFile == ""} {
          gomError "input file name (system) missing!"
	      return
    }

    import coordi pdb $InputFile append

# change to current directory
    lulChangeDirectory $InputFile

# move information to AutoDock namespace
    set lulAutoDock::lulAutoDockSystemFileName $InputFile

    return

  } elseif {[string match traj* $Action]} {
# trajectory
# check that there is an input file given...

    set InputFile [lindex $args 1]

    if { $InputFile == ""} {
          gomError "trajectory file name  missing!"
	      return
    }

# open widget
    lulAutoDock::MainControl
# get number of frames and set pointers
    lulAutoDock::GetNumberOfAutoDockFrames $InputFile

# move information to AutoDock namespace
    set lulAutoDock::lulAutoDockTrajFileName $InputFile

    return

  }
}
##########################################################################
# NAMESPACE AutoDock
#
namespace eval lulAutoDock {

     variable lulNumAtoms
     variable lulNumFrames
     variable lulAutoDockLigandFileName   ""
     variable lulAutoDockSystemFileName   ""
     variable lulAutoDockTrajFileName     ""
     variable lulAutoDockLogFileName      ""
     variable lulCurrentFrame
     variable lulLoopControl
     variable lulFramePointer
     variable lulFrameRetrieveMethod 1
     variable lulAutoDockHelp           "main_autodock_widget.html"
     variable lulCurrentAutoDockRun
     variable lulCurrentAutoDockCycle
     variable lulLigandAppend 0

###########################################################################
# PROC
proc GetNumberOfAutoDockFrames {FileName} {

     variable lulNumAtoms
     variable lulNumFrames
     variable lulAutoDockLigandFileName
     variable lulAutoDockSystemFileName
     variable lulAutoDockTrajFileName
     variable lulCurrentFrame
     variable lulFramePointer
     variable lulAnnealingRuns
     variable lulAnnealingRunArray
     variable lulAnnealingCycles
     variable lulAnnealingCycleArray
     variable lulAnnealingSteps
     variable lulCurrentAutoDockRun
     variable lulCurrentAutoDockCycle

     set i 0
     set lulNumAtoms [show numat 1]
     if {$lulNumAtoms < 1} {
         gomError "No structure available. Read in structure first"
         return 0
     }   
# open file ...
   set File_p [open $FileName r]
   if {$File_p == ""} {
       catch {gomError {ERROR: can't open AutoDock traj file!}} errRet
       if {$errRet != ""} {
          gomPrint "ERROR: can't open AutoDock traj file!"
          }
       error "ERROR: can't open AutoDock traj file!"
       return 0
   }

# change to current directory
   lulChangeDirectory $FileName

   set NumFrames 0
# runs
   set lulAnnealingRuns        1
   set lulAnnealingRunArray(1) 1
# cycles
   set lulAnnealingCycles(1)        1
   set lulAnnealingCycleArray(1,1) 1
#
   set lulAutoDockTrajFileName $FileName
   set RunsTrigger    0
   set CyclesTrigger  0
   set ww .gomautodock

   while {![eof $File_p]} {

    set lulFramePointer([expr $NumFrames + 1]) [tell $File_p]

    for {set i 1} {$i <= $lulNumAtoms} {incr i} {
     if {[eof $File_p]} {
          break
     }
     gets $File_p InputData
    }

     if {$InputData == ""} continue

     incr NumFrames
#
     set Dummy1 [lindex $InputData 12]
     set Dummy2 [lindex $InputData 11]

     if {$NumFrames > 1} {

         if {$RunsTrigger != $Dummy1} {
           incr lulAnnealingRuns
           set RunsTrigger   $Dummy1
           set CyclesTrigger $Dummy2
           set lulAnnealingRunArray($lulAnnealingRuns) $NumFrames
           set lulAnnealingCycles($lulAnnealingRuns)        1
           set lulAnnealingCycleArray($lulAnnealingRuns,1) $NumFrames
         } else {
           if {$CyclesTrigger != $Dummy2} {
             set CyclesTrigger $Dummy2
             incr lulAnnealingCycles($lulAnnealingRuns)
             set lulAnnealingCycleArray($lulAnnealingRuns,$lulAnnealingCycles($lulAnnealingRuns)) $NumFrames
           }
         }

     } else {
         set RunsTrigger   $Dummy1
         set CyclesTrigger $Dummy2
     }
#
     $ww.left.frame3.frame31.frames configure -state normal
     $ww.left.frame3.frame31.frames delete 0 end
     $ww.left.frame3.frame31.frames insert 0 $NumFrames
     $ww.left.frame3.frame31.frames configure -state disabled
# update widget
     update
   }

   puts "Number of frames found: $NumFrames"
   set lulNumFrames $NumFrames

# start from first frame
   set lulCurrentFrame  1
#
  set base .gomautodock.left.frame3.frame30.run
  $base configure -text "Run # 1"
  set base $base.menu
  destroy $base
  set m [menu $base]
  for {set i 1} {$i <= $lulAnnealingRuns} {incr i} {
   if {[expr $i % 20]} {
    $m add command -label "# $i"  -columnbreak 0 -command "lulAutoDock::ChooseAutoDockRun $i"
   } else {
    $m add command -label "# $i"  -columnbreak 1 -command "lulAutoDock::ChooseAutoDockRun $i"
   }
  }

  set base .gomautodock.left.frame3.frame30.cycle
  $base configure -text "Cycle # 1"
  set base $base.menu
  destroy $base
  set m [menu $base]
  for {set i 1} {$i <= $lulAnnealingCycles(1)} {incr i} {
   if {[expr $i % 20]} {
    $m add command -label "# $i"  -columnbreak 0 -command "lulAutoDock::ChooseAutoDockCycle $i"
   } else {
    $m add command -label "# $i"  -columnbreak 1 -command "lulAutoDock::ChooseAutoDockCycle $i"
   }
  }

# default is 1 , 1
     set lulCurrentAutoDockRun    1
     set lulCurrentAutoDockCycle  1

     GetFrameLimits
# done ...
   close $File_p
   return $NumFrames
} 
# end of proc

###########################################################################
# PROC
proc ChooseAutoDockRun { Which } {

     variable lulCurrentAutoDockRun
     variable lulCurrentAutoDockCycle

     set lulCurrentAutoDockRun    $Which
     set lulCurrentAutoDockCycle  1

     set base .gomautodock.left.frame3.frame30.run
     $base configure -text "Cycle # $Which"

     set base .gomautodock.left.frame3.frame30.cycle
     $base configure -text "Run # 1"

     GetFrameLimits
}

###########################################################################
# PROC
proc ChooseAutoDockCycle { Which } {

     variable lulCurrentAutoDockRun
     variable lulCurrentAutoDockCycle

     set lulCurrentAutoDockCycle $Which

     set base .gomautodock.left.frame3.frame30.cycle
     $base configure -text "Run # $Which"

     GetFrameLimits

}


###########################################################################
# PROC
proc GetAutoDockFrame {FileName FrameNumber} {

     global   gomADframeType
     variable lulNumAtoms
     variable lulNumFrames
     variable lulFramePointer
     variable lulFrameRetrieveMethod
     variable lulCurrentFrame

     if {$lulNumAtoms < 1} {
         gomError "No structure defined. read in structure first"
         return 1
     }   

     if {$FrameNumber > $lulNumFrames} {
         gomError "Frame index > total number of frames"
         return 1
     }   

# open file ...
   set File_p [open $FileName r]
   if {$File_p == ""} {
       catch {gomError {ERROR: can't open AutoDock traj file!}} errRet
       if {$errRet != ""} {
          gomPrint "ERROR: can't open AutoDock traj file!"
          }
       error "ERROR: can't open AutoDock traj file!"
       return 1
   }

# change to current directory
   lulChangeDirectory $FileName

# get translate array
   scan [show translation array] "%f %f %f" Xtrans Ytrans Ztrans

   if { !$lulFrameRetrieveMethod } {
     if {$FrameNumber > 1} {
       for {set i 1} {$i <= [expr $FrameNumber - 1]} {incr i} {
        for {set j 1} {$j <= $lulNumAtoms} {incr j} {
         gets $File_p InputData
        }
       }
     }
   } else {
     seek $File_p $lulFramePointer($FrameNumber) start
   }

   for {set i 1} {$i <= $lulNumAtoms} {incr i} {

          gets $File_p InputData

          set Xc       [expr [lindex $InputData 0] - $Xtrans]
          set Yc       [expr [lindex $InputData 1] - $Ytrans]
          set Zc       [expr [lindex $InputData 2] - $Ztrans]
          set code     [lindex $InputData 9]

          switch $gomADframeType {
            1 { if {$code < 3} {close $File_p;return}}
            2 { if {$code > 2} {close $File_p;return}}
          }

          define atom coordinates $Xc $Yc $Zc $i 1
#puts "atom coordinates $Xc $Yc $Zc $i 1"
   }
# plot stats...
   set energy [lindex $InputData 7]
   set temp   [lindex $InputData 8]
   set code   [lindex $InputData 9]
   plot -text
   if {$code < 3} {
      plot  text green 0.05 0.95 "Frame # $FrameNumber : E(T) = $energy : T = $temp"
   } else {
      plot  text red   0.05 0.95 "Frame # $FrameNumber : E(T) = $energy : T = $temp"
   }

   set lulCurrentFrame $FrameNumber

   close $File_p
   return
} 
# end of proc

###########################################################################
# PROC
proc Testing {FileName} {

     variable lulNumFrames

     GetNumberOfAutoDockFrames $FileName

     for {set i 1} {$i <= $lulNumFrames} {incr i} {

     GetAutoDockFrame $FileName $i
     display
     }
}
# end of proc
############################################################################
# PROC
proc MainControl {} {

     global gomHelpDir
     global gomHelpFile
     global gomTrajType
     global gomTrajFileName
     global gomControlFont
     global gomADframeType
     variable lulAutoDockLigandFileName
     variable lulAutoDockSystemFileName
     variable lulAutoDockTrajFileName
     variable lulAutoDockLogFileName
     variable lulCurrentFrame
     variable lulNumFrames
     variable lulFrameRetrieveMethod
     variable lulAutoDockHelp
     variable lulCurrentAutoDockRun    "*NULL*"
     variable lulCurrentAutoDockCycle  "*NULL*"
     variable lulLigandAppend

set w .gomautodock
catch {destroy $w}
toplevel $w 
wm title $w "AutoDock Main Control"
wm iconname $w "AutoDock Control"

#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
#button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
#        -command "lulDoTrajFileParams  $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
       "htmlShowHelp $lulAutoDockHelp"

pack   $w.buttons.dismiss $w.buttons.help -side left -expand 1
#    
	frame $w.left 

#	frame $w.right -borderwidth 2 -relief ridge 

	pack $w.left   -side top   -anchor n -fill both
#	pack $w.right  -side right  -anchor e

# Ligand file
    frame $w.left.frame11 -borderwidth 2 -relief raised
    set base $w.left.frame11
	label $base.filenamelabel \
		-text {Ligand file name:} -width 20

	entry $base.filename  -width 30

# make a <return> bind
    bind $base.filename <Return> "lulAutoDock::ImportFiles ligand $base.filename"

	button $base.filebrowse \
		-padx 10 \
		-pady 3 \
		-text Browse... -command "lulAutoDock::GetAutoDockFiles ligand $base.filename"

    pack $base -side top -anchor n
	pack $base.filenamelabel -side left
    pack $base.filename -side left -padx 3
	pack $base.filebrowse -side left

    if {$lulAutoDockLigandFileName != ""} {
         $base.filename delete 0 end
         $base.filename insert 0 $lulAutoDockLigandFileName
    }

    button $base.button -text "Import file" -command "lulAutoDock::ImportFiles ligand $base.filename;lulInvalidateDisplay" \
            -font {Arial 10 bold}
	pack $base.button -padx 4
    checkbutton $base.append   -text "Append structure" \
             -variable lulLigandAppend
    pack $base.append -side left -padx 3
	pack $base -side top -anchor n -fill x

# System
    frame $w.left.frame12 -borderwidth 2 -relief raised
    set base $w.left.frame12
	label $base.filenamelabel \
		-text {System file name:} -width 20

	entry $base.filename -width 30

# make a <return> bind
    bind $base.filename <Return> "lulAutoDock::ImportFiles system $base.filename"

	button $base.filebrowse \
		-padx 10 \
		-pady 3 \
		-text Browse... -command "lulAutoDock::GetAutoDockFiles system $base.filename"

    pack $base -side top -anchor n
	pack $base.filenamelabel -side left
    pack $base.filename -side left -padx 3
	pack $base.filebrowse -side left

    if {$lulAutoDockSystemFileName != ""} {
         $base.filename delete 0 end
         $base.filename insert 0 $lulAutoDockSystemFileName
    }

    button $base.button -text "Import file" -command "lulAutoDock::ImportFiles system $base.filename;lulInvalidateDisplay" \
            -font {Arial 10 bold}
	pack $base -side top -anchor n -fill x
	pack $base.button -padx 4

# trajectory
    frame $w.left.frame1 -borderwidth 2 -relief raised
    set base $w.left.frame1
	label $base.filenamelabel \
		-text {Traj file name:} -width 20

	entry $base.filename -width 30

# make a <return> bind
    bind $base.filename <Return> "lulAutoDock::ImportFiles traj $base.filename"

	button $base.filebrowse \
		-padx 10 \
		-pady 3 \
		-text Browse... -command "lulAutoDock::GetAutoDockFiles traj $base.filename"

    pack $base -side top -anchor n
	pack $base.filenamelabel -side left
    pack $base.filename -side left -padx 3
	pack $base.filebrowse -side left

    if {$lulAutoDockTrajFileName != ""} {
         $base.filename delete 0 end
         $base.filename insert 0 $lulAutoDockTrajFileName
    }

    button $base.button -text "Import file" -command "lulAutoDock::ImportFiles traj $base.filename" \
            -font {Arial 10 bold}
	pack $base -side top -anchor n -fill x
	pack $base.button -padx 4

###
# log file
    frame $w.left.frame13 -borderwidth 2 -relief raised
    set base $w.left.frame13
	label $base.filenamelabel \
		-text {Log file name:} -width 20

	entry $base.filename -width 30

# make a <return> bind
    bind $base.filename <Return> "lulAutoDock::ImportFiles log $base.filename"

	button $base.filebrowse \
		-padx 10 \
		-pady 3 \
		-text Browse... -command "lulAutoDock::GetAutoDockFiles log $base.filename"

    pack $base -side top -anchor n
	pack $base.filenamelabel -side left
    pack $base.filename -side left -padx 3
	pack $base.filebrowse -side left

    if {$lulAutoDockLogFileName != ""} {
         $base.filename delete 0 end
         $base.filename insert 0 $lulAutoDockLogFileName
    }

    button $base.button -text "Crunch!" -command {lulAutoDock::LogFileCruncher [.gomautodock.left.frame13.filename get]} \
            -font {Arial 10 bold} -width 9
	pack $base -side top -anchor n -fill x
	pack $base.button -padx 4

###
    frame $w.left.frame3 -borderwidth 2 -relief ridge
    pack  $w.left.frame3 -anchor w -side top -pady 4

    frame $w.left.frame3.frame30
    set base $w.left.frame3.frame30
	pack $base  -anchor w -side top
    menubutton $base.run -text "Run #: " -menu $base.run.menu \
                -relief raised
    pack   $base.run -side left -anchor w -pady 2
    set m [uplevel menu $base.run.menu]
    $m add command -label "*NULL*" -command {gomError "no trajectory available"}

	menubutton $base.cycle -text "Cycle #: " -menu $base.cycle.menu \
                -relief raised
    pack   $base.cycle -side left -anchor w -padx 6 -pady 2
    set m [uplevel menu $base.cycle.menu]
    $m add command -label "*NULL*" -command {gomError "no trajectory available"}

    frame $w.left.frame3.frame31
    set base $w.left.frame3.frame31
	label $base.frameslabel -text "Num of Frames:" -width 15
	entry $base.frames      -width 20 -state disabled
	pack $base  -anchor w -side top
	pack $base.frameslabel   -side left
	pack $base.frames        

    frame $w.left.frame3.frame32
    set base $w.left.frame3.frame32
    label $base.firstlabel  -text "First frame:" -width 15
	entry $base.first       -width 20
	pack  $base -anchor w -side top
	pack  $base.firstlabel   -side left
    pack  $base.first        

    frame $w.left.frame3.frame33
    set base $w.left.frame3.frame33
    label $base.lastlabel  -text "Last frame:" -width 15
	entry $base.last       -width 20
	pack  $base -anchor w -side top
	pack  $base.lastlabel $base.last  -side left

    frame $w.left.frame3.frame34
    set base $w.left.frame3.frame34
    label $base.steplabel  -text "Step frame(s):" -width 15
	entry $base.step       -width 20
	pack  $base -anchor w -side top
	pack  $base.steplabel $base.step  -side left

    frame $w.left.frame3.frame35
    set base $w.left.frame3.frame35
    label $base.currlabel  -text "Current frame:" -width 15
	entry $base.current       -width 20
	pack  $base -anchor w -side top
	pack  $base.currlabel $base.current  -side left

#    bind $base.current <Return> "lulAutoDock::JumpToFrame [$base.current get]"

    frame $w.left.frame2 -borderwidth 2 -relief ridge
    set base $w.left.frame2
    pack   $base -side left -anchor w

	button $base.first    -text "|<"   -command {lulAutoDock::GetFirstFrameInSet} \
            -font {Arial 12 bold}
	pack   $base.first    -side left  -padx 3 -pady 4 
	button $base.back     -text "<"    -command "lulAutoDock::TrajectoryControl -1" \
            -font {Arial 12 bold}
	pack   $base.back     -side left  -padx 3 -pady 4
	button $base.stop     -text "\[\]" -command "lulAutoDock::StopLooping"     \
            -font {Arial 12 bold}
   	pack   $base.stop     -side left  -padx 3 -pady 4
	button $base.play     -text "|>"   -command "lulAutoDock::TrajectoryControl 0"     \
            -font {Arial 12 bold}
	pack   $base.play     -side left  -padx 3 -pady 4
	button $base.forward  -text ">"    -command "lulAutoDock::TrajectoryControl +1"  \
            -font {Arial 12 bold}
	pack   $base.forward  -side left  -padx 3 -pady 4
	button $base.last     -text ">|"   -command {lulAutoDock::GetLastFrameInSet}     \
            -font {Arial 12 bold}
	pack   $base.last     -side left  -padx 3 -pady 4


      frame $w.left.frame4 
      pack  $w.left.frame4 -side top -anchor e -pady 5 -padx 10

#
      frame $w.left.frame4.which -borderwidth 2 -relief ridge
      pack  $w.left.frame4.which -side left -anchor w -pady 2 -padx 2 -fill x

      label   $w.left.frame4.which.label -text "Show frame type: "
      pack    $w.left.frame4.which.label -side top  -anchor w
      radiobutton $w.left.frame4.which.accepted -text "Accepted"  -value 2  -variable gomADframeType 
      radiobutton $w.left.frame4.which.rejected -text "Rejected"  -value 1  -variable gomADframeType 
      radiobutton $w.left.frame4.which.both     -text "Both"      -value 0  -variable gomADframeType 

      pack    $w.left.frame4.which.accepted $w.left.frame4.which.rejected $w.left.frame4.which.both     \
                  -side top -anchor w
#
      $w.left.frame4.which.both select

#
      frame $w.left.frame4.method -borderwidth 2 -relief ridge
      pack  $w.left.frame4.method -side top -anchor w -pady 2 -padx 2 -fill x

      label   $w.left.frame4.method.label -text "Trajectory retrieval: "
      pack    $w.left.frame4.method.label -side top  -anchor w
      radiobutton $w.left.frame4.method.fast  -text "Fast"  -value 1  -variable gomMethod \
                  -command  {set lulAutoDock::lulFrameRetrieveMethod 1}
      radiobutton $w.left.frame4.method.slow -text "Slow" -value 0    -variable gomMethod \
                  -command  {set lulAutoDock::lulFrameRetrieveMethod 0}
      pack    $w.left.frame4.method.fast $w.left.frame4.method.slow          \
                  -side top -anchor w
#
      if {$lulFrameRetrieveMethod} {
              $w.left.frame4.method.fast select
      } else {
              $w.left.frame4.method.slow select
      }

# write out section
# Ligand file
    frame $w.frame51 -borderwidth 2 -relief raised
    set base $w.frame51
	label $base.filenamelabel \
		-text {Ligand file name:} -width 20

	entry $base.filename  -width 30

# make a <return> bind
    bind $base.filename <Return> "lulAutoDock::MakePDBQCoordinates 1 $base.filename 1"

	button $base.filebrowse \
		-padx 10 \
		-pady 3 \
		-text Browse... -command "lulAutoDock::PutAutoDockFiles ligand $base.filename"

    pack $base -side top -anchor w
	pack $base.filenamelabel -side left
    pack $base.filename -side left -padx 3
	pack $base.filebrowse -side left

    button $base.button -text "Export file" -command "lulAutoDock::MakePDBQCoordinates 1 $base.filename 1" \
            -font {Arial 10 bold}
	pack $base -side top -anchor n -fill x
	pack $base.button -padx 4

# System
    frame $w.frame52 -borderwidth 2 -relief raised
    set base $w.frame52
	label $base.filenamelabel \
		-text {System file name:} -width 20

	entry $base.filename  -width 30

# make a <return> bind
    bind $base.filename <Return> "lulAutoDock::MakePDBQCoordinates 2 $base.filename 1"

	button $base.filebrowse \
		-padx 10 \
		-pady 3 \
		-text Browse... -command "lulAutoDock::PutAutoDockFiles system $base.filename"

    pack $base -side top -anchor w
	pack $base.filenamelabel -side left
    pack $base.filename -side left -padx 3
	pack $base.filebrowse -side left

    button $base.button -text "Export file" -command "lulAutoDock::MakePDBQCoordinates 2 $base.filename 1" \
            -font {Arial 10 bold}
	pack $base -side top -anchor n -fill x
	pack $base.button -padx 4

# end of proc
}

############################################################################
# PROC
proc GetFirstFrameInSet { } {

     variable lulAutoDockTrajFileName
     variable lulCurrentFrame
     variable lulNumFrames
     variable lulFrameRetrieveMethod
     variable lulAutoDockHelp
     variable lulCurrentAutoDockRun
     variable lulCurrentAutoDockCycle

     set Jump [string trim [.gomautodock.left.frame3.frame32.first get]]
     if {$Jump == ""} {
          gomError "'first frame' entry field is empty"
          return
     }

     set lulCurrentFrame $Jump
     GetAutoDockFrame $lulAutoDockTrajFileName $Jump
     UpdateCurrentFrame
     display
}
############################################################################
# PROC
proc GetLastFrameInSet { } {

     variable lulAutoDockTrajFileName
     variable lulCurrentFrame
     variable lulNumFrames
     variable lulFrameRetrieveMethod
     variable lulAutoDockHelp
     variable lulCurrentAutoDockRun
     variable lulCurrentAutoDockCycle

     set Jump [string trim [.gomautodock.left.frame3.frame33.last get]]
     if {$Jump == ""} {
          gomError "'last frame' entry field is empty"
          return
     }
 
     set lulCurrentFrame $Jump
     GetAutoDockFrame $lulAutoDockTrajFileName $Jump
     UpdateCurrentFrame
     display
}

############################################################################
# PROC
proc GetFrameLimits { } {

     variable lulCurrentAutoDockRun
     variable lulCurrentAutoDockCycle
     variable lulNumFrames
     variable lulAnnealingRuns
     variable lulAnnealingRunArray
     variable lulAnnealingCycles
     variable lulAnnealingCycleArray
     variable lulCurrentFrame

   if {$lulCurrentAutoDockRun == "*NULL*" || $lulCurrentAutoDockCycle == "*NULL*"} return
# first frame
   set FirstFrame $lulAnnealingCycleArray($lulCurrentAutoDockRun,$lulCurrentAutoDockCycle)
# last frame
   if {$lulAnnealingCycles($lulCurrentAutoDockRun) == $lulCurrentAutoDockCycle} {

       if {$lulCurrentAutoDockRun == $lulAnnealingRuns} {
         set LastFrame $lulNumFrames 
       } else {
         set LastFrame $lulAnnealingCycleArray([expr $lulCurrentAutoDockRun + 1],1) 
# correct last frame number...
         set LastFrame [expr $LastFrame - 1]
       }

   } else {
         set LastFrame $lulAnnealingCycleArray($lulCurrentAutoDockRun,[expr $lulCurrentAutoDockCycle + 1])
# correct last frame number...
         set LastFrame [expr $LastFrame - 1]
   }

   set eFrames  .gomautodock.left.frame3.frame31.frames
   $eFrames configure -state normal
   $eFrames delete 0 end
   $eFrames insert 0 [expr $LastFrame - $FirstFrame + 1]
   $eFrames configure -state disabled

   set eFirst   .gomautodock.left.frame3.frame32.first
   $eFirst configure -state normal
   $eFirst delete 0 end
   $eFirst insert 0 $FirstFrame
   $eFirst configure -state disabled

   set eLast    .gomautodock.left.frame3.frame33.last
   $eLast configure -state normal
   $eLast delete 0 end
   $eLast insert 0 $LastFrame
   $eLast configure -state disabled

   set eStep    .gomautodock.left.frame3.frame34.step
   $eStep configure -state normal
   $eStep delete 0 end
   $eStep insert 0 1


   set eCurrent .gomautodock.left.frame3.frame35.current
   $eCurrent configure -state normal
   $eCurrent delete 0 end
   $eCurrent insert 0 $FirstFrame
   $eCurrent configure -state disabled

   set lulCurrentFrame $FirstFrame

# get first frame and display it...
   GetFirstFrameInSet
}

############################################################################
# PROC
proc GetAutoDockFiles { what w } {

     variable lulAutoDockLigandFileName
     variable lulAutoDockSystemFileName
     variable lulAutoDockTrajFileName
     variable lulAutoDockLogFileName 


    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    
	if {$what == "ligand"} {
       set types {
         {"Lignad file"	".pdb .PDB .pdbq .PDBQ"	   TEXT}
         {"All files"		*}
       }
       set FileType 0
    } elseif {$what == "system"} {
       set types {
         {"System file"	".pdb .PDB .pdbq .PDBQ .pdbqs .PDBQS"	   TEXT}
         {"All files"		*}
       }
       set FileType 1
    } elseif {$what == "traj"} {
       set types {
         {"Traj file"	".trj .TRJ .tout .TOUT"			TEXT}
         {"All files"		*}
       }
       set FileType 2
    } elseif {$what == "log"} {
       set types {
         {"Traj file"	".dlg .DLG"			TEXT}
         {"All files"		*}
       }
       set FileType 3
    } 
    
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
        -title "AutoDock file"]]

    if [string compare $file ""] {
     $w delete 0 end
	 $w insert 0 $file

     switch FileType {
       0 {set lulAutoDockLigandFileName $file}
       1 {set lulAutoDockSystemFileName $file}
       2 {set lulAutoDockTrajFileName   $file}
       3 {set lulAutoDockLogFileName    $file}
     }
    }
# end of proc
}

############################################################################
# PROC
proc PutAutoDockFiles { what w } {

     variable lulAutoDockLigandFileName
     variable lulAutoDockSystemFileName
     variable lulAutoDockTrajFileName


    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    
	if {$what == "ligand"} {
       set types {
         {"Lignad file"	".pdb .PDB .pdbq .PDBQ"	   TEXT}
         {"All files"		*}
       }
    } else {
       set types {
         {"System file"	".pdb .PDB .pdbq .PDBQ"	   TEXT}
         {"All files"		*}
       }
    }
        
    set file [string trim [tk_getSaveFile -filetypes $types -parent $w \
        -title "AutoDock file"]]

    if [string compare $file ""] {
     $w delete 0 end
	 $w insert 0 $file
    }
# end of proc
}


############################################################################
# PROC
proc TrajectoryControl { What } {

     variable lulAutoDockTrajFileName
     variable lulCurrentFrame
     variable lulNumFrames
     variable lulLoopControl

# what to do?

# step back
  if {$What == "-1"} {
    set step [.gomautodock.left.frame3.frame34.step get]
    if {$step == ""} return

    set m1 [expr $lulCurrentFrame - $step]
    set Jump [string trim [.gomautodock.left.frame3.frame32.first get]]
    if {$Jump == ""} {
         gomError "'first frame' entry field is empty"
         return
    }

    if {$m1 < $Jump} {
       display
    } else {
       set lulCurrentFrame $m1
       eval lulAutoDock::GetAutoDockFrame $lulAutoDockTrajFileName $m1
       lulAutoDock::UpdateCurrentFrame;
       display
    }
# step forward
  } elseif {$What == "+1"} {
    set step [.gomautodock.left.frame3.frame34.step get]
    if {$step == ""} return

    set p1 [expr $lulCurrentFrame + $step]
    set Jump [string trim [.gomautodock.left.frame3.frame33.last get]]
    if {$Jump == ""} {
         gomError "'last frame' entry field is empty"
         return
    }

    if {$p1 > $Jump} {
       display
    } else {
       set lulCurrentFrame $p1
       eval lulAutoDock::GetAutoDockFrame $lulAutoDockTrajFileName $p1
       lulAutoDock::UpdateCurrentFrame;
       display
    }
# loop
  } elseif {$What == "0"} {

    set step [string trim [.gomautodock.left.frame3.frame34.step get]]
    if {$step == ""} return
    set Jump [string trim [.gomautodock.left.frame3.frame33.last get]]
    if {$Jump == ""} {
         gomError "'last frame' entry field is empty"
         return
    }

      set lulLoopControl 1

      while {$lulLoopControl} {

        set p1 [expr $lulCurrentFrame + $step]

        if {$p1 > $Jump} {
           set lulCurrentFrame [string trim [.gomautodock.left.frame3.frame32.first get]]
           set p1 $lulCurrentFrame
           eval lulAutoDock::GetAutoDockFrame $lulAutoDockTrajFileName $p1
           lulAutoDock::UpdateCurrentFrame;
           display
        } else {
           set lulCurrentFrame $p1
           eval lulAutoDock::GetAutoDockFrame $lulAutoDockTrajFileName $p1
           lulAutoDock::UpdateCurrentFrame;
           display
        }
        update idletasks
      }
  } 
# end of proc
}
############################################################################
# PROC
proc StopLooping { } {

     variable lulAutoDockTrajFileName
     variable lulCurrentFrame
     variable lulNumFrames
     variable lulLoopControl

     set lulLoopControl 0
# end of proc
}

############################################################################
# PROC
proc ImportFiles { Which w} {

     variable lulAutoDockLigandFileName
     variable lulAutoDockSystemFileName
     variable lulAutoDockTrajFileName
     variable lulCurrentFrame
     variable lulNumFrames
     variable lulLoopControl
     variable lulLigandAppend

     set FileName [$w get]

     if {$FileName == ""} {
        gomError "file name is not supplied"
        return
     }

     if {$Which == "ligand"} {

        if {$FileName != ""} {

        if {$lulLigandAppend} {
          eval import coord pdb $FileName append
		} else {
          eval import coord pdb $FileName
		}
        GetChargesFromPDBQfile $Which $FileName
        set lulAutoDockLigandFileName $FileName
        cd [file dirname $FileName]
        }

     } elseif {$Which == "system"} {

        if {$FileName != ""} {
        eval import coord pdb $FileName append
        GetChargesFromPDBQfile $Which $FileName
        set lulAutoDockSystemFileName $FileName
        cd [file dirname $FileName]
        }

     } elseif {$Which == "traj"} {

        if {$FileName != ""} {
        GetNumberOfAutoDockFrames $FileName
        set lulAutoDockTrajFileName $FileName
        cd [file dirname $FileName]

        set ww .gomautodock

        $ww.left.frame3.frame31.frames configure -state normal
        $ww.left.frame3.frame31.frames delete 0 end
        $ww.left.frame3.frame31.frames insert 0 $lulNumFrames
        $ww.left.frame3.frame31.frames configure -state disabled

        $ww.left.frame3.frame32.first delete 0  end
        $ww.left.frame3.frame32.first insert 0  1

        $ww.left.frame3.frame33.last delete 0  end
        $ww.left.frame3.frame33.last insert 0  $lulNumFrames

        $ww.left.frame3.frame34.step delete 0  end
        $ww.left.frame3.frame34.step insert 0  1

        $ww.left.frame3.frame35.current configure -state normal
        $ww.left.frame3.frame35.current delete 0  end
        $ww.left.frame3.frame35.current insert 0  1
        $ww.left.frame3.frame35.current configure -state disabled

        }
     }
# end of proc
}

############################################################################
# PROC
proc GetFrameCount {} {

     variable lulNumFrames

     return $lulNumFrames
#end of proc
}
############################################################################
# PROC
proc UpdateCurrentFrame { } {

     variable lulCurrentFrame

        set ww .gomautodock

        $ww.left.frame3.frame35.current configure -state normal
        $ww.left.frame3.frame35.current delete 0  end
        $ww.left.frame3.frame35.current insert 0  $lulCurrentFrame
        $ww.left.frame3.frame35.current configure -state disabled

# end of proc
}

############################################################################
# PROC
proc JumpToFrame { FrameNumber} {

     variable lulAutoDockTrajFileName

     GetAutoDockFrame $lulAutoDockTrajFileName $FrameNumber

#end of proc
}
###########################################################################
# PROC
proc GetChargesFromPDBQfile {Type FileName} {

     variable lulAutoDockLigandFileName
     variable lulAutoDockSystemFileName

   if {$Type == "ligand"} {
       set TotAtoms [show numatoms 1]
       set Struct  1
   } else {
       set TotAtoms [show numatoms 2]
       set Struct  2
   }
# open file ...
   set File_p [open $FileName r]
   if {$File_p == ""} {
       catch {gomError {ERROR: can't open AutoDock PDBQ file!}} errRet
       if {$errRet != ""} {
          gomPrint "ERROR: can't open AutoDock PDBQ file!"
          }
       error "ERROR: can't open AutoDock PDBQ file!"
       return 1
   }

# change to current directory
   lulChangeDirectory $FileName

   set Loop 1
   while {![eof $File_p]} {

     gets $File_p InputData
     
     if {[string match "ATOM*" $InputData]} {
#
       if {$Loop > $TotAtoms} {
          gomError "more atoms in file than defined for the structure"
          close $File_p
          return
       }
#
       set charge [string range $InputData 70 76]
       define atom charge $charge $Loop $Struct
       incr Loop
     }

   }

   close $File_p

# end of proc
}

###########################################################################
# PROC
proc MakePDBQCoordinates {Struct FileName Mask} {

   set def_title  "Default title for: PDBQ"

   if {[string index $FileName 0] == "."} {
       set FileName [$FileName get]
   }

   if {$FileName == ""} {
      gomError "file name is not supplied"
      return
   }

   set File_p [open $FileName w]

   if {$File_p == ""} {
       catch {lulErrorDialog {ERROR: can't open PDBQ output file for writing!}} errRet
       if {$errRet != ""} {
          gomPrint "ERROR: can't open PDBQ output file for writing!"
          }
       return
   }

# change to current directory
   lulChangeDirectory $FileName

   puts "Writing PDBQ coordinate file '$FileName' ..."

   set TotAtoms [show numatoms $Struct]

#   puts $File_p " $TotAtoms  $def_title"

   puts $File_p "REMARK ######################################"
   puts $File_p "REMARK      Default PDBQ file title"
   puts $File_p "REMARK ######################################"
   puts $File_p "ROOT"

   for {set i 1} {$i <= $TotAtoms } {incr i} {

# apply the display mask...
   if {[show atom displaystate $i $Struct] == 0 && $Mask} continue

   scan [show atom coordinates   $i $Struct] "%f %f %f" xc yc zc
   set  charge [show atom charge $i $Struct]
   set  vdw    [show atom vdw    $i $Struct]
   set  seg    [show atom segmen $i $Struct]
   set  resn   [show atom resnum $i $Struct]
   set  res    [show atom residu $i $Struct]
   set  atom   [show atom atomna $i $Struct]
   set  bval   0.0 ;#[show atom bvalue $i $Struct]

     switch [string index $atom 0] {
         {[1234567890]} {
         set OutPut [format "ATOM  %5d %-4.4s %-3.3s %1.1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f    %6.3f" \
           $i  $atom $res $seg $resn $xc $yc $zc "0.0" $bval $charge]
         puts $File_p $OutPut
         }
         default {
         set OutPut [format "ATOM  %5d  %-3.3s %-3.3s %1.1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f    %6.3f" \
           $i  $atom $res $seg $resn $xc $yc $zc "0.0" $bval $charge]
         puts $File_p $OutPut
         }
     }
   }

   puts $File_p "ENDROOT"
   close $File_p

   puts "Done!"

   return
# end of proc
}
###########################################################################
# PROC
proc LogFileCruncher { FileName } {

# check what version of AutoDock has been used
# open file ...
   if {[string trim $FileName] == ""} {return}

# input file
   set File_p [open $FileName r]
   if {$File_p == ""} {
       catch {gomError {ERROR: can't open AutoDock log file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open AutoDock log file!"
          }
       error "ERROR: can't open AutoDock log file!"
       return 0
   }

   set version 0.0
# What version of AutoDock
   while {![eof $File_p]} {

     gets $File_p InputData

     if {[string match "*AutoDock*" $InputData]} {

       set version [lindex $InputData 2]

	   break
     }
   }

   set type ""
# MONTE CARLO SIMULATED ANNEALING or GENETIC ALGORITHM-LOCAL SEARCH DOCKING
   seek $File_p 0 start
   while {![eof $File_p]} {

     gets $File_p InputData

     if {[string match "*BEGINNING MONTE CARLO SIMULATED ANNEALING*" $InputData]} {
       set type mc
       puts [string trim $InputData]
	   break
     } elseif {[string match "*BEGINNING GENETIC ALGORITHM-LOCAL SEARCH DOCKING*" $InputData]} {
       set type gn
       puts [string trim $InputData]
	   break
     }
   }

   close $File_p

   if {$version < 1.0} {
      gomError "Can't solve version number of used AutoDock"
      return
   }

   if {$version < 3.0} {
     gomPrint "This is an output from AutoDock version 2"
     LogFileCruncher2 "$FileName"
   } else {
     gomPrint "This is an output from AutoDock version 3"
     if {$type == "mc"} {
       LogFileCruncher3MC "$FileName"
     } elseif {$type == "gn"} {
       LogFileCruncher3 "$FileName"
     } else {
       gomError "unknown method in AutoDock log file (MC/Genetic)"
     }
   }

}
###########################################################################
# PROC
proc LogFileCruncher2 { FileName } {

# open file ...
   if {[string trim $FileName] == ""} {return}

# input file
   set File_p [open $FileName r]
   if {$File_p == ""} {
       catch {gomError {ERROR: can't open AutoDock log file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open AutoDock log file!"
          }
       error "ERROR: can't open AutoDock log file!"
       return 0
   }

# change to current directory
   lulChangeDirectory "$FileName"

   set OutFile1 "[file rootname $FileName]_INIT.xmol"

# output file (INIT)
   puts "Preparing INIT (xmol) trajectory file: $OutFile1"
   set File_po1 [open $OutFile1 w]
   if {$File_po1 == ""} {
       catch {gomError {ERROR: can't open INIT traj file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open INIT traj file!"
          }
       error "ERROR: can't open INIT traj file!"
       close $File_p
       return 0
   }

   set OutFile2 "[file rootname $FileName]_DOCKED.xmol"

# output file (DOCKED)
   gomPrint "Preparing DOCKED (xmol) trajectory file: $OutFile2"
   set File_po2 [open $OutFile2 w]
   if {$File_po2 == ""} {
       catch {gomError {ERROR: can't open DOCKED traj file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open DOCKED traj file!"
          }
       error "ERROR: can't open DOCKED traj file!"
       close $File_po1
       close $File_p
       return 0
   }

   set OutFile3 "[file rootname $FileName]_CLUSTER.xmol"

# output file (CLUSTER)
   gomPrint "Preparing CLUSTER (xmol) trajectory file: $OutFile3"
   set File_po3 [open $OutFile3 w]
   if {$File_po3 == ""} {
       catch {gomError {ERROR: can't open CLUSTER traj file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open CLUSTER traj file!"
          }
       error "ERROR: can't open CLUSTER traj file!"
       close $File_po2
       close $File_po1
       close $File_p
       return 0
   }

   set Runs 0

gomPrint "#################################################################"
gomPrint "#            A U T O D O C K (v. 2) log file cruncher           #"
gomPrint "#################################################################"
# scan first file for file info ...
   while {![eof $File_p]} {

     gets $File_p InputData

     if {$InputData == ""} continue
# file info ...
     if {[string match       "This file was created at:*"  $InputData]} {
        gomPrint $InputData
     } elseif {[string match "                   using:*"  $InputData]} {
        gomPrint $InputData
     } elseif {[string match "Number of atoms found in molecule =*" $InputData]} {
        gomPrint $InputData
        scan [string range $InputData [string length "Number of atoms found in molecule ="] end] "%d" NumAtoms
#        gomPrint "Number of atoms: $NumAtoms"
     } elseif {[string match "Total charge on molecule =*" $InputData]} {
        gomPrint $InputData
     } elseif {[string match "RUN *...*" $InputData]} {
        gomPrint "Found (INITIALIZING AUTOMATED DOCKING SIMULATION): $InputData"
        incr Runs
        
# open PDBQ file ...
        set OutFile "[file rootname $FileName]_INIT_$Runs.pdbq"
        set File_op [open $OutFile w]
          if {$File_op == ""} {
              catch {gomError {ERROR: can't open PDBQ output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDBQ output file!"
          }
        error "ERROR: can't open PDBQ output file!"
        return 0
        }
        
# open PDB file ...
        set OutFilePDB "[file rootname $FileName]_INIT_$Runs.pdb"
        set File_opPDB [open $OutFilePDB w]
          if {$File_opPDB == ""} {
              catch {gomError {ERROR: can't open PDB output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDB output file!"
          }
        error "ERROR: can't open PDB output file!"
        return 0
        }

# animations file
        puts $File_po1 " $NumAtoms"
        puts $File_po1 "INIT - Frame nr: $Runs"
         
# go and get the initialization coords (first some blanks!)
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gomPrint "Preparing: $OutFile..."
        while {![eof $File_p]} {

            gets $File_p InputData
            if {[string match "INITIAL-STATE:  *" $InputData]} {
              set Stuff    [string range $InputData 16 end]
              set StuffPDB [string range $InputData 16 70]
# temporary fix
              regsub {G[1-9][1-9]} $Stuff "C  " Stuff
              puts $File_op    $Stuff
              puts $File_opPDB $StuffPDB
              if {[string match "ATOM*" $Stuff]} {
                  set Atom   [lindex $Stuff 2]
                  set XC     [lindex $Stuff 5]
                  set YC     [lindex $Stuff 6]
                  set ZC     [lindex $Stuff 7]
                  set Charge [lindex $Stuff 10]
                  puts $File_po1 "$Atom $XC $YC $ZC $Charge"
              }
            } else {
              break
            }
        }

        close $File_op
        close $File_opPDB
        gomPrint "Done!"

     } elseif {[string match "Run Number *,*" $InputData]} {
        gomPrint "Found (FINAL DOCKED STATE): $InputData"
        
# open PDBQ file ...
        set OutFile "[file rootname $FileName]_DOCKED_$Runs.pdbq"
        set File_op [open $OutFile w]
          if {$File_op == ""} {
              catch {gomError {ERROR: can't open PDBQ output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDBQ output file!"
              }
           error "ERROR: can't open PDBQ output file!"
           return 0
          }

# open PDB file ...
        set OutFilePDB "[file rootname $FileName]_DOCKED_$Runs.pdb"
        set File_opPDB [open $OutFilePDB w]
          if {$File_opPDB == ""} {
              catch {gomError {ERROR: can't open PDB output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDB output file!"
              }
           error "ERROR: can't open PDB output file!"
           return 0
          }

# animations file
        puts $File_po2 " $NumAtoms"
        puts $File_po2 "DOCKED - Frame nr: $Runs"

# go and get the initialization coords (first some blanks!)
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gomPrint "Preparing: $OutFile..."
        while {![eof $File_p]} {

            gets $File_p InputData
            if {[string match "DOCKED: *" $InputData]} {
              set Stuff    [string range $InputData 8 end]
              set StuffPDB [string range $InputData 8 62]
# temporary fix
              regsub {G[1-9][1-9]} $Stuff "C  " Stuff
              puts $File_op    $Stuff
              puts $File_opPDB $StuffPDB
              if {[string match "ATOM*" $Stuff]} {
                  set Atom   [lindex $Stuff 2]
                  set XC     [lindex $Stuff 5]
                  set YC     [lindex $Stuff 6]
                  set ZC     [lindex $Stuff 7]
                  set Charge [lindex $Stuff 10]
                  puts $File_po2 "$Atom $XC $YC $ZC $Charge"
              }
            } else {
              break
            }
        }

        close $File_op
        close $File_opPDB
        gomPrint "Done!"

     } elseif {[string match "*LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER*" $InputData]} {
        gomPrint "Found (LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER)"

# go and get the initialization coords (first some blanks!)
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData

        set Clusters 1
        set OutFile "[file rootname $FileName]_CLUSTER_1.pdbq"        
        puts "Preparing file: $OutFile"
# open PDBQ file ...
        set File_op [open $OutFile w]
        if {$File_op == ""} {
            catch {gomError {ERROR: can't open PDBQ output file!}} errRet
            if {$errRet != ""} {
                puts "ERROR: can't open PDBQ output file!"
            }
            error "ERROR: can't open PDBQ output file!"
            return 0
        }

        set OutFilePDB "[file rootname $FileName]_CLUSTER_1.pdb"        
        puts "Preparing file: $OutFilePDB"
# open PDB file ...
        set File_opPDB [open $OutFilePDB w]
        if {$File_opPDB == ""} {
            catch {gomError {ERROR: can't open PDB output file!}} errRet
            if {$errRet != ""} {
                puts "ERROR: can't open PDB output file!"
            }
            error "ERROR: can't open PDB output file!"
            return 0
        }

        puts $File_po3 " $NumAtoms"
        puts $File_po3 "CLUSTER - nr: $Clusters"
        set FileIsOpen 1

        while {![eof $File_p]} {

            gets $File_p InputData

            if {[string match "END*" $InputData]} {

              puts  $File_op    $InputData
              puts  $File_opPDB $InputData
              close $File_op
              close $File_opPDB
              set FileIsOpen 0

            } elseif {[string trim $InputData] == ""} {
              break
            } else {

               if {!$FileIsOpen} {
                 incr Clusters
                 set OutFile "[file rootname $FileName]_CLUSTER_$Clusters.pdbq"
                 gomPrint "Preparing file: $OutFile"
# open PDBQ file ...
                 set File_op [open $OutFile w]
                 if {$File_op == ""} {
                     catch {gomError {ERROR: can't open PDBQ output file!}} errRet
                     if {$errRet != ""} {
                     puts "ERROR: can't open PDBQ output file!"
                     }
                     error "ERROR: can't open PDBQ output file!"
                     return 0
                 }

                 set OutFilePDB "[file rootname $FileName]_CLUSTER_$Clusters.pdb"
                 gomPrint "Preparing file: $OutFilePDB"
# open filePDB ...
                 set File_opPDB [open $OutFilePDB w]
                 if {$File_opPDB == ""} {
                     catch {gomError {ERROR: can't open PDB output file!}} errRet
                     if {$errRet != ""} {
                     puts "ERROR: can't open PDB output file!"
                     }
                     error "ERROR: can't open PDB output file!"
                     return 0
                 }

                 set FileIsOpen 1
                 puts $File_po3 " $NumAtoms"
                 puts $File_po3 "CLUSTER - nr: $Clusters"
               }
# temporary fix
               regsub {G[1-9][1-9]} $InputData "C  " InputData
               puts $File_op    $InputData
               puts $File_opPDB [string range $InputData 0 53]
               if {[string match "ATOM*" $InputData]} {
                   set Atom   [lindex $InputData 2]
                   set XC     [lindex $InputData 5]
                   set YC     [lindex $InputData 6]
                   set ZC     [lindex $InputData 7]
                   set Charge [lindex $InputData 10]
                   puts $File_po3 "$Atom $XC $YC $ZC $Charge"
               }
            }
        }     
     }
   }
# done ...

        close $File_po1
        close $File_po2
        close $File_po3
        close $File_p
#
        gomPrint "All is done! Number of INIT and DOCKED frames (*.pdbq files) is 2 * $Runs"
        gomPrint "All is done! Number of INIT and DOCKED frames (*.pdb  files) is 2 * $Runs"
        gomPrint "INIT    XMOL traj file: $OutFile1"
        gomPrint "DOCKED  XMOL traj file: $OutFile2"
        gomPrint "Cluster XMOL traj file: $OutFile3"
# end of proc
 }

###########################################################################
# PROC
proc LogFileCruncher3 { FileName } {

# open file ...
   if {[string trim $FileName] == ""} {return}

# input file
   set File_p [open $FileName r]
   if {$File_p == ""} {
       catch {gomError {ERROR: can't open AutoDock log file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open AutoDock log file!"
          }
       error "ERROR: can't open AutoDock log file!"
       return 0
   }

# change to current directory
   lulChangeDirectory "$FileName"

   set OutFile2 "[file rootname $FileName]_DOCKED.xmol"

# output file (DOCKED)
   gomPrint "Preparing DOCKED (xmol) trajectory file: $OutFile2"
   set File_po2 [open $OutFile2 w]
   if {$File_po2 == ""} {
       catch {gomError {ERROR: can't open DOCKED traj file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open DOCKED traj file!"
          }
       error "ERROR: can't open DOCKED traj file!"
       close $File_p
       return 0
   }

   set OutFile3 "[file rootname $FileName]_CLUSTER.xmol"

# output file (CLUSTER)
   gomPrint "Preparing CLUSTER (xmol) trajectory file: $OutFile3"
   set File_po3 [open $OutFile3 w]
   if {$File_po3 == ""} {
       catch {gomError {ERROR: can't open CLUSTER traj file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open CLUSTER traj file!"
          }
       error "ERROR: can't open CLUSTER traj file!"
       close $File_po2
       close $File_p
       return 0
   }

   set Runs 0

gomPrint "#################################################################"
gomPrint "#            A U T O D O C K (v.3) log file cruncher            #"
gomPrint "#################################################################"
# scan first file for file info ...
   while {![eof $File_p]} {

     gets $File_p InputData

     if {$InputData == ""} continue
# file info ...
     if {[string match       "This file was created at:*"  $InputData]} {
        gomPrint $InputData
     } elseif {[string match "                   using:*"  $InputData]} {
        gomPrint $InputData
     } elseif {[string match "Number of atoms found in molecule =*" $InputData]} {
        gomPrint $InputData
        scan [string range $InputData [string length "Number of atoms found in molecule ="] end] "%d" NumAtoms
#        gomPrint "Number of atoms: $NumAtoms"
     } elseif {[string match "Total charge on molecule =*" $InputData]} {
        gomPrint $InputData
     } elseif {[string match "INPUT PDBQ FILE:*" $InputData]} {

        gomPrint "Found $InputData"
        
# open PDBQ file ...
        set OutFile "[file rootname $FileName]_INIT.pdbq"
        set File_op [open $OutFile w]
          if {$File_op == ""} {
              catch {gomError {ERROR: can't open PDBQ output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDBQ output file!"
          }
        error "ERROR: can't open PDBQ output file!"
        return 0
        }
        
# open PDB file ...
        set OutFilePDB "[file rootname $FileName]_INIT.pdb"
        set File_opPDB [open $OutFilePDB w]
          if {$File_opPDB == ""} {
              catch {gomError {ERROR: can't open PDB output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDB output file!"
          }
        error "ERROR: can't open PDB output file!"
        return 0
        }
         
# go and get the initialization coords (first some blanks!)
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gomPrint "Preparing: $OutFile..."
        while {![eof $File_p]} {

            gets $File_p InputData
            if {[string match "INPUT-PDBQ:*" $InputData]} {
              set Stuff    [string range $InputData 12 end]
              set StuffPDB [string range $InputData 12 68]
# temporary fix
              regsub {G[1-9][1-9]} $Stuff "C  " Stuff
              puts $File_op    $Stuff
              puts $File_opPDB $StuffPDB
            } else {
              break
            }
        }

        close $File_op
        close $File_opPDB
        gomPrint "Done!"

     } elseif {[string match "Run:*" $InputData]} {

	    set Run1 [lindex $InputData 1]
		set Run2 [lindex $InputData 3]

        gomPrint "$InputData"
        
# open PDBQ file ...
        set OutFile "[file rootname $FileName]_DOCKED_$Run1-$Run2.pdbq"
        set File_op [open $OutFile w]
          if {$File_op == ""} {
              catch {gomError {ERROR: can't open PDBQ output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDBQ output file!"
              }
           error "ERROR: can't open PDBQ output file!"
           return 0
          }

# open PDB file ...
        set OutFilePDB "[file rootname $FileName]_DOCKED_$Run1-$Run2.pdb"
        set File_opPDB [open $OutFilePDB w]
          if {$File_opPDB == ""} {
              catch {gomError {ERROR: can't open PDB output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDB output file!"
              }
           error "ERROR: can't open PDB output file!"
           return 0
          }

# animations file
        puts $File_po2 " $NumAtoms"
        puts $File_po2 "DOCKED - Frame nr: $Run1 / $Run2"

# go and get the initialization coords (first some blanks!)

        while {![eof $File_p]} {

            gets $File_p InputData
            if {[string match "DOCKED: MODEL*" $InputData]} {
			  set thisrun [lindex $InputData 2]
			  if {$thisrun != $Run1} {
                gomPrint "Run: $Run1 != Model: $thisrun"
			  }
			  break
			}
        }

        gomPrint "Preparing: $OutFile..."

        while {![eof $File_p]} {
           gets $File_p InputData
           if {[string match "*x * y * z * vdW * Elec * q*" $InputData]} {
             break
		   }
		}

        gets $File_p InputDat

        while {![eof $File_p]} {

            gets $File_p InputData
            if {[string match "DOCKED: ATOM *" $InputData]} {
              set Stuff    [string range $InputData 8 end]
              set StuffPDB [string range $InputData 8 62]
# temporary fix
              regsub {G[1-9][1-9]} $Stuff "C  " Stuff
              puts $File_op    $Stuff
              puts $File_opPDB $StuffPDB
              if {[string match "ATOM*" $Stuff]} {
                  set Atom   [lindex $Stuff 2]
                  set XC     [lindex $Stuff 5]
                  set YC     [lindex $Stuff 6]
                  set ZC     [lindex $Stuff 7]
                  set Charge [lindex $Stuff 10]
                  puts $File_po2 "$Atom $XC $YC $ZC $Charge"
              }
            } else {
              break
            }
        }

        close $File_op
        close $File_opPDB
        gomPrint "Done!"

     } elseif {[string match "*LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER*" $InputData]} {
		gomPrint "Found (LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER)"

            gets $File_p InputData
            gets $File_p InputData
            gets $File_p InputData
            gets $File_p InputData
            gets $File_p InputData
            gets $File_p InputData

        while {![eof $File_p]} {

            gets $File_p InputData

            if {[string match "MODEL*" $InputData]} {

			    set model [lindex $InputData 1]
                gets $File_p InputData
			    set  run [lindex $InputData 3]
                gets $File_p InputData
			    set clusterrank [lindex $InputData 4]

                set OutFile "[file rootname $FileName]_CLUSTER_$model-$run-$clusterrank.pdbq"        
                gomPrint "Preparing file: $OutFile"
# open PDBQ file ...
                set File_op [open $OutFile w]
                if {$File_op == ""} {
                  catch {gomError {ERROR: can't open PDBQ output file!}} errRet
                  if {$errRet != ""} {
                     puts "ERROR: can't open PDBQ output file!"
                  }
                  error "ERROR: can't open PDBQ output file!"
                  return 0
                }

                set OutFilePDB "[file rootname $FileName]_CLUSTER_$model-$run-$clusterrank.pdb"        
                gomPrint "Preparing file: $OutFilePDB"
# open PDB file ...
                set File_opPDB [open $OutFilePDB w]
                if {$File_opPDB == ""} {
                  catch {gomError {ERROR: can't open PDB output file!}} errRet
                  if {$errRet != ""} {
                     puts "ERROR: can't open PDB output file!"
                  }
                  error "ERROR: can't open PDB output file!"
                  return 0
                }

                puts $File_po3 " $NumAtoms"
                puts $File_po3 "CLUSTER - $model / $run / $clusterrank"

				continue


            } elseif {[string match "ENDMDL*" $InputData]} {

              puts  $File_op    $InputData
              puts  $File_opPDB $InputData
              close $File_op
              close $File_opPDB

            } elseif {[string trim $InputData] == ""} {
              break
            } elseif {[string match "ATOM*" $InputData]} {


# temporary fix
               regsub {G[1-9][1-9]} $InputData "C  " InputData
               puts $File_op    $InputData
               puts $File_opPDB [string range $InputData 0 53]

                   set Atom   [lindex $InputData 2]
                   set XC     [lindex $InputData 5]
                   set YC     [lindex $InputData 6]
                   set ZC     [lindex $InputData 7]
                   set Charge [lindex $InputData 10]
                   puts $File_po3 "$Atom $XC $YC $ZC $Charge"
            }
        }     
     }
   }
# done ...

        close $File_po2
        close $File_po3
        close $File_p
#
        gomPrint "DOCKED  XMOL traj file: $OutFile2"
        gomPrint "Cluster XMOL traj file: $OutFile3"
# end of proc
 }
###########################################################################
# PROC
proc LogFileCruncher3MC { FileName } {

# open file ...
   if {[string trim $FileName] == ""} {return}

# input file
   set File_p [open $FileName r]
   if {$File_p == ""} {
       catch {gomError {ERROR: can't open AutoDock log file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open AutoDock log file!"
          }
       error "ERROR: can't open AutoDock log file!"
       return 0
   }

# change to current directory
   lulChangeDirectory "$FileName"

   set OutFile1 "[file rootname $FileName]_INIT.xmol"

# output file (INIT)
   gomPrint "Preparing INIT (xmol) trajectory file: $OutFile1"
   set File_po1 [open $OutFile1 w]
   if {$File_po1 == ""} {
       catch {gomError {ERROR: can't open INIT traj file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open INIT traj file!"
          }
       error "ERROR: can't open INIT traj file!"
       close $File_p
       return 0
   }

   set OutFile2 "[file rootname $FileName]_DOCKED.xmol"

# output file (DOCKED)
   gomPrint "Preparing DOCKED (xmol) trajectory file: $OutFile2"
   set File_po2 [open $OutFile2 w]
   if {$File_po2 == ""} {
       catch {gomError {ERROR: can't open DOCKED traj file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open DOCKED traj file!"
          }
       error "ERROR: can't open DOCKED traj file!"
       close $File_po1
       close $File_p
       return 0
   }

   set OutFile3 "[file rootname $FileName]_CLUSTER.xmol"

# output file (CLUSTER)
   gomPrint "Preparing CLUSTER (xmol) trajectory file: $OutFile3"
   set File_po3 [open $OutFile3 w]
   if {$File_po3 == ""} {
       catch {gomError {ERROR: can't open CLUSTER traj file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open CLUSTER traj file!"
          }
       error "ERROR: can't open CLUSTER traj file!"
       close $File_po2
       close $File_po1
       close $File_p
       return 0
   }

   set Runs 0

gomPrint "#################################################################"
gomPrint "#            A U T O D O C K (v. 3) log file cruncher           #"
gomPrint "#################################################################"
# scan first file for file info ...
   while {![eof $File_p]} {

     gets $File_p InputData

     if {$InputData == ""} continue
# file info ...
     if {[string match       "This file was created at:*"  $InputData]} {
        gomPrint $InputData
     } elseif {[string match "                   using:*"  $InputData]} {
        gomPrint $InputData
     } elseif {[string match "Number of atoms found in molecule =*" $InputData]} {
        gomPrint $InputData
        scan [string range $InputData [string length "Number of atoms found in molecule ="] end] "%d" NumAtoms
#        gomPrint "Number of atoms: $NumAtoms"
     } elseif {[string match "Total charge on molecule =*" $InputData]} {
        gomPrint $InputData
     } elseif {[string match "RUN *...*" $InputData]} {
        gomPrint "Found (INITIALIZING AUTOMATED DOCKING SIMULATION): $InputData"
        incr Runs
        
# open PDBQ file ...
        set OutFile "[file rootname $FileName]_INIT_$Runs.pdbq"
        set File_op [open $OutFile w]
          if {$File_op == ""} {
              catch {gomError {ERROR: can't open PDBQ output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDBQ output file!"
          }
        error "ERROR: can't open PDBQ output file!"
        return 0
        }
        
# open PDB file ...
        set OutFilePDB "[file rootname $FileName]_INIT_$Runs.pdb"
        set File_opPDB [open $OutFilePDB w]
          if {$File_opPDB == ""} {
              catch {gomError {ERROR: can't open PDB output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDB output file!"
          }
        error "ERROR: can't open PDB output file!"
        return 0
        }

# animations file
        puts $File_po1 " $NumAtoms"
        puts $File_po1 "INIT - Frame nr: $Runs"
         
# go and get the initialization coords (first some blanks!)
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gomPrint "Preparing: $OutFile..."
        while {![eof $File_p]} {

            gets $File_p InputData
            if {[string match "INITIAL-STATE:  *" $InputData]} {
              set Stuff    [string range $InputData 16 end]
              set StuffPDB [string range $InputData 16 70]
# temporary fix
              regsub {G[1-9][1-9]} $Stuff "C  " Stuff
              puts $File_op    $Stuff
              puts $File_opPDB $StuffPDB
              if {[string match "ATOM*" $Stuff]} {
                  set Atom   [lindex $Stuff 2]
                  set XC     [lindex $Stuff 5]
                  set YC     [lindex $Stuff 6]
                  set ZC     [lindex $Stuff 7]
                  set Charge [lindex $Stuff 10]
                  puts $File_po1 "$Atom $XC $YC $ZC $Charge"
              }
            } else {
              break
            }
        }

        close $File_op
        close $File_opPDB
        gomPrint "Done!"

     } elseif {[string match "Run Number *,*" $InputData]} {
        gomPrint "Found (FINAL DOCKED STATE): $InputData"
        
# open PDBQ file ...
        set OutFile "[file rootname $FileName]_DOCKED_$Runs.pdbq"
        set File_op [open $OutFile w]
          if {$File_op == ""} {
              catch {gomError {ERROR: can't open PDBQ output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDBQ output file!"
              }
           error "ERROR: can't open PDBQ output file!"
           return 0
          }

# open PDB file ...
        set OutFilePDB "[file rootname $FileName]_DOCKED_$Runs.pdb"
        set File_opPDB [open $OutFilePDB w]
          if {$File_opPDB == ""} {
              catch {gomError {ERROR: can't open PDB output file!}} errRet
              if {$errRet != ""} {
              puts "ERROR: can't open PDB output file!"
              }
           error "ERROR: can't open PDB output file!"
           return 0
          }

# animations file
        puts $File_po2 " $NumAtoms"
        puts $File_po2 "DOCKED - Frame nr: $Runs"

# go and get the initialization coords (first some blanks!)
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gets $File_p InputData
        gomPrint "Preparing: $OutFile..."
        while {![eof $File_p]} {

            gets $File_p InputData
            if {[string match "DOCKED: *" $InputData]} {
              set Stuff    [string range $InputData 8 end]
              set StuffPDB [string range $InputData 8 62]
# temporary fix
              regsub {G[1-9][1-9]} $Stuff "C  " Stuff
              puts $File_op    $Stuff
              puts $File_opPDB $StuffPDB
              if {[string match "ATOM*" $Stuff]} {
                  set Atom   [lindex $Stuff 2]
                  set XC     [lindex $Stuff 5]
                  set YC     [lindex $Stuff 6]
                  set ZC     [lindex $Stuff 7]
                  set Charge [lindex $Stuff 10]
                  puts $File_po2 "$Atom $XC $YC $ZC $Charge"
              }
            } else {
              break
            }
        }

        close $File_op
        close $File_opPDB
        gomPrint "Done!"

     } elseif {[string match "*LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER*" $InputData]} {
        gomPrint "Found (LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER)"

# go and get the initialization coords (first some blanks!)
            gets $File_p InputData
            gets $File_p InputData
            gets $File_p InputData
            gets $File_p InputData
            gets $File_p InputData
            gets $File_p InputData

        while {![eof $File_p]} {

            gets $File_p InputData

            if {[string match "MODEL*" $InputData]} {

			    set model [lindex $InputData 1]
                gets $File_p InputData
			    set  run [lindex $InputData 3]
                gets $File_p InputData
			    set clusterrank [lindex $InputData 4]

                set OutFile "[file rootname $FileName]_CLUSTER_$model-$run-$clusterrank.pdbq"        
                gomPrint "Preparing file: $OutFile"
# open PDBQ file ...
                set File_op [open $OutFile w]
                if {$File_op == ""} {
                  catch {gomError {ERROR: can't open PDBQ output file!}} errRet
                  if {$errRet != ""} {
                     puts "ERROR: can't open PDBQ output file!"
                  }
                  error "ERROR: can't open PDBQ output file!"
                  return 0
                }

                set OutFilePDB "[file rootname $FileName]_CLUSTER_$model-$run-$clusterrank.pdb"        
                gomPrint "Preparing file: $OutFilePDB"
# open PDB file ...
                set File_opPDB [open $OutFilePDB w]
                if {$File_opPDB == ""} {
                  catch {gomError {ERROR: can't open PDB output file!}} errRet
                  if {$errRet != ""} {
                     puts "ERROR: can't open PDB output file!"
                  }
                  error "ERROR: can't open PDB output file!"
                  return 0
                }

                puts $File_po3 " $NumAtoms"
                puts $File_po3 "CLUSTER - $model / $run / $clusterrank"

				continue


            } elseif {[string match "ENDMDL*" $InputData]} {

              puts  $File_op    $InputData
              puts  $File_opPDB $InputData
              close $File_op
              close $File_opPDB

            } elseif {[string trim $InputData] == ""} {
              break
            } elseif {[string match "ATOM*" $InputData]} {


# temporary fix
               regsub {G[1-9][1-9]} $InputData "C  " InputData
               puts $File_op    $InputData
               puts $File_opPDB [string range $InputData 0 53]

                   set Atom   [lindex $InputData 2]
                   set XC     [lindex $InputData 5]
                   set YC     [lindex $InputData 6]
                   set ZC     [lindex $InputData 7]
                   set Charge [lindex $InputData 10]
                   puts $File_po3 "$Atom $XC $YC $ZC $Charge"
            }
        }     
     }
   }
# done ...

        close $File_po1
        close $File_po2
        close $File_po3
        close $File_p
#
        gomPrint "All is done! Number of INIT and DOCKED frames (*.pdbq files) is 2 * $Runs"
        gomPrint "All is done! Number of INIT and DOCKED frames (*.pdb  files) is 2 * $Runs"
        gomPrint "INIT    XMOL traj file: $OutFile1"
        gomPrint "DOCKED  XMOL traj file: $OutFile2"
        gomPrint "Cluster XMOL traj file: $OutFile3"
# end of proc
 }

}


