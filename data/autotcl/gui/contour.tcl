#####################################
# contour display defaults
#####################################
#                         solid/mesh
set gomDefaultContourType solid
#                         value between 1.0 and 0.0
set gomDefaultContourOpaque 1.0
#                          1/0 (on/off)
set gomDefaultContourSmooth 0
#
# array gomContourColour(1...MAX)
set gomMaxContourLevels 10
for {set i 1} { $i <= $gomMaxContourLevels} {incr i} { 
     set gomContourColour($i)     #ffffff
}
##.......................##

set gomContourID            1
set gomContourExt          ".plt .PLT"
set gomContourFileName     ""
set gomContourSmoothState   0
set gomHelpFile(contour)             "contour_widget.html"
############################################################################
# PROC
proc lulContourControlStart {} {

  global lulContourControlStartup

  if {$lulContourControlStartup} { 
    set w .gomcontour
    if {[winfo exists $w]} {

       lulContourControl
    }
  }
}




############################################################################
# PROC
proc lulShowAlphaValue { w value } {

$w delete 0 end
$w insert 0 [format "%.3f" $value]
}
############################################################################
# PROC
proc lulApplyContourDetails { w Contour Level } {

     global gomContourDisplayState
	 global gomContourSmoothState
	 global gomContourType
     global gomContourCullFaceState

set  ContourName [show contour name $Contour]

# get opacity
set Alpha [$w.alpha.entry get]
contour alpha   $ContourName $Alpha    $Level
# if the value in input field is different from scale correct it!
$w.alpha.scale set $Alpha

contour display  $ContourName $gomContourDisplayState  $Level
contour smooth   $ContourName $gomContourSmoothState   $Level
contour type     $ContourName $gomContourType          $Level
contour cullface $ContourName $gomContourCullFaceState $Level

  lulInvalidateDisplay
}

############################################################################
proc lulApplyContourCommand { w } {

   global gomContourColour
   global gomMaxContourLevels
   global gomContourID
   global gomMappingIndex
   global gomContourDisplayMethod

   set CatString ""

# Fixed 1998-07-28/LUL to return right contour mapping name
# set MappingContourName [show contour name [lindex $gomMappingIndex 0]]
# New Fix 1999-01-10/LUL. Introduced new command "contour mapping Name1 Name2
# where the values from Name2 will be colour mapped on Name1

  set MappingContourName [lindex [show contour mapping $gomContourID] 0]

  if {$gomContourID == $MappingContourName} {

      for {set i 1} {$i <= $gomMaxContourLevels} {incr i} {

          set value [lindex [$w.left.textframe.text.frame$i.entry$i get] 0] 
          if {[string trim $value] == ""} break
          set colour [lulColourHex2Float $gomContourColour($i)]
          set CatString [concat $CatString $value "{" $colour "}"]
      }

# nothing supplied, reset the level counter ... 
      if {$i == 1 && $value == ""} {
          eval "contour plot $gomContourID"
	      lulInvalidateDisplay
	      return
      }

      if {[show contour defined]  > 0} {
          eval "contour plot $gomContourID $CatString"
          lulSelectContourDisplayMethod $gomContourDisplayMethod
          lulInvalidateDisplay
      }
# apply contour mapping...
  } else { 
      for {set i 1} {$i <= $gomMaxContourLevels} {incr i} {

          set tvalue [$w.left.textframe.text.frame$i.entry$i get] 
          if {[string trim $tvalue] == ""} break
          set value [lindex $tvalue 0]
          set min   [lindex $tvalue 1]
          set max   [lindex $tvalue 2]
          set CatString [concat $CatString $value "\{$min\}" "\{$max\}" ]
	  }

# nothing supplied, reset the level counter ... 
      if {$i == 1 && $tvalue == ""} {
          eval "contour plot $gomContourID"
          lulInvalidateDisplay
          return
      }

      if {[show contour defined]  > 0} {
          eval "contour plot { $gomContourID $MappingContourName } $CatString"
          lulSelectContourDisplayMethod $gomContourDisplayMethod
          lulInvalidateDisplay
      }
  }
}

############################################################################
# PROC
proc lulGetContourFile { w } {


    global gomContourExt
	global gomContourFileName
    global gomContourID

# reset contour file name 
    set gomContourFileName ""

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    
       set types {
         {"Contour files"		".plt .PLT"		       TEXT}
         {"All files"		*}
       }
    set file [string trim [tk_getOpenFile -filetypes $types -parent $w \
         -defaultextension $gomContourExt -title "Contour file"]]

    if [string compare $file ""] {
    $w delete 0 end
	$w insert 0 $file

    set gomContourFileName $file
    }
}

############################################################################
# PROC
proc lulReadContourFile { w } {

   global gomContourColour
   global gomContourFileName
   global gomContourID
   global gomMaxContourLevels
   global gomMappingIndex
   global gomCombineIndex
   global lulTraceTransformationValue
   global gomStructureNames

   set gomContourFileName [$w.left.frame1.filename get]
    
   if {$gomContourFileName != ""} {

      set  NumContours [show contour defined]
      incr NumContours 
      contour file "$gomContourFileName" $NumContours

      set Name [file tail $gomContourFileName]
	  set Min  [format "%.3f" [show contour minimum $NumContours]]
	  set Max  [format "%.3f" [show contour maximum $NumContours]]
      set base $w.left.contourframe.sub$NumContours

      frame       $base
      button      $base.contmap$NumContours -text "Mapping" \
                   -command "gomContourMappingControl $NumContours"
	  radiobutton $base.contour$NumContours               \
         -text "($NumContours) $Name Min: $Min Max: $Max" \
		 -value $NumContours                              \
		 -variable gomContourID                           \
		 -command "lulUpdateContourLevels $w"

      button      $base.combine$NumContours -text "Glue ==> $gomStructureNames(1)" \
                  -command "gomContourCombineControl $NumContours $base.combine$NumContours"

      set gomCombineIndex($NumContours)   1

      pack $base                      -side top  -anchor w -fill x
      pack $base.contmap$NumContours  -side left -anchor w
      pack $base.contour$NumContours  -side left -anchor w
      pack $base.combine$NumContours  -side right -anchor e

#      if {[show transformation type] == "global"} {
#          $base.combine$NumContours configure -state disabled
#      }

      $base.contour$NumContours select

      set gomMappingIndex "$NumContours : No mapping"

      for {set i 1} {$i <= $gomMaxContourLevels} {incr i} {

	 set gomContourColour($i) #ffffff
         $w.left.textframe.text.frame$i.button$i configure -bg white -fg black
         $w.left.textframe.text.frame$i.entry$i  delete 0 end

      }
         $w.left.textframe.text.frame1.entry1  icursor 0

# change to current directory
         lulChangeDirectory "$gomContourFileName"
   }

}

############################################################################
# PROC
proc lulContourSmoothing { w entry } {

     global gomContourSmoothState

     if {![show contour defined]} {
        $w deselect
        lulErrorDialog {ERROR: no contour(s) defined so far. Read first at least one contour file!}
        return
	 }

     set Levels [show contour levels $entry]

	 if {$Levels} {
	    if {$gomContourSmoothState} {
            contour smooth  $entry on
	    } else {
            contour smooth  $entry off
		}
	 } else {
        $w deselect
        lulErrorDialog {ERROR: no contour levels defined so far. Press the 'Apply' button first to define contour levels!}
	 }
}

############################################################################
# PROC
proc lulContourDisplayState { w entry } {

     global gomContourDisplayState
     global gomContourID

     if {![show contour defined]} {
        $w deselect
        lulErrorDialog {ERROR: no contour(s) defined so far. Read first at least one contour file!}
        return
	 }

     set Levels [show contour levels $entry]

	 if {$Levels} {
	   if {[show contour display $state]} {
           contour display  $entry on
	   } else {
           contour display  $entry off
	   }
	    
	 } else {
        $w deselect
        lulErrorDialog {ERROR: no contour levels defined so far. Press the 'Apply' button first to define contour levels!}
	 }
}

############################################################################
# PROC
proc lulFillContourWidget { w } {

  global gomContourColour
  global gomContourSmoothState
  global gomContourID
  global gomMappingIndex
  global gomStructureNames

  set NumContours [show contour defined]
  set gomMappingIndex [lindex [show contour mapping $gomContourID] 0]

  if {$NumContours  > 0} {

      for { set i 0} { $i < $NumContours} {incr i} {

        set j    [expr $i + 1]
        set Name [file tail [show contour filename $j]]
        set Min  [show contour minimum  $j]
	    set Max  [show contour maximum  $j]
        set base $w.left.contourframe.sub$j

      frame       $base
      button      $base.contmap$j -text "Mapping" \
                   -command "gomContourMappingControl $j"
	  radiobutton $base.contour$j \
         -text     "($j) $Name Min: $Min Max: $Max" \
		 -value    [show contour name $j]           \
		 -variable  gomContourID                    \
		 -command  "lulUpdateContourLevels $w"

      set Glue [show contour structure $j]
      button      $base.combine$j -text "Glue ==> $gomStructureNames($Glue)" \
                   -command "gomContourCombineControl $j $base.combine$j"

      pack $base            -side top  -anchor w -fill x
      pack $base.contmap$j  -side left -anchor w
      pack $base.contour$j  -side left -anchor w
      pack $base.combine$j  -side right -anchor e

#      if {[show transformation type] == "global"} {
#          $base.combine$j configure -state disabled
#      }

   }

# select always the first one ...
      set gomContourID [show contour name 1]
#      $base.contour2 select
      set CValue [show contour mapping [show contour name 1]]

      if {$gomContourID == [lindex $CValue 0]} {

       set Levels [show contour levels $gomContourID]

       for {set i 1} {$i <= $Levels} {incr i} {
        set ContVal [show contour values $gomContourID $i]
        scan $ContVal "%e %f %f %f" Value Red Green Blue

        $w.left.textframe.text.frame$i.entry$i  delete 0 end
	    $w.left.textframe.text.frame$i.entry$i  insert 0 $Value
        $w.left.textframe.text.frame$i.button$i configure \
	       -bg [lulColourFloat2Hex $Red $Green $Blue]
# corrected colour 1998-07-28/LUL to update correctly when clicked
        set gomContourColour($i) [lulColourFloat2Hex $Red $Green $Blue]
	    }

      } else {

        set Levels [show contour levels $gomContourID]

        set Red   1.0
        set Green 1.0
        set Blue  1.0

        for {set i 1} {$i <= $Levels} {incr i} {
          set ContVal [show contour values $gomContourID $i]
	      scan $ContVal "%e" Value
          set min [format "%.4e" [lindex $CValue [expr 2 * $i    ]]]
          set max [format "%.4e" [lindex $CValue [expr 2 * $i + 1]]]
          set Value "$Value  $min  $max"

          $w.left.textframe.text.frame$i.entry$i  delete 0 end
	      $w.left.textframe.text.frame$i.entry$i  insert 0 $Value
          $w.left.textframe.text.frame$i.button$i configure \
	         -bg [lulColourFloat2Hex $Red $Green $Blue]
# corrected colour 1998-07-28/LUL to update correctly when clicked
        set gomContourColour($i) [lulColourFloat2Hex $Red $Green $Blue]
	    }
      }

      $w.left.textframe.text.frame1.entry1  icursor 0

  }
}
############################################################################


############################################################################
# PROC
proc lulUpdateContourLevels { w } {

     global gomContourID
     global gomMaxContourLevels

# blank input field and reset button colour

     for {set i 1} {$i <= $gomMaxContourLevels} {incr i} {

       $w.left.textframe.text.frame$i.button$i configure -bg white -fg black
       $w.left.textframe.text.frame$i.entry$i  delete 0 end

     }

    set CValue [show contour mapping $gomContourID]

    if {$gomContourID == [lindex $CValue 0]} {

      set Levels [show contour levels $gomContourID]

      for {set i 1} {$i <= $Levels} {incr i} {
        set ContVal [show contour values $gomContourID $i]
	    scan $ContVal "%e %f %f %f" Value Red Green Blue

        $w.left.textframe.text.frame$i.entry$i  delete 0 end
	    $w.left.textframe.text.frame$i.entry$i  insert 0 $Value
        $w.left.textframe.text.frame$i.button$i configure \
	       -bg [lulColourFloat2Hex $Red $Green $Blue]
	  }

    } else {
    
      set Levels [show contour levels $gomContourID]

      set Red   1.0
      set Green 1.0
      set Blue  1.0

      for {set i 1} {$i <= $Levels} {incr i} {
        set ContVal [show contour values $gomContourID $i]
	    scan $ContVal "%e" Value
        set min [format "%.4e" [lindex $CValue [expr 2 * $i    ]]]
        set max [format "%.4e" [lindex $CValue [expr 2 * $i + 1]]]
        set Value "$Value  $min  $max"

        $w.left.textframe.text.frame$i.entry$i  delete 0 end
	    $w.left.textframe.text.frame$i.entry$i  insert 0 $Value
        $w.left.textframe.text.frame$i.button$i configure \
	       -bg [lulColourFloat2Hex $Red $Green $Blue] \
	       -fg [lulGetVisibleForegroundColour [
	           lulColourFloat2Hex $Red $Green $Blue]]		
	  }
    } 
        $w.left.textframe.text.frame1.entry1  icursor 0

}
############################################################################
# PROC
proc gomContourCombineControl { Which ww} {

     global gomControlFont
     global gomHelpFile
     global gomCombineIndex
     global gomStructureNames

# return if no contours defined
    if {[show contour defined]   < 1} {
        lulErrorDialog {ERROR: no contour available. Read a contour first!}
        return
	}

      if {[show transformation type] == "global"} {
          lulErrorDialog {This button has no meaning for the 'global' transformation mode}
          return
      }

set w .gomcontourcombine
catch {destroy $w}
toplevel $w 
wm title $w "Contour Structure Combine"
wm iconname $w "Contour Combine"
#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss  -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply    -text Apply    -font "$gomControlFont" \
        -command "lulCombineContour2Structure $Which $ww"
button $w.buttons.help     -text Help     -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(contglueing)"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Contour glue structure:" -font "$gomControlFont"
pack   $w.label -side top -anchor w

set CurrentStructure [show contour structure $Which]

set temp [file tail [show contour filename $Which]]
frame  $w.frame0 
pack   $w.frame0 -side top -anchor w
label  $w.frame0.text -font {Arial 8 bold} -text "Glue contour file '$temp' to structure:"
pack   $w.frame0.text -side top -anchor w -padx 4 -pady 4

frame  $w.frame 
pack   $w.frame -side top -anchor w

set SumString ""

for {set i 1} {$i <= [show molstructures]} {incr i} {
#  set temp [file tail [show contour filename $i]]
  lappend SumString "$i : $gomStructureNames($i)"
}

eval tk_optionMenu $w.frame.options gomCombineIndex($Which) $SumString
pack   $w.frame.options -side left -anchor w
}

############################################################################
# PROC
proc lulCombineContour2Structure { Which w } {

     global gomCombineIndex
     global gomStructureNames

     contour combine $Which [lindex $gomCombineIndex($Which) 0]

     $w configure -text "Glue ==> $gomStructureNames([lindex $gomCombineIndex($Which) 0])"
}

############################################################################
# PROC
proc gomContourMappingControl { Which } {

     global gomControlFont
     global gomHelpFile
     global gomMappingIndex
     global gomDefaultMappingIndex

# return if no contours defined
    if {[show contour defined]   < 1} {
        lulErrorDialog {ERROR: no contour available. Read a contour first!}
        return
	}

set w .gomcontourmapping
catch {destroy $w}
toplevel $w 
wm title $w "Contour Mapping Control"
wm iconname $w "Contour Mapping Control"
#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "gomSetDefaultMappingIndex;destroy $w"
button $w.buttons.apply   -text Accept    -font "$gomControlFont" \
        -command "gomApplyMapping $Which;destroy $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(contmapping)"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -text "Contour mapping:" -font "$gomControlFont"
pack   $w.label -side top -anchor w

set FileTail [file tail [show contour filename $Which]]

frame  $w.frame0 
pack   $w.frame0 -side top -anchor w
label  $w.frame0.text -font {Arial 8 bold} -text "Mapping selected surface (below) on surface: \"$FileTail\""
pack   $w.frame0.text -side top -anchor w -padx 4 -pady 4

frame  $w.frame 
pack   $w.frame -side top -anchor w

set SumString ""

for {set i 1} {$i <= [show contour defined]} {incr i} {
 if {$Which == $i} {
  lappend SumString "$i : No mapping"
 } else {
  set temp [file tail [show contour filename $i]]
  lappend SumString "$i : $temp"
 }
}

  if {[show contour name $Which] == [lindex [show contour mapping $Which] 0]} {
      set gomMappingIndex "$Which : No mapping"
  } else {
# corrected 1998-07-26/LUL to display correct mapping
      set gomMappingIndex "$Which : [file tail [show contour filename [lindex [show contour mapping $Which] 0]]]"
  } 

  set gomDefaultMappingIndex $gomMappingIndex

eval tk_optionMenu $w.frame.options gomMappingIndex $SumString
pack   $w.frame.options -side left -anchor w
}

############################################################################
# PROC
proc gomSetDefaultMappingIndex { } {
     global gomDefaultMappingIndex
     global gomMappingIndex

     set gomMappingIndex $gomDefaultMappingIndex
}
############################################################################
# PROC
proc gomApplyMapping { Which } {
     global gomMappingIndex

contour mapping $Which [lindex $gomMappingIndex 0]

}
############################################################################
# PROC
proc lulFillMappingDefault {w Which} {

     global gomMappingIndex

     if {$Which == [lindex $gomMappingIndex 0]} {
         gomError "Can't map on itself"
         return
         }

     set Min [show contour minimum $Which]
     set Max [show contour maximum $Which]

$w.frame.min delete 0 end
$w.frame.min insert 0 $Min

$w.frame.max delete 0 end
$w.frame.max insert 0 $Max

}
############################################################################
# PROC
proc lulContourPlaneControl {} {

     global gomControlFont
     global gomContourColour
     global gomMaxContourLevels
     global gomContourSmoothState
     global gomHelpFile
     global gomContourID
     global gomCutX
     global gomCutY
     global gomCutZ
     global gomXYZ1
     global gomXYZ2
     global gomXYZ3

# return if no contours defined
    if {[show contour defined]   < 1} {
        lulErrorDialog {ERROR: no contour available. Read a contour first!}
        return
	}

set w .gomcontourplane
catch {destroy $w}
toplevel $w 
wm title $w "Contour Plane Control"
wm iconname $w "Contour Plane Control"
# if this window is closed the ".gomanimatecontplane" window
# has to be closed as well!
wm protocol .gomcontourplane WM_DELETE_WINDOW  {catch {destroy .gomcontourplane; destroy .gomanimatecontplane; destroy .gomcontourplanexyz}}

#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w;catch {destroy .gomanimatecontplane}"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyContourPlaneCommand $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(contourplane)"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

#
    scan [show contour cube $gomContourID] "%f %f %f %f %f %f %d %d %d" \
	      Minx Maxx Miny Maxy Minz Maxz Xdim Ydim Zdim

    frame $w.frame0
	pack  $w.frame0 -side top -anchor w

    frame $w.frame0.left -borderwidth 2 -relief ridge
    label $w.frame0.left.label0 -text "Box\[cube\] dimensions: "
    label $w.frame0.left.label1 -text "Min x: [format %6.3e $Minx], \
            max x: [format %6.3e $Maxx]"
    label $w.frame0.left.label2 -text "Min y: [format %6.3e $Miny], \
            max y: [format %6.3e $Maxy]"
    label $w.frame0.left.label3 -text "Min z: [format %6.3e $Minz], \
            max z: [format %6.3e $Maxz]"
   	pack  $w.frame0.left -side left       -anchor w
	pack  $w.frame0.left.label0 -side top -anchor w
	pack  $w.frame0.left.label1 -side top
	pack  $w.frame0.left.label2 -side top
	pack  $w.frame0.left.label3 -side top
	 
    frame  $w.frame0.right -borderwidth 2 -relief sunken
	button $w.frame0.right.delete -text "Delete all contours" \
	       -command "contour delete;destroy $w;lulInvalidateDisplay"
	pack   $w.frame0.right -side left -padx 5
	pack   $w.frame0.right.delete -side left 

    frame $w.box -borderwidth 2 -relief raised

	pack $w.box   -side top   -anchor n -fill both

    frame       $w.box.frame1 
	label       $w.box.frame1.labelplane   -text "Plane" -width 10
	label       $w.box.frame1.labelcoord   -text "\[X/Y/Z\]-coord"
	label       $w.box.frame1.labelmin     -text "Min" -width 10
	label       $w.box.frame1.labelmax     -text "Max" -width 10
    label       $w.box.frame1.labelanimate -text "Animate" -width 18
	label       $w.box.frame1.labelsmooth  -text "Smoothing"
	label       $w.box.frame1.displaytype  -text "Plot type" -width 12

	pack        $w.box.frame1              -side top    -anchor w
	pack        $w.box.frame1.labelplane   -side left   -anchor w
	pack        $w.box.frame1.labelcoord   -side left   -anchor w
	pack        $w.box.frame1.labelmin     -side left   -anchor w
	pack        $w.box.frame1.labelmax     -side left   -anchor w
    pack        $w.box.frame1.labelanimate -side left   -anchor w
	pack        $w.box.frame1.labelsmooth  -side left   -anchor w
	pack        $w.box.frame1.displaytype  -side left   -anchor w

# 1
	frame       $w.box.frame2
	checkbutton $w.box.frame2.yz -text "YZ-plane" -width 10 \
	            -variable gomYZ -command {plot -cutplane x;lulInvalidateDisplay}
	entry       $w.box.frame2.x         -width 10
	entry       $w.box.frame2.min       -width 10
	entry       $w.box.frame2.max       -width 10
    button      $w.box.frame2.animatex  -width 10 \
                 -text "X - ---> +" -command "lulAnimateCutplane x "
	menubutton  $w.box.frame2.smooth -text "Smooth type" \
	            -menu $w.box.frame2.smooth.menu          \
				-borderwidth 2 -relief raised
	menu  $w.box.frame2.smooth.menu
	set m $w.box.frame2.smooth.menu
	$m add radiobutton -label "No smoothing" -variable gomCutX -value 0 
    $m add radiobutton -label "Sqrt"         -variable gomCutX -value 1 
    $m add radiobutton -label "Log10"        -variable gomCutX -value 2 
    button      $w.box.frame2.plottype       \
                 -text "Type \[2D/3D\]" -command "lulSetCutplanePlotType  x $w"

	pack        $w.box.frame2          -side top
	pack        $w.box.frame2.yz       -side left
	pack        $w.box.frame2.x        -side left
	pack        $w.box.frame2.min      -side left
	pack        $w.box.frame2.max      -side left
    pack        $w.box.frame2.animatex -side left -padx 4
	pack        $w.box.frame2.smooth   -side left
    pack        $w.box.frame2.plottype -side left -padx 4

	$w.box.frame2.x delete 0 end
	$w.box.frame2.x inser  0 [format "%.5f" [expr ($Maxx + $Minx)/2.0]]

    $w.box.frame2.smooth.menu invoke 1

    $w.box.frame2.min delete 0 end
    $w.box.frame2.min inser  0 [show contour minimum $gomContourID]

    $w.box.frame2.max delete 0 end
    $w.box.frame2.max inser  0 [show contour maximum $gomContourID]

# 2
	frame       $w.box.frame3
	checkbutton $w.box.frame3.xz -text "XZ-plane" -width 10 \
	            -variable gomXZ -command {plot -cutplane y;lulInvalidateDisplay}
	entry       $w.box.frame3.y  -width 10
	entry       $w.box.frame3.min -width 10
	entry       $w.box.frame3.max -width 10
    button      $w.box.frame3.animatey  -width 10        \
                 -text "Y - ---> +" -command "lulAnimateCutplane y "
	menubutton  $w.box.frame3.smooth -text "Smooth type" \
	            -menu $w.box.frame3.smooth.menu          \
				-borderwidth 2 -relief raised
	menu  $w.box.frame3.smooth.menu
	set m $w.box.frame3.smooth.menu
	$m add radiobutton -label "No smoothing" -variable gomCutY -value 0
    $m add radiobutton -label "Sqrt"         -variable gomCutY -value 1
    $m add radiobutton -label "Log10"        -variable gomCutY -value 2
    button      $w.box.frame3.plottype       \
                 -text "Type \[2D/3D\]" -command "lulSetCutplanePlotType  y $w"

	pack        $w.box.frame3           -side top
	pack        $w.box.frame3.xz        -side left
	pack        $w.box.frame3.y         -side left
	pack        $w.box.frame3.min       -side left
	pack        $w.box.frame3.max       -side left
    pack        $w.box.frame3.animatey  -side left -padx 4
    pack        $w.box.frame3.smooth    -side left
    pack        $w.box.frame3.plottype  -side left -padx 4

	$w.box.frame3.y delete 0 end
	$w.box.frame3.y inser  0 [format "%.5f" [expr ($Maxy + $Miny)/2.0]]

    $w.box.frame3.smooth.menu invoke 1

    $w.box.frame3.min delete 0 end
    $w.box.frame3.min inser  0 [show contour minimum $gomContourID]

    $w.box.frame3.max delete 0 end
    $w.box.frame3.max inser  0 [show contour maximum $gomContourID]

# 3
	frame       $w.box.frame4
	checkbutton $w.box.frame4.xy -text "XY-plane" -width 10 \
	            -variable gomXY -command {plot -cutplane z;lulInvalidateDisplay}
	entry       $w.box.frame4.z         -width 10
	entry       $w.box.frame4.min       -width 10
	entry       $w.box.frame4.max       -width 10
    button      $w.box.frame4.animatez  -width 10        \
                 -text "Z - ---> +" -command "lulAnimateCutplane z "
	menubutton  $w.box.frame4.smooth -text "Smooth type" \
	            -menu $w.box.frame4.smooth.menu          \
				-borderwidth 2 -relief raised
	menu  $w.box.frame4.smooth.menu
	set m $w.box.frame4.smooth.menu
	$m add radiobutton -label "No smoothing"   -variable gomCutZ -value 0
    $m add radiobutton -label "Sqrt"           -variable gomCutZ -value 1
    $m add radiobutton -label "Log10"          -variable gomCutZ -value 2
    button      $w.box.frame4.plottype         \
                 -text "Type \[2D/3D\]" -command "lulSetCutplanePlotType  z $w"

	pack        $w.box.frame4 -side top
	pack        $w.box.frame4.xy       -side left
	pack        $w.box.frame4.z        -side left
	pack        $w.box.frame4.min      -side left
	pack        $w.box.frame4.max      -side left
    pack        $w.box.frame4.animatez -side left -padx 4
    pack        $w.box.frame4.smooth   -side left
    pack        $w.box.frame4.plottype -side left -padx 4

	$w.box.frame4.z delete 0 end
	$w.box.frame4.z inser  0 [format "%.5f" [expr ($Maxz + $Minz)/2.0]]

    $w.box.frame4.smooth.menu invoke 1

    $w.box.frame4.min delete 0 end
    $w.box.frame4.min inser  0 [show contour minimum $gomContourID]

    $w.box.frame4.max delete 0 end
    $w.box.frame4.max inser  0 [show contour maximum $gomContourID]

    frame  $w.plane   -borderwidth 2 -relief raised
	pack   $w.plane   -side top   -anchor w -pady 5

    label  $w.plane.label1  -text "Define XYZ plane #1:"
    button $w.plane.button1 -text "Click here" \
            -command "lulPlotContourPlaneXYZ 1 $w"
    pack   $w.plane.label1  -side left -anchor w
    pack   $w.plane.button1 -side left -anchor w -padx 5

    label  $w.plane.label2  -text "Define XYZ plane #2:"
    button $w.plane.button2 -text "Click here" \
            -command "lulPlotContourPlaneXYZ 2 $w"
    pack   $w.plane.label2  -side left -anchor w
    pack   $w.plane.button2 -side left -anchor w -padx 5

    label  $w.plane.label3  -text "Define XYZ plane #3:"
    button $w.plane.button3 -text "Click here" \
            -command "lulPlotContourPlaneXYZ 3 $w"
    pack   $w.plane.label3  -side left -anchor w
    pack   $w.plane.button3 -side left -anchor w -padx 5

    frame $w.profile   -borderwidth 2 -relief raised
	pack  $w.profile   -side top   -anchor w -pady 5

    label  $w.profile.label -text "Peek contour profile:"
    button $w.profile.button -text "Click here" \
            -command "lulPlotContourProfile $w"
    label  $w.profile.label1 -text "Number of bins:"
    entry  $w.profile.nbins  -width 6
    label  $w.profile.label2 -text "Cut value: "
    entry  $w.profile.cutv   -width 8

    pack   $w.profile.label  -side left -anchor w
    pack   $w.profile.button -side left -anchor w -padx 5
    pack   $w.profile.label1 -side left -anchor w
    pack   $w.profile.nbins  -side left -anchor w -pady 5
    pack   $w.profile.label2 -side left -anchor w 
    pack   $w.profile.cutv   -side left -anchor w -pady 5

    $w.profile.nbins insert 0  "41"
    $w.profile.cutv  insert 0  "1.0"

bind    $w.profile.nbins <Return> {lulPlotContourProfile .gomcontourplane}
bind    $w.profile.cutv  <Return> {lulPlotContourProfile .gomcontourplane}
    
}

############################################################################
# PROC
proc lulPlotContourPlaneXYZ { Which w } {

   global gomXYZ1
   global gomXYZ2
   global gomXYZ3
   global gomContourID
   global gomControlFont
   global gomContourColour
   global gomHelpFile

set w .gomcontourplanexyz
catch {destroy $w}
toplevel $w 
wm title $w "Contour XYZ Plane Control"
wm iconname $w "Contour XYZ Plane Control"
#
#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w;catch {destroy .gomanimatecontplane}"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyCutplaneXYZ $Which $w;lulInvalidateDisplay"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(contourxyzplane)"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Define Contour Plane #$Which:" -font $gomControlFont
pack   $w.label  -side top -anchor w -pady 4

frame  $w.inplabel
pack   $w.inplabel -side top
label  $w.inplabel.text  -text " (ON/OFF)                    (x, y, z)                     (Min/Max)"
pack   $w.inplabel.text  -side left

frame  $w.frame1 -borderwidth 2 -relief ridge
pack   $w.frame1 -side top -anchor w

	checkbutton $w.frame1.xyz -text "XYZ$Which-plane" -width 10 \
	            -variable gomXYZ$Which
	pack        $w.frame1.xyz                  -side left

    frame       $w.frame1.coord
	pack        $w.frame1.coord                -side left

    frame       $w.frame1.scale
	pack        $w.frame1.scale                -side left

	entry       $w.frame1.coord.p1xyz          -width 25
	entry       $w.frame1.coord.p2xyz          -width 25
	entry       $w.frame1.coord.p3xyz          -width 25
	pack        $w.frame1.coord.p1xyz          -side top
	pack        $w.frame1.coord.p2xyz          -side top
	pack        $w.frame1.coord.p3xyz          -side top

	entry       $w.frame1.scale.scmin          -width 10
	entry       $w.frame1.scale.scmax          -width 10
	pack        $w.frame1.scale.scmin          -side top
	pack        $w.frame1.scale.scmax          -side top

}

############################################################################
# PROC
proc lulApplyCutplaneXYZ { Which w } {

   global gomXYZ1
   global gomXYZ2
   global gomXYZ3
   global gomContourID

          set XYZ1 [string trim [$w.frame1.coord.p1xyz get]]
          if {[string index $XYZ1 0] == "\["} {
              eval set XYZ1 $XYZ1
          }
          set XYZ2 [string trim [$w.frame1.coord.p2xyz get]]
          if {[string index $XYZ2 0] == "\["} {
              eval set XYZ2 $XYZ2
          }
          set XYZ3 [string trim [$w.frame1.coord.p3xyz get]]
          if {[string index $XYZ3 0] == "\["} {
              eval set XYZ3 $XYZ3
          }

          set ScaleMin [string trim [$w.frame1.scale.scmin get]]
          set ScaleMax [string trim [$w.frame1.scale.scmax get]]

   switch $Which {

     1 { 
        if {$gomXYZ1} {
          eval plot cutplane xyz1 $gomContourID $XYZ1 $XYZ2 $XYZ3 $ScaleMin $ScaleMax
        } else {
          plot -cutplane xyz1
        }
     }
     2 { 
        if {$gomXYZ2} {
          eval plot cutplane xyz2 $gomContourID $XYZ1 $XYZ2 $XYZ3 $ScaleMin $ScaleMax
        } else {
          plot -cutplane xyz2
        }
     }
     3 { 
        if {$gomXYZ3} {
          eval plot cutplane xyz3 $gomContourID $XYZ1 $XYZ2 $XYZ3 $ScaleMin $ScaleMax
        } else {
          plot -cutplane xyz3
        }
     }
   }

}
############################################################################
# PROC
proc lulSetCutplanePlotType { Axis w} {

     global gomControlFont
     global gomContourColour
     global gomCutplaneType
     global gomHelpFile

set w .gomsetcutplaneplottype
catch {destroy $w}
toplevel $w 
wm title $w "Cutplane plot type"
wm iconname $w "Cutplane plot type"

#
frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
        -command "lulApplyCutplaneType $Axis $w"
button $w.buttons.help    -text Help     -font "$gomControlFont" \
        -command \
       "htmlShowHelp $gomHelpFile(cutplane_plottype)"

pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label  -text "Define cutplane plot type and damping factor:" -font $gomControlFont
pack   $w.label  -side top -anchor w -pady 4

frame  $w.frame1 -borderwidth 2 -relief ridge
pack   $w.frame1 -side top -anchor w

# display type
    set type  [show cutplane type]

if {$Axis == "x"} {

label  $w.frame1.label -text "Define plot type for X-direction \[2D/3D\]"
pack   $w.frame1.label -side left -padx 2

} elseif {$Axis == "y"} {

label  $w.frame1.label -text "Define plot type for Y-direction \[2D/3D\]"
pack   $w.frame1.label -side left -padx 2

} elseif {$Axis == "z"} {

label  $w.frame1.label -text "Define plot type for Z-direction \[2D/3D\]"
pack   $w.frame1.label -side left -padx 2

} else {

    gomError {unknown axis direction '$Axis'}
    return
}

    radiobutton $w.frame1.type2d -text "2D-type" -value 0 -variable gomCutplaneType \
                -command {}
    pack  $w.frame1.type2d -side left -anchor w
    radiobutton $w.frame1.type3d -text "3D-type" -value 1 -variable gomCutplaneType \
                -command {}
    pack  $w.frame1.type3d -side left -anchor w
    label $w.frame1.label2 -text "Damping: "
    pack  $w.frame1.label2 -side left -anchor w
    entry $w.frame1.damping -width 10
    pack  $w.frame1.damping -side left -anchor w

# place in damping factor
    set damping [show cutplane damping]
    $w.frame1.damping delete 0 end
    $w.frame1.damping insert 0 $damping

    switch $Axis {
      x {set type [lindex $type 0]}    
      y {set type [lindex $type 1]}    
      z {set type [lindex $type 2]}    
    }

    if {[string match "2d*" "$type"]} {
        set gomCutplaneType 0
        $w.frame1.type2d select
    } else {
        set gomCutplaneType 1
        $w.frame1.type3d select
    }

}

############################################################################
# PROC
proc lulApplyCutplaneType { Axis w} {

     global gomCutplaneType

  set damping [$w.frame1.damping get]
  define cutplane damping $damping

if { $Axis == "x" } {
  
  if {$gomCutplaneType} {
       plot cutplane 3dx
  } else {
       plot cutplane 2dx
  }

  lulInvalidateDisplay
} elseif { $Axis == "y" } {
  
  if {$gomCutplaneType} {
       plot cutplane 3dy
  } else {
       plot cutplane 2dy
  }

  lulInvalidateDisplay
} elseif { $Axis == "z" } {
  
  if {$gomCutplaneType} {
       plot cutplane 3dz
  } else {
       plot cutplane 2dz
  }

  lulInvalidateDisplay
} else {

    gomError "unknown axis direction '$Axis'"
    return
}

}


############################################################################
# PROC
proc lulAnimateCutplane { Axis } {

    global gomControlFont
    global gomContourColour
    global gomMaxContourLevels
    global gomContourSmoothState
    global gomCutplaneLoopState
    global gomHelpFile
    global gomContourID
    global gomCutX
    global gomCutY
    global gomCutZ

# return if no contours defined
    if {[show contour defined]   < 1} {
        lulErrorDialog {ERROR: no contour available. Read a contour first!}
        return
    }
    
    set CapitalAxis [string toupper $Axis]
    
    set w .gomanimatecontplane$Axis
    catch {destroy $w}
    toplevel $w 
    wm title $w "Animate Contour Plane Control $CapitalAxis"
    wm iconname $w "Animate Contour Plane Control $CapitalAxis"

#
    frame  $w.buttons -borderwidth 2 -relief raised
    pack   $w.buttons -side bottom -fill x -pady 2m
    button $w.buttons.dismiss -text Dismiss  -font "$gomControlFont" \
	-command "destroy $w"
    button $w.buttons.apply   -text Apply    -font "$gomControlFont" \
	-command "lulApplyContourPlaneAnimate $w" -state disabled
    button $w.buttons.help    -text Help     -font "$gomControlFont" \
	-command \
	"htmlShowHelp $gomHelpFile(contplane_animate)"

    pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1


#
    set  ww .gomcontourplane
    scan [show contour cube $gomContourID] "%f %f %f %f %f %f %d %d %d" \
	Min(x) Max(x) Min(y) Max(y) Min(z) Max(z) Dim(x) Dim(y) Dim(z)

    frame $w.settings -borderwidth 2 -relief ridge
    pack  $w.settings -side top -anchor w
    label $w.settings.label -text "$CapitalAxis axis:"
    pack  $w.settings.label -side top -anchor w
    frame $w.settings.type
    pack  $w.settings.type -side top -anchor w
    frame $w.settings.range
    pack  $w.settings.range -side top -anchor w


    frame $w.settings.range.min
    pack  $w.settings.range.min -side left
    label $w.settings.range.min.label -text "Min:"
    entry $w.settings.range.min.input -width 15
    pack  $w.settings.range.min.label $w.settings.range.min.input -side left
    $w.settings.range.min.input insert 0 [format %6.3e $Min($Axis)]

    frame $w.settings.range.max
    pack  $w.settings.range.max -side left
    label $w.settings.range.max.label -text "Max:"
    entry $w.settings.range.max.input -width 15
    pack  $w.settings.range.max.label $w.settings.range.max.input -side left
    $w.settings.range.max.input insert 0 [format %6.3e $Max($Axis)]

    frame $w.settings.range.steps
    pack  $w.settings.range.steps -side left
    label $w.settings.range.steps.label -text "Steps:"
    entry $w.settings.range.steps.input -width 5
    pack  $w.settings.range.steps.label $w.settings.range.steps.input -side left
    $w.settings.range.steps.input insert 0 $Dim($Axis)

    frame $w.settings.type.3d
    pack  $w.settings.type.3d -side left
    label $w.settings.type.3d.label -text "3D: "
    pack  $w.settings.type.3d.label -side left
    radiobutton $w.settings.type.3d.on -text "on" -value 1 -variable gomCutplaneType$CapitalAxis \
	-command "plot cutplane 3d$Axis;lulInvalidateDisplay"
    radiobutton $w.settings.type.3d.off -text "off" -value 0 -variable gomCutplaneType$CapitalAxis \
	-command "plot cutplane 2d$Axis;lulInvalidateDisplay"
    pack  $w.settings.type.3d.on $w.settings.type.3d.off -side left

    frame $w.settings.type.damping
    pack  $w.settings.type.damping -side left
    label $w.settings.type.damping.label -text "Damping: "
    entry $w.settings.type.damping.input -width 10
    button $w.settings.type.damping.apply -text "Apply" \
	-command "lulSetCutPlaneDamping $w.settings.type.damping.input"
    pack  $w.settings.type.damping.label $w.settings.type.damping.input \
	 $w.settings.type.damping.apply -side left
    $w.settings.type.damping.input insert 0 [show cutplane damping]

    scan [show cutplane type] "%s %s %s" Type(x) Type(y) Type(z)

    if {$Type($Axis) == "3dx"} {
	$w.settings.type.3d.on  select
    } else {
	$w.settings.type.3d.off select
    }

#	set   Delta [expr ($Max - $Min)/($Xdim - 1)]

    switch [expr \$gomCut$CapitalAxis] {
	0 { set CutSmooth ""}
	1 { set CutSmooth "sqrt"}
	2 { set CutSmooth "log10"}
    }

    set Command  "plot cutplane $Axis $gomContourID "
    set TickMark "$ww.box.frame2.[regsub $Axis "xyz" ""] select"
    set NewValue "$ww.box.frame2.$Axis"

    frame $w.control -borderwidth 2 -relief ridge
    pack  $w.control -side top -pady 5 -padx 5

    set gomCutplaneLoopState 0

    button $w.control.first  -text "|<" \
	-command "lulAnimatePlaneControl $w $ww $Axis start {$CutSmooth}" \
	-font {Arial 12 bold}
    pack   $w.control.first    -side left  -padx 3 -pady 4 
    button $w.control.back     -text "<"    \
	-command "lulAnimatePlaneControl $w $ww $Axis -1  {$CutSmooth}" \
	-font {Arial 12 bold}
    pack   $w.control.back     -side left  -padx 3 -pady 4
    button $w.control.reverse   -text "<|"   \
	-command "lulAnimatePlaneControl $w $ww $Axis reverse {$CutSmooth}" \
	-font {Arial 12 bold}
    pack   $w.control.reverse  -side left  -padx 3 -pady 4
    button $w.control.stop     -text "\[\]" \
	-command "set gomCutplaneLoopState 0" \
	-font {Arial 12 bold}
    pack   $w.control.stop     -side left  -padx 3 -pady 4
    button $w.control.play     -text "|>"   \
	-command "lulAnimatePlaneControl $w $ww $Axis play {$CutSmooth}" \
	-font {Arial 12 bold}
    pack   $w.control.play     -side left  -padx 3 -pady 4
    button $w.control.forward  -text ">"    \
	-command "lulAnimatePlaneControl $w $ww $Axis 1  {$CutSmooth}"  \
	-font {Arial 12 bold}
    pack   $w.control.forward  -side left  -padx 3 -pady 4
    button $w.control.last     -text ">|"                   \
	-command "lulAnimatePlaneControl $w $ww $Axis end {$CutSmooth}" \
	-font {Arial 12 bold}
    pack   $w.control.last     -side left  -padx 3 -pady 4
}


############################################################################
# PROC
proc lulSetCutPlaneDamping { w } {
   define cutplane damping [$w get]
   lulInvalidateDisplay
}

############################################################################
# PROC
proc lulAnimatePlaneControl { w ww Axis Action Smooth } {

    global gomControlFont
    global gomContourColour
    global gomMaxContourLevels
    global gomCutplaneLoopState
    global gomContourSmoothState
    global gomHelpFile
    global gomContourID
    global gomCutX
    global gomCutY
    global gomCutZ

    scan [show contour cube $gomContourID] "%f %f %f %f %f %f %d %d %d" \
	Min(x) Max(x) Min(y) Max(y) Min(z) Max(z) Dim(x) Dim(y) Dim(z)

    array set boxframes {x frame2 y frame3 z frame4}

    set MinV     [string trim [$w.settings.range.min.input get]]
    set MaxV     [string trim [$w.settings.range.max.input get]]
    set Steps    [string trim [$w.settings.range.steps.input get]]
    if {$MinV < $Min($Axis)} {
	set MinV $Min($Axis)
	$w.settings.range.min.input delete 0 end
	$w.settings.range.min.input insert 0 [format %6.3e $MinV]
    }
    if {$MaxV > $Max($Axis)} {
	set MaxV $Max($Axis)
	$w.settings.range.max.input delete 0 end
	$w.settings.range.max.input insert 0 [format %6.3e $MaxV]
    }

    set Delta    [expr ($MaxV-$MinV)/($Steps - 1.0)]
    set CurrentV [string trim [$ww.box.$boxframes($Axis).$Axis get]]
    set Fmin     [string trim [$ww.box.$boxframes($Axis).min get]]
    set Fmax     [string trim [$ww.box.$boxframes($Axis).max get]]
    set frame    [expr round((1.0*$CurrentV-$MinV)/(($MaxV-$MinV)/($Steps-1.0)))]

    set TickMark "$ww.box.$boxframes($Axis).[regsub $Axis "xyz" ""] select"
    set NewValue "$ww.box.$boxframes($Axis).$Axis"

    switch -- $Action {
      start {
	eval plot cutplane $Axis $gomContourID $MinV $Fmin $Fmax $Smooth
	display
	eval $TickMark
	eval $NewValue delete 0 end
	eval $NewValue insert 0 $MinV
      }
      reverse -
      play {
	if {$gomCutplaneLoopState} {return};
	set gomCutplaneLoopState 1;

	if {$Action == "play"} {set Op +} else {set Op -}

	if {$frame <  0     } {set frame 0}
	if {$frame >= $Steps} {set frame [expr $Steps - 1]}
	set frame [expr $frame $Op 1];
	if {$frame >= 0 && $frame < $Steps} {eval $TickMark}
	while {$frame >= 0 && $frame < $Steps} {
	    update idletasks;
	    if {!$gomCutplaneLoopState} {break};

	    if {$frame == $Steps - 1} {
		set CurrentV $MaxV
	    } else {
		set CurrentV [expr $MinV + $frame * $Delta]
	    }
	    eval plot cutplane $Axis $gomContourID $CurrentV $Fmin $Fmax $Smooth
	    display
	    eval $NewValue delete 0 end
	    eval $NewValue insert 0 $CurrentV
	    set frame [expr $frame $Op 1]
	}
	set gomCutplaneLoopState 0;
      }
      end {
	eval plot cutplane $Axis $gomContourID $MaxV $Fmin $Fmax $Smooth
	display
	eval $TickMark
	eval $NewValue delete 0 end
	eval $NewValue insert 0 $MaxV
      }
      default {
	set frame [expr $frame + $Action];
	if {$frame >= $Steps - 1} {
	    set CurrentV [expr $MaxV + ($frame - $Steps + 1) * $Delta]
	} else {
	    set CurrentV [expr $MinV + $frame * $Delta]
	}
	if {$CurrentV < $Min($Axis) || $CurrentV > $Max($Axis)} {
	    lulMessageDialog "Value: $CurrentV is outside of range $Min($Axis):$Max($Axis)"
	    return
	}
	eval plot cutplane $Axis $gomContourID $CurrentV $Fmin $Fmax $Smooth
	display
	eval $TickMark
	eval $NewValue delete 0 end
	eval $NewValue insert 0 $CurrentV
      }
    } 
}

############################################################################
# PROC
proc lulPlotContourProfile { w } {

	 global gomContourID
     global gomYZ
	 global gomXZ
	 global gomXY

# return if no contours defined
    if {[show contour defined]   < 1} {
        return
	}

    if { $gomYZ } {
	   set Command "plot cutplane profile $gomContourID  x  \
                                   [$w.profile.nbins   get]  \
	                               [$w.box.frame2.min  get]  \
								   [$w.box.frame2.max  get]  \
                                   [$w.profile.cutv get]"
	   eval $Command
	}  elseif { $gomXZ } {
	   set Command "plot cutplane profile $gomContourID  y  \
                                   [$w.profile.nbins   get]  \
                                   [$w.box.frame3.min  get]  \
                                   [$w.box.frame3.max  get]  \
                                   [$w.profile.cutv get]"
	   eval $Command
    } elseif { $gomXY } {
	   set Command "plot cutplane profile $gomContourID  z  \
                                   [$w.profile.nbins   get]  \
	                               [$w.box.frame4.min  get]  \
								   [$w.box.frame4.max  get]  \
                                   [$w.profile.cutv get]"
	   eval $Command
	}
} 

############################################################################
proc lulApplyContourPlaneCommand { w } {

     global gomContourID
     global gomYZ
	 global gomXZ
	 global gomXY
     global gomCutX
     global gomCutY
     global gomCutZ

# return if no contours defined
    if {[show contour defined]   < 1} {
        return
	}

    if { $gomYZ } {
       set CutSmooth ""
	   switch $gomCutX {
	     1 { set CutSmooth "sqrt"}
		 2 { set CutSmooth "log10"}
		}
	   set Command "plot cutplane x $gomContourID [$w.box.frame2.x get]    \
	                                              [$w.box.frame2.min get] \
									              [$w.box.frame2.max get] \
												  $CutSmooth"
	   eval $Command
    } else {
	   plot -cutplane x
	}

    if { $gomXZ } {
       set CutSmooth ""
	   switch $gomCutY {
	     1 { set CutSmooth "sqrt"}
		 2 { set CutSmooth "log10"}
		}
	   set Command "plot cutplane y $gomContourID [$w.box.frame3.y get]    \
	                                              [$w.box.frame3.min get] \
									              [$w.box.frame3.max get] \
												  $CutSmooth"
	   eval $Command
    }  else {
	   plot -cutplane y
    }

    if { $gomXY } {
       set CutSmooth ""
	   switch $gomCutZ {
	     1 { set CutSmooth "sqrt"}
		 2 { set CutSmooth "log10"}
		}
	   set Command "plot cutplane z $gomContourID [$w.box.frame4.z get]    \
	                                              [$w.box.frame4.min get] \
									              [$w.box.frame4.max get] \
												  $CutSmooth"
	   eval $Command
    } else {
	   plot -cutplane z
	}

    lulInvalidateDisplay
}
