############################################################################
# This is a first approximation for a script to prepare
# animations (movies) from a trajectory file.
#
# Copyright Leif Laaksonen/CSC 1997 Version 0.0
#
proc lulMakeMovie { } {

     global gomControlFont
     global gomPictureType
     global gomPictureExt
     global gomHelpDir
     global gomHelpFile
     global env

# return if no molecular systems defined
     if {[show trajectory frames] < 1} {
         lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
        return
        }

set w .makemovie
catch {destroy $w}
toplevel $w
wm title $w "Make Movie"
wm iconname $w "Make Movie"

frame  $w.buttons -borderwidth 2 -relief raised
pack   $w.buttons -side bottom -fill x -pady 2m
button $w.buttons.dismiss -text Dismiss -font "$gomControlFont" \
        -command "destroy $w"
button $w.buttons.apply   -text Apply   -font "$gomControlFont" \
        -command "lulApplyMakeMovie $w"
button $w.buttons.help    -text Help    -font "$gomControlFont" \
        -command \
          "htmlShowHelp $gomHelpFile(makemovie)"
pack   $w.buttons.dismiss $w.buttons.apply $w.buttons.help -side left -expand 1

label  $w.label -font $gomControlFont -text "Make movie:"
pack   $w.label -side top -anchor w

frame        $w.frame1 -borderwidth 2 -relief ridge
label        $w.frame1.label -text "File formats: "
radiobutton  $w.frame1.bmp -text "BMP"        -variable gomPictureType -value bmp
radiobutton  $w.frame1.rgb -text "RGB"        -variable gomPictureType -value rgb
radiobutton  $w.frame1.tga -text "TGA"        -variable gomPictureType -value tga
radiobutton  $w.frame1.xwd -text "XWD"        -variable gomPictureType -value xwd

pack   $w.frame1 -side top -anchor w
pack   $w.frame1.label $w.frame1.bmp $w.frame1.rgb $w.frame1.tga $w.frame1.xwd -side left

frame  $w.frame2
label  $w.frame2.label     -text "Base file name: "
entry  $w.frame2.filename  -width 40
button $w.frame2.browse    -text "Browse..." -command "lulDoMakeMovieAction $w"

pack   $w.frame2 -side top -pady 6
pack   $w.frame2.label $w.frame2.filename -side left
pack   $w.frame2.browse -padx 3 -side left
bind   $w.frame2.filename <Return> "eval lulApplyMakeMovie $w"

$w.frame1.tga select
set gomPictureType "tga"


}


########################################################################
proc lulDoMakeMovieAction { w } {

     global gomPictureType
     global gomPictureExt

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    if {$gomPictureType == "bmp"} {
       set types {
         {"BMP files"		".bmp .BMP"		       TEXT}
         {"All files"		*}
       }
    } elseif {$gomPictureType == "rgb"} {
       set types {
         {"SGI RGB files"	".rgb .RGB"         	TEXT}
         {"All files"		*}
       }
    } elseif {$gomPictureType == "tga"} {
       set types {
         {"Targa files"	".tga .TGA"         	TEXT}
         {"All files"		*}
       }
    } elseif {$gomPictureType == "xwd"} {
       set types {
         {"XWD files"	".xwd .XWD"         	TEXT}
         {"All files"		*}
       }
    }

    set file [tk_getSaveFile -filetypes $types -parent $w   \
         -defaultextension *.$gomPictureType                \
		 -title "Make Picture File for a Movie"]

    if [string compare $file ""] {
       regsub -all {[0-9]} $file "" file
       set file [file rootname $file]
       $w.frame2.filename delete 0 end
       $w.frame2.filename insert 0 $file
    }

}
###########################################################################
# PROC
proc lulApplyMakeMovie { w } {

     global gomPictureType
     global gomEnv
     global tcl_platform
     global env

     set File [$w.frame2.filename get]
     regsub -all {[0-9]} $File "" File

# check window size...
     set WinInfo [show window paramete]
     set Width   [lindex $WinInfo 2]
     set Heigth  [lindex $WinInfo 3]

     set Size [expr $Width * $Heigth]
     if {$Size > 160000} {
       gomError "your graphics window is too big ( > 400 * 400 pixels)"
       return
     }
# ready to go if file name is supplied

     if {$File != ""} {

	 if {[file isdirectory $File]} {
      gomError "your file name is a directory"
      return
     }
#
# some string manipulation ...
#
         set RootName     [file rootname $File]
         set DirName      [file dirname  $File]
         set FileNameFull [file tail     $File]
         set FileNameRoot [file root     $FileNameFull]

         if {[string length $FileNameRoot] > 5} {
           gomError "This feature is using an old DOS\nprogram that can't handle long file names.\nPlease use short names ( < 6 chars long)" 
           return
         }

# if no directory defined direct it to gOpenMol temp
         if {$DirName == "" || $DirName == "."} {
           set DirName "$gomEnv(GOM_TEMP)"
           set RootName "[file join $gomEnv(GOM_TEMP) $RootName]"
         }

# return if no molecular systems defined
     set NumFrames [show trajectory frames]
	 if {$NumFrames < 1} {
            lulErrorDialog {ERROR: no trajectory is defined. Define a trajectory first!}
            return
	 }

# loop the frames
         set Loop 0
         for {set i 1} { $i <= $NumFrames} { incr i} {
           import coord frame $i
           display 1
#           eval "hardcopy 1 $gomPictureType $RootName$Loop.$gomPictureType"
#           incr Loop
           eval "hardcopy 1 $gomPictureType $RootName$Loop.$gomPictureType"
           incr Loop
#
# Put in here the exec to convert your file into a form/type that can
# be processed by your mpeg "gluer"
#

	 }

# glue the targa files into a mpeg file

if {$tcl_platform(platform) == "windows"} {

    if {$gomPictureType != "tga"} {
        gomError "Only Targa (tga) files can be made into mpeg!"
        return
        }

      puts "Will now glue the files into file '$RootName.mpg'"

# save first old location
      set CurrDir [pwd]
# Change first to the desired directory
      cd $DirName
# Copy the program to the directory with the Targa files ...
      file copy -force [file join $gomEnv(GOM_ROOT) CMpeg1.0 cmpeg.exe] $DirName
# Copy the control file ...
      file copy -force [file join $gomEnv(GOM_ROOT) CMpeg1.0 ipb.ctl]   $DirName

# then run the program

       set Name [file rootname [file tail $File]]

       set Command "cmpeg.exe -v1 -t0 ipb.ctl $Name%d.tga $Name.mpg"

       puts "Starting the job. Please stand by!"

       set f [open "|$Command" r]

       set i 1

       while {![eof $f]} { 
          gets $f Text
          puts $Text
          incr i
          puts "Please stand by...!"
          update
       }

#close pipe
       close $f

# go back to original directory
      cd $CurrDir

       } else {
         puts "The making of mpg files is so far not supported on UNIX!"
       }
      puts "Done! File '$RootName.mpg' is ready"
     }
}
###########################################################################
# PROC
proc lulApplyMakeMovieRotation { File axis angle steps } {

     global gomPictureType
     global gomEnv
     global tcl_platform
     global env

     puts "Starting making movie..."
     
     if {[string trim $File] == ""} {
       set File "Zm"
     }
     regsub -all {[0-9]} $File "" File

# check window size...
     set WinInfo [show window paramete]
     set Width   [lindex $WinInfo 2]
     set Heigth  [lindex $WinInfo 3]
     set gomPictureType "tga"
     set angle [format "%f" $angle]
     set steps [format "%f" $steps]

     puts "Window size is $Width x $Heigth"
     
     set Size [expr $Width * $Heigth]
     if {$Size > 160000} {
       gomError "your graphics window is too big ( > 400 * 400 pixels)"
       return
     }
# ready to go if file name is supplied

     puts "Input root file name is $File"
     
     if {$File != ""} {

	 if {[file isdirectory $File]} {
      gomError "your file name is a directory"
      return
     }
#
# some string manipulation ...
#
         set RootName [file rootname $File]
         set DirName  [file dirname  $File]
         set FileNameFull [file tail     $File]
         set FileNameRoot [file root     $FileNameFull]

         if {[string length $FileNameRoot] > 5} {
           gomError "This feature is using an old DOS\nprogram that can't handle long file names.\nPlease use short names ( < 6 chars long)" 
           return
         }

# if no directory defined direct it to gOpenMol temp
         if {$DirName == "" || $DirName == "."} {
           set DirName "$gomEnv(GOM_TEMP)"
           regsub -all {\\} $DirName "/" DirName
           set RootName "[file join $gomEnv(GOM_TEMP) $RootName]"
         }

         puts "Starting the loop with the command 'hardcopy 1 $gomPictureType $RootName$Loop.$gomPictureType'"
         
# loop the frames
         set Loop 0
         set sum 0.0
#           eval "hardcopy 1 $gomPictureType $RootName$Loop.$gomPictureType"
#           incr Loop
           eval "hardcopy 1 $gomPictureType $RootName$Loop.$gomPictureType"
           incr Loop
         set dstep [expr $angle/$steps]
         for {set i 1} { $i <= $steps} { incr i} {
           if {$axis == "x"} {
             rotate display $dstep 0.0 0.0
           } elseif {$axis == "y"} {
             rotate display 0.0 $dstep 0.0
           } elseif {$axis == "z"} {
             rotate display 0.0 0.0 $dstep 
             set sum [expr $sum + $dstep]
             puts "sum: $sum , dstep: $dstep"
           } else {
             gomError "unknown angle parameter '$axis'. Has to be x, y or z"
             return
           }
           display 1
           eval "hardcopy 1 $gomPictureType $RootName$Loop.$gomPictureType"
           incr Loop
           eval "hardcopy 1 $gomPictureType $RootName$Loop.$gomPictureType"
           incr Loop
#
# Put in here the exec to convert your file into a form/type that can
# be processed by your mpeg "gluer"
#

	 }

# glue the targa files into a mpeg file

if {$tcl_platform(platform) == "windows"} {

    if {$gomPictureType != "tga"} {
        gomError "Only Targa (tga) files can be made into mpeg!"
        return
        }

      puts "Will now glue the files into file '$RootName.mpg'"

# save first old location
      set CurrDir [pwd]
# Change first to the desired directory
      cd $DirName
# Copy the program to the directory with the Targa files ...
      file copy -force [file join $gomEnv(GOM_ROOT) CMpeg1.0 cmpeg.exe] $DirName
# Copy the control file ...
      file copy -force [file join $gomEnv(GOM_ROOT) CMpeg1.0 ipb.ctl]   $DirName

# then run the program

       set Name [file rootname [file tail $File]]

       set Command "cmpeg.exe -v1 -t0 ipb.ctl $Name%d.tga $Name.mpg"

       puts "Starting the job. Please stand by!"

       set f [open "|$Command" r]

       set i 1

       while {![eof $f]} { 
          gets $f Text
          puts $Text
          incr i
          puts "Please stand by...!"
          update
       }

#close pipe
       close $f
# delete temp files
       gomPrint "Deleting files '$Name*.tga'" 
       eval file delete [glob $Name*.tga]
# go back to original directory
       cd $CurrDir

       } else {
         puts "The making of mpg files is so far not supported on UNIX!"
       }
      puts "Done! File '$RootName.mpg' is ready"
     }
}


