
package require http 2.0


##########################################################################
# NAMESPACE lulUtility
#
namespace eval lulUtility {

   variable gomURLtoken

#########################################################
# PROC
# use the http protocol to retrieve an URL
proc gomHttpGetURL { url } {

   variable gomURLtoken

set token [::http::geturl $url]
set gomURLtoken $token

 if {[::http::status $token] != "ok"} {
    gomError "Can't retrieve URL: $url"
    return ""
 }

 return [::http::data $token]

}

#########################################################
# PROC
# use the http protocol to retrieve a coordinate file as an URL
proc gomGetCoordinatesURL { url } {

     global gomEnv
     global gomURLFileName
     variable gomURLtoken

 set CoordinateFile [gomHttpGetURL $url]  
 set ErrorCode [lindex [::http::code $gomURLtoken] 1]
 if {$ErrorCode != 200} {
	  set gomURLFileName ""
	  gomError "Can't get the URL '$url'. Error code '$ErrorCode'"
      ::http::cleanup $gomURLtoken
	  return
	  }

 set FileName [file tail $url]

 if {$FileName == ""} {
     gomError "File name is missing in URL '$url'"
     ::http::cleanup $gomURLtoken
	 return
	 }

 set gomURLFileName "$gomEnv(GOM_TEMP)/$FileName"
 set  File [open $gomURLFileName w]

 if {$File == ""} {
    gomError "Can't open file '$gomURLFileName' for writing"
    ::http::cleanup $gomURLtoken
	return
	}

 puts  $File $CoordinateFile
 close $File

 ::http::cleanup $gomURLtoken

 return
}
#########################################################
# PROC
# get current gom version from http://www.csc.fi/gopenmol/current-version.html
#
proc gomGetCurrentVersion { } {

return [gomHttpGetURL http://www.csc.fi/gopenmol/current-version.html]

}

#########################################################
# PROC
proc gomPeekCurrentProgramVersion {} {

set gomRunningVersion "[format "%.2f" [expr [string trim [show version]]/100.]]"

set gomCurrentVersion [gomGetCurrentVersion]
set length   [llength  $gomCurrentVersion]

set  AppString ""
for {set i 1} {$i < $length} {incr i} {
     lappend AppString [lindex $gomCurrentVersion $i]
     }

gomPrint "Your current version is: $gomRunningVersion"
gomPrint "Latest available version is: [lindex $gomCurrentVersion 0]"
gomPrint "Available for the platforms: $AppString" 

lulMessageDialog "Your version is: $gomRunningVersion\nLatest available version is: [lindex $gomCurrentVersion 0]\nAvailable for the platforms: $AppString"
}
#########################################################
# PROC
# calculate the centroid for a surface
#
proc gomCalculateSurfaceCentroid { } {

  global tcl_platform

  if {![show contour defined]} {
    gomError "No contours defined to get the mesh data"
    return
  }

  if {[show contour polygon method] != "save"} {
     gomError "You have to turn the 'save' option on\nwith the command 'contour method save'"
     return
  }

  set Entries [show contour polygon entries]

  if {$Entries < 1} {
     gomError "No contour polygons available"
     return
  }

     set sumx 0.0
     set sumy 0.0
     set sumz 0.0
 
  for {set i 1} {$i <= $Entries} {incr i} {

     set Text [show contour polygon data $i] 
# puts "Done with: $i ($Entries)"
     set x1 [lindex $Text 0]
     set y1 [lindex $Text 1]
     set z1 [lindex $Text 2]

     set x2 [lindex $Text 6]
     set y2 [lindex $Text 7]
     set z2 [lindex $Text 8]

     set x3 [lindex $Text 12]
     set y3 [lindex $Text 13]
     set z3 [lindex $Text 14]

     set sumx [expr $sumx + $x1 + $x2 + $x3]
     set sumy [expr $sumy + $y1 + $y2 + $y3]
     set sumz [expr $sumz + $z1 + $z2 + $z3]
  }

     set sumx [expr $sumx / (3 * $Entries)]
     set sumy [expr $sumy / (3 * $Entries)]
     set sumz [expr $sumz / (3 * $Entries)]

     gomPrint "Surface centroid at (x,y,z): $sumx , $sumy , $sumz"

# save to clipboard on Windows
     if {$tcl_platform(platform) == "windows"} {
       copy text "$sumx  $sumy  $sumz" 
       gomPrint "The following text is also put into the 'Paste' buffer:\n$sumx  $sumy  $sumz"
     }

     return "$sumx  $sumy  $sumz"
}

#
# small tcl script to shuffle the triangular mesh data
# from the contour display over to the "plot triangle"
# command
#
proc gomContourData2Polygon {} {

   set Method    [show contour polygon method]
   if {$Method == "direct"} {
       gomError "Rerun the plot but turn first the\nsave option on with command:\ncontour method save"
       return
   }
   set Triangles [show contour polygon entries] 

   if {$Triangles < 1} {
       gomError "No contour polygon data available"
       return
   }
gomPrint "Converting now $Triangles triangles"

   for {set i 1} {$i <= $Triangles} {incr i} {

     set Fvalue [expr $i./$Triangles.]
     lulTimeFliesWrapper $Fvalue

     set Text [show contour polygon data $i]
       set x1  [lindex $Text 0]
       set y1  [lindex $Text 1]
       set z1  [lindex $Text 2]
       set nx1 [lindex $Text 3]
       set ny1 [lindex $Text 4]
       set nz1 [lindex $Text 5]
       set x2  [lindex $Text 6]
       set y2  [lindex $Text 7]
       set z2  [lindex $Text 8]
       set nx2 [lindex $Text 9]
       set ny2 [lindex $Text 10]
       set nz2 [lindex $Text 11]
       set x3  [lindex $Text 12]
       set y3  [lindex $Text 13]
       set z3  [lindex $Text 14]
       set nx3 [lindex $Text 15]
       set ny3 [lindex $Text 16]
       set nz3 [lindex $Text 17]
       set c1  [lindex $Text 18]
       set c2  [lindex $Text 19]
       set c3  [lindex $Text 20]

     plot triangle $x1 $y1 $z1 $nx1 $ny1 $nz1 \
                   $x2 $y2 $z2 $nx2 $ny2 $nz2 \
                   $x3 $y3 $z3 $nx3 $ny3 $nz3 \
                   $c1 $c2 $c3 1.0 append 
   }
# reset display
   lulTimeFliesWrapper -1.0
}

#
# Small script to merge all the currently read structures
# into one. At this stage the model viewing matrix is also
# applied on the atoms. 
#
proc MergeStructures {} {

set  NumStruct [show molstructures]

     if {$NumStruct < 2 } {
        gomError "only one structure is available"
        return
     }
     set TotNumAtoms 0
     set UpToAtoms   0
     set k 1
     set resnumloop 0
             
# save first data from current structures 
     for {set i 1} {$i <= $NumStruct} {incr i} {

       set TotAtoms [show numatoms $i]
       set AtomsInStructure($i) $TotAtoms
       set TotNumAtoms [expr $TotNumAtoms + $TotAtoms]

       for {set j  1} {$j <= $TotAtoms} {incr j} {

       scan [show atom coordinates   $j $i] "%f %f %f" xc($k) yc($k) zc($k)
       set  charge($k) [show atom charge $j $i]
       set  vdw($k)    [show atom vdw    $j $i]
       set  resn($k)   [expr $resnumloop + [show atom resnum $j $i]]
       set  seg($k)    [show atom segmen $j $i]
       set  res($k)    [show atom residu $j $i]
       set  atom($k)   [show atom atomna $j $i]
# connectivity
       set  text       [show atom connect $j $i]
       set  connections($k,1) [lindex $text 0]
       if {$connections($k,1) > 0} {
          for {set kk 1} {$kk <= $connections($k,1)} {incr kk} {
            set kkk [expr $kk + 1]
            set connections($k,$kkk) [expr $UpToAtoms + [lindex $text $kk]]
          }
       }
       incr k
       }

       set resnumloop [expr $resnumloop + [show atom resnum $TotAtoms $i]]
       set UpToAtoms [expr $UpToAtoms + $TotAtoms]
     }

     gomPrint "Total number of atoms: $TotNumAtoms"

# delete old data structure and set up the new one
     reset gopenmol
     define structure "Merged structure" $TotNumAtoms new

     for {set i 1} {$i <= $TotNumAtoms} {incr i} {
     
# put the values into the structure ...
         define atom    label    $atom($i)               $i 1
         unset atom($i)
         define residue label    $res($i)                $i 1
         unset res($i)
         define segment label    $seg($i)                $i 1
         unset seg($i)
         define atom coordinates $xc($i) $yc($i) $zc($i) $i 1
         unset xc($i)
         unset yc($i)
         unset zc($i)
         define atom charge      $charge($i)             $i 1
         unset charge($i)
         define atom resnumber   $resn($i)               $i 1
         unset resn($i)     
     } 

     import postproc 1 {Merged structure}

# break all new bonds created at the post-processing stage 
     edit bond break * * *
# now take the bonds saved from the original set

     for {set i 1} {$i <= $TotNumAtoms} {incr i} {
        if {$connections($i,1)} {
          for {set j 1} {$j <= $connections($i,1)} {incr j} {
            set jj [expr $j + 1]
            edit bond create * * $i * * $connections($i,$jj)
          }
        } 
     }
}
#
# Script to be executed at a gOpenMol reset
#
proc gOpenMolReset { } {

   global gomEnv
   global gomSecondaryStructureRecord
   global gomSecondaryStructureRecords

   # default transformation mode is global
   define transform global
   # reset all plumber data
   namespace inscope :: {source [file join $gomEnv(GOM_DATA) autotcl general plumber.tcl]}
}
#################################################################
#
# Joining two gamess irc-files.
# Both files start with the same transtion state, one going
# forwards and the other backwards.
# The order of the coordinate sets is reversed in the first file
# so the sequence in the result file is 
#
# reactant(s) - transition state - product(s)
#
# or in reversed order depending on the order in whih the files
# are given.
#
# G Höjer, 2001-05-15
#
#################################################################
# 
# User edits next three statements with filenames including
# paths.
#
proc join-gamess-irc {IrcFile1 IrcFile2 IrcFile3} {

# first irc file
#set IrcFile1 c:\\science\\gamess\\clase\\hcn-tsb.irc

# second irc file
#set IrcFile2 c:\\science\\gamess\\clase\\hcn-tsf.irc

# resulting irc file
#set IrcFile3 c:\\science\\gamess\\clase\\hcnxxx.irc

# open files ...

set File_p1 [open $IrcFile1 r]
if {$File_p1 == ""} {
  catch {lulErrorDialog {ERROR: can't open GAMESS IRC file 1 for input!}} errRet
  if {$errRet != ""} {
    puts "ERROR: can't open GAMESS IRC file '$File_p1' for input"
  }
  error "ERROR: can't open GAMESS IRC file 1 for input!"
  return
}

set File_p2 [open $IrcFile2 r]
if {$File_p2 == ""} {
  catch {lulErrorDialog {ERROR: can't open GAMESS IRC file 2 for input!}} errRet
  if {$errRet != ""} {
    puts "ERROR: can't open GAMESS IRC file '$File_p2' for input"
  }
  error "ERROR: can't open GAMESS IRC file 2 for input!"
  return
}

set File_p3 [open $IrcFile3 w]
if {$File_p3 == ""} {
  catch {lulErrorDialog {ERROR: can't open GAMESS IRC file 3 for output!}} errRet
  if {$errRet != ""} {
    puts "ERROR: can't open GAMESS IRC file '$File_p3' for output"
  }
  error "ERROR: can't open GAMESS IRC file 3 for output!"
  return
}

# do the first part of the job ...
# find number of coord_sets and number of atoms in first set
# there must be the same number in each set

gomPrint "Reading GAMESS IRC coordinate file '$IrcFile1'..."

set coord_sets  0
set Trigger     0

while {![eof $File_p1]} {
  gets $File_p1 Text

  if {[string match "CARTESIAN COORDINATES (BOHR)" $Text]} {
    incr coord_sets
    set natoms 0

    while {![eof $File_p1]} {
      gets $File_p1 Text
      if {[lindex $Text 0] == "MASS-WEIGHTED"} {
        break
      }
      incr natoms
    }

    if {$Trigger ==0} {
      set natoms1 $natoms
    }
    set Trigger 1
    if {($Trigger ==1) && ($natoms !=$natoms1)} {
      catch {lulErrorDialog {ERROR: Wrong number of atoms in a set}} errRet
      if {$errRet != ""} {
        gomPrint "ERROR: Wrong number of atoms in set '$coord_sets'"
      }
      error "ERROR: Wrong number of atoms in a set"
      return
    }
  }
}

gomPrint "Found '$coord_sets' coordinate sets and '$natoms' atoms in each set in first file"

if {$natoms < 1} {
  gomError "can't find atoms in the irc file 1"
  return
}

gomPrint "Joining '$IrcFile1' and '$IrcFile2' to form '$IrcFile3'"

# do the second part of the job reordering the first file ...
set nlines [expr (2 * $natoms + 4)]
gomPrint "Lines per set: $nlines"

for {set i 1} {$i <= $coord_sets} {incr i} {
  set j [expr (($coord_sets - $i) * $nlines)]
  seek $File_p1 0 start

# position  file ..

  for {set k 1} {$k <= $j} {incr k} {
    gets $File_p1 Text
  }  

  for {set k 1} {$k <= $nlines} {incr k} {
    gets $File_p1 Text
# add line to output file
    puts $File_p3 "$Text"
  }  
}

# do the last part of the job appending the second file to the output file ... 
while {![eof $File_p2]} {
  gets $File_p2 Text
  puts $File_p3 "$Text"
}

# close file and return ...
close $File_p1
close $File_p2
close $File_p3
puts Done!
}
#
# Procedure to write out time series as a html page and read
# them into either the internal browser or an external one
#
proc TimeSeries2HTML { } {
#
# open the temp file in the temp directory
#
set HTMLfile [file join $gomEnv(GOM_TEMP) temp-file.html]

}
# end
}
