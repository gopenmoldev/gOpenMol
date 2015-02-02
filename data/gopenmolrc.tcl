############################################################################
#                       g O p e n M o l                                    # 
# This is the system startup file general for all non-graphics parts       #
#                                                                          #
# Copyright Leif Laaksonen 2005                                            #
#                                                                          #
############################################################################
#
# Check that you have at least version tcl/tk 8.4

set version [info tclversion]
if {$version < 8.4} {
   puts "\n\n"
   puts "****************************************************"
   puts "You have to be running Tcl/Tk version 8.4 or higher!" 
   puts "****************************************************"
   exit 1
}
   
# Default parser (tcl/python)
#
set gomParser "tcl"
set gomAllowConnection {127.0.0.1}
#
# Gnuplot
#
if {[regexp -nocase "windows" $tcl_platform(platform)]} {
    set gnuplotexe wgnupl32.exe
} else {
}
#
# Define the elements ...
#
set gomAtom_SN(1)   "H    Hydrogen" 
set gomAtom_SN(2)   "He   Helium" 
set gomAtom_SN(3)   "Li   Lithium" 
set gomAtom_SN(4)   "Be   Beryllium" 
set gomAtom_SN(5)   "B    Boron"
set gomAtom_SN(6)   "C    Carbon"
set gomAtom_SN(7)   "N    Nitrogen" 
set gomAtom_SN(8)   "O    Oxygen" 
set gomAtom_SN(9)   "F    Fluorine" 
set gomAtom_SN(10)  "Ne   Neon" 
set gomAtom_SN(11)  "Na   Sodium" 
set gomAtom_SN(12)  "Mg   Magnesium" 
set gomAtom_SN(13)  "Al   Aluminum" 
set gomAtom_SN(14)  "Si   Silicon" 
set gomAtom_SN(15)  "P    Phosphorus" 
set gomAtom_SN(16)  "S    Sulfur" 
set gomAtom_SN(17)  "Cl   Chlorine" 
set gomAtom_SN(18)  "Ar   Argon" 
set gomAtom_SN(19)  "K    Potassium" 
set gomAtom_SN(20)  "Ca   Calcium" 
set gomAtom_SN(21)  "Sc   Scandium" 
set gomAtom_SN(22)  "Ti   Titanium" 
set gomAtom_SN(23)  "V    Vanadium" 
set gomAtom_SN(24)  "Cr   Chromium" 
set gomAtom_SN(25)  "Mn   Manganese" 
set gomAtom_SN(26)  "Fe   Iron" 
set gomAtom_SN(27)  "Co   Cobalt" 
set gomAtom_SN(28)  "Ni   Nickel" 
set gomAtom_SN(29)  "Cu   Copper" 
set gomAtom_SN(30)  "Zn   Zinc" 
set gomAtom_SN(31)  "Ga   Gallium" 
set gomAtom_SN(32)  "Ge   Germanium" 
set gomAtom_SN(33)  "As   Arsenic" 
set gomAtom_SN(34)  "Se   Selenium" 
set gomAtom_SN(35)  "Br   Bromine" 
set gomAtom_SN(36)  "Kr   Krypton" 
set gomAtom_SN(37)  "Rb   Rubidium" 
set gomAtom_SN(38)  "Sr   Strontium" 
set gomAtom_SN(39)  "Y    Yttrium" 
set gomAtom_SN(40)  "Zr   Zirconium" 
set gomAtom_SN(41)  "Nb   Niobium" 
set gomAtom_SN(42)  "Mo   Molybdenum" 
set gomAtom_SN(43)  "Tc   Technetium" 
set gomAtom_SN(44)  "Ru   Ruthenium" 
set gomAtom_SN(45)  "Rh   Rhodium" 
set gomAtom_SN(46)  "Pd   Palladium" 
set gomAtom_SN(47)  "Ag   Silver" 
set gomAtom_SN(48)  "Cd   Cadmium" 
set gomAtom_SN(49)  "In   Indium" 
set gomAtom_SN(50)  "Sn   Tin" 
set gomAtom_SN(51)  "Sb   Antimony" 
set gomAtom_SN(52)  "Te   Tellurium" 
set gomAtom_SN(53)  "I   Iodine" 
set gomAtom_SN(54)  "Xe   Xenon" 
set gomAtom_SN(55)  "Cs   Cesium" 
set gomAtom_SN(56)  "Ba   Barium" 
# Lanthanide series
set gomAtom_SN(57)  "La   Lanthanum" 
set gomAtom_SN(58)  "Ce   Cerium" 
set gomAtom_SN(59)  "Pr   Praseodymium" 
set gomAtom_SN(60)  "Nd   Neodymium" 
set gomAtom_SN(61)  "Pm   Promethium" 
set gomAtom_SN(62)  "Sm   Samarium" 
set gomAtom_SN(63)  "Eu   Europium" 
set gomAtom_SN(64)  "Gd   Gadolinium" 
set gomAtom_SN(65)  "Tb   Terbium" 
set gomAtom_SN(66)  "Dy   Dysprosium" 
set gomAtom_SN(67)  "Ho   Holmium" 
set gomAtom_SN(68)  "Er   Erbium" 
set gomAtom_SN(69)  "Tm   Thulium" 
set gomAtom_SN(70)  "Yb   Ytterbium" 
set gomAtom_SN(71)  "Lu   Lutetium" 
set gomAtom_SN(72)  "Hf   Hafnium" 
set gomAtom_SN(73)  "Ta   Tantalum" 
set gomAtom_SN(74)  "W    Tungsten" 
set gomAtom_SN(75)  "Re   Rhenium" 
set gomAtom_SN(76)  "Os   Osmium" 
set gomAtom_SN(77)  "Ir   Iridium" 
set gomAtom_SN(78)  "Pt   Platinum" 
set gomAtom_SN(79)  "Au   Gold" 
set gomAtom_SN(80)  "Hg   Mercury" 
set gomAtom_SN(81)  "Tl   Thallium" 
set gomAtom_SN(82)  "Pb   Lead" 
set gomAtom_SN(83)  "Bi   Bismuth" 
set gomAtom_SN(84)  "Po   Polonium" 
set gomAtom_SN(85)  "At   Astatine" 
set gomAtom_SN(86)  "Rn   Radon" 
set gomAtom_SN(87)  "Fr   Francium" 
set gomAtom_SN(88)  "Ra   Radium"
# Actinide series
set gomAtom_SN(89)  "Ac   Actinium" 
set gomAtom_SN(90)  "Th   Thorium" 
set gomAtom_SN(91)  "Pa   Protactinium" 
set gomAtom_SN(92)  "U    Uranium" 
set gomAtom_SN(93)  "Np   Neptunium" 
set gomAtom_SN(94)  "Pu   Plutonium" 
set gomAtom_SN(95)  "Am   Americium" 
set gomAtom_SN(96)  "Cm   Curium" 
set gomAtom_SN(97)  "Bk   Berkelium" 
set gomAtom_SN(98)  "Cf   Californium" 
set gomAtom_SN(99)  "Es   Einsteinium" 
set gomAtom_SN(100) "Fm   Fermium" 
set gomAtom_SN(101) "Md   Mendelevium" 
set gomAtom_SN(102) "No   Nobelium" 
set gomAtom_SN(103) "Lr   Lawrencium" 

# Start the Tcl  definitions .........
gomPrint " "
gomPrint "#################################"
gomPrint "Processing System Startup File..."
gomPrint "#################################"
gomPrint " "

# set gOpenMol version
set gomVersion "[format "%.2f" [expr [string trim [show version]]/100.]]"

set gOpenMolInfo "*** gOpenMol version $gomVersion ([show release]) ***
This is the gOpenMol program for the display and analysis of molecular structures.\n 
Copyright Leif Laaksonen (1997 - 2005)\nScientific discovery and engineering: \nKevin Boyd, Eero Häkkinen, Leif Laaksonen"

###################################################################################
#
# This controls the atom select routine:
# == 0 , no selection list will be built
# != 0 , a list is built and can be found in the variable "gomAtomHitList"
# The selected atoms are collected as running numbers in the "gomAtomHitList"
# variable. Example: {1 2 3 5 6 8 9}
#
set gomAtomHitListActive 0

gomPrint "loading macro 'm_distance x1 y1 z1 x2 y2 z2' ..."
#
# m_distance x1 y1 z1 x2 y2 z2
#
proc m_distance {x1 y1 z1 x2 y2 z2 } {

set dx [expr $x1 - $x2]
set dy [expr $y1 - $y2]
set dz [expr $z1 - $z2]  

set dr [expr sqrt($dx*$dx + $dy*$dy + $dz*$dz)]

#distance is ...
return $dr
}

gomPrint "loading macro 'm_angle x1 y1 z1 x2 y2 z2 x3 y3 z3' ..."
#
# m_angle x1 y1 z1 x2 y2 z2 x3 y3 z3
#
proc m_angle {x1 y1 z1 x2 y2 z2 x3 y3 z3} {

set  tmp1  [expr $x1 - $x2]
set  tmp2  [expr $y1 - $y2]
set  tmp3  [expr $z1 - $z2]

set  d212  [expr $tmp1 * $tmp1 + $tmp2 * $tmp2 + $tmp3 * $tmp3]

set  tmp1  [expr $x2 - $x3]
set  tmp2  [expr $y2 - $y3]
set  tmp3  [expr $z2 - $z3]

set  d223  [expr $tmp1 * $tmp1 + $tmp2 * $tmp2 + $tmp3 * $tmp3]

set tmp1    [expr $x1 - $x3]
set tmp2    [expr $y1 - $y3]
set tmp3    [expr $z1 - $z3]

set  d213  [expr $tmp1 * $tmp1 + $tmp2 * $tmp2 + $tmp3 * $tmp3]

set  temp  [expr 0.5 * ($d212 + $d223 - $d213)/sqrt($d212 * $d223)]

if {$temp >  1.0} {set temp  1.0}
if {$temp < -1.0} {set temp -1.0}

set angijk  [expr acos( $temp )]

# angle is ...
return $angijk
}

gomPrint "loading macro 'm_torsion x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4' ..."
#
# m_torsion x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4
#
proc m_torsion {iax iay iaz ibx iby ibz icx icy icz idx idy idz} {

#        bond lengths and unit vectors    

set   eabx [expr $ibx - $iax]
set   eaby [expr $iby - $iay]
set   eabz [expr $ibz - $iaz]

set   rab  [expr sqrt ($eabx * $eabx + $eaby * $eaby + $eabz * $eabz)]
set   xrab [expr 1.0 / $rab]

set   eabx [expr $eabx * $xrab]
set   eaby [expr $eaby * $xrab]
set   eabz [expr $eabz * $xrab]

set   ebcx [expr $icx - $ibx]
set   ebcy [expr $icy - $iby]
set   ebcz [expr $icz - $ibz]

set   rbc  [expr sqrt ($ebcx * $ebcx + $ebcy * $ebcy + $ebcz * $ebcz)]
set   xrbc [expr 1.0 / $rbc]

set   ebcx [expr $ebcx * $xrbc]
set   ebcy [expr $ebcy * $xrbc]
set   ebcz [expr $ebcz * $xrbc]

set   ecdx [expr $idx - $icx]
set   ecdy [expr $idy - $icy]
set   ecdz [expr $idz - $icz]

set   rcd  [expr sqrt ($ecdx * $ecdx + $ecdy * $ecdy + $ecdz * $ecdz)]
set   xrcd [expr 1.0 / $rcd]

set   ecdx [expr $ecdx * $xrcd]
set   ecdy [expr $ecdy * $xrcd]
set   ecdz [expr $ecdy * $xrcd]

#        cross and dot products between unit vectors, and bond (!)
#        angles

set   abbcx [expr $eaby * $ebcz - $eabz * $ebcy]
set   abbcy [expr $eabz * $ebcx - $eabx * $ebcz]
set   abbcz [expr $eabx * $ebcy - $eaby * $ebcx]

set   cosb  [expr -($eabx * $ebcx + $eaby * $ebcy + $eabz * $ebcz)]
set   phib  [expr acos($cosb)]
set   xsinb [expr 1.0 / sin($phib)]

set   dccbx [expr $ecdy * $ebcz - $ecdz * $ebcy]
set   dccby [expr $ecdz * $ebcx - $ecdx * $ebcz]
set   dccbz [expr $ecdx * $ebcy - $ecdy * $ebcx]

set   cosc  [expr -($ecdx * $ebcx + $ecdy * $ebcy + $ecdz * $ebcz)]
set   phic  [expr  acos($cosc)]
set   xsinc [expr 1.0 / sin($phic)]

#        torsional angle   

set   abcdx [expr -( $abbcy * $dccbz - $abbcz * $dccby )]
set   abcdy [expr -( $abbcz * $dccbx - $abbcx * $dccbz )]
set   abcdz [expr -( $abbcx * $dccby - $abbcy * $dccbx )]

set   temp  [expr $abcdx * $ebcx + $abcdy * $ebcy + $abcdz * $ebcz]

set   signum  1.0
if { $temp < 0.0 } {set signum  -1.0}

set   cosdel [expr -($abbcx * $dccbx + $abbcy * $dccby + $abbcz * $dccbz) *  $xsinb * $xsinc]

    if { $cosdel < -1.0} {set cosdel  -1.0}
    if { $cosdel >  1.0} {set cosdel   1.0}

set angle [expr $signum * acos ($cosdel)]

    return $angle
}

gomPrint "loading macro 'p_rmsd seg res atom scale fit/nofit' ..."
#
# this procedure will calculate the root mean square displacement
# for a set of atoms and plot them as ellipsoids
#
#  Leif Laaksonen 1996
#
#  Usage:
#
#          p_rmsd Segment Residue Atom Scale fit/nofit
#
#
proc p_rmsd {seg res atm scale fit} {
puts "Root mean square fluctuation for: $seg:$res:$atm $fit"
set numat [show numatoms]
if {$scale == ""} {set scale 1.0}
puts "Will use scale factor: $scale"
#check atoms
if {$numat < 1} {
   puts "no atoms defined"
   return 1
}
puts "Number of atoms: $numat"
#check trajectory
set frames [show numframes]
if {$frames < 1} {
   puts "no trajectory defined"
   return 1
}
puts "Number of frames: $frames"

# calculate RMSD
eval calculate rmsf $seg $res $atm $fit
#
set rmsflu [show rmsfl length]
# reset the sphere plot
plot -sphere
for {set i 1} {$i <= $rmsflu} {incr i} {
set temp   [show rmsflu $i]
scan $temp "%f %f %f %f" rad2 xf2 yf2 zf2
set aindex [show rmsfl atomin $i]
set  temp1 [show atom colour $aindex]
scan $temp1 "%f %f %f" rc gc bc
set atomc  [show atom coordinate $aindex]
scan $atomc "%f %f %f" xc yc zc
#

plot sphere $xc $yc $zc [expr $scale * sqrt($rad2)] "$rc $gc $bc" \
            [expr sqrt($xf2)] [expr sqrt($yf2)] [expr sqrt($zf2)] append

}
puts "Job is done..."
return 0
}
#
# Make the conversion from coordinates to angles
#
proc gomVector2Angles {vx vy vz px py pz} {

     set x  [expr $px - $vx];
     set y  [expr $py - $vy];
     set z  [expr $pz - $vz];
     set r  [expr sqrt($x * $x + $y * $y + $z * $z)];
     set x  [expr $x / $r];
     set y  [expr $y / $r];
     set z  [expr $z / $r];

     set beta  [expr atan2($y,$x)];
     set r     [expr $y * sin($beta) + $x * cos($beta)];
     set alpha [expr atan2($r,$z)];

	 return "$alpha $beta"
}


#
# convert a colour given as {float float float} to hex
#
proc lulColourFloat2Hex {red green blue} {

     set Ired   [expr round($red   * 255.0)]
     set Igreen [expr round($green * 255.0)]
     set Iblue  [expr round($blue  * 255.0)]

     return [format "#%02x%02x%02x" $Ired $Igreen $Iblue]
}
#
# convert a colour from hex to {float float float}
#
proc lulColourHex2Float {colour} {
    set colour [string trim $colour]
    set fmt    "%[expr ([string length $colour] - 1)/3]s"
    scan $colour "#$fmt$fmt$fmt" red green blue
    set factor "1.0/0x[regsub -all . $red f]"
    return [list \
		[expr 0x$red   * $factor] \
		[expr 0x$green * $factor] \
		[expr 0x$blue  * $factor] \
		]
}

#
# Optimized version of [llength [lulExpandIndexRangeList $list]]
#
proc lulCalculateIndexRangeListLength { list } {
    set length 0

    foreach range [split $list ,] {
	set iStart [lindex [split $range -]   0]
	set iStop  [lindex [split $range -] end]
	
	set length [expr $length+($iStop-$iStart+1)]
    }

    return $length
}
#
# convert from 1-3,5-7 to {1 2 3 5 6 7}
#
proc lulExpandIndexRangeList { list } {
    set indexlist {}

    foreach range [split $list ,] {
	set iStart [lindex [split $range -]   0]
	set iStop  [lindex [split $range -] end]
	for {set i $iStart} {i<=$iStop} {incr i} {lappend indexlist $i}
    }

    return $indexlist
}
#
# convert from {1-3,*O*,5-7 3-6} to {1-7,*O*}
#
proc lulJoinAtomList { lists } {
    set indexlist {}
    
    foreach list $lists {
	foreach part [split $list ,] {
	    switch -regexp $part {
	    {^[0-9]+(-[0-9]+)?$} {lappend indexlist   $part}
	    default              {set     patterns($part) 1}
	    }
	}
    }

    set list {}
    set iStart 0
    set iStop  -1
    foreach range [lsort -dictionary $indexlist] {
	puts $range
	set rStart [lindex [split $range -]   0]
	set rStop  [lindex [split $range -] end]
	if {$rStart > $iStop + 1} {
	    if {$iStart > 0} {
		if {$iStop > $iStart} {lappend list $iStart-$iStop} else {lappend list $iStart}
	    }
	    set iStart $rStart
	}
	if {$rStop  > $iStop} {set iStop $rStop}
    }
    if {$iStart > 0} {
	if {$iStop > $iStart} {lappend list $iStart-$iStop} else {lappend list $iStart}
    }

    return [join [concat $list [array names patterns]] ,]
}

############################# run #################################
#
# run command
#
# PROC
proc lulRun { Input } {

     global env

set Params    [string trim $Input]
set Elements  [llength $Params]

if {$Elements < 1} {
   puts "name of program to run is missing"
   return
   }

set Program   [lindex $Params 0]

# run probesurf program
if {$Program == "probesurf"} {

   set InputFile   ""
   set OutputFile  ""
   set Extra       ""

   for {set i 1} {$i < $Elements} {incr i} {

        switch $i {
		1 {set InputFile   [string trim [lindex $Params $i]]}
		2 {set OutputFile  [string trim [lindex $Params $i]]}
		3 {set Extra       [string trim [lindex $Params $i]]}
		}
   }  

   if {[string trim $OutputFile] != ""} {
       if {$Extra != ""} {
       set Program    "probsurf -o$OutputFile -l$Extra \< $InputFile "
       } else {
       set Program    "probsurf -o$OutputFile          \< $InputFile "
       }
   } else {
       if {$Extra != ""} {
       set Program    "probsurf -l$Extra \< $InputFile "
       } else {
       set Program    "probsurf          \< $InputFile "
       }
   }

   set Command [file join $gomEnv(GOM_BIN) $Program]

   puts "Executing command: $Command"

   set f [open "|$Command" r]

   while {![eof $f]} { 
     gets $f Text
     puts $Text
   }
# run ICON8 program
} elseif {$Program == "icon8"} {
} else {
  puts "Unknown program name: $Program"
  return 1
}

return

}

#################################################################
# PROC make GAMESS input file
proc lulMakeGAMESSinput { Which File } {

   if {$Which < 1 || $Which > [show molstruc]} {
      gomError "structure index out of range"
      return
   }

   set def_title  "Default title for: GAMESS"

   set File_p [open $File w]
   if {$File_p == ""} {
       catch {lulErrorDialog {ERROR: can't open GAMESS input file for writing!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open GAMESS input file for writing!"
          }
       return
   }


# calculate first number of hydrogens and heavy atoms 

   set nh 0
   set na 0

   set TotAtoms [show numatoms $Which]

   for {set i  1} {$i <= $TotAtoms} {incr i} {

# apply the display mask...
      if {[show atom displaystate $i $Which] == 0} continue

      set AtomName [show atom atomname $i $Which]

      if {[string match H* [show atom atomname $i $Which]]} {
        incr nh
      } else {
        incr na}
   }

   puts "\n=> Writing RAW GAMESS input to: $File "
   puts "\n\n **** GAMESS input generator"
   puts " Number of heavy atoms    : $na"
   puts " Number of hydrogen atoms : $nh"

# start writing input file 
   puts $File_p " \$CONTRL SCFTYP=RHF RUNTYP=ENERGY    \$END"
   puts $File_p " \$BASIS GBASIS=STO NGAUSS=3 \$END"
   puts $File_p " \$DATA"
   puts $File_p "$def_title"
   puts $File_p "C1"

   for {set i 1} {$i <= $TotAtoms } {incr i} {

# apply the display mask...
   if {[show atom displaystate $i 1] == 0} continue

   scan [show atom coordinates   $i $Which] "%f %f %f" xc yc zc

   puts $File_p "[show atom atomname      $i $Which] \
                 [show atom nuclearcharge $i $Which] \
                 [format "%.5f" $xc] [format "%.5f" $yc] [format "%.5f" $zc]"
   }

   puts $File_p " \$END"

   close $File_p

   return


}

#################################################################
# PROC make USER input file
proc lulMakeUSERinput { Which File } {


   if {$Which < 1 || $Which > [show molstruc]} {
      gomError "structure index out of range"
      return
   }

   set File_p [open $File w]
   if {$File_p == ""} {
       catch {lulErrorDialog {ERROR: can't open USER input file for writing!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open USER input file for writing!"
          }
       return
   }

#####################################################################################
# put your user code here
#####################################################################################

   gomPrint "USER own output generator"

   return


}

#################################################################
# PROC make MOPAC6 input file
proc lulMakeMOPACinput { Which File } {


   if {$Which < 1 || $Which > [show molstruc]} {
      gomError "structure index out of range"
      return
   }

   set def_title  "Default title for: MOPAC6"

   set File_p [open $File w]
   if {$File_p == ""} {
       catch {lulErrorDialog {ERROR: can't open MOPAC6 input file for writing!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open MOPAC6 input file for writing!"
          }
       return
   }


# calculate first number of hydrogens and heavy atoms 

   set nh 0
   set na 0

   set TotAtoms [show numatoms $Which]

   for {set i  1} {$i <= $TotAtoms} {incr i} {

# apply the display mask...
      if {[show atom displaystate $i $Which] == 0} continue

      set AtomName [show atom atomname $i $Which]

      if {[string match H* [show atom atomname $i $Which]]} {
        incr nh
      } else {
        incr na}
   }

   puts "\n=> Writing RAW MOPAC6 input to: $File "
   puts "\n\n **** MOPAC6 input generator"
   puts " Number of heavy atoms    : $na"
   puts " Number of hydrogen atoms : $nh"

# start writing input file 
   puts $File_p "AM1 GRAPH"
   puts $File_p "$def_title"
   puts $File_p "  "

   for {set i 1} {$i <= $TotAtoms } {incr i} {

# apply the display mask...
   if {[show atom displaystate $i $Which] == 0} continue

   scan [show atom coordinates   $i $Which] "%f %f %f" xc yc zc

   puts $File_p "[show atom atomname      $i $Which]  \
                 [format "%.5f" $xc] 0 [format "%.5f" $yc] 0 [format "%.5f" $zc] 0"
   }

   close $File_p

   return


}

############################# Change directory ##################
# PROC
# Change directory command
proc lulChangeDirectory { Input } {

  set Directory [file dirname $Input]

  if {$Directory != ""} {
       cd $Directory
  }
}

#################################################################
# PROC
# Input file name:                          FileName
# Read new structure or append to old list: Action (= 0 new  != 0 append)
proc lulReadUSERCoordinates {FileName Action} {

#
# Read USER (PDB look-a-like) coordinate file 
#
# open file ...
   set File_p [open $FileName r]
   if {$File_p == ""} {
       catch {lulErrorDialog {ERROR: can't open USER file!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open USER file!"
          }
       error "ERROR: can't open USER file!"
       return
       }

# do the job ...
   puts "Reading USER coordinate file '$FileName'..."

   set NumAtoms [lulGetUSERAtoms $File_p]

   if {!$Action} {
    define structure "$FileName" $NumAtoms new
   } else {
    define structure "$FileName" $NumAtoms append
   }

   set StrucNum [show molstructures]

   set i 0

  while {![eof $File_p]} { 

    gets $File_p Text

    set Option [lindex $Text 0]

    if {([string trim $Text] == "") || ($Option != "ATOM")} {
         continue}

        incr i

        set Dummy  [lindex $Text 1]
        set Atm    [lindex $Text 2]
        set Res    [lindex $Text 3]
        set ResN   [lindex $Text 4]
        set Seg    "USER"
        set Xc     [lindex $Text 5]
        set Yc     [lindex $Text 6]
        set Zc     [lindex $Text 7]
        set Ch     "0.0"

# put the values into the structure ...
     define atom    label    $Atm        $i $StrucNum
     define residue label    $Res        $i $StrucNum
     define segment label    $Seg        $i $StrucNum
     define atom coordinates $Xc $Yc $Zc $i $StrucNum
     define atom charge      $Ch         $i $StrucNum
     define atom resnumber   $ResN       $i $StrucNum

    }
# close file and return ...
     puts Done!
     close $File_p

}

#################################################################
proc lulGetUSERAtoms { File_p } {

  set i 0

# rewind it to be able to start reading from beginning...
  seek $File_p 0 start

  while {![eof $File_p]} { 

    gets $File_p Text

    set Option [lindex $Text 0]

    if {([string trim $Text] == "") || ($Option != "ATOM")} {
         continue}

    incr i
  }

# rewind it to be able to start reading from beginning...
  seek $File_p 0 start
  return $i
}

#################################################################
# PROC
proc lulMakeUserCoordinates {FileName} {

     error "ERROR: The code has to be defined!"
}

#################################################################
#   NEW COMMANDS
# * run command
proc run args {

    global env
    global gomEnv

if {[lindex $args 0] == ""} {
	    lulErrorDialog "ERROR - please define program name to run!"
		return
}

set program(probesurf) prob*
set program(xvibs)     xvib*

set Program [lindex $args 0]

# RUN Probesurf program
if {[string match $program(probesurf) $Program]} {

# chack that there is an input file given...

  set InputFile [lindex $args 1]
  if { $InputFile == ""} {
        lulErrorDialog "ERROR - Probesurf input file name missing!"
		return
  }

# chack that the file really exists
  if {![file exists $InputFile]} {
        lulErrorDialog "ERROR - input file '$InputFile' does not exist!"
		return
  }

# check if output grid file name is given ...
  set GridFile [lindex $args 2]

  if {$GridFile == ""} {
    set GridFile "[file join [file dirname $InputFile] probesurf.plt]"
  }
    set Runit    "probsurf -o$GridFile   \< $InputFile "

  set Command [file join $gomEnv(GOM_BIN) $Runit]
# RUN Xvibs
} elseif {[string match $program(xvibs) $Program]} {

# chack that there is an input file given...

  set InputFile [lindex $args 1]
  if { $InputFile == ""} {
        lulErrorDialog "ERROR - Xvibs input file name missing!"
		return
  }

# chack that the file really exists
  if {![file exists $InputFile]} {
        lulErrorDialog "ERROR - input file '$InputFile' does not exist!"
		return
  }

# check modes to save ...
  set Modes [lindex $args 2]
  if { $Modes == ""} {
        lulErrorDialog "ERROR - modes to save are missing \{all \| 1..(3*N-6)\}"
		return
  }

# palindrome on?
  set Palindrome [lindex $args 3]

  if {$Palindrome != "" && [string match $Palindrome* palindrome]} {
    set Runit "xvibs $InputFile $Modes palindrome"
  } else {
    set Runit "xvibs $InputFile $Modes"
  } 

  set Command [file join $gomEnv(GOM_BIN) $Runit]
# DEFAULT = ERROR option!
} else {
        lulErrorDialog "ERROR - can't resolve command '$args'"
		return

}

# finally run the program

puts $Command

set f [open "|$Command" r]

  set i 1

  while {![eof $f]} { 
    gets $f Text
    update idletasks
    puts $Text
    incr i
  }

catch {close $f}
puts "Job Dome!"
}

#################################################################
# PROC
# Load a plugin
proc lulLoadPlugin { Name } {
    global tcl_platform
    global gomEnv

    set libext [info sharedlibextension]
    if {[file exists [file join $gomEnv(GOM_PLUGINS) ${Name}${libext}]]} {
	load [file join $gomEnv(GOM_PLUGINS) ${Name}${libext}]
    } else {
	load [file join $gomEnv(GOM_PLUGINS) lib${Name}${libext}]
    }
}

#################################################################
# PROC
# Send commands to the gnuplot program
#
proc lulSend2GnuPlot {String} {

   global gomEnv
   global gnuplotexe

   set Command [file join $gomEnv(GOM_BIN) $gnuplotexe]
#
# open temp file into the temp directory
#
#   set filename [file join $gomEnv(GOM_TEMP) tmp[clock seconds].gpl]

#   set fi [open $filename w]
#   puts $fi $String
#   puts $fi "raw_input(\"press ENTER to continue\")"
#   close $fi

# done
    
   puts "Executing command: $Command"

   set f [open "| $Command " w]
   fconfigure $f -buffering line

   puts $f $String
 
   close $f

}
#################################################################
# PROC
# exe this proc after reading a frame in the trajectory control if
# variable 'lulPostReadFrameTrigger' is set to '1'
#
proc lulPostReadFrame { i } {
}
#################################################################

########################### openmol #############################
#if {[file exists [file join $gomEnv(GOM_DATA) openmol.tcl]]} {
#puts "loading 'openmol.tcl'..."
#source [file join $gomEnv(GOM_DATA) openmol.tcl]
#}
#################################################################

########################### AutoDock ############################
#if {[file exists [file join $gomEnv(GOM_DATA) autodock.tcl]]} {
#puts "loading 'autodock.tcl'..."
#source [file join $gomEnv(GOM_DATA) autodock.tcl]
#}

########################### Gaussian ############################
#if {[file exists [file join $gomEnv(GOM_DATA) gaussianxx.tcl]]} {
#puts "loading 'gaussianxx.tcl'..."
#source [file join $gomEnv(GOM_DATA) gaussianxx.tcl]
#}

########################### Utility ############################
#if {[file exists [file join $gomEnv(GOM_DATA) utility.tcl]]} {
#puts "loading 'utility.tcl'..."
#source [file join $gomEnv(GOM_DATA) utility.tcl]
#}
# use tcllib
#   lappend auto_path [file join $gomEnv(GOM_ROOT) lib tcllib1.4]
#   package require tcllib

########################### DDE #################################
if {[regexp -nocase "windows" $tcl_platform(platform)]} {
	package require dde 1.2
	dde servername gopenmol
}

#############################################################
#
# source automatically all files in the {auto,init}tcl directories
#
puts "Sourcing necessary tcl files in the '[file join $gomEnv(GOM_DATA) pkgtcl]' directory ..."
set i [file join $gomEnv(GOM_DATA) pkgtcl]
foreach i [glob -directory $i -nocomplain -type d * */* */*/* */*/*/*] {
    puts "    $i"
    lappend auto_path $i
}
package require gom::general

puts "Sourcing all tcl files in the '[file join $gomEnv(GOM_DATA) autotcl]' directory ..."
foreach i [expr [show graphics] ? {"general gui"} : {"general"}] {
    set i [file join $gomEnv(GOM_DATA) autotcl $i]
    foreach i [glob -directory $i -nocomplain *.tcl] {
	puts "    $i"
	if { [catch {source $i} errRet] } {
	    gomError "Parse error in file '$i':\n$errRet"
	}
    }
}

########################### GUI #################################
#
# read the GUI part now ...
#
set lulGraphicsAvailable [show graphics]

if {$lulGraphicsAvailable} {
   puts "    Graphics is available..."
# use BWidget
#   lappend auto_path [file join $gomEnv(GOM_ROOT) lib BWidget-1.6.0]
    package require BWidget
    source [file join $gomEnv(GOM_DATA) gopenmol_guirc.tcl]
    package require gom::gui
    package require gom::Plugins
} else {
   puts "    Graphics is not available..."
}
#################################################################

#
# Change to the last visited directory
#
set tPath [file split $env(HOME)]
set Path ""
foreach i "$tPath" {
  set Path [file join $Path $i]
}
if {[file exists [file join $Path .gopenmolrc.tcl]]} {

   catch {source [file join $env(HOME)/.gopenmolrc.tcl]}
   puts "Sourcing file: [file join $env(HOME)]/.gopenmolrc.tcl"
}
#
##################################################################
#
proc lulShutdownProcedure { } {
#
   global env

set tPath [file split $env(HOME)]
set Path ""
foreach i "$tPath" {
  set Path [file join $Path $i]
}

set File [open [file join $Path .gopenmolrc.tcl] w]
puts $File "# [clock format [clock seconds] -format {%d/%m/%y at %H:%I}]"
set tPath [file split [pwd]]
set Path ""
foreach i "$tPath" {
  set Path [file join $Path $i]
}
puts $File "cd \"$Path\""
puts $File "# End"
close $File
}
# this is the "trigger" string!
set lulProgramShutdownString "lulShutdownProcedure ;#Running shutdown script..."

# .... daiga daiga duu

puts "Done!"
