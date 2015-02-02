#
# Take the molecular 3D coordinates and transform them into 2D
# 
# Leif Laaksonen 1997
#
#  xp(i) =  y(i) - x(i)
#  yp(i) = (x(i) + y(i) - 2 * z(i))/2
#
#
proc lulMol3Dto2D { } {

  if {[show molstructures] < 1} {
      lulErrorDialog {ERROR: no structure available. Read a structure first!}
      return
      }

  set NumAtoms [show numatoms 1]

puts "Calculating the new x and y coordinates ..."

      for {set i  1} { $i <= $NumAtoms} { incr i} {

      scan [show atom coord $i 1] "%f %f %f" xc yc zc

      set xp($i) [expr $yc - $xc]
      set yp($i) [expr ($xc + $yc - 2.0 * $zc)/2.0]

      }

puts "Calculating the translation to center system ..."

# sum up and center the system...
      set sumx 0.0
      set sumy 0.0

      for {set i 1} { $i <= $NumAtoms} {incr i} {

      set sumx [expr $sumx + $xp($i)] 
      set sumy [expr $sumy + $yp($i)] 

      }

puts "Making the translation ..."

      set sumx [expr $sumx / $NumAtoms]
      set sumy [expr $sumy / $NumAtoms]
 
      for {set i 1} { $i <= $NumAtoms} {incr i} {

      define atom coord [expr $xp($i) - $sumx] [expr $yp($i) - $sumy] 0.0 $i 1

      }
puts "Job Done!"
}