load d:\\queens\\surfarea\\debug\\surfarea.dll
#
# save polygon entries for external program
#
set File c:/temp/polygons.txt
   set File_p [open $File w]
   if {$File_p == ""} {
    exit
   }

contour method save
display
set numPoly [show contour polygon entries]
puts $numPoly
set contcrap [show contour cell 1]
set cxmin [lindex $contcrap 0]
set cxmax [lindex $contcrap 1]
set cymin [lindex $contcrap 2]
set cymax [lindex $contcrap 3]
set czmin [lindex $contcrap 4]
set czmax [lindex $contcrap 5]
set ncx [lindex $contcrap 6]
set ncy [lindex $contcrap 7]
set ncz [lindex $contcrap 8]
set contptr [show address contour 1]
sacontour 1 $cxmin $cxmax $cymin $cymax $czmin $czmax $ncx $ncy $ncz $contptr
set sArea 0
set pArea 0
set dArea 0
#
# loop
#
#set numPoly 3
for {set i 1} {$i <= $numPoly} {incr i} {

  set Text [show contour polygon data $i]
# 1
  set x1   [lindex $Text 0]
  set y1   [lindex $Text 1]
  set z1   [lindex $Text 2]
  set u1   [lindex $Text 3]
  set v1   [lindex $Text 4]
  set w1   [lindex $Text 5]
# 2
  set x2   [lindex $Text 6]
  set y2   [lindex $Text 7]
  set z2   [lindex $Text 8]
  set u2   [lindex $Text 9]
  set v2   [lindex $Text 10]
  set w2   [lindex $Text 11]
# 3
  set x3   [lindex $Text 12]
  set y3   [lindex $Text 13]
  set z3   [lindex $Text 14]
  set u3   [lindex $Text 15]
  set v3   [lindex $Text 16]
  set w3   [lindex $Text 17]

puts $File_p "smooth_triangle \{  \n\
<  $x1,   $y1, $z1 > <    $u1,    $v1,  $w1 > \n\
<  $x2,   $y2, $z2 > <    $u2,    $v2,  $w2 > \n\
<  $x3,   $y3, $z3 > <    $u3,    $v2,  $w3 > \n\
\}"
  set qq [triangles $x1 $y1 $z1 $x2 $y2 $z2 $x3 $y3 $z3 $u1 $v1 $w1 $u2 $v2 $w2 $u3 $v3 $w3 0] 
  puts $qq
  puts $i
  set temp [lindex $qq 0]
  set pArea [expr {$pArea+$temp}]
  set temp [lindex $qq 1]
  set sArea [expr {$sArea+$temp}]
  set temp [lindex $qq 3]
  set dArea [expr {$dArea+$temp*$temp}]
}
set dArea [expr {sqrt($dArea)}]
puts $pArea
puts $sArea
puts $dArea
close $File_p
contour method direct
display

puts "Job done!"

