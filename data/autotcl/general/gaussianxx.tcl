##########################################################################
# NAMESPACE GaussianXX
#
namespace eval lulGaussian {

     variable lulNumAtoms
     variable lulGaussianLogFileName

############################################################################
# PROC
proc ExtractXMOLfileFromGaussianLogfile { FileName } {

   global gomAtom_SN

   if {$FileName == ""} {
      gomError "file name is not supplied"
      return
   }

   set f [open $FileName r]

   if {$f == ""} {
       catch {lulErrorDialog {ERROR: can't open Gaussian log file for reading!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open Gaussian log file for reading!"
          }
       return
   }

   set FileNameOut "[file rootname $FileName].xmol"

   set fo [open $FileNameOut w]

   if {$fo == ""} {
       catch {lulErrorDialog {ERROR: can't open XMOL output file for writing!}} errRet
       if {$errRet != ""} {
          puts "ERROR: can't open XMOL output file for writing!"
          }
       return
   }


# change to current directory
   lulChangeDirectory $FileName

   puts "Writing XMOL coordinate file '$FileNameOut' ..."

# scan the Gaussian log file to get the atoms ...

    set j 0
    while {![eof $f]} { 
      gets $f Text
      update idletasks
# look for "Standard orientation:"
      if {[string match "*Standard orientation:*" $Text]} {
# 4 text lines first
          gets $f Text
          gets $f Text
          gets $f Text
          gets $f Text
          set  i 0
          incr j
          while {![eof $f]} { 
           gets $f Text
           set Text [string trim $Text]
           if {($Text == "") || ([string match "*----------*" $Text])} break
           incr i
           set atom($i) [lindex $gomAtom_SN([lindex $Text 1]) 0]
           set xc($i)   [lindex $Text 2]
           set yc($i)   [lindex $Text 3]
           set zc($i)   [lindex $Text 4]
          }
          puts $fo $i
          puts $fo "Coordinate set nr: $j"
          for {set k 1} {$k <= $i} {incr k} {
            puts $fo "$atom($k) $xc($k) $yc($k) $zc($k)"  
          }
      }
    }


   close $f
   close $fo

   puts "Done!"

   return
# end of proc
}

}
