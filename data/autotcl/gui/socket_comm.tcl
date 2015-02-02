##########################################################################
# Socket communication code
# Leif Laaksonen 2002 
# Based on code from the Practical Programming in Tcl/Tk
# by Brent B. Welch
#

##########################################################################
# NAMESPACE SocketComm
#
namespace eval lulSocketComm {

   variable socket 2309
   variable socketport 0
   variable echo
   variable Commands "file coordinates quit"
   variable status   0
   variable filename
   variable filetype
   variable fileappend
   variable file_p

#######################################################################
proc Server {port} {
#######################################################################
	variable echo
    variable socket
    variable socketport

    if {($port == $socketport) || ($port == 0)} {
     close $echo(main)
     unset echo(main)
     set socketport 0
    }
	set echo(main) [socket -server lulSocketComm::Accept $port]
	set socketport $port
}
#######################################################################
proc Accept {sock addr port} {
#######################################################################
	variable echo
    variable socket
    global   gomAllowConnection

    set allow   0
    set socket $sock
    foreach i $gomAllowConnection {
      if {$i == $addr} {
        puts "Connection: $sock from '$addr' port '$port'"
        set allow 1
        break
      }
    }
    if {!$allow} {
      puts "No connection allowed from '$addr'"
      return
    }
#    set socket $sock
	set echo(addr,$sock) [list $addr $port]
	fconfigure $sock -buffering line
	fileevent $sock readable [list lulSocketComm::Process $sock]
}
#######################################################################
proc Process {sock} {
#######################################################################
	variable echo
    variable status
    variable filename
    variable filetype
    variable fileappend
    variable file_p
    global gomEnv

	if {[eof $sock] || [catch {gets $sock line}]} {
		# end of file or abnormal connection drop
		close $sock
		puts "Close $echo(addr,$sock)"
		unset echo(addr,$sock)
        set status 0
	} else {
# quit: end socket connection
		if {[string compare $line "quit"] == 0} {
            puts $sock "OK \{$line\}"
			# Prevent new connections.
			# Existing connections stay open.
			close $echo(main)
# A file was received
            if {$status == 1} {
              close $file_p
              import coord $filetype $filename $fileappend
              display
            } elseif {$status == 2} {
              display
            }
            set status 0
            return
# file: receiving a file
		} elseif {[string compare [lindex $line 0] "file"] == 0} {
            # A full coordinate file is coming
            set status 1
            set filetype    [lindex $line 1]
            set fileappend  [lindex $line 2]
            set filename "[file join $gomEnv(GOM_TEMP) junk.tmp]"
            set file_p [open $filename w]
            puts $sock "OK \{$line\}"
            return
# coordinates: receiving just the coords
        } elseif {[string compare [lindex $line 0] "coords"] == 0} {
          set status 2
          set atomnr [lindex $line 1]
          set xc     [lindex $line 2]
          set yc     [lindex $line 3]
          set zc     [lindex $line 4]
          set system [lindex $line 5]
          if {$system == ""} {
            set system 1
          }

          define atom coord $xc $yc $zc $atomnr $system
#          display 
        }
		puts $sock "OK:$line"
        puts "Local data: $line"

        if {$status == 0} {
        } elseif {$status == 1} {
# write data to file
          puts $file_p $line
        }
        
	}
}

#######################################################################
proc Client {host port} {
#######################################################################
	set s [socket $host $port]
	fconfigure $s -buffering line
	return $s
}
#######################################################################
proc ClientSend {host stext} {
#######################################################################

    variable socket

    set s [Client $host $socket]
# send file
    if {[lindex $stext 0] == "sendfile"} {
       
       set filename [lindex $stext 1]
       set filetype [lindex $stext 2]
# open file to read
       set f [open $filename r]
       if {$f == ""} {
         puts "Can't open file '$filename'"
         return
       }
       set i 1
# send first definition
       puts $s "file $filetype" 
       while {![eof $f]} { 
          gets $f Text
          if {[string trim $Text] == ""} {
            break
          }
          puts $s $Text
          puts "Sending:$Text"
          incr i
       }
#close pipe at end of file
       puts "Number of lines sent: $i"
       close $f
# break connection
       puts $s "quit"
    } else {
     puts $s $stext
#    return [gets $s]
     return
    }
}
}

