#! /bin/sh
# <<<<< Tcl comment begins <<<<< \
case $0 in *[\\/]*) \
    dir=`echo "$0" | sed 's,[^\\\\/]*$,,'`; \
    cd "$dir" \
    ;; \
esac; \
case $0 in *.tcl) exec ./install; exit ;; esac; \
@SHLIBPATHVAR@=`pwd`/lib${@SHLIBPATHVAR@+@PATH_SEPARATOR@$@SHLIBPATHVAR@}; \
export @SHLIBPATHVAR@; \
set x bin/tclsh*; exec "$2" ./install; exit; \
# <<<<< Tcl comment ends <<<<<
#
# Make the installation of gOpenMol using Tcl
#

puts " ####################################################"
puts " #                                                  #"
puts " #  This is the gOpenMol installation program       #"
puts " #  Leif Laaksonen 2005                             #"
puts " #  Scarecrow computing                             #"
puts " #                                                  #"
puts " ####################################################"

puts " "

puts "Will inspect your system now ..."

set machine   $tcl_platform(machine)
puts "Your machine is:   $machine"

set os        $tcl_platform(os)
puts "Your OS is:        $os"

set osversion $tcl_platform(osVersion)
puts "Your OSversion is: $osversion"

set platform  $tcl_platform(platform)
puts "Your platform is:  $platform"

# this does not exist on Windows 95
#set systemroot $env(SYSTEMROOT)
#puts "SYSTEMROOT:        $systemroot"
#
#if {$systemroot == ""} {
#     puts "?ERROR - can't get 'SYSTEMROOT'"
#     exit
#}

# check that a TEMP/temp directory is defined ...

  if {[string match -nocase "*windows*" $os]} {

    if { ! [info exists env(TEMP)] && ! [info exists env(temp)] } {
        puts "?ERROR - there is no 'TEMP' or 'temp' directory define. Can't continue"
        exit
    }
  }

# Ok, here we go ...

  if {[string match -nocase "*windows*" $os]} {
       set gopenmol_run "rungOpenMol.bat"
  } else {
       set gopenmol_run "rungOpenMol"
  }

# if this is a WIN95 computer the OpenGL stuff is not included ...
#if { $os == "Windows 95" || $os == "Windows 98"} {
if { $os == "Windows 95" } {
      file copy -force bin/OpenGL/GLU32.DLL     bin/
      file copy -force bin/OpenGL/OPENGL32.DLL  bin/
}

##########################################################################
# Ok, here we go ... start with the run script ...
# write the bat file to run gOpenMol
  puts "Will now create the file 'bin/$gopenmol_run' to run gOpenMol..."

  set id [open bin/$gopenmol_run "w"]
  if {$id == 0} {
       puts "?ERROR - can't open 'bin/$gopenmol_run' for writing"
       exit
  }

  if {[info exist env(PWD)]} {
    # PWD contains a logical path to start up directory of the process.
    set currdir $env(PWD)
  } else {
    # Tcl pwd command gives current physical path.
    set currdir [pwd]
  }

  if {[regexp {[[:space:]]} $currdir]} {
    if {[string match -nocase "*windows*" $os]} {
      set currdir [file attribute $currdir -shortname]
    } else {
      puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      puts "!!!!!  gOpenMol can not be installed in a directory  !!!!!"
      puts "!!!!!  with a name containing space(s)               !!!!!"
      puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      exit
    }
  }

  if {[string match -nocase "*windows*" $os]} {
    set currdir [file native $currdir]

    puts $id "REM bat file to run gOpenMol"
    puts $id "REM contains only essential assignments"
    puts $id "REM the rest will be placed in the 'environment.txt' file"
    puts $id "REM Leif Laaksonen 2003"
    puts $id "set GOM_ROOT=$currdir"
    puts $id "set TCL_LIBRARY=$currdir\\lib\\tcl8.4"
    puts $id "set TK_LIBRARY=$currdir\\lib\\tk8.4"

    if {$os == "Windows NT"} {
      puts $id "if %ERRORLEVEL% LSS 1 goto okay"
    } else {
      puts $id "if ERRORLEVEL 0 goto okay"
    }
	puts $id "echo Hit an error! Please check the environment space"
	puts $id "echo One source of error is a too small environment space"
	puts $id "echo on the Windows 95/98 platform"

	puts $id "pause"
	puts $id "exit"
	puts $id ":okay"
	
    puts $id "\"%GOM_ROOT%\\bin\\gopenmol.exe\" %1 %2 %3 %4 %5 %6 %7 %8 %9"
    puts $id "pause"

  } else {
    puts $id "#! /bin/sh"
    puts $id "# script to run gOpenMol"
    puts $id "# contains only essential assignments"
    puts $id "# the rest will be placed in the 'environment.txt' file"
    puts $id "# Leif Laaksonen 2004"
    puts $id "GOM_ROOT='$currdir'"
    puts $id {# GOM_ROOT=${GOM_ROOT:-$(dirname "$(dirname "$(readlink -f "$(which "$0")" || ( cd "$(dirname "$(which "$0")")" && echo "$(pwd)/dummy" ))")")}}
    puts $id {TCL_LIBRARY=$GOM_ROOT/lib/tcl8.4}
    puts $id {TK_LIBRARY=$GOM_ROOT/lib/tk8.4}
    puts $id {@SHLIBPATHVAR@=$GOM_ROOT/lib${@SHLIBPATHVAR@+@PATH_SEPARATOR@$@SHLIBPATHVAR@}}
    puts $id {export GOM_ROOT TCL_LIBRARY TK_LIBRARY @SHLIBPATHVAR@}
    puts $id {LC_ALL=C xdpyinfo 2>/dev/null |}
    puts $id {	sed '/^number of extensions:/,/:/!d' |}
    puts $id {	grep '^ *GLX$' >/dev/null || use_mesa=yes}
    puts $id {( LC_ALL=C ldd -r "$GOM_ROOT/bin/gopenmol" 2>&1 || echo error: $? ) |}
    puts $id {	grep :. >/dev/null && use_mesa=yes}
    puts $id {if test 'yes' = "$use_mesa"}
    puts $id {then}
    puts $id {    @SHLIBPATHVAR@=$GOM_ROOT/lib/Mesa@PATH_SEPARATOR@$@SHLIBPATHVAR@}
    puts $id {    export @SHLIBPATHVAR@}
    puts $id {fi}
    puts $id {case $# in}
    puts $id {0) exec "$GOM_ROOT/bin/gopenmol"      ;;}
    puts $id {*) exec "$GOM_ROOT/bin/gopenmol" "$@" ;;}
    puts $id {esac}
    file attribute bin/$gopenmol_run -permission 0755
  }

  close $id

  puts "File $gopenmol_run ... Ready!"

puts " "
puts "*****************************"
puts "The installation is complete!"
puts "*****************************"
puts " "
puts "To run gopenmol please run the '$gopenmol_run' script "
puts "in the bin directory"
puts " "
