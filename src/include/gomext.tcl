#! /bin/sh
# The next line is part of Tcl comment but is a Posix shell command. \
exec tclsh "$0" "$@"

proc getapifile { header dstdir } {
    # Our header file in public gomlib directory must begin with "gom"
    # to avoid header file collision.
    switch -glob $header {
	gom*    {return [file join $dstdir gomlib $header]}
	default {return [file join $dstdir gomlib "gom$header"]}
    }
}

proc getextfile { header dstdir } {
    return [file join $dstdir gomext $header]
}

proc getmtime { file } {
    set mtime 0
    catch { set mtime [file mtime $file] }
    return $mtime
}

proc process_header_file { header srcdir dstdir } {
    set FNsrc [file join $srcdir $header]
    set FNapi [getapifile $header $dstdir]
    set FNext [getextfile $header $dstdir]
    set srcmtime [getmtime $FNsrc]
    if { [getmtime $FNapi] > $srcmtime &&
	 [getmtime $FNext] > $srcmtime } {
	puts "$header is up to date"
	return
    }
    regsub -all {[^A-Za-z0-9]} "INC_GOPENMOL_[string toupper $header]" \
	    {_} incldef
    set precmd   {
	puts "processing $header ..."
	set Fapi [open $FNapi w]
	set Fext [open $FNext w]
	puts $Fapi "#ifndef $incldef"
	puts $Fapi "#define $incldef"
	puts $Fapi ""
	puts $Fapi "#include <gomdefs.h>"
	puts $Fapi ""
    }
    set postcmd  {
	file delete $FNapi $FNext
    }
    set type ""
    set input [open $FNsrc r]
    while {[gets $input line] >= 0} {
	switch -glob -- ${type}${line} {
	    {* GOMAPI END *} {
		puts $Fapi ""
		set type ""
	    }
	    "\nGOMAPI\n*" {
		# Create #define macros which change "gomp_" functions
		# to "gom_" functions.
		if {[regexp {^extern.*gomp_([A-Za-z0-9_]*)} $line dummy func]} {
		    puts $Fext "#define gomp_$func gom_$func"
		}
		# Public entries are prefixed with "gom_" instead of "gomp_".
		set gomapi $line
		regsub -all {gomp_} $gomapi {gom_} gomapi
		# Change extern to GOPENMOLAPI for Windows compatibility.
		regsub {^extern } $gomapi {GOPENMOLAPI } gomapi
		# Use <> style include directives.
		if {[regexp "(\#\[ \t\]*include\[ \t\]*)\"(.*)\"" \
			$line dummy prefix file]} {
		    set gomapi "$prefix<$file>"
		}
		puts $Fapi $gomapi;
	    }
	    {* GOMAPI BEGIN *} {
		eval $precmd
		set precmd ""
		set postcmd {
		    puts  $Fapi "#endif /* $incldef */"
		    puts  $Fext "extern void dummy_func(void);"
		    close $Fapi
		    close $Fext
		    puts "    $FNapi is updated"
		    puts "    $FNext is updated"
		}
		set type "\nGOMAPI\n"
	    }
	}
    }
    eval $postcmd
}

set srcdir [file dirname $argv0]

if { "-o" == [lindex $argv 0] } {
    set dstdir [lindex $argv 1]
    set argv [lrange $argv 2 end]
} else {
    set dstdir $srcdir
}

if { [llength $argv] == 0 } {
    set cwd [pwd]
    cd $srcdir
    set argv [glob *.h]
    cd $cwd
}

foreach header $argv {
    process_header_file $header $srcdir $dstdir
}
