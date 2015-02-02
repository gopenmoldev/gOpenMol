##############################################################################
#                          Copyright (c) 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded 2004 by: Eero Häkkinen
##############################################################################

package provide gom::PackageUtilities 1.0

namespace eval gom {
namespace eval PackageUtilities {

proc SourcePackageFiles { subdir } {
    set code ""
    foreach file [ListPackageFiles $subdir] {
	append code "source [list $file]\n"
    }
    namespace eval :: $code
}

proc _GetPackageSubDirectory { subdir } {
    global gomEnv
    return [eval file join \
		[file split $gomEnv(GOM_DATA)] \
		pkgtcl \
		[file split $subdir]]
}

proc ListPackageFiles { subdir } {
    set files [list]
    foreach file [glob \
		      -directory [_GetPackageSubDirectory $subdir] \
		      -nocomplain \
		      *.tcl] {
	switch [string tolower [file tail $file]] {
	    init.tcl     {}
	    pkgindex.tcl {}
	    default      { lappend files $file }
	}
    }
    return $files
}

namespace export {[A-Z]*}

}

namespace import PackageUtilities::*

}
