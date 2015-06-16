##############################################################################
#                            Copyright (c) 2003 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

top_builddir = /Users/shkljo14/src/gOpenMol-3.00-osx/src
builddir     = /Users/shkljo14/src/gOpenMol-3.00-osx/src/readwrite

set_subdir   = subdir=`echo '$(builddir)/' | sed 's:$(top_builddir)/::'`

# Rule to update existing files.
*: FORCE
	@$(set_subdir); \
	case '$@' in \
	Makefile) ;; \
	*) cd '$(top_builddir)' && $(MAKE) $(MFLAGS) "$${subdir}$@" ;; \
	esac

delegate     = $(set_subdir); target='$@'; \
               case "$$target" in *.*) target="$${subdir}$$target" ;; esac; \
               cd '$(top_builddir)' && $(MAKE) $(MFLAGS) "$$target"

# Rules to create files.
.SUFFIXES:

%:
	@$(delegate)

.DEFAULT:
	@$(delegate)

FORCE:

.PHONY: FORCE
