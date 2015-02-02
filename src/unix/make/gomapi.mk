##############################################################################
#                           Copyright (c) 2003 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

# Create API header file depedencies.
gomapi: FORCE
gomapi graphics/gomapi.mk:
	{ cd '$(incdir)' && ls *.h | \
	    sed -e 's,.*,include/gomext/& include/gomlib/gom&: $$(incdir)/&,' \
		-e 's, \(include/gomlib/gom\)gom, \1,' \
		-e 'p' \
		-e 's/:.*//' \
		-e 's/ /: /'; } > graphics/gomapi.mk

# Create API headers. Do not recreate anything which depend on this file.
apiheaders: FORCE
apiheaders graphics/stamp-gomapi:
	$(MKDIR_P) include/gomext
	$(MKDIR_P) include/gomlib
	$(TCLSH) $(incdir)/gomext.tcl -o include
	echo "created" >graphics/stamp-gomapi
	$(MAKEOLD) graphics/stamp-gomapi

# graphics/stamp-gomapi target changes the same files so we need
# the prerequisite.
include/gomlib/*.h: graphics/stamp-gomapi
	file=`echo '$@' | sed 's,.*[\\\\/]gom,,'`; \
	test -f '$(incdir)'/"$$file" || file="gom$$file"; \
	$(TCLSH) $(incdir)/gomext.tcl -o include "$$file"

# graphics/stamp-gomapi target changes the same files so we need
# the prerequisite.
include/gomext/*.h: graphics/stamp-gomapi
	file=`echo '$@' | sed 's,.*[\\\\/],,'`; \
	case $$file in "*.h") exit ;; gom*) ;; *) file="gom$$file" ;; esac; \
	sed -n 's/^GOPENMOLAPI.*gom_\([^ ()]*\)(.*/#define gomp_\1 gom_\1/p' \
		< "include/gomlib/$$file" > $@

# List symbols to be exported.
# Note: it's OK to have prerequisite "include/gomext/*.h"
# (see case statement above).
graphics/gomapi.exp: include/gomext/*.h
	{ cat include/gomext/*.h | sed 's/^.*gom_/gom_/'; echo gom_main; } >$@

clean: clean-gomapi clean-stamps FORCE

distclean: clean-mks FORCE

clean-gomapi: FORCE
	$(RM) include/gomext/*.h
	$(RM) include/gomlib/gom*.h
	$(RM) graphics/gomapi.exp

clean-mks: clean-mk-gomapi FORCE

clean-mk-gomapi: FORCE
	$(RM) graphics/gomapi.mk

clean-stamps: FORCE
	$(RM) graphics/stamp-*

FORCE:

.PHONY: FORCE
