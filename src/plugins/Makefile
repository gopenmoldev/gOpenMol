##############################################################################
#                        Copyright (c) 2002 - 2003 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

build-all: build FORCE ;

all install install-local uninstall: FORCE
	@for dir in *; do \
	    if test -f "$$dir/Makefile"; then \
		( cd "$$dir" && $(MAKE) $(MFLAGS) $@ \
			gomroot='$(gomroot)' \
			gomdataroot='$(gomdataroot)' ) || exit 1;\
	    fi ;\
	done

build compile clean clean-libs clean-objs: FORCE
	-@for dir in *; do \
	    if test -f "$$dir/Makefile"; then \
		( cd "$$dir" && $(MAKE) $(MFLAGS) $@ ) ;\
	    fi ;\
	done

distclean maintainer-clean: FORCE
	-@for dir in *; do \
	    if test -f "$$dir/Makefile"; then \
		( cd "$$dir" && $(MAKE) $(MFLAGS) $@ ) ;\
	    fi ;\
	done
	rm -f config.mk

FORCE:

.PHONY: FORCE
