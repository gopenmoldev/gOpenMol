##############################################################################
#                       Copyright (c) 2003 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

compile: $(OBJECTS) ;

# List all libtool objects.
objects: FORCE
objects graphics/objects.mk:
	{ \
	cd '$(top_srcdir)' && \
	for dir in $(SUBDIRS); do \
	    name=`echo $$dir | sed 's,[\\\\/],_,g'`_OBJECTS; \
	    echo "$$name = \\"; \
	    for file in "$$dir"/*.c; do \
		echo "	$$file \\"; \
	    done | sed 's/\.c \\$$/.lo \\/'; \
	    echo '	$$(end_of_list)'; \
	done; \
	echo "OBJECTS = \\"; \
	echo $(SUBDIRS) | sed 's,/,_,g;s/\([^ ][^ ]*\)/$$(\1_OBJECTS)/g'; \
	echo "OBJLISTS = \\"; \
	echo $(SUBDIRS) | sed 's,\([^ ][^ ]*\),\1/all_lo,g'; \
	} > graphics/objects.mk

# Find all dependency files, create missings ones and include all.
depends: FORCE
depends graphics/depends.mk:
	set -e; \
	{ \
	bld=`pwd`; \
	cd '$(top_srcdir)' && \
	for dir in $(SUBDIRS); do \
	    objs=`echo $$dir | sed 's,[\\\\/],_,g'`_OBJECTS; \
	    echo "$$dir/all_lo: \$$($$objs)"; \
	    for file in "$$dir"/*.c; do \
		file=`echo "\$$file" | sed 's/\\.c$$/.d/'`; \
		test -f "$$bld/$$file" || echo "# dummy" >"$$bld/$$file"; \
		echo "include $$file"; \
	    done; \
	done; \
	} > graphics/depends.mk

# Rule to compile sources to objects.
# API headers are needed during compiling.
$(OBJECTS): graphics/stamp-gomapi
	file=`echo '$@' | sed 's/\\.lo\$$//'`; \
	source='$(top_srcdir)'/"$$file.c"; \
	source="$$source" object="$$file.lo" depfile="$$file.d" libtool=yes \
	$(CCDEPCOMP) \
	$(COMPILE) '$@' "$$source" $(COMPILEFLAGS)

# Libtool doesn't create libtool setup files for reloadable objects.
# Do it here.
# Link all objects in a single directory to a single reloadable object.
# We use file list files, so we can continue even if the linkage fails.
$(OBJLISTS):
	$(RM) $@ $@.*
	ls `echo $@ | sed 's,[^\\\\/]*$$,,'`*.lo >$@
	trap '$(RM) $@.c' 0 1 2 15; \
	echo 'int i = 0;' >$@.c; \
	$(COMPILE) $@.lo $@.c $(COMPILEFLAGS)
	( $(set_show); link='$(LINK)'; \
	dir=`echo $@ | sed 's,[^\\\\/]*$$,,'`; \
	case $$dir in "") . ./$@.lo ;; *) . $@.lo ;; esac; \
	case $$non_pic_object in none) ;; *) \
	set $$link $${dir}$${non_pic_object} `sed 's/.lo$$/.o/' <$@`; \
	$$show "$$@" && "$$@" || exit 1 ;; \
	esac; \
	case $$pic_object in none) ;; *) \
	set $$link $${dir}$${pic_object} -objectlist $$list; \
	$$show "$$@" && "$$@" || exit 1 ;; \
	esac ) && echo $@.lo >$@ || \
	echo "not using reloadable objects"

clean: clean-objs clean-objlists clean-libtool-dirs FORCE

distclean: clean-deps clean-mks FORCE

# $(REMOVE) *.lo may remove both the .lo and .o files.
# Do it before removing .o files for speed.
clean-objs: FORCE
	@for dir in $(SUBDIRS); do \
	    $(REMOVE) $$dir/*.lo; \
	    $(REMOVE) $$dir/*.o; \
	done

clean-objlists: FORCE
	@for dir in $(SUBDIRS); do \
	    $(REMOVE) $$dir/all_lo $$dir/all_lo.*; \
	done

clean-deps: FORCE
	@$(set_show); for dir in $(SUBDIRS); do \
	    $$show $(RM) $$dir/\*.d; \
	    $(RM) $$dir/*.d; \
	done

clean-mks: clean-mk-compile FORCE

clean-mk-compile: FORCE
	$(RM) graphics/depends.mk
	$(RM) graphics/objects.mk

clean-libtool-dirs: FORCE
	@$(set_show); for dir in $(SUBDIRS); do \
	    $$show $(RM_R) "$$dir/.libs" "$$dir/_libs"; \
	    $(RM_R) "$$dir/.libs" "$$dir/_libs"; \
	done

FORCE:

.PHONY: FORCE
