##############################################################################
#                       Copyright (c) 2003 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

LIBRARY    = lib$(PLUGINNAME).la
objdeps    = FORCE
libdeps    = FORCE
depfile    = $(gomsrcdir)/plugins/nodeps.d

all: build FORCE ;

include $(depfile)

build: $(LIBRARY) FORCE ;

compile: $(OBJECTS) FORCE

# Rule to link.
# Top-level make:
#  * libdeps is FORCE (FORCEd run).
#  * Concatenate existing dependency files, set depencies.
#  * Run sub make using the same target.
# Sub make:
#  * libdeps is the list of objects.
#  * Make checks whether the library and the objects are up to date.
#  * Run sub sub make using the target link-$(LANG).
$(LIBRARY): $(libdeps)
	@objects='$(OBJECTS)'; \
	case '$(libdeps)' in \
	FORCE) \
	    cat `echo "\$$objects " | \
		sed 's/\.[^.]* /.d /g'` 2>/dev/null >alldeps.dT || :; \
	    exec $(MAKE) $(MFLAGS) '$@' \
		libdeps="$$objects" objdeps= depfile=alldeps.dT ;;\
	*) exec $(MAKE) $(MFLAGS) link-$(LANG) ;; \
	esac
	@$(RM) alldeps.dT

link-CC: FORCE
	$(GOM_LINK_CC) $(PLUG_CFLAGS) -o $(LIBRARY) $(OBJECTS) \
	    $(PLUG_LDFLAGS) $(PLUG_LIBS) $(GOM_LINKFLAGS)

link-CXX: FORCE
	$(GOM_LINK_CXX) $(PLUG_CXXFLAGS) -o $(LIBRARY) $(OBJECTS) \
	    $(PLUG_LDFLAGS) $(PLUG_LIBS) $(GOM_LINKFLAGS)

# Rule to compile the source files.
# First make:
#  * objdeps is FORCE (FORCEd run).
#  * Check if dependency file exists. Set depencies.
#  * Run sub make using the same target.
# Sub make:
#  * objdeps is empty (only dependencies in dependency file are used).
#  * Make checks whether the object is up to date.
#  * Find the source file.
#  * Run sub sub make using the target compile-<source_file_lang>.
$(OBJECTS): $(objdeps)
	@file="`echo '$@' | sed 's:\.$(OBJSUFFIX)::'`"; \
	case '$(objdeps)' in \
	FORCE)	test ! -f "$$file.d" || depfile="depfile=$$file.d"; \
		exec $(MAKE) $(MFLAGS) '$@' objdeps= $$depfile ;; \
	*) \
	    for src in $$file.* '${gomsrcdir}/plugins/${DIRNAME}'/$$file.*; do\
		case $$src in \
		*.c) lang=CC ;; \
		*.cc | *.cpp) lang=CXX ;; \
		*) continue ;; \
		esac; \
		exec $(MAKE) $(MFLAGS) compile-$$lang obj='$@' src="$$src"; \
	    done; echo "no source file for $@" >&2; exit 1 ;; \
	esac

# Compile source file written in C and optionally generate dependency file.
compile-CC: FORCE
	source='$(src)'; \
	object='$(obj)'; \
	file=`echo "\$$object" | sed 's/\.[^.]*\$$//'`; \
	source="$$source" object="$$object" depfile="$$file.d" libtool=yes \
	$(CCDEPCOMP) \
	$(GOM_COMPILE_CC) $(PLUG_CPPFLAGS) $(PLUG_CFLAGS) $(PLUG_DEFS) \
	    -c -o "$$object" "$$source"

# Compile source file written in C++ and optionally generate dependency file.
compile-CXX: FORCE
	source='$(src)'; \
	object='$(obj)'; \
	file=`echo "\$$object" | sed 's/\.[^.]*\$$//'`; \
	source="$$source" object="$$object" depfile="$$file.d" libtool=yes \
	$(CXXDEPCOMP) \
	$(GOM_COMPILE_CXX) $(PLUG_CPPFLAGS) $(PLUG_CXXFLAGS) $(PLUG_DEFS) \
	    -c -o "$$object" "$$source"

# Rules to install.
install: install-lib install-tcl FORCE ;

install-local: install-lib FORCE ;

install-lib: FORCE
	$(MKDIR_P) '$(gompluginsdir)' || test -d '$(gompluginsdir)'
	$(GOM_PREINSTALL)  $(gompluginsdir)
	test ! -f $(LIBRARY) || \
	    $(GOM_INSTALL) $(LIBRARY) $(gompluginsdir)
	$(GOM_POSTINSTALL) $(gompluginsdir)

install-tcl: FORCE
	$(MKDIR_P) '$(gompluginstcldir)' || test -d '$(gompluginstcldir)'
	for dir in . '${gomsrcdir}/plugins/${DIRNAME}'; do \
	    file="$$dir/$(PLUGINNAME).tcl"; \
	    test ! -f "$$file" || \
		$(GOM_INSTALL_DATA) "$$file" $(gompluginstcldir); \
	done

# Rule to uninstall.
uninstall:
	$(GOM_REMOVE) $(gompluginsdir)/lib$(PLUGINNAME).*
	$(RM) $(gompluginstcldir)/$(PLUGINNAME).tcl

# Rules to clean up.
clean: clean-libs clean-objs FORCE ;
	$(RM_R) .libs _libs

clean-libs: FORCE
	$(GOM_REMOVE) *.a *.la

clean-objs: FORCE
	$(GOM_REMOVE) $(OBJECTS)
	$(GOM_REMOVE) *.o *.lo
	$(RM) `echo '$(OBJECTS) ' | sed 's/\.[^.]* /.d /g'`
	$(RM) *.d *.dT

distclean: clean FORCE ;

maintainer-clean: clean FORCE ;

FORCE:

.PHONY: FORCE
