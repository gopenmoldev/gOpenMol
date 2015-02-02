##############################################################################
#                      Copyright (c) 2003 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

install: install-local-local install-plugins install-shared FORCE ;

install-local: install-local-local install-plugins-local FORCE ;

# Install in place to `$(top_srcdir)/..'.
altinstall: FORCE
	@rootdir=`cd '$(rootdir)' && pwd`; \
	$(MAKE) $(MFLAGS) \
	 	install-local-generated \
	 	install-plugins \
	 	gomroot="$$rootdir" gomdataroot="$$rootdir"

patch = ( eval "`$(LIBTOOL) --config | grep shlibpath_var`"; \
	sed \
	-e 's/^.INS *//' \
	-e 's%^TCL_LIBRARY=.*%TCL_LIBRARY='"'"'$(TCL_LIBRARY)'"'%" \
	-e 's%^TK_LIBRARY=.*%TK_LIBRARY='"'"'$(TK_LIBRARY)'"'%" \
	-e 's%^GOM_ROOT=.*%GOM_ROOT='"'"'$(gomroot)'"'%" \
	-e 's%@PATH_SEPARATOR@%$(PATH_SEPARATOR)%' \
	-e "s%@SHLIBPATHVAR@%$${shlibpath_var}%g" )

mkinstalldirs: mkinstalldirs-local mkinstalldirs-shared FORCE ;

# The first MKDIR_P may fail in case $(gomroot) and $(gomdataroot) are
# equal and -j is given.
mkinstalldirs-local: FORCE
	$(MKDIR_P) '$(gomroot)' || test -d '$(gomroot)'
	@$(set_show); set -e; \
	for dir in bin src/plugins src/unix; do \
	    dir='$(gomroot)'/$$dir; \
	    $$show $(MKDIR_P) "$$dir"; \
	    $(MKDIR_P) "$$dir" || test -d "$$dir"; \
	done

mkinstalldirs-shared: FORCE
	$(MKDIR_P) '$(gomdataroot)' || test -d '$(gomdataroot)'
	@$(set_show); set -e; \
	for dir in utility; do \
	    dir='$(gomdataroot)'/$$dir; \
	    $$show $(MKDIR_P) "$$dir"; \
	    $(MKDIR_P) "$$dir" || test -d "$$dir"; \
	done

install-local-generated: install-local-bin \
                         install-local-lib \
                         install-local-install \
                         install-local-utility \
                         FORCE ;

# Install best of gom-{dynamic,shared,noshared,static} to bin/.
install-local-bin: mkinstalldirs-local FORCE
	for file in gom-dynamic gom-shared gom-noshared gom-static; do \
	    if test -f graphics/$$file; then \
		exec $(LTINSTALL) \
		    graphics/$$file '$(gomroot)/bin/gopenmol'; \
	    fi; \
	done; exit 1

# Install libgopenmmol.la if it exists.
install-local-lib: mkinstalldirs-local FORCE
	$(LTPREINSTALL) '$(gomroot)/lib'
	test ! -f graphics/libgopenmol.la || \
	    $(LTINSTALL) graphics/libgopenmol.la '$(gomroot)/lib'
	$(LTPOSTINSTALL) '$(gomroot)/lib'

# Install GOMROOT/install.
install-local-install: mkinstalldirs-local FORCE
	$(INSTALL_SCRIPT) '$(rootdir)/install.tcl' '$(gomroot)/install'
	$(patch) <'$(rootdir)/install.tcl' >'$(gomroot)/install'

# Install utilities.
install-local-utility: FORCE
	cd utility && $(MAKE) $(MFLAGS) install \
		gomroot='$(gomroot)' gomdataroot='$(gomdataroot)'

install-local-local: install-local-generated \
                     install-local-environment \
                     install-local-rungopenmol \
                     install-local-data \
                     FORCE ;

# Copy GOMROOT/environment.txt.
install-local-environment: mkinstalldirs-local FORCE
	$(INSTALL_DATA) \
		'$(rootdir)/environment.txt' \
		'$(gomroot)/environment.txt'

# Install GOMROOT/bin/rungOpenMol.
install-local-rungopenmol: mkinstalldirs-local FORCE
	$(MKDIR_P) '$(exec_prefix)/bin'
	rootdir=`cd '$(rootdir)' && pwd` && \
	cd '$(exec_prefix)' && $(TCLSH) "$$rootdir/install.tcl"
	$(RM) graphics/rungOpenMol
	sed 's/x/x/' <'$(exec_prefix)/bin/rungOpenMol' >graphics/rungOpenMol
	( \
	echo '#! /bin/sh'; \
	echo "dataroot='"'$(gomdataroot)'"'"; \
	echo 'GOM_DATA="$$dataroot/data"'; \
	echo 'GOM_HELP="$$dataroot/help"'; \
	echo 'GOM_DEMO="$$dataroot/demo"'; \
	echo 'export GOM_DATA GOM_HELP GOM_DEMO'; \
	$(patch) <graphics/rungOpenMol ) >'$(exec_prefix)/bin/rungOpenMol'
	$(RM) graphics/rungOpenMol

install-local-data: install-local-data-include \
                    install-local-data-src \
                    FORCE ;

# Install header files:
#     src/include/gomlib     copy from builddir and from srcdir
install-local-data-include: FORCE
	$(MAKE) $(MFLAGS) install-data-dir \
		srcdir='$(top_srcdir)/include/gomlib' \
		dstdir='$(gomincdir)'
	$(MAKE) $(MFLAGS) install-data-dir \
		srcdir='$(top_builddir)/include/gomlib' \
		dstdir='$(gomincdir)'

# Install to system specific directories:
#     src/plugins/nodeps.d   copy from srcdir
#     src/plugins/rules.mk   copy from srcdir
#     src/plugins/config.mk  patch
#     src/unix/depcomp       copy from srcdir
#     src/unix/install-sh    copy from srcdir
#     src/libtool            copy from builddir
install-local-data-src: mkinstalldirs-local FORCE
	@$(set_show); $$show 'Installing $(gomroot)/src/plugins/config.mk'
	@$(INSTALL_DATA) \
		'$(top_builddir)/plugins/config.mk' \
		'$(gomroot)/src/plugins/config.mk'
	@$(patch) <'$(top_builddir)/plugins/config.mk' \
		>'$(gomroot)/src/plugins/config.mk'
	$(INSTALL_DATA) \
		'$(top_srcdir)/plugins/nodeps.d' \
		'$(gomroot)/src/plugins/nodeps.d'
	$(INSTALL_DATA) \
		'$(top_srcdir)/plugins/rules.mk' \
		'$(gomroot)/src/plugins/rules.mk'
	$(INSTALL_SCRIPT) \
		'$(rootdir)/src/unix/depcomp' \
		'$(gomroot)/src/unix/depcomp'
	$(INSTALL_SCRIPT) \
		'$(rootdir)/src/unix/install-sh' \
		'$(gomroot)/src/unix/install-sh'
	$(INSTALL_SCRIPT) '$(top_builddir)/libtool' '$(gomroot)/src/libtool'

# Install plugins.
install-plugins: FORCE
	cd plugins && $(MAKE) $(MFLAGS) install \
		gomroot='$(gomroot)' gomdataroot='$(gomdataroot)'

install-plugins-local: FORCE
	cd plugins && $(MAKE) $(MFLAGS) install-local \
		gomroot='$(gomroot)' gomdataroot='$(gomdataroot)'

install-shared: install-shared-data FORCE ;

# Install COPYRIGHT and README files and data, demo, docs and help
# directories to shared data directory.
install-shared-data: mkinstalldirs-shared FORCE
	@$(set_show); set -e; \
	for file in `(cd '$(rootdir)' && \
	    echo COPYRIGHT  README* utility/*.html)`; do \
	    $$show $(INSTALL_DATA) '$(rootdir)'/"$$file" '$(gomdataroot)'/"$$file"; \
	    $(INSTALL_DATA) '$(rootdir)'/"$$file" '$(gomdataroot)'/"$$file"; \
	done
	@set -e; for dir in data demo docs help; do \
	    $(MAKE) $(MFLAGS) install-data-dir \
		srcdir='$(rootdir)'/"$$dir" dstdir='$(gomdataroot)'/"$$dir";\
	done

# Install data files recursively.
# .so files are treated as libraries.
skip = ////
install-data-dir: FORCE
	$(MKDIR_P) '$(dstdir)'
	@$(set_show); set -e; \
	cwd=`pwd`; dirs='.'; src='$(srcdir)'; dst='$(dstdir)'; \
	while test -n "$$dirs"; do \
	    IFS=':'; set x $$dirs; shift; dir="$$1"; shift; dirs="$$*"; \
	    IFS=' 	'; test -n "$$dir" || continue; \
	    cd "$$src/$$dir"; set x * .*; shift; cd "$$cwd"; \
	    $$show $(MKDIR_P) "$$dst/$$dir"; \
	    $(MKDIR_P) "$$dst/$$dir"; \
	    for file in "*" "$$@"; do \
		case "$$dir/$$file" in \
		*/CVS*|*/.|*/..|*/"*"|*/".*"|*/'$$@'|$(skip)) continue ;; \
		*.so|*.so.*) set $(INSTALL_PROGRAM) ;; \
		*) \
		    if test -d "$$src/$$dir/$$file"; then \
			dirs="$$dir/$$file:$$dirs"; continue; \
		    fi; \
		    set $(INSTALL_DATA) ;; \
		esac; \
		$$show $$* "$$src/$$dir/$$file" "$$dst/$$dir/$$file"; \
		eval $$* '"$$src/$$dir/$$file" "$$dst/$$dir/$$file"'; \
	    done; \
	done

uninstall: uninstall-local FORCE
	$(RM_R) '$(gomdataroot)'

uninstall-local: FORCE
	$(RM) '$(bindir)/rungOpenMol'
	$(RM_R) '$(gomroot)'

FORCE:

.PHONY: FORCE
