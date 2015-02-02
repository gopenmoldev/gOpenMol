##############################################################################
#                      Copyright (c) 2003 - 2004 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

# Extend all.
all: build-all

build-all:        build build-plugins build-utilities FORCE
build-all-bin:    build-all-$(BUILDTARGET) FORCE ;
build-all-dynamic build-all-static build-all-noshared: \
			build-dynamic build-static build-noshared FORCE ;
build-all-shared: build-shared FORCE ;

build:          build-$(BUILDTARGET) FORCE
build-dynamic:  graphics/gom-dynamic$(EXEEXT)  FORCE ;
build-static:   graphics/gom-static$(EXEEXT)   FORCE ;
build-noshared: graphics/gom-noshared$(EXEEXT) FORCE ;
build-shared:   graphics/gom-shared$(EXEEXT)   FORCE ;

build-plugins:  build-plugins-$(BUILDTARGET) FORCE ;
# Plugins need header files.
build-plugins-dynamic build-plugins-noshared build-plugins-shared \
plugins: graphics/stamp-gomapi FORCE
	cd plugins && $(MAKE) $(MFLAGS) build

build-utilities utilities: FORCE
	cd utility && $(MAKE) $(MFLAGS) build

# Build a single binary which exports API entries for plugins.
graphics/gom-dynamic$(EXEEXT):  $(top_srcdir)/graphics/main/main.c \
                                graphics/libgom.a \
                                graphics/gomapi.exp
	$(LINK) $@ $(top_srcdir)/graphics/main/main.c graphics/libgom.a \
		$(LINKFLAGS) -export-dynamic \
		-export-symbols graphics/gomapi.exp

# Build a single binary which doesn't exports API entries for plugins.
graphics/gom-static$(EXEEXT):   $(top_srcdir)/graphics/main/main.c \
                                graphics/libgom.a
	$(LINK) $@ $(top_srcdir)/graphics/main/main.c \
		graphics/libgom.a $(LINKFLAGS)

# Don't use libraries or reloadable objects at all.
graphics/gom-noshared$(EXEEXT): $(top_srcdir)/graphics/main/main.c \
                                $(OBJECTS)
	$(LINK) $@ $(top_srcdir)/graphics/main/main.c \
		$(OBJECTS) $(LINKFLAGS)

# Use a shareable object. Executable contains only the main function.
graphics/gom-shared$(EXEEXT):   $(top_srcdir)/graphics/main/main.c \
                                graphics/libgopenmol.la
	$(LINK) $@ $(top_srcdir)/graphics/main/main.c \
		graphics/libgopenmol.la $(LINKFLAGS)

# Create shareable object from object files.
graphics/libgopenmol.la: $(OBJLISTS) graphics/gomapi.exp
	$(LINK_SO) $@ `echo $(OBJLISTS) | \
		sed 's/[^ ][^ ]*/-objectlist &/g'` \
		$(LINKFLAGS_SO) -export-symbols graphics/gomapi.exp

# Create static library from object files.
graphics/libgom.a: $(OBJLISTS)
	$(RM) $@
	$(LINK_A) $@ `echo $(OBJLISTS) | \
		sed 's/[^ ][^ ]*/-objectlist &/g'` $(LINKFLAGS_A)

# Extend clean.
clean: clean-bins clean-libs FORCE

clean-bins: FORCE
	$(REMOVE) graphics/gom-*
	cd utility && $(MAKE) $(MFLAGS) clean

clean-libs: FORCE ;
	$(REMOVE) graphics/libgopenmol.*
	$(REMOVE) graphics/libgom.*
	cd plugins && $(MAKE) $(MFLAGS) clean

FORCE:

.PHONY: FORCE
