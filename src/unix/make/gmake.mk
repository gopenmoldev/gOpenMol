##############################################################################
#                           Copyright (c) 2003 by:
#          Leif Laaksonen, CSC - Scientific Computing Ltd, ESPOO, FINLAND 
#                      Confidential unpublished property of 
#                              Leif Laaksonen   
#                            All rights reserved
#
#          Coded by: Eero Häkkinen
##############################################################################

# Retrieve object files in all gOpenMol directories.
OBJECTS  := \
	$(patsubst $(top_srcdir)/%.c,%.lo,\
		$(foreach dir,\
			$(SUBDIRS),\
			$(wildcard $(top_srcdir)/$(dir)/*.c)))
OBJLISTS := $(patsubst %,%/all_lo,$(SUBDIRS))
DEPENDS  := $(patsubst %.lo,%.d,$(OBJECTS))

# Include dependencies of object list files.
%/all_lo.mk:
	echo '$*/all_lo: $$(patsubst $$(top_srcdir)/%.c,%.lo,$$(wildcard $$(top_srcdir)/$*/*.c))' >$@

-include $(patsubst %,%/all_lo.mk,$(SUBDIRS))

# Include dependency files.
-include $(DEPENDS)
-include graphics/gomapi.mk

# Update API dependencies then new header files are created.
graphics/gomapi.mk: $(incdir)

# Make libtool follow the silence of make.
ifneq (,$(findstring s,$(MAKEFLAGS)))
LIBTOOL += --silent
endif
