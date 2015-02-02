###############################################################################
# MAKEFILE
#	Make file for the SGI  version of icon8

FFLAGS =  -O -Olimit 800 -static 
DEFS =  
CFLAGS = -O

OBJS = icon8.o \
       output.o \
       fmalloc.o \
       rsp.o

icon8 : $(OBJS) 
	f77 $(FFLAGS) $(DEFS) -o ../bin/icon8  $(OBJS) 

icon8.o : icon8.f 
	f77 $(DEFS) $(FFLAGS) -c  icon8.f

output.o : output.f
	f77 $(DEFS) $(FFLAGS) -c  output.f

rsp.o : rsp.f
	f77 $(DEFS) $(FFLAGS) -c  rsp.f

fmalloc.o : fmalloc.c
	cc $(DEFS) $(CFLAGS) -c   fmalloc.c

