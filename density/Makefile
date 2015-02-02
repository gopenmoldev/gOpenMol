###############################################################################
# MAKEFILE
#	Make file for the SGI  version density

CFLAGS =  -O2  
DEFS = 

OBJS = density.o


density : $(OBJS) 
	cc $(DEFS) -o ../bin/density  $(OBJS) -lm 

density.o : density.c
	cc $(DEFS) $(CFLAGS) -c  density.c

