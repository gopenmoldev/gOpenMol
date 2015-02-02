###############################################################################
# MAKEFILE
#	Make file for the SGI  version vss

CFLAGS =  -O2  
DEFS = 

OBJS = probsurf.o


probsurf : $(OBJS) 
	cc $(DEFS) -o ../bin/probsurf  $(OBJS) -lm 

probsurf.o : probsurf.c
	cc $(DEFS) $(CFLAGS) -c  probsurf.c

