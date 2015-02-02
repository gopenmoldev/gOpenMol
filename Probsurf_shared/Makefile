###############################################################################
# MAKEFILE
#	Make file for probsurf

CFLAGS =  -O2  
DEFS = 

OBJS = probsurf2.0.o


probsurf : $(OBJS) 
	cc $(DEFS) -o ../bin/probsurf  $(OBJS) -lm 

probsurf.o : probsurf.c
	cc $(DEFS) $(CFLAGS) -c  probsurf.c

