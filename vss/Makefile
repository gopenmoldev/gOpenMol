###############################################################################
# MAKEFILE
#	Make file for vss

CFLAGS =  -O2  
DEFS = 

OBJS = vssmod.o


vss : $(OBJS) 
	cc $(DEFS) -o ../bin/vss  $(OBJS) -lm 

vssmod.o : vssmod.c
	cc $(DEFS) $(CFLAGS) -c  vssmod.c

