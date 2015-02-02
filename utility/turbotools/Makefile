#
#  Makefile for tmol2plt
#
#  Written by Jonas Juselius, 1998
#
#

CC=gcc
CFLAGS=-O2 

.PHONY: all 
all: tmole2plt

tmole2plt: tmole2plt.o

install: 
	-cp tmole2plt $(HOME)/bin/$(HOSTTYPE)

clean:
	-rm -f tmole2plt tmole2plt.o core 
