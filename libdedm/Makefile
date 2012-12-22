#####################################################################################
# Makefile for libdedm.a , a GADGET-2 library used to enable simulations of coupled 
# dark energy - dark matter simulations						    
#										    
# Edoardo Carlesi 2010 - 2012							    
# edoardo.carlesi@uam.es							    
#####################################################################################
include ../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS} ${DEFINEFLAGS} $(EXTRA_LIBS) $(EXTRA_INCL)

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS    = dedmvars.o integrate.o interpolate.o mass.o read_tables.o
OBJS_DE = dedmvars.o interpolate.o read_tables.o

all: libdedm.a

libde.a: ${OBJS_DE}
	${AR} -r libde.a ${OBJS_DE}

libdedm.a: ${OBJS} 
	${AR} -r libdedm.a ${OBJS}

clean: 
	rm -rf *.o libdedm.a libde.a
