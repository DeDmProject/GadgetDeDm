# Makefile for libvde.a , a GADGET-2 library used to enable simulations of Vector Dark Energy models
include ../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS} ${DEFINEFLAGS} $(EXTRA_LIBS) $(EXTRA_INCL)

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS    = vdevars.o interpolate.o read_tables.o

all: libvde.a

libvde.a: ${OBJS}
	${AR} -r libvde.a ${OBJS}

clean: 
	rm -rf *.o libvde.a
