# Makefile for libscott.a , a GADGET-2 library used to implement correct reading of input 
# scalar field values
include ../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS} ${DEFINEFLAGS}

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS    = spline.o read_scott_tables.o splineimpl_gsl.o subspline.o

all: libscott.a

libscott.a: ${OBJS} 
	${AR} -r libscott.a ${OBJS}

clean: 
	rm -rf *.o libscott.a
