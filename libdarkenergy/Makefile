include ../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS} ${DEFINEFLAGS} $(EXTRA_LIBS) $(EXTRA_INCL)

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS   = darkenergy.o

all: libdarkenergy.a

libdarkenergy.a: ${OBJS}
		${AR} -r libdarkenergy.a ${OBJS}

clean: 
	rm -rf *.o *.a
