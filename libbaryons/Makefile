include ../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS} ${DEFINEFLAGS}

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS   = kinfb_mhm.o sfr_mhm.o fof.o  sfr_eff.o  \
	 cooling.o conduction.o \
	 bsmooth.o  \
	 read_ic_cluster_gas.o \
         read_ic_cluster.o

INCL   = ../proto.h cooling.h domain.h  c_metals.h cosmic_rays.h chemistry.h 

all: libbaryons.a

libbaryons.a: ${OBJS}
		${AR} -r libbaryons.a ${OBJS}

clean: 
	rm -rf *.o *.a
