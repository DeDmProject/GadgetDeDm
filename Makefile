include Makefile.config

OPTIONS =  $(OPTIMIZE) $(OPT)

OBJS   = main.o run.o  predict.o begrun.o endrun.o global.o  \
	 timestep.o  init.o restart.o  io.o  \
	 accel.o   read_ic.o  ngb.o \
	 system.o  allocate.o  density.o  \
	 gravtree.o hydra.o  driftfac.o  \
	 domain.o  allvars.o potential.o \
         forcetree.o  peano.o gravtree_forcetest.o \
	 pm_periodic.o pm_nonperiodic.o longrange.o 

INCL   = allvars.h  proto.h  forcetree.h  domain.h 

LIBS   = $(HDF5LIB) -g $(MPICHLIB)  $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIB) -lm 

ifeq (${DEDM},"on")
INCL    += libdedm/dedmvars.h libdedm/interpolate.h libdedm/read_tables.h libdedm/mass.h libdedm/integrate.h
INCL    += libscott/spline.h libscott/splineimpl.h libscott/read_scott_tables.h libscott/subspline.h
MLIBS   += libdedm/libdedm.a
MLIBS   += libscott/libscott.a
endif

ifeq (${DEONLY},"on")
INCL    += libdedm/dedmvars.h libdedm/interpolate.h libdedm/read_tables.h
MLIBS   += libdedm/libde.a
endif

ifeq (${VDE},"on")
INCL    += libvde/interpolate.h libvde/read_tables.h libvde/vdevars.h
MLIBS   += libvde/libvde.a
endif

ifeq (${EXTRA_BARYONS}, "on")
INCL	+= libbaryons/cooling.h libbaryons/c_metals.h libbaryons/cosmic_rays.h libbaryons/chemistry.h
MLIBS   += libbaryons/libbaryons.a
endif

ifeq (${DARK_ENERGY}, "on")
INCL	+= libdarkenergy/darkenergy.h
MLIBS   += libdarkenergy/libdarkenergy.a
endif

LIBS   += $(MLIBS)
SOURCES = $(OBJS:.o = .c)

INCL	+= Makefile 

CFLAGS =   $(OPTIONS)  $(GSL_INCL) $(FFTW_INCL)  $(HDF5INCL)

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))    # fftw installed with type prefix?
  FFTW_LIB = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif

$(EXEC): $(SOURCES) $(MLIBS)
	$(CC) $(OBJS) $(LIBS) -o  $(EXEC) $(MLIBS) $(MLIBS)

$(OBJS): $(INCL) 

pm_periodic.o:   pm_hpm.c

clean:
	rm -f $(OBJS) $(EXEC) *~
	cd libbaryons ; ${MAKE} clean
	cd libdarkenergy ; ${MAKE} clean
	cd libdedm ; ${MAKE} clean
	cd libscott ; ${MAKE} clean
	cd libvde ; ${MAKE} clean

libbaryons/libbaryons.a:
	cd libbaryons/ ; ${MAKE} libbaryons.a

libdarkenergy/libdarkenergy.a: 
	cd libdarkenergy/ ; ${MAKE} libdarkenergy.a

libdedm/libdedm.a: 
	cd libdedm/ ; ${MAKE} libdedm.a

libdedm/libde.a: 
	cd libdedm/ ; ${MAKE} libde.a

libscott/libscott.a: 
	cd libscott/ ; ${MAKE} libscott.a

libvde/libvde.a: 
	cd libvde/ ; ${MAKE} libvde.a
