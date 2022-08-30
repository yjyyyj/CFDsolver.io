MARK = 4
FC  = gfortran
#FC=ifort
# DEBUG = -Wall -fbounds-check -Wuninitialized -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace -g
DEBUG = -Wall -fbounds-check -Wuninitialized -ffpe-trap=invalid,zero,denormal -fbacktrace -g  # removed (overflow,underflow)
OPT = -O3
OPT = -O ${DEBUG}
# FFLAGS =-fconvert=little-endian  -frecord-marker=${MARK} ${OPT}
FFLAGS =-fconvert=big-endian  -frecord-marker=${MARK} ${OPT}
# FFLAGS2=${FFLAGS}

OBJS = module.o flux.o scheme.o main.o init.o outf.o bc.o sumdf.o muscl.o step.o visflux.o
PROGRAM = ../job.out

all:	$(OBJS)
	# $(FC) -o $(PROGRAM) $(FFLAGS) $(OBJS) -lblas -llapack
	$(FC) -o $(PROGRAM) $(FFLAGS) $(OBJS)
%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	rm -f *.out *.o *.mod
