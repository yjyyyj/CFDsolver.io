MARK = 4
FC  = gfortran
#FC=ifort
DEBUG = -Wall -fbounds-check -Wuninitialized -ffpe-trap=invalid,zero,denormal -fbacktrace -g 
OPT = -O3
# OPT = -O ${DEBUG}
LDFLAGS = -fopenmp
FFLAGS =-fconvert=little-endian  -frecord-marker=${MARK} ${LDFLAGS} ${OPT}

OBJS = param.o flux.o scheme.o main.o outf.o init.o  bc.o sumdf.o muscl.o step.o visflux.o read_data.o
PROGRAM = ../job.out

all:	$(OBJS)
	# $(FC) -o $(PROGRAM) $(FFLAGS) $(OBJS) -lblas -llapack
	$(FC) -o $(PROGRAM) $(FFLAGS) $(OBJS)
%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	rm -f *.out *.o *.mod
