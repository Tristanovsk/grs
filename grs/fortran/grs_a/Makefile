
FC     = gfortran
FFLAGS = -m64 -fpic -lm -g -O5 -Wall
VPATH = interp:minpack:


all: main_algo binding

#---------------------------------------
MINPACK= enorm.o  lmdif.o  qrsolv.o  dpmpar.o fdjac2.o  lmpar.o  qrfac.o
INTERP = rgrd1.o rgrd2.o rgrd3.o interp4d.o
module = common.o lsq.o solver.o
SRC = main_algo.o

main_algo: $(MINPACK) $(INTERP) $(module) $(SRC) 
	$(FC) -c $(MINPACK) $(INTERP) $(module) $(SRC) 

binding:
	f2py -c -m main_algo main_algo.f95 $(MINPACK) $(INTERP) $(module) --opt='$(FFLAGS)'

clean :
	rm -rf *.o *.mod


%.o:%.f95
	$(FC) $(FFLAGS) -c $< 
%.o:%.f90
	$(FC) $(FFLAGS) -c $< 
%.o:%.f
	$(FC) $(FFLAGS) -c $<

