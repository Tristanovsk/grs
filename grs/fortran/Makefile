
FC     = gfortran
FFLAGS = -m64 -fpic -lm -g -O5 -Wall
VPATH = interp:minpack:


MINPACK= enorm.o  lmdif.o  qrsolv.o  dpmpar.o fdjac2.o  lmpar.o  qrfac.o
INTERP = rgrd1.o rgrd2.o rgrd3.o interp4d.o
SRC = main_algo.o

all: $(MINPACK) $(INTERP) $(SRC) 
	$(FC)  -c $(MINPACK) $(INTERP) $(SRC) 
	f2py3 -c -m main_algo main_algo.f95 $(MINPACK) $(INTERP)  --opt='$(FFLAGS)'

clean :
	rm -rf *.o *.mod


%.o:%.f95
	$(FC) $(FFLAGS) -c $^ -o $@

%.o:%.f
	        $(FC) $(FFLAGS) -c $^ -o $@



