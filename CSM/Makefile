mkllib=/opt/intel/mkl/lib/intel64
mklinc=/opt/intel/mkl/include
TAR = ~/programs/CSM/CSM
SRCS = $(wildcard *.f90) #obtain every .f90 file in the present directory
OBJ = $(patsubst %.f90,%.o,$(SRCS)) #replace .f90 with .o and obtain every .o file
FCCFLAG= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread
FCC = ifx

$(TAR): $(OBJ)
	$(FCC) -o $(TAR) -I$(mklinc) -L$(mkllib) $(FCCFLAG) -g $^ #$^ means "every dependent file" -O0

%.o:%.f90 %.mod
	$(FCC) -c *.f90 -g -O0

%.mod:%.f90
	$(FCC) -c *.f90 -g -O0

clean:
	rm -rf *.o *.mod
