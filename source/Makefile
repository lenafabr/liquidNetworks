OBJS = minnetutil.o keys.o readkey.o genutil.o mt19937.o sort2int.o inputparams.o browndyn.o stackutil.o

FC = gfortran
FLAGS = -O3 
# DEBUGFLAGS = -march=native -ffast-math -funroll-loops
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,invalid -fcheck=all -fbounds-check
LIBS =  -lblas -llapack

# executable name
EXE = ../minnetBD.exe
TESTEXE = ../testrun.exe

.SUFFIXES: .o .f90 

# general suffix rules
%.o: %.f90 
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

%.o: tools/%.f*
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

clean: 
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
minnetutil.o: keys.o mt19937.o genutil.o inputparams.o stackutil.o
readkey.o: inputparams.o keys.o genutil.o mt19937.o
genutil.o: mt19937.o 
browndyn.o: mt19937.o minnetutil.o
