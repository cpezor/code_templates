## Compiler
FC=gfortran

##	Flags for fortran
FFLAGS= -Wall -Wextra -std=f2008

##	Modules
SRC=

# SRC=tools.f90 bidiagonal.f90 main.f90

##	Obj, compilated
OBJ=${SRC:.f90=.o}

##	Any .f90, transform to .o
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

##	Linking the program
main: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

#	Clean up the file (explicit)
clean:
	@rm -f *.o *.mod main
