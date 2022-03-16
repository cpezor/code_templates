## Compiler
CC=gcc

##	Flags for C
CFLAGS=-I. -lm

##	Modules (.c)
SRC=<++>.c

##	Dependencys (.h)
DEPS = 

##	Obj, compilated
OBJ=${SRC:.c=.o}

##	Any .f90, transform to .o
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

##	Linking the program
main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

#	Clean up the file (explicit)
clean:
	@rm -f *.o *.mod main
