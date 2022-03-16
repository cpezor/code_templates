## Compiler
CC=g++

##	Flags for C++
CFLAGS=-I. -lm

##	Modules (.cpp)
SRC=<++>.cpp

##	Dependencys (.hpp)
DEPS=<++>.hpp

##	Obj, compilated
OBJ=${SRC:.cpp=.o}

##	Any .cpp, transform to .o
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

##	Linking the program
main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

#	Clean up the file (explicit)
clean:
	@rm -f *.o *.mod main
